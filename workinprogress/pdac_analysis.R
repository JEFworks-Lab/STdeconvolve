## PDAC analysis building off of Brendan's

pdacAst1Path <- "~/Desktop/GSE111672_RAW/GSM3036911_PDAC-A-ST1-filtered.txt"
pdacAst1 <- read.table(pdacAst1Path, header = TRUE, sep = "\t")
head(pdacAst1)

spotIDs <- unlist(lapply(colnames(pdacAst1)[2:ncol(pdacAst1)], function(i) {
  ix <- strsplit(i, "X")[[1]][2]
  ix
}))
pdacAst1.Pos <- do.call(rbind, lapply(spotIDs, function(ix) {
  coords <- as.numeric(strsplit(ix, "x")[[1]])
  coords
}))
colnames(pdacAst1.Pos) <- c("x", "y")
rownames(pdacAst1.Pos) <- spotIDs

pdacAst1.t <- t(pdacAst1)
colnames(pdacAst1.t) <- pdacAst1$Genes
pdacAst1.t <- pdacAst1.t[2:ncol(pdacAst1),]
pdacAst1 <- apply(pdacAst1.t, 2,FUN = as.numeric)
rownames(pdacAst1) <- rownames(pdacAst1.Pos)
plot(pdacAst1.Pos[,"x"], pdacAst1.Pos[,"y"], col='red', pch="+")

############## Recapitulate previous analysis

## Remove poor genes (seems like datasets already cleaned)
pdacAst1.Clean <- MERINGUE::cleanCounts(counts = t(pdacAst1), # genes x spots mtx
                                        min.reads = 100,
                                        min.lib.size = 0,
                                        verbose=TRUE)

# CPM normalize
dim(pdacAst1.Clean)
pdacAst1.CPM <- MERINGUE::normalizeCounts(counts = pdacAst1.Clean,
                                          verbose=TRUE, log=FALSE)

## Get overdispsersed genes
counts <- pdacAst1.Clean
odGenes <- MERINGUE::getOverdispersedGenes(counts, plot = TRUE, details = TRUE, alpha = 0.01)
countsFilt <- counts[odGenes$ods,]

## remove genes that are present in more than X% of spots
## if we don't remove perplexity poor
## also remove genes that are too rare
library(Matrix)
vi <- rowSums(countsFilt > 0) > ncol(countsFilt)*0.1 & rowSums(countsFilt > 0) < ncol(countsFilt)*0.9
table(vi)
countsFiltnoUnifGenes <- countsFilt[vi,]

dim(countsFiltnoUnifGenes)
hist(log10(colSums(countsFiltnoUnifGenes)+1))
hist(log10(rowSums(countsFiltnoUnifGenes)+1))

# Dimensionality reduction by PCA on log10 CPM expression values
pcs.info <- prcomp(t(as.matrix(pdacAst1.CPM[rownames(countsFiltnoUnifGenes),])), center=TRUE)
plot(pcs.info$sdev[1:30])
nPcs <- 15
pcs <- pcs.info$x[,1:nPcs]
# 2D embedding by tSNE
emb <- Rtsne::Rtsne(pcs,
                    is_distance=FALSE,
                    perplexity=30,
                    num_threads=1,
                    verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)
# Graph-based cluster detection
k <- 30
com <- MERINGUE::getClusters(pcs, k, method=igraph::cluster_louvain)

par(mfrow=c(1,2), mar=rep(1,4))
MUDAN::plotEmbedding(emb, groups=com,
              show.legend=TRUE, xlab=NA, ylab=NA,
              verbose=FALSE)
MUDAN::plotEmbedding(pdacAst1.Pos, groups=com,
              cex=1, xlab=NA, ylab=NA,
              verbose=FALSE)


############################### Topic Modeling
pdac_corpus <- t(as.matrix(countsFiltnoUnifGenes))
pdac_corpus_slamMtx <- slam::as.simple_triplet_matrix(pdac_corpus)

# optimal k?
#Ks = seq(10,100,by=10)
#Ks
#pdac_lda <- fitLDA(pdac_corpus_slamMtx, Ks=Ks)
#par(mfrow=c(1,1), mar=rep(5,4))
#plot(x=Ks, pdac_lda$perplexities)

# pick k
#kopt = Ks[which(pdac_lda$perplexities == min(pdac_lda$perplexities))]
kopt = 30
pdac_lda.k <- topicmodels::LDA(pdac_corpus_slamMtx, k=kopt)
topicmodels::perplexity(pdac_lda.k, pdac_corpus_slamMtx)

pdac_lda.tmResult<- topicmodels::posterior(pdac_lda.k)
pdac_lda.theta <- pdac_lda.tmResult$topics
pdac_lda.beta <- pdac_lda.tmResult$terms
pdac_lda.topicFreqsOverall <- colSums(pdac_lda.theta) / nrow(pdac_corpus_slamMtx)

library(dynamicTreeCut)
pdac_lda.clust <- clusterTopics(beta = pdac_lda.beta,
                                distance = "euclidean",
                                clustering = "complete",
                                dynamic = "hybrid",
                                deepSplit = 4,
                                plotDendro = TRUE)

clusterColorsk50 <- pdac_lda.clust$clusters
levels(clusterColorsk50) <- gg_color_hue(length(levels(clusterColorsk50)))

m <- t(pdac_lda.beta*1e6)
range(m)
m <- scale(m)
m <- t(scale(t(m)))
range(m)
m[m < -2] <- -2
m[m > 2] <- 2
rc <- hclust(as.dist(1-cor(t(m))), method='ward.D2')
par(mfrow=c(1,1), mar=c(8,8,3,2))
gplots::heatmap.2(m,
                  Colv = pdac_lda.clust$dendro,
                  Rowv = as.dendrogram(rc),
                  density.info = "none",
                  trace = "none",
                  ColSideColors = as.vector(clusterColorsk50),
                  col = colorRampPalette(c('blue', 'white', 'red'))(100),
                  cexRow=0.5,cexCol=0.7,margins=c(6,3),
                  main = "pdac_lda.beta",
                  lhei = c(1,5),
                  key.title = NA,
                  scale = 'none')

pdac_lda.betaCombined <- combineTopics(mtx = pdac_lda.beta, clusters = clusterColorsk50, mtxType = "b")

#m <- t(log10(pdac_lda.betaCombined*1e6+1))
m <- t(pdac_lda.betaCombined*1e6)
range(m)
m <- scale(m)
m <- t(scale(t(m)))
range(m)
m[m < -2] <- -2
m[m > 2] <- 2
par(mfrow=c(1,1), mar=c(8,8,3,2))
gplots::heatmap.2(m,
                  #Colv = pdac_lda.clust$dendro,
                  density.info = "none",
                  trace = "none",
                  #ColSideColors = as.vector(clusterColorsk50),
                  col = colorRampPalette(c('blue', 'white', 'red'))(100),
                  cexRow=0.5,cexCol=0.7,margins=c(6,3),
                  main = "pdac_lda.betaCombined",
                  lhei = c(1,5),
                  key.title = NA)

## compare
betaProxy <- do.call(rbind, lapply(levels(com), function(cell_layer){
  layerCounts <- colSums(pdac_corpus[names(com[which(com == cell_layer)]),])
  layerProportions <- layerCounts / sum(layerCounts)
  layerProportions
}))
rownames(betaProxy) <- levels(com)
pdac_lda.betaGTcorMtx <- correlationBetweenBetas(beta1 = pdac_lda.beta,
                                                 beta2 = betaProxy,
                                                 thresh = NULL)
par(mfrow=c(1,1), mar=c(8,8,3,2))
gplots::heatmap.2(pdac_lda.betaGTcorMtx,
          Rowv = pdac_lda.clust$dendro,
          density.info = "none",
          trace = "none",
          RowSideColors = as.vector(clusterColorsk50),
          # ColSideColors = as.vector(classColors),
          col = correlation_palette,
          breaks = correlation_breaks,
          cexRow=0.5,cexCol=0.7,margins=c(6,3),
          main = "pdac_lda.betaGTcorMtx",
          lhei = c(1,5),
          key.xlab = "Correlation",
          key.title = NA)

pdac_lda.thetaCombined <- combineTopics(mtx = pdac_lda.theta, clusters = clusterColorsk50, mtxType = "t")
library(ggplot2)
library(scatterpie)
vizAllTopics(theta = pdac_lda.thetaCombined,
             pos = pdacAst1.Pos,
             topicOrder = seq_len(length(colnames(pdac_lda.thetaCombined))),
             cluster_cols = colnames(pdac_lda.thetaCombined),
             groups = NA,
             group_cols = NA,
             r = 0.4,
             lwd = 0.01)

# plot
par(mfrow=c(5,6), mar=rep(1,4))
sapply(1:ncol(pdac_lda.theta ), function(i) {
  MERINGUE::plotEmbedding(pdacAst1.Pos[rownames(pdac_lda.theta ),], colors=pdac_lda.theta[,i], cex=1, main=i)
})

######### Look at cell-types
m <- t(pdac_lda.beta*1e6)
range(m)
m <- scale(m)
m <- t(scale(t(m)))

i <- 12
head(sort(m[,i], decreasing=TRUE), n=20)
#barplot(sort(pdac_lda.betaCombined[i,], decreasing=TRUE))

######### Looking at specific genes
g <- 'KRT19'
g %in% rownames(pdacAst1.CPM)
MERINGUE::plotEmbedding(pdacAst1.Pos, col=pdacAst1.CPM[g,], main=g, cex=1)

gs <- rownames(pdacAst1.CPM)[grepl('^IL', rownames(pdacAst1.CPM))]
sapply(gs, function(g) {
  g %in% rownames(pdacAst1.CPM)
  MERINGUE::plotEmbedding(pdacAst1.Pos, col=pdacAst1.CPM[g,], main=g, cex=1)
})

gs <- rownames(pdacAst1.CPM)[grepl('^CD', rownames(pdacAst1.CPM))]
gs
sapply(gs, function(g) {
  g %in% rownames(pdacAst1.CPM)
  MERINGUE::plotEmbedding(pdacAst1.Pos, col=pdacAst1.CPM[g,], main=g, cex=1)
})


