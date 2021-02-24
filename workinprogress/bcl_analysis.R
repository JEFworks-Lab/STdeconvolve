library(MERINGUE)
data(BCL)
names(BCL)
head(BCL$pos)
head(BCL$counts)

pos <- BCL$pos[,1:2]
slice <- BCL$pos[,3]
counts <- BCL$counts
names(slice) <- rownames(pos)

########## Clustering

# clean
bcl.Clean <- MERINGUE::cleanCounts(counts = counts, # genes x spots mtx
                                        min.reads = 100,
                                        min.lib.size = 100,
                                        verbose=TRUE)

# CPM normalize
dim(bcl.Clean)
bcl.CPM <- MERINGUE::normalizeCounts(counts = bcl.Clean,
                                          verbose=TRUE, log=FALSE)

## Get overdispsersed genes
counts <- bcl.Clean
odGenes <- Reduce(intersect, lapply(1:4, function(i) {
  MERINGUE::getOverdispersedGenes(bcl.CPM[, slice[colnames(counts)]==i], plot = TRUE, details = TRUE, alpha = 0.2)$ods
}))
length(odGenes)
countsFilt <- counts[odGenes,]

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
pcs.info <- prcomp(t(as.matrix(bcl.CPM[rownames(countsFiltnoUnifGenes),])), center=TRUE)
plot(pcs.info$sdev[1:30])
nPcs <- 10
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

par(mfrow=c(2,2), mar=rep(1,4))
MUDAN::plotEmbedding(emb, groups=com,
                     show.legend=TRUE, xlab=NA, ylab=NA,
                     verbose=FALSE)
MUDAN::plotEmbedding(emb, groups=slice,
                     show.legend=TRUE, xlab=NA, ylab=NA,
                     verbose=FALSE)
MUDAN::plotEmbedding(pos, groups=com,
                     cex=1, xlab=NA, ylab=NA,
                     verbose=FALSE)
MUDAN::plotEmbedding(pos, groups=slice,
                     cex=1, xlab=NA, ylab=NA,
                     verbose=FALSE)

############################### Topic Modeling
bcl_corpus <- t(as.matrix(countsFiltnoUnifGenes))
bcl_corpus_slamMtx <- slam::as.simple_triplet_matrix(bcl_corpus)

# optimal k?
#Ks = seq(5,30,by=5)
#Ks
#bcl_lda <- fitLDA(bcl_corpus_slamMtx, Ks=Ks)
#par(mfrow=c(1,1), mar=rep(5,4))
#plot(x=Ks, bcl_lda$perplexities)

# pick k
#kopt = Ks[which(bcl_lda$perplexities == min(bcl_lda$perplexities))]
#kopt
#kopt = 10 ## manual
kopt = 6
bcl_lda.k <- topicmodels::LDA(bcl_corpus_slamMtx, k=kopt)
topicmodels::perplexity(bcl_lda.k, bcl_corpus_slamMtx)

bcl_lda.tmResult<- topicmodels::posterior(bcl_lda.k)
bcl_lda.theta <- bcl_lda.tmResult$topics
bcl_lda.beta <- bcl_lda.tmResult$terms
#bcl_lda.topicFreqsOverall <- colSums(bcl_lda.theta) / nrow(bcl_corpus_slamMtx)

library(dynamicTreeCut)
bcl_lda.clust <- clusterTopics(beta = bcl_lda.beta,
                                distance = "euclidean",
                                clustering = "ward.D2",
                                dynamic = "hybrid",
                                deepSplit = 4,
                                plotDendro = TRUE)
clusterColorsk50 <- bcl_lda.clust$clusters
levels(clusterColorsk50) <- gg_color_hue(length(levels(clusterColorsk50)))

## remove topics with too few words
#m <- scale(bcl_lda.beta)
#dim(m)
#heatmap(m, scale='none')
#hist(m[1,], breaks=100)
#table(m[1,]>0)
#hist(rowSums(m>0.1), breaks=20)
#good.topics <- rowSums(m>0)> 16
#good.topics <- names(which(good.topics))

#m <- t(bcl_lda.beta[good.topics,]*1e6)
#m <- log10(t(bcl_lda.beta*1e6)+1)
m <- t(bcl_lda.beta)
range(m)
m <- scale(m)
#m <- t(scale(t(m)))
range(m)
m[m < -2] <- -2
m[m > 2] <- 2
rc <- hclust(as.dist(1-cor(t(m))), method='ward.D2')
par(mfrow=c(1,1), mar=c(8,8,3,2))
gplots::heatmap.2(m,
                  Colv = bcl_lda.clust$dendro,
                  Rowv = as.dendrogram(rc),
                  density.info = "none",
                  trace = "none",
                  ColSideColors = as.vector(clusterColorsk50),
                  col = colorRampPalette(c('blue', 'white', 'red'))(100),
                  cexRow=0.5,cexCol=0.7,margins=c(6,3),
                  main = "bcl_lda.beta",
                  lhei = c(1,5),
                  key.title = NA,
                  scale = 'none')

bcl_lda.betaCombined <- combineTopics(mtx = bcl_lda.beta, clusters = clusterColorsk50, mtxType = "b")

#m <- t(log10(bcl_lda.betaCombined*1e6+1))
m <- t(bcl_lda.betaCombined*1e6)
### do not combine
#m <- t(bcl_lda.beta[good.topics,]*1e6)
range(m)
m <- scale(m)
#m <- t(scale(t(m)))
range(m)
m[m < -2] <- -2
m[m > 2] <- 2
par(mfrow=c(1,1), mar=c(8,8,3,2))
gplots::heatmap.2(m,
                  #Colv = bcl_lda.clust$dendro,
                  density.info = "none",
                  trace = "none",
                  #ColSideColors = as.vector(clusterColorsk50),
                  col = colorRampPalette(c('blue', 'white', 'red'))(100),
                  cexRow=0.5,cexCol=0.7,margins=c(6,3),
                  main = "bcl_lda.betaCombined",
                  lhei = c(1,5),
                  key.title = NA)

## compare
betaProxy <- do.call(rbind, lapply(levels(com), function(cell_layer){
  layerCounts <- colSums(bcl_corpus[names(com[which(com == cell_layer)]),])
  layerProportions <- layerCounts / sum(layerCounts)
  layerProportions
}))
rownames(betaProxy) <- levels(com)
bcl_lda.betaGTcorMtx <- correlationBetweenBetas(beta1 = bcl_lda.beta,
                                                 beta2 = betaProxy,
                                                 thresh = NULL)
par(mfrow=c(1,1), mar=c(8,8,3,2))
gplots::heatmap.2(bcl_lda.betaGTcorMtx,
                  Rowv = bcl_lda.clust$dendro,
                  density.info = "none",
                  trace = "none",
                  RowSideColors = as.vector(clusterColorsk50),
                  # ColSideColors = as.vector(classColors),
                  col = correlation_palette,
                  breaks = correlation_breaks,
                  cexRow=0.5,cexCol=0.7,margins=c(6,3),
                  main = "bcl_lda.betaGTcorMtx",
                  lhei = c(1,5),
                  key.xlab = "Correlation",
                  key.title = NA)

bcl_lda.thetaCombined <- combineTopics(mtx = bcl_lda.theta, clusters = clusterColorsk50, mtxType = "t")
#bcl_lda.thetaCombined <- bcl_lda.theta[, good.topics]
#bcl_lda.thetaCombined <- bcl_lda.theta
library(ggplot2)
library(scatterpie)
vizAllTopics(theta = bcl_lda.thetaCombined,
             pos = pos[slice==1,],
             topicOrder = seq_len(length(colnames(bcl_lda.thetaCombined))),
             cluster_cols = rainbow(ncol(bcl_lda.thetaCombined)),
             groups = NA,
             group_cols = NA,
             r = 0.4,
             lwd = 0.01)

vizAllTopics(theta = bcl_lda.thetaCombined,
             pos = pos[slice==2,],
             topicOrder = seq_len(length(colnames(bcl_lda.thetaCombined))),
             cluster_cols = rainbow(ncol(bcl_lda.thetaCombined)),
             groups = NA,
             group_cols = NA,
             r = 0.4,
             lwd = 0.01)

vizAllTopics(theta = bcl_lda.thetaCombined,
             pos = pos[slice==3,],
             topicOrder = seq_len(length(colnames(bcl_lda.thetaCombined))),
             cluster_cols = rainbow(ncol(bcl_lda.thetaCombined)),
             groups = NA,
             group_cols = NA,
             r = 0.4,
             lwd = 0.01)

vizAllTopics(theta = bcl_lda.thetaCombined,
             pos = pos[slice==4,],
             topicOrder = seq_len(length(colnames(bcl_lda.thetaCombined))),
             cluster_cols = rainbow(ncol(bcl_lda.thetaCombined)),
             groups = NA,
             group_cols = NA,
             r = 0.4,
             lwd = 0.01)



