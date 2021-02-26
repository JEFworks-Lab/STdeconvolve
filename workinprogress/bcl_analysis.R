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
                                        min.reads = 10,
                                        min.lib.size = 10,
                                        verbose=TRUE)

'CCR7' %in% rownames(bcl.Clean)

# CPM normalize
dim(bcl.Clean)
bcl.CPM <- MERINGUE::normalizeCounts(counts = bcl.Clean,
                                          verbose=TRUE, log=FALSE)

## Get overdispsersed genes
counts <- bcl.Clean
odGenes <- Reduce(intersect, lapply(1:4, function(i) {
  MERINGUE::getOverdispersedGenes(counts[, slice[colnames(counts)]==i], plot = TRUE, details = TRUE, alpha = 0.2)$ods
}))
length(odGenes)
countsFilt <- counts[odGenes,]

## remove genes that are present in more than X% of spots
## if we don't remove perplexity poor
## also remove genes that are too rare
#library(Matrix)
#vi <- rowSums(countsFilt > 0) > ncol(countsFilt)*0.1 & rowSums(countsFilt > 0) < ncol(countsFilt)*0.9
#table(vi)
#countsFiltnoUnifGenes <- countsFilt[vi,]

countsFiltnoUnifGenes <- countsFilt
#dim(countsFiltnoUnifGenes)
hist(log10(colSums(countsFiltnoUnifGenes)+1))
hist(log10(rowSums(countsFiltnoUnifGenes)+1))

# Dimensionality reduction by PCA on log10 CPM expression values
pcs.info <- prcomp(t(as.matrix(bcl.CPM[rownames(countsFiltnoUnifGenes),])), center=TRUE)
plot(pcs.info$sdev[1:30])
nPcs <- 5
pcs <- pcs.info$x[,1:nPcs]
# 2D embedding by tSNE
emb <- Rtsne::Rtsne(pcs,
                    is_distance=FALSE,
                    perplexity=30,
                    num_threads=1,
                    verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)
# Graph-based cluster detection
k <- 300
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

############################### Overlay ontop image
library(jpeg)
he1 <- readJPEG('~/Desktop/bcl/HE_layer1_BC.jpg')
#rasterImage(he1, 0, 0, dim(he1)[1], dim(he1)[2])
dim(he1)

posScale <- pos*32
posScale[,1] = posScale[,1] - min(posScale[,1])
posScale[,2] = posScale[,2] - min(posScale[,2])
range(posScale)

par(mfrow=c(1,1), mar=rep(5,4))
plot.new()
plot(posScale[slice==1,"x"], posScale[slice==1,"y"], col='red', pch="+")
rasterImage(he1,
             -10, -5,
             dim(he1)[2], dim(he1)[1])
points(posScale[slice==1,"x"], posScale[slice==1,"y"], col='black', cex=3)
points(posScale[slice==1,"x"], posScale[slice==1,"y"], col=MERINGUE:::fac2col(com), pch=16)

############################### Topic Modeling
bcl_corpus <- t(as.matrix(countsFiltnoUnifGenes))
bcl_corpus_slamMtx <- slam::as.simple_triplet_matrix(bcl_corpus)

# optimal k?
Ks = seq(2,10,by=1)
Ks
bcl_lda <- fitLDA(bcl_corpus_slamMtx, Ks=Ks)
par(mfrow=c(1,1), mar=rep(5,4))
plot(x=Ks, bcl_lda$perplexities)

# pick k
#kopt = Ks[which(bcl_lda$perplexities == min(bcl_lda$perplexities))]
#kopt
#kopt = 9 ## manual
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
                                clustering = "ward.D",
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

#bcl_lda.thetaCombined <- combineTopics(mtx = bcl_lda.theta, clusters = clusterColorsk50, mtxType = "t")
#bcl_lda.thetaCombined <- bcl_lda.theta[, good.topics]
bcl_lda.thetaCombined <- bcl_lda.theta
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

####### plot
# plot
par(mfrow=c(kopt,4), mar=rep(1,4))
sapply(1:ncol(bcl_lda.theta ), function(i) {
  cc = intersect(names(which(slice==1)), rownames(bcl_lda.theta))
  MERINGUE::plotEmbedding(pos[cc,], colors=bcl_lda.theta[cc,i], cex=1, main=i)
  cc = intersect(names(which(slice==2)), rownames(bcl_lda.theta))
  MERINGUE::plotEmbedding(pos[cc,], colors=bcl_lda.theta[cc,i], cex=1, main=i)
  cc = intersect(names(which(slice==3)), rownames(bcl_lda.theta))
  MERINGUE::plotEmbedding(pos[cc,], colors=bcl_lda.theta[cc,i], cex=1, main=i)
  cc = intersect(names(which(slice==4)), rownames(bcl_lda.theta))
  MERINGUE::plotEmbedding(pos[cc,], colors=bcl_lda.theta[cc,i], cex=1, main=i)
})

par(mfrow=c(3,4), mar=rep(1,4))
sapply(1:ncol(bcl_lda.thetaCombined ), function(i) {
  cc = intersect(names(which(slice==1)), rownames(bcl_lda.theta))
  MERINGUE::plotEmbedding(pos[cc,], colors=bcl_lda.thetaCombined[cc,i], cex=2, main=i)
  cc = intersect(names(which(slice==2)), rownames(bcl_lda.theta))
  MERINGUE::plotEmbedding(pos[cc,], colors=bcl_lda.thetaCombined[cc,i], cex=2, main=i)
  cc = intersect(names(which(slice==3)), rownames(bcl_lda.theta))
  MERINGUE::plotEmbedding(pos[cc,], colors=bcl_lda.thetaCombined[cc,i], cex=2, main=i)
  cc = intersect(names(which(slice==4)), rownames(bcl_lda.theta))
  MERINGUE::plotEmbedding(pos[cc,], colors=bcl_lda.thetaCombined[cc,i], cex=2, main=i)
})

### overlay onto images
par(mfrow=c(1,1), mar=rep(5,4))
cc = intersect(names(which(slice==1)), rownames(bcl_lda.theta))
plot.new()
plot(posScale[cc,"x"], posScale[cc,"y"], col='red', pch="+")
rasterImage(he1,
            -10, -5,
            dim(he1)[2], dim(he1)[1])
#points(posScale[cc,"x"], posScale[cc,"y"], col='black', cex=3)
points(posScale[cc,"x"], posScale[cc,"y"], col=MERINGUE:::map2col(bcl_lda.thetaCombined[cc,3]), pch=16, cex=3)
points(posScale[cc,"x"], posScale[cc,"y"], col=MERINGUE:::map2col(bcl_lda.thetaCombined[cc,2]), pch=16, cex=3)
points(posScale[cc,"x"], posScale[cc,"y"], col=MERINGUE:::map2col(bcl_lda.thetaCombined[cc,4]), pch=16, cex=3)


#### proportional co-localization
head(bcl_lda.theta)
rowSums(bcl_lda.theta)

pcor <- do.call(rbind, lapply(1:kopt, function(i) {
  p1 <- bcl_lda.theta[,i]
  unlist(lapply(1:kopt, function(j) {
    p2 <- bcl_lda.theta[,j]
    cor(p1, p2)
  }))
}))
diag(pcor) <- NA
range(pcor)
heatmap(pcor, col=correlation_palette)
which(pcor == max(pcor), arr.ind=TRUE)

tcor <- do.call(rbind, lapply(1:kopt, function(i) {
  p1 <- bcl_lda.beta[i,]
  unlist(lapply(1:kopt, function(j) {
    p2 <- bcl_lda.beta[j,]
    cor(p1, p2)
  }))
}))
diag(tcor) <- NA
range(tcor)
heatmap(tcor, col=correlation_palette, scale='none')
which(tcor == max(tcor), arr.ind=TRUE)

foo <- pcor - tcor
range(foo)
diag(foo) <- NA
heatmap(foo, col=correlation_palette, scale='none')
which(foo == max(foo, na.rm=TRUE), arr.ind=TRUE)


## look at genes
m <- t(bcl_lda.beta*1e6)
range(m)
#m <- scale(m)
#m <- t(scale(t(m)))

i <- 5
head(sort(m[,i], decreasing=TRUE), n=9)
#barplot(sort(bcl_lda.betaCombined[i,], decreasing=TRUE))

par(mfrow=c(5,6), mar=rep(1,4))
gs <- names(head(sort(m[,i], decreasing=TRUE), n=30))
sapply(gs, function(g) {
  g %in% rownames(bcl.CPM)
  gexp <- scale(bcl.CPM[g,])[,1]
  gexp[gexp > 2] <- 2
  gexp[gexp < -2] <- -2
  MERINGUE::plotEmbedding(pos[slice[rownames(pos)]==1, ], col=gexp, main=g, cex=1)
})

## specific
par(mfrow=c(4,5), mar=rep(1,4))

## slice
lapply(1:4, function(s) {
  ## topic
  i = 4
  cc = intersect(names(which(slice==s)), rownames(bcl_lda.theta))
  MERINGUE::plotEmbedding(pos[cc,], colors=bcl_lda.thetaCombined[cc,i], cex=2, main=paste0('slice:', s, '+ topic:', i))
  g <- 'CCL21'
  g %in% rownames(bcl.CPM)
  gexp <- scale(bcl.CPM[g,cc])[,1]
  gexp[gexp > 2] <- 2
  gexp[gexp < -2] <- -2
  MERINGUE::plotEmbedding(pos[cc, ], col=gexp, main=g, cex=2)
  g <- 'IGLL5'
  g %in% rownames(bcl.CPM)
  gexp <- scale(bcl.CPM[g,cc])[,1]
  gexp[gexp > 2] <- 2
  gexp[gexp < -2] <- -2
  MERINGUE::plotEmbedding(pos[cc, ], col=gexp, main=g, cex=2)
  g <- 'JCHAIN'
  g %in% rownames(bcl.CPM)
  gexp <- scale(bcl.CPM[g,cc])[,1]
  gexp[gexp > 2] <- 2
  gexp[gexp < -2] <- -2
  MERINGUE::plotEmbedding(pos[cc, ], col=gexp, main=g, cex=2)
  g <- 'CD7'
  g %in% rownames(bcl.CPM)
  gexp <- scale(bcl.CPM[g,cc])[,1]
  gexp[gexp > 2] <- 2
  gexp[gexp < -2] <- -2
  MERINGUE::plotEmbedding(pos[cc, ], col=gexp, main=g, cex=2)
})

lapply(1:4, function(s) {
  ## topic
  i = 5
  cc = intersect(names(which(slice==s)), rownames(bcl_lda.theta))
  MERINGUE::plotEmbedding(pos[cc,], colors=bcl_lda.thetaCombined[cc,i], cex=2, main=paste0('slice:', s, '+ topic:', i))
  g <- 'PGR'
  g %in% rownames(bcl.CPM)
  gexp <- scale(bcl.CPM[g,cc])[,1]
  gexp[gexp > 2] <- 2
  gexp[gexp < -2] <- -2
  MERINGUE::plotEmbedding(pos[cc, ], col=gexp, main=g, cex=2)
  g <- 'FN1'
  g %in% rownames(bcl.CPM)
  gexp <- scale(bcl.CPM[g,cc])[,1]
  gexp[gexp > 2] <- 2
  gexp[gexp < -2] <- -2
  MERINGUE::plotEmbedding(pos[cc, ], col=gexp, main=g, cex=2)
  g <- 'CXCL10'
  g %in% rownames(bcl.CPM)
  gexp <- scale(bcl.CPM[g,cc])[,1]
  gexp[gexp > 2] <- 2
  gexp[gexp < -2] <- -2
  MERINGUE::plotEmbedding(pos[cc, ], col=gexp, main=g, cex=2)
  g <- 'CXCL9'
  g %in% rownames(bcl.CPM)
  gexp <- scale(bcl.CPM[g,cc])[,1]
  gexp[gexp > 2] <- 2
  gexp[gexp < -2] <- -2
  MERINGUE::plotEmbedding(pos[cc, ], col=gexp, main=g, cex=2)
})


## groups of genes
gs <- rownames(bcl.CPM)[grepl('^CD', rownames(bcl.CPM))]
sapply(gs, function(g) {
  g %in% rownames(bcl.CPM)
  gexp <- scale(bcl.CPM[g,])[,1]
  gexp[gexp > 2] <- 2
  gexp[gexp < -2] <- -2
  MERINGUE::plotEmbedding(pos[slice[rownames(pos)]==1, ], col=gexp, main=g, cex=2)
})
