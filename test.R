library(MUDAN)

########### Create single cell reference
data(pbmcA)

## filter out poor genes and cells
cd <- cleanCounts(pbmcA,
                  min.reads = 10,
                  min.detected = 10,
                  verbose=FALSE)
## CPM normalization
mat <- normalizeCounts(cd,
                       verbose=FALSE)
## variance normalize, identify overdispersed genes
matnorm.info <- normalizeVariance(mat,
                                  details=TRUE,
                                  verbose=FALSE)
## log transform
matnorm <- log10(matnorm.info$mat+1)
## 30 PCs on overdispersed genes
pcs <- getPcs(matnorm[matnorm.info$ods,],
              nGenes=length(matnorm.info$ods),
              nPcs=30,
              verbose=FALSE)
## get tSNE embedding on PCs
emb <- Rtsne::Rtsne(pcs,
                    is_distance=FALSE,
                    perplexity=30,
                    num_threads=parallel::detectCores(),
                    verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)

## get clusters
set.seed(0)
com <- getComMembership(pcs, k=15, method=igraph::cluster_louvain)
table(com)

par(mfrow=c(3,3), mar=c(0.5,0.5,2,0.5))
plotEmbedding(emb, groups=com, mark.clusters = TRUE)

## plot expression of marker genes
marker.genes <- c('MS4A1', 'CD14', 'FCGR3A', 'GZMA', 'CD8A', 'CD4', 'CD3E')
invisible(lapply(marker.genes, function(g) {
  gcol <- scale(mat[g,])[,1] ## plot a gene
  plotEmbedding(emb,
                color=gcol,
                mark.clusters=TRUE,
                main=g, xlab=NA, ylab=NA,
                verbose=FALSE, alpha=0.1)
}))

annot <- as.character(com)
names(annot) <- names(com)
#annot[com == 5] <- 'Nonclassical Monocytes'
#annot[com == 10] <- 'Classical Monocytes'
annot[com == 4] <- 'B Cells'
annot[com == 3] <- 'NK Cells'
#annot[com %in% c(1,2,6,7,8,9)] <- 'T Cells'
#annot[com == 7] <- 'CD8+ Effector T Cells'
#annot[com == 9] <- 'CD8+ Memory T Cells'
#annot[com == 8] <- 'CD8+ Naive T Cells'
#annot[com == 6] <- 'CD4+ Memory T Cells 1'
#annot[com == 1] <- 'CD4+ Memory T Cells 2'
#annot[com == 2] <- 'CD4+ Memory T Cells 3'
annot[com %in% c(1,2,6,8,9,5,10,7)] <- NA
annot <- factor(annot)
table(annot)

par(mfrow=c(1,1))
plotEmbedding(emb, groups=annot, mark.clusters = TRUE, mark.cluster.cex = 1)
annot <- na.omit(annot)

## summarize cell-types
mm <- model.matrix(~ 0 + annot)
colnames(mm) <- levels(annot)
## take alook
dim(mm)
annot.summary <- cd[, names(annot)] %*% mm
annot.summary <- log10(normalizeCounts(annot.summary)+1)
head(annot.summary)

## restrict to od genes
#od.genes <- rownames(matnorm)[matnorm.info$ods]
#annot.summary <- annot.summary[od.genes,]
#dim(annot.summary)
heatmap(as.matrix(annot.summary), scale='row')

##### Create mixtures
N <- 50 ## number samples
#split <- sample(paste0('spot', seq_len(N)), length(annot), replace=TRUE)
#split <- factor(split)
#ground.truth <- table(split, annot)
#head(ground.truth)
ground.truth <- cbind(1:N, N-(1:N))
colnames(ground.truth) <- levels(annot)
rownames(ground.truth) <- paste0('spot', seq_len(N))
heatmap(ground.truth, scale='none', Rowv = NA, Colv=NA)

## make model matrix
#mm <- model.matrix(~ 0 + split)
#colnames(mm) <- paste0('spot', seq_len(N))
### take alook
#dim(mm)
#cd.summary <- cd[od.genes,] %*% mm
#head(cd.summary)

set.seed(0)
cd.summary <- do.call(cbind, lapply(1:N, function(i) {
  cells.sample <- c(
    sample(names(annot)[annot==levels(annot)[1]], ground.truth[i,1], replace=TRUE),
    sample(names(annot)[annot==levels(annot)[2]], ground.truth[i,2], replace=TRUE)
  )
  #rowSums(cd[od.genes, cells.sample])
  rowSums(cd[,cells.sample])
}))
colnames(cd.summary) <- rownames(ground.truth)
dim(cd.summary)
head(cd.summary)

foo <- normalizeVariance(cd.summary, details=TRUE)
od.genes <- rownames(cd.summary)[foo$ods]
annot.summary <- annot.summary[od.genes,]
cd.summary <- cd.summary[od.genes,]
heatmap(as.matrix(annot.summary), scale='none', Colv=NA)

##### LDA
library(topicmodels)
library(slam)

ap_lda <- LDA(slam::as.simple_triplet_matrix(t(cd.summary)), k=length(levels(annot)), control = list(seed = 0))
ap_lda

## compare discovered topics to real cell-types
library(tidytext)
ap_topics <- tidy(ap_lda, matrix = "beta")
ap_topics
ap_topics.mat <- t(cast_sparse(ap_topics, topic, term, beta))
dim(ap_topics.mat)
head(ap_topics.mat)
#par(mfrow=c(ncol(ap_topics.mat), ncol(annot.summary)))
cor.mat <- do.call(rbind, lapply(colnames(ap_topics.mat), function(x) {
  sapply(colnames(annot.summary), function(y) {
    #plot(ap_topics.mat[,x], annot.summary[,y],
    #     main=paste0('Topic:', x, '\nAnnot:', y))
    cor(ap_topics.mat[,x], annot.summary[,y],)
  })
}))
dim(cor.mat)
rownames(cor.mat) <- colnames(ap_topics.mat)
head(cor.mat)
heatmap(cor.mat, scale='none')

## compare discovered proportions to real proportions
ap_documents <- tidy(ap_lda, matrix = "gamma")
ap_documents
ap_documents.mat <- cast_sparse(ap_documents, document, topic, gamma)
dim(ap_documents.mat)
head(ap_documents.mat)
dim(ground.truth)

cor.mat <- do.call(rbind, lapply(rownames(ap_documents.mat), function(x) {
  sapply(rownames(ground.truth), function(y) {
    cor(ap_documents.mat[x,], ground.truth[y,],)
  })
}))
dim(cor.mat)
head(cor.mat)
rownames(cor.mat) <- rownames(ap_documents.mat)
head(cor.mat)
cor.mat[is.na(cor.mat)] <- 0
heatmap(cor.mat, scale='none')
