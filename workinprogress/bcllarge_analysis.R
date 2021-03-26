library(Seurat)
cd <- Read10X_h5("~/Desktop/BCL/Parent_Visium_Human_BreastCancer_filtered_feature_bc_matrix.h5")
pos.info <- read.csv("~/Desktop/BCL/tissue_positions_list.csv", header=FALSE)
pos <- pos.info[,5:6]
rownames(pos) <- pos.info[,1]
head(pos)

counts <- cleanCounts(cd, min.lib.size = 0, min.reads = 200, verbose=TRUE)
pos <- pos[colnames(counts),]
plotEmbedding(pos[colnames(counts),], col=colSums(counts))
plotEmbedding(pos[colnames(counts),], col=counts[,])

## clustering
mat <- MERINGUE::normalizeCounts(counts, log=FALSE)
ods <- MERINGUE::getOverdispersedGenes(mat, plot=TRUE, details = TRUE, alpha=0.01)
pcs <- prcomp(t(mat[ods$ods,]))
plot(pcs$sdev[1:10], type="l")
emb <- Rtsne::Rtsne(pcs$x[,1:5])$Y
rownames(emb) <- colnames(counts)
com <- MERINGUE::getClusters(pcs$x[,1:5], k=30)
table(com)
par(mfrow=c(1,2))
plotEmbedding(emb, groups=com)
plotEmbedding(pos, groups=com)

## LDA
corpus <- restrictCorpus(counts, alpha = 0.01, t = 1)
#lda <- fitLDA(corpus, Ks = seq(5, 20, by = 5), plot=TRUE)
lda <- fitLDA(corpus, Ks = 10, plot=TRUE)
#lda$kOpt1
#lda$kOpt2
results <- getBetaTheta(lda$models[["10"]])
library(dynamicTreeCut)
order <- clusterTopics(results$beta, plot=TRUE)
library(ggplot2)
library(scatterpie)
vizAllTopics(theta = results$theta,
             pos = pos,
             topicOrder = order$order)

par(mfrow=c(4,4))
lapply(1:10, function(i) {
  plotEmbedding(pos, col=results$theta[,i], main=i)
})

head(sort(results$beta[1,], decreasing=TRUE))
head(sort(results$beta[2,], decreasing=TRUE))
head(sort(results$beta[3,], decreasing=TRUE))
head(sort(results$beta[4,], decreasing=TRUE))
head(sort(results$beta[5,], decreasing=TRUE))
head(sort(results$beta[6,], decreasing=TRUE))
head(sort(results$beta[7,], decreasing=TRUE))
head(sort(results$beta[8,], decreasing=TRUE))
head(sort(results$beta[9,], decreasing=TRUE))
head(sort(results$beta[10,], decreasing=TRUE))

## GSEA on each to interpret
## enrichment
library(liger)
data(org.Hs.GO2Symbol.list) ## load built in GO gene sets
go.env <- org.Hs.GO2Symbol.list

## annotate GO terms with names
library(GO.db)
library(AnnotationDbi)
desc <- AnnotationDbi::select(
  GO.db,
  keys = names(go.env),
  columns = c("TERM"),
  multiVals = 'CharacterList'
)
names(go.env) <- paste(names(go.env), desc$TERM)

vals <- sort(results$beta[8,], decreasing=TRUE)
gsea.results <- iterative.bulk.gsea(values=vals, set.list=go.env, rank=TRUE)

gsea.sig.up <- gsea.results[gsea.results$q.val < 0.05 & gsea.results$sscore > 0 & gsea.results$edge > 0,]
head(gsea.sig.up)

## visualize top hits
print(rownames(gsea.sig.up)[1])
gsea(values=vals, geneset=go.env[[rownames(gsea.sig.up)[1]]], plot=TRUE, rank=TRUE)

## plot gene sets
par(mfrow=c(2,2))
desc[grepl('CD8-positive', desc$TERM),]
desc[grepl('ribosome', desc$TERM),]
i <- 380 ## angiogenesisi
i <- 728 ## immune
i <- 8459 ## cd8
i <- 2257
gexp <- scale(colSums(mat[intersect(go.env[[i]], rownames(mat)),]))[,1]
#gexp <- scale(colSums(mat))[,1]
range(gexp)
gexp[gexp > 2] <- 2
plotEmbedding(emb, col=gexp, main=names(go.env)[i])
plotEmbedding(pos, col=gexp)
