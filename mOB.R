library(MERINGUE)
data(mOB)
names(mOB)

cd.summary <- cleanCounts(mOB$counts, min.lib.size = 100, min.detected = 100)
dim(cd.summary)
annot <- mOB$annot[colnames(cd.summary)]
levels(annot)
pos <- mOB$pos
par(mfrow=c(1,1))
plotEmbedding(pos, annot, cex=3)

foo <- normalizeVariance(cd.summary, details=TRUE, alpha = 0.2)
od.genes <- rownames(cd.summary)[foo$ods]

#K = length(levels(annot))
K = 9
ap_lda <- LDA(slam::as.simple_triplet_matrix(t(cd.summary[od.genes,])), k=K, control = list(seed = 0))
ap_lda

library(tidytext)
ap_topics <- tidy(ap_lda, matrix = "beta")

ap_topics.mat <- t(cast_sparse(ap_topics, topic, term, beta))
dim(ap_topics.mat)
heatmap(as.matrix(ap_topics.mat), scale='row')

ap_documents <- tidy(ap_lda, matrix = "gamma")
ap_documents
ap_documents.mat <- cast_sparse(ap_documents, document, topic, gamma)
dim(ap_documents.mat)
head(ap_documents.mat)

## test each topic with fgsea
## if share pathways, then presumably not so independent
library(fgsea)
library(liger)
data("org.Hs.GO2Symbol.list")

fgseaRes.all <- do.call(rbind, lapply(1:ncol(ap_topics.mat), function(i) {
  print(i)
  ranks <- ap_topics.mat[,i]
  ranks <- sort(ranks, decreasing=TRUE)
  names(ranks) <- toupper(names(ranks))
  #barplot(ranks)

  fgseaRes <- fgsea(pathways = org.Hs.GO2Symbol.list,
                    stats    = ranks,
                    eps      = 0.0,
                    minSize  = 0,
                    maxSize  = Inf,
                    scoreType = 'pos')
  fgseaRes <- data.frame(fgseaRes)
  rownames(fgseaRes) <- fgseaRes$pathway
  fgseaRes <- fgseaRes[names(org.Hs.GO2Symbol.list),]

  out <- -log10(fgseaRes$padj)
  is.na(out) <- 0
  return(out)
}))
colnames(fgseaRes.all) <- names(org.Hs.GO2Symbol.list)
fgseaRes.all[1:5,1:5]

## at least one significant
vi <- names(which(colSums(fgseaRes.all >= -log10(0.2)) >= 1))
length(vi)
head(vi)
fgsea.sig <- fgseaRes.all[, vi]
fgsea.sig[1:5,1:5]

heatmap(fgsea.sig, scale='none')

## plot
par(mfrow=c(3,3))
sapply(1:ncol(ap_documents.mat), function(i){
  plotEmbedding(pos[rownames(ap_documents.mat),], col=ap_documents.mat[,i], cex=1, main=paste0('Proportion of Cell-Type:', i))
})

## is there a higher resolution mOB dataset we can use?
