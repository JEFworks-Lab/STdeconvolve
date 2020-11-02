## Compare latent direchlet allocation
## with correlated topics modeling

################### Simulation
set.seed(0)
base = abs(round(rnorm(10, 10)))
names(base) <- paste0('gene', 1:10)
head(base)

## cell-type A upregulates genes 1 to 4
x1 = base
x1[1:4] = base[1:4] + 10

## cell-type B upregulates genes 7 to 10
x2 = base
x2[7:10] = base[7:10] + 10

## cell-type C upregulates genes 5 to 10 slightly
x3 = base
x3[5:10] = base[5:10] + 5

## visualize transcriptional distinctness
mat <- rbind(x1, x2, x3)
heatmap(mat, scale='none')
## some cell-types are more correlated than others
cor(t(mat))

## training data is differing proportions of
## cell-type A and cell-type B with random cell-type C
set.seed(0)
train.data <- do.call(rbind, lapply(seq(1,80, by=5), function(i) {
  rand = round(runif(1)*20)
  ground.truth <- c(i, 80-i, rand)
  names(ground.truth) <- c('ctA', 'ctB', 'ctC')
  y = ground.truth[1]*x1 + ground.truth[2]*x2 + rand*x3
  c(y, ground.truth/sum(ground.truth))
}))
rownames(train.data) <- paste0('sim', 1:nrow(train.data))
head(train.data)

## remove lables
train.data.sub <- as.data.frame(train.data[, names(x1)])
head(train.data.sub)

heatmap(as.matrix(train.data.sub), scale='none', Rowv = NA, Colv = NA)

############### LDA
library(topicmodels)
library(slam)
library(Matrix)
K = 3
ap_lda <- topicmodels::LDA(slam::as.simple_triplet_matrix(as.matrix(train.data.sub)),
                           k=K, control = list(seed = 0))
ap_lda

## association of each gene with each cell-type
library(tidytext)
ap_topics <- tidy(ap_lda, matrix = "beta")
ap_topics

ap_topics.mat <- t(cast_sparse(ap_topics, topic, term, beta))
dim(ap_topics.mat)
heatmap(as.matrix(ap_topics.mat), scale='none', Rowv=NA, Colv=NA)

## proportional representation of each topic in each document
## or each cell-type in each spot
ap_documents <- tidy(ap_lda, matrix = "gamma")
ap_documents

ap_documents.mat <- cast_sparse(ap_documents, document, topic, gamma)
dim(ap_documents.mat)
head(ap_documents.mat)

plot(ap_documents.mat[,1], train.data[, 'ctB'])
cor(ap_documents.mat[,1], train.data[, 'ctB'])
plot(ap_documents.mat[,2], train.data[, 'ctA'])
cor(ap_documents.mat[,2], train.data[, 'ctA'])
plot(ap_documents.mat[,3], train.data[, 'ctC'])
cor(ap_documents.mat[,3], train.data[, 'ctC'])

################ CTM
K = 3
ap_ctm <- topicmodels::CTM(slam::as.simple_triplet_matrix(as.matrix(train.data.sub)),
                           k=K, control = list(seed = 0))
ap_ctm

## association of each gene with each cell-type
library(tidytext)
ap_topics <- tidy(ap_ctm, matrix = "beta")
ap_topics

ap_topics.mat <- t(cast_sparse(ap_topics, topic, term, beta))
dim(ap_topics.mat)
heatmap(as.matrix(ap_topics.mat), scale='row', Rowv=NA, Colv=NA)

## proportional representation of each topic in each document
## or each cell-type in each spot
ap_documents <- tidy(ap_lda, matrix = "gamma")
ap_documents

ap_documents.mat <- cast_sparse(ap_documents, document, topic, gamma)
dim(ap_documents.mat)
head(ap_documents.mat)

plot(ap_documents.mat[,1], train.data[, 'ctB'])
cor(ap_documents.mat[,1], train.data[, 'ctB'])
plot(ap_documents.mat[,2], train.data[, 'ctA'])
cor(ap_documents.mat[,2], train.data[, 'ctA'])
plot(ap_documents.mat[,3], train.data[, 'ctC'])
cor(ap_documents.mat[,3], train.data[, 'ctC'])

################ Apply to real data
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
ap_lda <- CTM(slam::as.simple_triplet_matrix(t(cd.summary[od.genes,])), k=K, control = list(seed = 0))
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

## plot
par(mfrow=c(3,3))
sapply(1:ncol(ap_documents.mat), function(i){
  plotEmbedding(pos[rownames(ap_documents.mat),], col=ap_documents.mat[,i], cex=1, main=paste0('Proportion of Cell-Type:', i))
})

## Conclusion:
## LDA seems better
## more independent spatial patterns
## suggesting identifying more independent cell-types
## We should be able to assume that cell-types, while can be correlated,
## should also have distinct marker genes (else questionable whether it is an independent cell-type)
