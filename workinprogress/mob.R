library(MERINGUE)
library(Matrix)
data(mOB)
pos <- mOB$pos
cd <- mOB$counts
annot <- mOB$annot

## Remove poor datasets and genes
par(mfrow=c(1,1), mar=rep(5,4))
counts <- MERINGUE::cleanCounts(counts = cd,
                                min.reads = 100,
                                min.lib.size = 100,
                                plot=TRUE,
                                verbose=TRUE)

## overdispersed genes only
odGenes <- MERINGUE::getOverdispersedGenes(counts, plot = TRUE, details = TRUE, alpha = 0.05)
countsFilt <- counts[odGenes$ods,]

## remove genes that are present in more than X% of spots
## if we don't remove perplexity poor
vi <- rowSums(countsFilt > 0) > ncol(countsFilt)*0.9
table(vi)
countsFiltnoUnifGenes <- countsFilt[!vi,]

## check clusters
pcs.info <- prcomp(t(as.matrix(MERINGUE::normalizeCounts(countsFiltnoUnifGenes, log=TRUE))), center=TRUE, scale=TRUE)
plot(pcs.info$sdev[1:30])
nPcs <- 5
pcs <- pcs.info$x[,1:nPcs]

## 2D embedding by tSNE
emb <- Rtsne::Rtsne(pcs,
                    is_distance=FALSE,
                    perplexity=30,
                    num_threads=1,
                    verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)

# Graph-based cluster detection
k <- 30
com2 <- getClusters(pcs, k, weight=TRUE)

par(mfrow=c(2,2), mar=rep(1,4))
plotEmbedding(emb, groups=annot,
              show.legend=TRUE, xlab=NA, ylab=NA,
              verbose=FALSE)
plotEmbedding(pos, groups=annot,
              cex=1, xlab=NA, ylab=NA,
              verbose=FALSE)
plotEmbedding(emb, groups=com2,
              show.legend=TRUE, xlab=NA, ylab=NA,
              verbose=FALSE)
plotEmbedding(pos, groups=com2,
              cex=1, xlab=NA, ylab=NA,
              verbose=FALSE)

# topic modeling
mob_corpus <- t(as.matrix(countsFiltnoUnifGenes))
mob_corpus_slamMtx <- slam::as.simple_triplet_matrix(mob_corpus)

# optimal k
Ks = seq(2,20,by=5)
Ks
mob_lda <- fitLDA(mob_corpus_slamMtx, Ks=Ks)
par(mfrow=c(1,1), mar=rep(5,4))
plot(x=Ks, mob_lda$perplexities)

# pick k
kopt = Ks[which(mob_lda$perplexities == min(mob_lda$perplexities))]
kopt = 18
mob_lda.k <- topicmodels::LDA(mob_corpus_slamMtx, k=kopt)
topicmodels::perplexity(mob_lda.k, mob_corpus_slamMtx)

mob_lda.tmResult<- topicmodels::posterior(mob_lda.k)
mob_lda.theta <- mob_lda.tmResult$topics
mob_lda.beta <- mob_lda.tmResult$terms
mob_lda.topicFreqsOverall <- colSums(mob_lda.theta) / nrow(mob_corpus_slamMtx)

# plot
par(mfrow=c(4,3), mar=rep(1,4))
sapply(1:ncol(mob_lda.theta), function(i) {
  plotEmbedding(pos[rownames(mob_lda.theta),], colors=mob_lda.theta[,i], cex=3, main=i)
})

# compare with ground truth annotations
cpm.deconvolved <- t(mob_lda.beta * 1e6)

heatmap(cpm.deconvolved, scale='row')
d <- as.dist(1-cor(cpm.deconvolved ))
#d <- dist(t(cpm.deconvolved))
hc <- hclust(d, method='ward.D2')
groups <- cutree(hc, 4)
groups
groups <- factor(groups)


d2 <- as.dist(1-cor(t(cpm.deconvolved )))
hc2 <- hclust(d2, method='ward.D2')

m <- log10(cpm.deconvolved+1)
m <- t(scale(t(m)))
#m <- scale(m)
range(m)
m[m > 2] <- 2
m[m < -2] <- -2
heatmap(m,
        scale='none',
        ColSideColors = rainbow(length(unique(groups)))[groups],
        Colv = as.dendrogram(hc),
        Rowv = as.dendrogram(hc2),
        col=colorRampPalette(c('blue', 'white', 'red'))(100))

mm <- model.matrix(~ 0 + annot)
colnames(mm) <- levels(annot)
true <- countsFiltnoUnifGenes %*% mm
head(true)
cpm.true <- MERINGUE::normalizeCounts(true, log=FALSE)

head(cpm.deconvolved)
head(cpm.true)

# rmse between each
rmse.results <- do.call(rbind, lapply(1:ncol(cpm.true), function(i) {
  sapply(1:ncol(cpm.deconvolved), function(j) {
    #Metrics::rmse(cpm.true[,i], cpm.deconvolved[,j])
    cor(cpm.true[,i], cpm.deconvolved[,j]) ## correlation instead?
  })
}))
rownames(rmse.results) <- colnames(cpm.true)
colnames(rmse.results) <- colnames(cpm.deconvolved)
head(rmse.results)

heatmap(t(rmse.results), scale='none')

## use hungarian sort algorithm
## to get optimal matching
library(clue)
compare <- rmse.results
order <- clue::solve_LSAP(compare-min(compare), maximum = TRUE)
#heatmap(t(compare[seq_along(order), order]), Rowv=NA, Colv=NA, scale='none')
## also visualize things without matches
no.pairs <-setdiff(as.numeric(colnames(rmse.results)), as.numeric(order))
no.pairs
heatmap(t(compare[, c(no.pairs, order)]), Rowv=NA, Colv=NA, scale='row', mar=c(20,5))

