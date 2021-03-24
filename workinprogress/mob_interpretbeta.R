## get single cell MOB data

cd <- read.csv('~/Downloads/GSE121891_OB_6_runs.raw.dge.csv.gz')
head(cd)
colnames(cd)

meta <- read.csv('~/Downloads/GSE121891_OB_metaData_seurat.csv.gz', row.names = 1)
head(meta)

wt <- rownames(meta)[which(meta$orig.ident %in% c("WT1", "WT2"))]
head(wt)

table(wt %in% colnames(cd))
cd.wt <- cd[, wt]
cd.wt[1:5,1:5]

annot.wt <- meta[wt,]$ClusterName
names(annot.wt) <- wt
annot.wt <- factor(annot.wt)
head(annot.wt)

## sum expression within cell-types
mm <- model.matrix(~ 0 + annot.wt)
colnames(mm) <- levels(annot.wt)
counts.summary.mm <- as.matrix(cd.wt) %*% mm
head(counts.summary.mm)
dim(counts.summary.mm)


################ LDA model
library(MERINGUE)
library(Matrix)
data(mOB)
pos <- mOB$pos
cd <- mOB$counts
annot <- mOB$annot

## Remove poor datasets and genes
par(mfrow=c(1,1), mar=rep(5,4))
counts <- MERINGUE::cleanCounts(counts = cd,
                                min.reads = 200,
                                min.lib.size = 100,
                                plot=TRUE,
                                verbose=TRUE)

## overdispersed genes only
odGenes <- MERINGUE::getOverdispersedGenes(counts, plot = TRUE, details = TRUE, alpha = 0.1)
countsFilt <- counts[odGenes$ods,]

## remove genes that are present in more than X% of spots
## if we don't remove perplexity poor
vi <- rowSums(countsFilt > 0) > ncol(countsFilt)*0.95
table(vi)
countsFiltnoUnifGenes <- countsFilt[!vi,]
dim(countsFiltnoUnifGenes)

## topic modeling
mob_corpus <- as.matrix(countsFiltnoUnifGenes)
mob_corpus_slamMtx <- slam::as.simple_triplet_matrix(mob_corpus)

## optimal k (use known for know)
Ks = length(levels(annot.wt))
Ks
mob_lda <- fitLDA(mob_corpus_slamMtx, Ks=Ks)
mob_lda.tmResult<- topicmodels::posterior(mob_lda$models[[1]])
mob_lda.theta <- mob_lda.tmResult$topics
mob_lda.beta <- mob_lda.tmResult$terms

## scale reference to million then log transform
good.genes <- intersect(colnames(mob_lda.beta), rownames(counts.summary.mm))
mat.summary <- MERINGUE::normalizeCounts(counts.summary.mm[good.genes,], log=TRUE)
## scale to million and log transform
beta.cpm <- MERINGUE::normalizeCounts(t(mob_lda.beta[,good.genes]), log=TRUE)
beta.cpm <- as.matrix(beta.cpm)
range(mat.summary)
range(beta.cpm)

## compare beta to average gene expression
## take one well characterize cell-type
levels(annot.wt)
ct <- "OEC1"
ct <- "OPC"
ct <- "N1"
ct.gexp <- mat.summary[,ct]
## which topic is best correlated
ct.cor <- cor(ct.gexp, beta.cpm)[1,]
ct.order <- ct.cor[order(ct.cor, decreasing=TRUE)]
topic.test <- names(ct.order)[1]
topic.test
## look at correspondence
topic.test.cpm <- beta.cpm[, topic.test]
par(mfrow=c(1,1), mar=rep(5,4))
plot(ct.gexp, topic.test.cpm, main=paste0(ct, '\n', ct.order[topic.test]))

# plot
plotEmbedding(pos[rownames(mob_lda.theta),], colors=mob_lda.theta[,topic.test], cex=3, main=topic.test)

