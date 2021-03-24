## Interpreting beta values

## load old model
load('merfish_209topics.RData')
## 209 topics
## 125 genes
dim(ldamodel.beta)

## Get one animal
load("~/Desktop/mpoa_merfish_clean.RData")
head(annot.table)
head(celltype)
head(features)
head(counts)
colnames(counts)
## remove blanks
counts <- counts[, !grepl('Blank', colnames(counts))]
colnames(counts)
dim(counts)
## get one animal
table(features$dataset_name)
vi <- features$dataset_name %in% c('171021_FN7_2_M22_M26', '171023_FN7_1_M22_M26')
selected_cells <- rownames(features)[vi]
annot.table.sub <- annot.table[selected_cells,]
good_cells <- !is.na(annot.table.sub$Cell_class)
annot.table.sub <- annot.table.sub[good_cells,]
features.sub <- features[rownames(annot.table.sub),]
## get fine cell types
celltype.fine <- paste0(annot.table.sub$Cell_class, ':', annot.table.sub$Neuron_cluster_ID)
names(celltype.fine) <- rownames(annot.table.sub)
celltype.fine <- factor(celltype.fine)
levels(celltype.fine)
head(celltype.fine)

############ Compare beta value to average gene expression for each cell-type

## sum expression within cell-types
mm <- model.matrix(~ 0 + celltype.fine)
colnames(mm) <- levels(celltype.fine)
counts.summary.mm <- Matrix::t(counts)[, names(celltype.fine)] %*% mm
head(counts.summary.mm)
dim(counts.summary.mm)
## scale to million then log transform
mat.summary <- MERINGUE::normalizeCounts(counts.summary.mm, log=TRUE)
range(mat.summary)
sum(mat.summary[,1])

## beta value are probabilities of gene expression for each topic
## check that it sums to 1
sum(ldamodel.beta[1,])
## scale to million and log transform
beta.cpm <- log10(ldamodel.beta*1e6+1)

## compare beta to average gene expression
## take one well characterize cell-type
levels(celltype.fine)
ct <- 'Astrocyte:'
ct <- "Inhibitory:I-2"
ct <- "OD Immature 1:"
ct <- "Pericytes:"
ct <- "Ependymal:"
ct.gexp <- mat.summary[,ct]
## which topic is best correlated
ct.cor <- cor(ct.gexp, t(beta.cpm))[1,]
ct.order <- ct.cor[order(ct.cor, decreasing=TRUE)]
topic.test <- names(ct.order)[1]
topic.test
## look at correspondence
topic.test.cpm <- beta.cpm[topic.test, ]
plot(ct.gexp, topic.test.cpm, main=paste0(ct, '\n', ct.order[topic.test]))


######### Compare with differential expression z-scores and p-values
counts.sub <- Matrix::t(counts)[, names(celltype.fine)] ## just this animal
## remove rare cell-types
bad.ct <- names(which(table(celltype.fine) < 30))
celltype.test <- celltype.fine[!(celltype.fine %in% bad.ct)]
dg <- MERINGUE::getDifferentialGenes(counts.sub[, names(celltype.test)], factor(celltype.test))

## compare beta to Z-score
levels(celltype.fine)
ct <- "Ependymal:"
ct.Z <- dg[[ct]]$Z

## which topic is best correlated to raw beta values
ct.cor <- cor(ct.Z, t(ldamodel.beta), method='spearman')[1,]
ct.order <- ct.cor[order(ct.cor, decreasing=TRUE)]
topic.test <- names(ct.order)[1]
topic.test
## look at correspondence
topic.test.beta <- ldamodel.beta[topic.test, ]
plot(ct.Z, topic.test.beta, main=paste0(ct, '\n', ct.order[topic.test]))

## compare beta to -log10(p-value)
levels(celltype.fine)
ct <- 'Astrocyte:'
ct.pv <- -log10(dg[[ct]]$p.adj)

## which topic is best correlated to raw beta values
ct.cor <- cor(ct.Z, t(ldamodel.beta), method='pearson')[1,]
ct.order <- ct.cor[order(ct.cor, decreasing=TRUE)]
topic.test <- names(ct.order)[1]
topic.test
## look at correspondence
topic.test.beta <- ldamodel.beta[topic.test, ]
plot(ct.pv, topic.test.beta, main=paste0(ct, '\n', ct.order[topic.test]))
