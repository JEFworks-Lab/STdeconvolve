## evaluating merfish gene properties to potentially improve marker selection

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
## all
selected_cells <- rownames(counts)
annot.table.sub <- annot.table[selected_cells,]
good_cells <- !is.na(annot.table.sub$Cell_class)
annot.table.sub <- annot.table.sub[good_cells,]
features.sub <- features[rownames(annot.table.sub),]
## get fine cell types
celltype.fine <- paste0(annot.table.sub$Cell_class, ':', annot.table.sub$Neuron_cluster_ID)
names(celltype.fine) <- rownames(annot.table.sub)
celltype.fine <- factor(celltype.fine)
levels(celltype.fine)
table(celltype.fine)

## how many marker genes per cell-type
dg <- MERINGUE::getDifferentialGenes(t(counts)[,good_cells], celltype.fine)

markers <- lapply(dg, function(x) {
  ## significant only
  x <- x[x$Z > 10,]
  #x <- x[x$highest,]
  x <- x[x$fe > 0.5,]
  #x <- x[order(x$Z, decreasing=TRUE),]
  rownames(x)
})
lapply(markers,length)
## ~ 10

## what is the variance of these marker genes
## what is the magnitude of these marker genes
cd <- t(counts)[,good_cells]
vg <- apply(cd, 1, var)
mg <- apply(cd, 1, mean)
par(mfrow=c(1,1))
plot(log10(mg+1), log10(vg+1))


######################### Compare with MOB
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
