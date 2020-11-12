load('~/Desktop/mpoa_merfish_clean.RData')
## contains counts, annot.table, and features

## focus on one animal
table(features$dataset_name)
#vi <- features$dataset_name %in% c('171023_FN7_1_M22_M26', '171021_FN7_2_M22_M26') ## whole animal
vi <- features$dataset_name %in% c('171021_FN7_2_M22_M26') ## one dataset (half of animal)

cells.good <- rownames(features)[vi]
length(cells.good)
table(features[cells.good,]$brain_pos)

## focus on one tissue section (faster testing)
#vi <- features[cells.good,]$brain_pos == 10
#table(vi)
#cells.good <- cells.good[vi]
#table(features[cells.good,]$brain_pos)

head(annot.table)
#spatial_position <- features[cells.good,c('centroid_1', 'centroid_2', 'brain_pos')]
spatial_position <- annot.table[cells.good,c('Centroid_X', 'Centroid_Y', 'Bregma')]
annot <- paste0(annot.table[cells.good, 'Cell_class'], ':', annot.table[cells.good, 'Neuron_cluster_ID'])
annot.main <- as.character(annot.table[cells.good, 'Cell_class'])
annot.main[grepl('OD', annot.main)] <- 'OD'
annot.main[grepl('Endothelia', annot.main)] <- 'Endothelia'
table(annot.main)
names(annot.main) <- names(annot) <- cells.good
annot <- factor(annot)
annot.main <- factor(annot.main)

par(mfrow=c(1,1))
plot(spatial_position[,1:2], col=MERINGUE:::fac2col(annot, v=0.8, s=0.8), pch=16, cex=0.5)

################## Create spatial pools of data
cells.have <- intersect(cells.good, rownames(counts))
length(cells.have)
length(cells.good)

annot <- annot[cells.have]
annot.main <- annot.main[cells.have]
pos = spatial_position[cells.have, 1:2]
ct = annot[cells.have]
plot(pos, col=MERINGUE:::fac2col(ct, v=0.8, s=0.8), pch=16, cex=0.5)

fov = paste0(features[cells.have,]$dataset_name, '-', features[cells.have,]$primary_fovID)
names(fov) <- cells.have
fov <- as.factor(fov)

## remove blanks
good.genes <- colnames(counts)[!grepl('Blank', colnames(counts))]
good.genes
cd <- counts[cells.have, good.genes]

library(MUDAN)
plotEmbedding(pos, groups=fov, shuffle.colors = TRUE, cex=0.5)
length(table(fov))

## Sum up expression by FOV to simulate spots
## https://jef.works/blog/2020/04/06/quickly-creating-pseudobulks/
## TODO: create function that spits cells into NxN grid
## so that we can sum up expression by NxN grid to
## simulate spots with more or fewer cells
## and see if that impacts performance
## (like for larger human samples with fewer cells)
mm <- model.matrix(~ 0 + fov)
colnames(mm) <- levels(fov)
dim(cd)
dim(mm)
spots <-  t(t(cd) %*% mm)
head(spots)
dim(spots)

## average cell position within fov as spot positions
mm <- model.matrix(~ 0 + fov)
colnames(mm) <- levels(fov)
spots.pos <- t(t(pos) %*% mm)/colSums(mm)
head(spots.pos)

## Run hierarchical dirichlet process
library(hdp)
quick_hdp <- hdp_quick_init(as.matrix(spots))
quick_chain <- hdp_posterior(quick_hdp, burnin=100, n=100, space=10, seed=1234) ## SSLOOOWWW, how to speed up?

## check chain convergence
par(mfrow=c(1,3))
plot_lik(quick_chain, bty="L")
plot_numcluster(quick_chain, bty="L")
plot_data_assigned(quick_chain, bty="L")

#################### Evaluate performance
## extract results for cell-type proportions
quick_chain <- hdp_extract_components(quick_chain)
results <- comp_dp_distn(quick_chain)$mean
dim(results)
dim(spots)
head(results)
results <- results[-1,] ## bug it seems
rownames(results) <- rownames(spots)
colnames(results) <- paste0('ct', 1:ncol(results))
head(results)
dim(results)

## proportion of cell-types per spot
#ground.truth <- table(fov,annot[names(fov)])
ground.truth <- table(fov,annot.main[names(fov)]) ## count
ground.truth = ground.truth/rowSums(ground.truth) ## proportion
head(ground.truth)
head(results)

## compare via correlation for each predicted cell-type
## to each real cell-type in terms of proportions
## across all spots
compare <- do.call(rbind, lapply(1:ncol(ground.truth), function(i) {
  sapply(1:ncol(results), function(j) {
    cor(ground.truth[,i], results[,j])
    #dist(rbind(ground.truth[,i], results[,j]))
  })
}))
table(is.na(compare))
#compare[is.na(compare)] <- 0
rownames(compare) <- levels(annot.main)
#rownames(compare) <- levels(annot)
colnames(compare) <- colnames(results)
heatmap(compare)

## focus only on common cell-types
#good.ct <- sort(table(annot.main[names(fov)]), decreasing=TRUE)
#good.ct <- good.ct[!(names(good.ct) %in% c('Ambiguous', 'Unstable'))]
#good.ct <- names(good.ct)
#good.ct <- names(good.ct[1:ncol(results)])
#good.ct
#heatmap(compare[good.ct,])

## use hungarian sort algorithm
## to get optimal matching
library(clue)
#compare <- t(compare)
order <- clue::solve_LSAP(compare-min(compare), maximum = TRUE)
heatmap(compare[seq_along(order), order], Rowv=NA, Colv=NA, scale='row')
## also visualize things without matches
no.pairs <-setdiff(1:length(levels(annot.main)), as.numeric(order))
heatmap(compare[, c(order, no.pairs)], Rowv=NA, Colv=NA, scale='row')

######## Compare in terms of expression similarity
output <- comp_categ_distn(quick_chain)$mean
dim(spots)
dim(output)
rownames(output) <- paste0('ct', 1:nrow(output))
colnames(output) <- colnames(spots)
head(output)

## ground truth
mm <- model.matrix(~ 0 + annot.main)
colnames(mm) <- levels(annot.main)
dim(cd)
dim(mm)
ct.gexp <-  t(t(cd) %*% mm)
ct.gexp <- ct.gexp/rowSums(ct.gexp)
head(ct.gexp)

compare2 <- do.call(rbind, lapply(1:nrow(ct.gexp), function(i) {
  sapply(1:nrow(output), function(j) {
    cor(ct.gexp[i,], output[j,])
  })
}))
rownames(compare2) <- rownames(ct.gexp)
colnames(compare2) <- rownames(output)
head(compare2)
heatmap(compare2)

## optimal match
#compare2 <- t(compare2)
order <- clue::solve_LSAP(compare2-min(compare2), maximum = TRUE)
heatmap(compare2[seq_along(order), order], Rowv=NA, Colv=NA, scale='row')
no.pairs <-setdiff(1:length(levels(annot.main)), as.numeric(order))
heatmap(compare[, c(order, no.pairs)], Rowv=NA, Colv=NA, scale='row')

