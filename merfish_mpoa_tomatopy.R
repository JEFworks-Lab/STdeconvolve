library(reticulate)

######### Use Tomatopy
## previously in bash
## > conda create -n tomatopy python=3.8
## > conda activate tomatopy
## > pip install tomotopy
use_condaenv("tomatopy", required = TRUE)
tp = import("tomotopy")
print(tp$isa)

######### Load MERFISH
load('~/Desktop/mpoa_merfish_clean.RData')
## contains counts, annot.table, and features

## focus on one animal
table(features$dataset_name)
vi <- features$dataset_name %in% c('171023_FN7_1_M22_M26', '171021_FN7_2_M22_M26') ## whole animal
cells.good <- rownames(features)[vi]
spatial_position <- annot.table[cells.good,c('Centroid_X', 'Centroid_Y', 'Bregma')]
annot <- paste0(annot.table[cells.good, 'Cell_class'], ':', annot.table[cells.good, 'Neuron_cluster_ID'])
names(annot) <- cells.good

cells.have <- intersect(cells.good, rownames(counts))
annot <- annot[cells.have]
annot <- as.factor(annot)
pos = spatial_position[cells.have, 1:2]
ct = annot[cells.have]

slice = paste0(features[cells.have,]$dataset_name, '-', features[cells.have,]$brain_pos)
slice <- factor(slice)
par(mfrow=c(4,3), mar=rep(1,4))
lapply(levels(slice), function(s) {
  vi <- slice == s
  plot(pos[vi,], col=MERINGUE:::fac2col(ct, v=0.8, s=0.8), pch=16, cex=0.3,
       axes=FALSE, main=s); box()
})

## remove blanks
good.genes <- colnames(counts)[!grepl('Blank', colnames(counts))]
good.genes
cd <- counts[cells.have, good.genes]

########### Simulate spots by grouping by FOV
fov = paste0(features[cells.have,]$dataset_name, '-', features[cells.have,]$primary_fovID)
names(fov) <- cells.have
fov <- as.factor(fov)
library(MUDAN)
par(mfrow=c(4,3))
lapply(levels(slice), function(s) {
  vi <- slice == s
  plotEmbedding(pos[vi,], groups=fov, shuffle.colors = TRUE, cex=0.3, main=s)
})

## Sum up expression by FOV to simulate spots
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

############################# Topic modeling in tomatopy
#mdl = tp$HDPModel()
mdl = tp$LDAModel(k=as.integer(length(levels(annot))))
## add document to model
sapply(1:nrow(spots), function(d) {
  #print(d)
  words = spots[d,]
  doc = unlist(lapply(1:length(words), function(i) {
    #print(i)
    rep(names(words)[i], words[i])
  }))
  #table(doc)
  mdl$add_doc(doc)
})
mdl$train(100L) ## does running longer lead to better performance?

###### Expression similarity
dim(spots)
#num_topics <- mdl$live_k ## for HDP
num_topics <- length(levels(annot))-1
num_topics
## python uses 0 index
output <- do.call(cbind, lapply(0:num_topics, function(i) {
  print(i)
  out <- do.call(rbind, lapply(mdl$get_topic_words(as.integer(i), top_n=125L), function(x) unlist(x)))
  rownames(out) <- out[,1]
  out <- out[,-1]
  as.numeric(out[colnames(spots)])
}))
rownames(output) <- colnames(spots)
colnames(output) <- paste0('ct', 1:ncol(output))
output <- t(output)
head(output)

## ground truth
mm <- model.matrix(~ 0 + annot)
colnames(mm) <- levels(annot)
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
compare2[is.na(compare2)] <- 0
heatmap(compare2, scale='none')

## optimal match
#compare2 <- t(compare2)
order <- clue::solve_LSAP(compare2-min(compare2), maximum = TRUE)
heatmap(compare2[seq_along(order), order], Rowv=NA, Colv=NA, scale='none')
no.pairs <-setdiff(1:length(levels(annot.main)), as.numeric(order))
heatmap(compare2[, c(order, no.pairs)], Rowv=NA, Colv=NA, scale='none')

## look at specific genes
g <- 'Gad1'
par(mfrow=c(1,1), mar=rep(5,4))
plot(ct.gexp[, g], output[as.numeric(order),g], main=g)

####### Proportion
## proportion across spots (do we need to recompute again?)
results2 <- do.call(rbind, lapply(1:nrow(spots), function(d) {
  words = spots[d,]
  doc = unlist(lapply(1:length(words), function(i) {
    #print(i)
    rep(names(words)[i], words[i])
  }))
  #table(doc)
  doc_inst = mdl$make_doc(doc)
  results = mdl$infer(doc_inst)
  results[[1]]
}))
dim(results2)
dim(spots)
rownames(results2) <- rownames(spots)
colnames(results2) <- paste0('ct', 1:ncol(results2))
results2

ground.truth <- table(fov,annot[names(fov)])
ground.truth = ground.truth/rowSums(ground.truth) ## proportion
head(ground.truth)
ct.prop <- ground.truth

## compare with ground truth
compare <- do.call(rbind, lapply(1:ncol(ct.prop), function(i) {
  sapply(1:ncol(results2), function(j) {
    cor(ct.prop[,i], results2[,j])
  })
}))
dim(compare)
rownames(compare) <- colnames(ct.prop)
colnames(compare) <- colnames(results2)
head(compare)
heatmap(compare, scale='none')
order <- clue::solve_LSAP(compare-min(compare), maximum = TRUE)
heatmap(compare[seq_along(order), order], Rowv=NA, Colv=NA, scale='none')

