############################################### MERFISH
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

## remove all genes with detection worse than blanks
#sort(Matrix::colSums(counts), decreasing=TRUE)
#blanks <- counts[,grepl('Blank', colnames(counts))]
#counts <- counts[, !grepl('Blank', colnames(counts))]
#which(Matrix::colSums(counts) < mean(Matrix::colSums(blanks)))
#counts <- counts[, Matrix::colSums(counts) > mean(Matrix::colSums(blanks))]
#colnames(counts)
#dim(counts)

######### Get one animal
table(features$dataset_name)
vi <- features$dataset_name %in% c('171021_FN7_2_M22_M26', '171023_FN7_1_M22_M26')
selected_cells <- rownames(features)[vi]

annot.table.sub <- annot.table[selected_cells,]
good_cells <- !is.na(annot.table.sub$Cell_class)
annot.table.sub <- annot.table.sub[good_cells,]
features.sub <- features[rownames(annot.table.sub),]

## fine cell types
celltype.fine <- paste0(annot.table.sub$Cell_class, ':', annot.table.sub$Neuron_cluster_ID)
names(celltype.fine) <- rownames(annot.table.sub)
celltype.fine <- factor(celltype.fine)
levels(celltype.fine)

## coarse cell types
celltype <- annot.table.sub$Cell_class
names(celltype) <- rownames(annot.table.sub)
table(celltype)
## combine glial cell types
celltype[grepl('OD', celltype)] <- 'OD'
celltype[grepl('Endo', celltype)] <- 'Endo'
table(celltype)
celltype <- factor(celltype)
levels(celltype)

spatial_position_and_class <- data.frame(annot.table.sub[, c('Centroid_X', 'Centroid_Y', 'Bregma')], Cell_class = celltype)
head(spatial_position_and_class)
dim(spatial_position_and_class)
spatial_position_and_class <- na.omit(spatial_position_and_class) # remove rows with NA
dim(spatial_position_and_class)

######### Plot
library(ggplot2)
ggplot(spatial_position_and_class, aes(x=Centroid_X, y=Centroid_Y, color=Cell_class)) +
  geom_point(size=0.01) +
  facet_wrap(vars(Bregma), scales = 'free') +
  theme_classic() +
  theme(legend.position="none")

## Simulate spots
table(spatial_position_and_class$Bregma)
bids <- sort(unique(spatial_position_and_class$Bregma))
bids
simulatedMerfishSpots <- lapply(bids, function(bid) {
  print(bid)

  selected_bregma <- spatial_position_and_class[which(spatial_position_and_class$Bregma == bid),]
  #ggplot(selected_bregma, aes(x=Centroid_X, y=Centroid_Y, color=Cell_class)) + geom_point()
  foo = data.frame(selected_bregma$Centroid_X, selected_bregma$Centroid_Y)
  bar = selected_bregma$Cell_class
  rownames(foo) <- names(bar) <- rownames(selected_bregma)
  MERINGUE::plotEmbedding(foo, groups=bar, shuffle.colors = TRUE)

  # Sequence of X-coord positions for left edge of each patch:
  x_edges <- seq(min(selected_bregma$Centroid_X), max(selected_bregma$Centroid_X), 100)
  # drop first and last to avoid any issues with the edges of the whole region
  inner_x_edges <- x_edges[2:(length(x_edges)-1)]
  # Sequence of Y-coord positions for bottom edge of each patch:
  y_edges <- seq(min(selected_bregma$Centroid_Y), max(selected_bregma$Centroid_Y), 100)
  inner_y_edges <- y_edges[2:(length(y_edges)-1)]

  selected_bregma$patch_id <- character(length(rownames(selected_bregma)))
  for (x in inner_x_edges) {
    for (y in inner_y_edges) {
      patch_id <- paste0(as.character(x), "_", as.character(y))
      patch <- selected_bregma[which( (selected_bregma$Centroid_X > x) & (selected_bregma$Centroid_X < x+100) & (selected_bregma$Centroid_Y > y) & (selected_bregma$Centroid_Y < y+100) ),]
      selected_bregma[rownames(patch),]$patch_id <- patch_id
    }
  }

  #par(mfrow=c(1,1))
  #ggplot(selected_bregma, aes(x=Centroid_X, y=Centroid_Y, color=patch_id)) +
  #  geom_point() +
  #  theme_classic() +
  #  theme(legend.position="none")

  foo = data.frame(selected_bregma$Centroid_X, selected_bregma$Centroid_Y)
  bar = selected_bregma$patch_id
  rownames(foo) <- names(bar) <- rownames(selected_bregma)
  MERINGUE::plotEmbedding(foo, groups=bar, shuffle.colors = TRUE)

  selected_bregma_patches <- selected_bregma[which(selected_bregma$patch_id != ""),c("patch_id", "Cell_class")]
  selected_bregma_patch_cells <- table(selected_bregma_patches[])
  par(mfrow=c(1,1), mar=rep(5,4))
  hist(rowSums(selected_bregma_patch_cells), breaks = 20)

  unique_types_per_patch <- sapply(seq_len(length(rownames(selected_bregma_patch_cells))), function(i) {
    cell_type_count <- length(which(selected_bregma_patch_cells[i,] != 0))
    cell_type_count
  })
  hist(unique_types_per_patch, breaks=10)

  ## and ground truth cell-type counts
  gtDocTopics <- do.call(rbind, lapply(unique(selected_bregma_patches$patch_id), function(pid) {
    vi <- selected_bregma_patches$patch_id == pid
    ## count cell types in spot
    cc <- selected_bregma_patches[vi, 'Cell_class']
    cct <- table(cc)
    ## need to convert to numeric
    cct2 <- as.numeric(cct)
    names(cct2) <- names(cct)
    ## fill in unrepresented celltypes with 0s
    cct3 <- cct2[levels(celltype)]
    names(cct3) <- levels(celltype)
    cct3[is.na(cct3)] <- 0
    return(cct3)
  }))
  rownames(gtDocTopics) <- unique(selected_bregma_patches$patch_id)
  dim(gtDocTopics)
  head(gtDocTopics)

  ## save gene expression per patch
  gtTopicWords <- do.call(rbind, lapply(unique(selected_bregma_patches$patch_id), function(pid) {
    vi <- selected_bregma_patches$patch_id == pid
    ## cell names in each spot
    cc <- rownames(selected_bregma_patches[vi, ])
    ## sum gexp
    if(length(cc) > 1) {
      totgexp <- Matrix::colSums(counts[cc, ])
    } else if(length(cc) == 1) {
      totgexp <- counts[cc, ]
    } else {
      totgexp <- rep(0, ncol(counts))
    }
    return(totgexp)
  }))
  rownames(gtTopicWords) <- unique(selected_bregma_patches$patch_id)
  dim(gtTopicWords)
  head(gtTopicWords)

  ## save positions
  patchPos <- do.call(rbind, lapply(unique(selected_bregma_patches$patch_id), function(x) as.numeric(strsplit(x, '_')[[1]])))
  rownames(patchPos) <- unique(selected_bregma_patches$patch_id)
  colnames(patchPos) <- c('X', 'Y')
  head(patchPos)

  ## plot some examples
  #MERINGUE::plotEmbedding(patchPos, col=gtTopicWords[,'Gad1'], cex=3)
  #MERINGUE::plotEmbedding(patchPos, col=gtDocTopics[,'Ependymal:'], cex=3)

  out <- list(gtTopicWords, gtDocTopics, data.frame(patchPos, bid))
  return(out)
})
names(simulatedMerfishSpots) <- bids
length(simulatedMerfishSpots)
length(simulatedMerfishSpots[[1]])
save(simulatedMerfishSpots, file='simulatedMerfishSpots100umMajorCt.RData')

########## Test LDA
#load('simulatedMerfishSpots100um.RData')
load('simulatedMerfishSpots100umMajorCt.RData')
dim(simulatedMerfishSpots[[1]][[1]]) ## words
dim(simulatedMerfishSpots[[1]][[2]]) ## ground truth topics
dim(simulatedMerfishSpots[[1]][[3]]) ## pos

## concat together
head(simulatedMerfishSpots[[1]][[1]])
length(simulatedMerfishSpots)
spotGexp <- do.call(rbind, lapply(simulatedMerfishSpots, function(x) x[[1]]))
spotPos <- do.call(rbind, lapply(simulatedMerfishSpots, function(x) x[[3]]))
rownames(spotPos) <- rownames(spotGexp)

## test out on one first
#spotGexp <- simulatedMerfishSpots[[6]][[1]]
#spotPos <- simulatedMerfishSpots[[6]][[3]]
#dim(spotGexp)

## make corpus
corpus <- as.matrix(spotGexp)
corpus_slamMtx <- slam::as.simple_triplet_matrix(corpus)
head(corpus_slamMtx)

## search
#perp <- recursiveFitLDA(corpus_slamMtx, interval = c(5, 10))
Ks <- seq(10, 200, by=10)
perp <- fitLDA(corpus_slamMtx, Ks)
par(mfrow=c(1,1), mar=rep(5,4))
plot(Ks, perp$perplexities, type='l')
abline(v=45,col='red')

#kopt = perp[,1][which(perp[,2]==min(perp[,2]))]
kopt = 45
#kopt = length(levels(celltype.fine)) - 1
kopt
ldamodel <- topicmodels::LDA(corpus_slamMtx, k=kopt) ## is there some progress bar?
topicmodels::perplexity(ldamodel, corpus_slamMtx)

ldamodel.tmResult<- topicmodels::posterior(ldamodel)
ldamodel.theta <- ldamodel.tmResult$topics
ldamodel.beta <- ldamodel.tmResult$terms

## plot a few examples
par(mfrow=c(5,5), mar=rep(1,4))
sapply(1:ncol(ldamodel.theta), function(i) {
  MERINGUE::plotEmbedding(spotPos[rownames(ldamodel.theta),1:2], colors=ldamodel.theta[,i], cex=0.1, main=i)
})

############## Cluster
library(dynamicTreeCut)
ldamodel.clust <- clusterTopics(beta = ldamodel.beta,
                               distance = "euclidean",
                               clustering = "ward.D",
                               dynamic = "hybrid",
                               deepSplit = 4,
                               plotDendro = TRUE)
clusterColorsk <- ldamodel.clust$clusters
levels(clusterColorsk) <- gg_color_hue(length(levels(clusterColorsk)))

m <- t(ldamodel.beta)
range(m)
m <- scale(m)
#m <- t(scale(t(m)))
range(m)
m[m < -2] <- -2
m[m > 2] <- 2
rc <- hclust(as.dist(1-cor(t(m))), method='ward.D2')
par(mfrow=c(1,1), mar=c(8,8,3,2))
gplots::heatmap.2(m,
                  Colv = ldamodel.clust$dendro,
                  Rowv = as.dendrogram(rc),
                  density.info = "none",
                  trace = "none",
                  ColSideColors = as.vector(clusterColorsk),
                  col = colorRampPalette(c('blue', 'white', 'red'))(100),
                  cexRow=0.5,cexCol=0.7,margins=c(6,3),
                  main = "ldamodel.beta",
                  lhei = c(1,5),
                  key.title = NA,
                  scale = 'none')

ldamodel.betaCombined <- combineTopics(mtx = ldamodel.beta, clusters = clusterColorsk, mtxType = "b")

cpmCombined.deconvolved <- t(ldamodel.betaCombined * 1e6)

mm <- model.matrix(~ 0 + celltype)
colnames(mm) <- levels(celltype)
true <- t(as.matrix(counts[names(celltype),])) %*% mm
head(true)
cpm.true <- MERINGUE::normalizeCounts(true, log=FALSE)
cpm.true <- cpm.true[, which(colnames(cpm.true) != 'Ambiguous')]
dim(cpm.true)

hc <- hclust(as.dist(1-cor(as.matrix(cpm.true))), method='ward')
plot(hc)
cpm.true <- cpm.true[,rev(hc$order)]

gexp.results <- do.call(rbind, lapply(1:ncol(cpm.true), function(i) {
  sapply(1:ncol(cpmCombined.deconvolved), function(j) {
    #Metrics::rmse(cpm.true[,i], cpm.deconvolved[,j])
    cor(cpm.true[,i], cpmCombined.deconvolved[,j]) ## correlation instead?
  })
}))
rownames(gexp.results) <- colnames(cpm.true)
colnames(gexp.results) <- colnames(cpmCombined.deconvolved)
head(gexp.results)

heatmap(gexp.results, scale='none')

############## All identified topics
## compare with ground truth cell type gene expressions
#cpm.deconvolved <- t(ldamodel.beta * 1e6)
cpm.deconvolved <- t(ldamodel.beta)

mm <- model.matrix(~ 0 + celltype.fine)
colnames(mm) <- levels(celltype.fine)
true <- t(as.matrix(counts[names(celltype.fine),])) %*% mm
head(true)
#cpm.true <- MERINGUE::normalizeCounts(true, log=FALSE)
cpm.true <- true ## no normalization
cpm.true <- cpm.true[, which(colnames(cpm.true) != 'Ambiguous:')]
dim(cpm.true)

hc <- hclust(as.dist(1-cor(as.matrix(cpm.true))), method='ward.D')
plot(hc)
cpm.true <- cpm.true[,rev(hc$order)]

gexp.results <- do.call(rbind, lapply(1:ncol(cpm.true), function(i) {
  sapply(1:ncol(cpm.deconvolved), function(j) {
    #Metrics::rmse(cpm.true[,i], cpm.deconvolved[,j])
    cor(cpm.true[,i], cpm.deconvolved[,j]) ## correlation instead?
  })
}))
rownames(gexp.results) <- colnames(cpm.true)
colnames(gexp.results) <- colnames(cpm.deconvolved)
head(gexp.results)

heatmap(gexp.results, scale='none')

library(clue)
compare <- t(gexp.results)
order <- clue::solve_LSAP(compare-min(compare), maximum = TRUE)
#order2 <- unlist(lapply(rev(hc$order), function(i) which(order==i)))
#m <- t(gexp.results)[order2,rev(hc$order)]
m <- t(gexp.results)[, order]
#m <- scale(m)
#m <- t(scale(t(m)))
range(m)
m[m < 0] <- 0
#m[m < -3] <- -3
#m[m > 3] <- 3
heatmap(t(m), Rowv=NA, Colv=NA, scale='none', mar=c(10,5), col=colorRampPalette(c('blue', 'white', 'red'))(100))

library(gplots)
rc = hclust(as.dist(1-cor(t(m))), method='ward.D2')
hc = hclust(as.dist(1-cor(m)), method='ward.D2')
gplots::heatmap.2(t(m), Colv=as.dendrogram(rc), Rowv=as.dendrogram(hc),
                  trace='none', scale='none', mar=c(10,5), col=colorRampPalette(c('blue', 'white', 'red'))(100))

## reorder within groups?
groups1 <- cutree(rc, 3)
groups2 <- cutree(hc, 3)

od1 <- rc$order
od2 <- hc$order

## compare with ground truth cell type proportions
gtProp <- do.call(rbind, lapply(simulatedMerfishSpots, function(x) x[[2]]))
head(gtProp)
gtProp <- gtProp/rowSums(gtProp) ## divide per spot
rowSums(gtProp, na.rm=TRUE) ## double check adds up to 100%
gtProp[is.nan(gtProp)] <- 0
head(gtProp)
rowSums(ldamodel.theta)
head(ldamodel.theta)
## remove ambiguous
gtProp <- gtProp[, which(colnames(gtProp) != 'Ambiguous:')]

## correlation across spots
prop.results <- do.call(rbind, lapply(1:ncol(gtProp), function(i) {
  sapply(1:ncol(ldamodel.theta), function(j) {
    cor(gtProp[,i], ldamodel.theta[,j])
  })
}))
rownames(prop.results) <- colnames(gtProp)
colnames(prop.results) <- colnames(ldamodel.theta)
head(prop.results)

prop.results[is.na(prop.results)] <- 0

heatmap(t(prop.results), scale='none')

m2 <- prop.results[rownames(m), colnames(m)]
m2 <- scale(m2)
m2 <- t(scale(t(m2)))
range(m2)
m2[m2 < -2.5] <- -2.5
m2[m2 > 2.5] <- 2.5
heatmap(t(m2), Rowv=NA, Colv=NA, scale='none', mar=c(10,5),
        col=colorRampPalette(c('blue', 'white', 'red'))(100))


#save.image('merfish_45topics.RData')
#load('merfish_209topics.RData')
