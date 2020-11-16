library(reticulate)

## previously in bash
## > conda create -n tomatopy python=3.8
## > conda activate tomatopy
## > pip install tomotopy
use_condaenv("tomatopy", required = TRUE)
tp = import("tomotopy")
print(tp$isa)

## baseline expression
set.seed(0)
base = abs(round(rnorm(10, 10)))
names(base) <- paste0('gene', 1:10)
head(base)
## cell-type A upregulates genes 1 to 5
x1 = base
x1[1:5] = base[1:5] + 10
## cell-type B upregulates genes 6 to 10
x2 = base
x2[6:10] = base[6:10] + 10
## visualize transcriptional distinctness
ct.gexp <- cbind(x1, x2)
ct.gexp
heatmap(ct.gexp, scale='none')
## create test
data <- do.call(rbind, lapply(seq(1,100, by=5), function(i) {
  ground.truth <- c(i, 100-i)
  names(ground.truth) <- c('ctA', 'ctB')
  y = ground.truth[1]*x1 + ground.truth[2]*x2
  c(y, ground.truth)
}))
rownames(data) <- paste0('sim', 1:nrow(data))
ct.prop <- data[,11:12]
ct.prop <- ct.prop/rowSums(ct.prop)
spots <- data[, 1:10]
heatmap(spots, scale='none', Rowv=NA, Colv=NA)

## LDA in tomatopy
#mdl = tp$LDAModel(k=2L) ## integer
mdl = tp$HDPModel()
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
mdl$train(100L)

## TODO: automate by number of topics
## first topic
r1 = do.call(rbind, lapply(mdl$get_topic_words(0L, top_n=10L), function(x) unlist(x)))
## second topic
r2 = do.call(rbind, lapply(mdl$get_topic_words(1L, top_n=10L), function(x) unlist(x)))
rownames(r1) <- r1[,1]; r1 = r1[,-1]
rownames(r2) <- r2[,1]; r2 = r2[,-1]
r1
r2

## expression of cell-types
results <- data.matrix(cbind(
  as.numeric(r1[colnames(data)[1:10]]),
  as.numeric(r2[colnames(data)[1:10]])
  ))
rownames(results) <- colnames(data)[1:10]
colnames(results) <- paste0('ct', 1:2)
results

heatmap(results, scale='none')

## compare with ground truth
compare <- do.call(rbind, lapply(1:ncol(ct.gexp), function(i) {
  sapply(1:ncol(results), function(j) {
    cor(ct.gexp[,i], results[,j])
  })
}))
rownames(compare) <- colnames(ct.gexp)
colnames(compare) <- colnames(results)
head(compare)
heatmap(compare, scale='none')
order <- clue::solve_LSAP(compare-min(compare), maximum = TRUE)
heatmap(compare[seq_along(order), order], Rowv=NA, Colv=NA, scale='row')

## proportion across spots
results2 <- do.call(rbind, lapply(1:nrow(data), function(d) {
  words = spots[d,]
  doc = unlist(lapply(1:length(words), function(i) {
    print(i)
    rep(names(words)[i], words[i])
  }))
  table(doc)
  doc_inst = mdl$make_doc(doc)
  results = mdl$infer(doc_inst)
  results[[1]]
}))
rownames(results2) <- rownames(spots)
colnames(results2) <- paste0('ct', 1:2)
results2

## compare with ground truth
compare <- do.call(rbind, lapply(1:ncol(ct.prop), function(i) {
  sapply(1:ncol(results2), function(j) {
    cor(ct.prop[,i], results2[,j])
  })
}))
rownames(compare) <- colnames(ct.gexp)
colnames(compare) <- colnames(results2)
head(compare)
heatmap(compare, scale='none')
order <- clue::solve_LSAP(compare-min(compare), maximum = TRUE)
heatmap(compare[seq_along(order), order], Rowv=NA, Colv=NA, scale='row')


