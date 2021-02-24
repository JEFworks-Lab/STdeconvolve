### Test LDA
load('simulatedMerfishSpots100umMajorCt.RData')
dim(simulatedMerfishSpots[[1]][[1]]) ## words
dim(simulatedMerfishSpots[[1]][[2]]) ## ground truth topics
dim(simulatedMerfishSpots[[1]][[3]]) ## pos

## test out on one first
spotGexp <- simulatedMerfishSpots[[6]][[1]]
spotPos <- simulatedMerfishSpots[[6]][[3]]
dim(spotGexp)

## make corpus
corpus <- as.matrix(spotGexp)
corpus_slamMtx <- slam::as.simple_triplet_matrix(corpus)
head(corpus_slamMtx)

## search for optimal k (needs to be improved)
recursiveFitLDA <- function(corpus, interval, lowerperplexity=NULL, upperperplexity=NULL, seed = 0) {

  controls <- list(seed = seed,
                   verbose = 1, keep = 1,
                   alpha = 1, estimate.alpha = TRUE)

  lower = min(interval);
  mid = round((min(interval) + max(interval))/2);
  upper = max(interval);

  f <- function(k) {
    fitted_model <- topicmodels::LDA(corpus, k=k, control = controls)
    p <- topicmodels::perplexity(fitted_model, corpus)
    return(p)
  }

  if(is.null(lowerperplexity)) {
    fl <- f(lower)
  }
  if(is.null(upperperplexity)) {
    fu <- f(upper)
  }
  fm <- f(mid)

  if(lower == mid | mid == upper) {
    ## optimum found
    return (c(mid, f(mid)))
  }

  if(is.null(lowerperplexity) || is.null(upperperplexity)) {
    out <- rbind(
      c(lower, fl),
      c(mid, fm),
      c(upper, fu)
    )
  } else {
    out <- c(mid, fm)
  }

  if (fl < fu) {
    print('searching lower half')
    return(
      rbind(out, recursiveFitLDA(corpus, interval = c(lower, mid), lowerperplexity = fl, upperperplexity = fm))
    )
  } else {
    print('searching upper half')
    return(
      rbind(out, recursiveFitLDA(corpus, interval = c(mid, upper), lowerperplexity = fm, upperperplexity = fu))
    )
  }

}

## search
perp <- recursiveFitLDA(corpus_slamMtx, interval = c(5, 10))
par(mfrow=c(1,1), mar=rep(5,4))
plot(perp[order(perp[,1]),1], perp[order(perp[,1]),2], type='l')

#kopt = perp[,1][which(perp[,2]==min(perp[,2]))]
kopt = 84
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

############## Performance evaluation
## compare with ground truth cell type proportions
gtProp <- do.call(rbind, lapply(simulatedMerfishSpots, function(x) x[[2]]))
head(gtProp)
gtProp <- gtProp/rowSums(gtProp) ## divide per spot
rowSums(gtProp, na.rm=TRUE) ## double check adds up to 100%
gtProp[is.nan(gtProp)] <- 0
head(gtProp)
rowSums(ldamodel.theta)
head(ldamodel.theta)

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

library(clue)
## transpose as needed
if(nrow(prop.results) > ncol(prop.results)) {
  compare <- t(prop.results)
} else {
  compare <- prop.results
}
order <- clue::solve_LSAP(compare-min(compare), maximum = TRUE)
no.pairs <-setdiff(as.numeric(colnames(compare)), as.numeric(order))
no.pairs
heatmap(t(compare[, c(no.pairs, order)]), Rowv=NA, Colv=NA, scale='none', mar=c(10,5))
heatmap(t(compare[, c(order)]), Rowv=NA, Colv=NA, scale='none', mar=c(10,5))

## compare with ground truth cell type gene expressions
cpm.deconvolved <- t(ldamodel.beta * 1e6)

mm <- model.matrix(~ 0 + celltype)
colnames(mm) <- levels(celltype)
true <- t(as.matrix(counts[names(celltype),])) %*% mm
head(true)
cpm.true <- MERINGUE::normalizeCounts(true, log=FALSE)

gexp.results <- do.call(rbind, lapply(1:ncol(cpm.true), function(i) {
  sapply(1:ncol(cpm.deconvolved), function(j) {
    #Metrics::rmse(cpm.true[,i], cpm.deconvolved[,j])
    cor(cpm.true[,i], cpm.deconvolved[,j]) ## correlation instead?
  })
}))
rownames(gexp.results) <- colnames(cpm.true)
colnames(gexp.results) <- colnames(cpm.deconvolved)
head(gexp.results)

heatmap(t(gexp.results), scale='none')
## use same ordering as previous
heatmap(t(gexp.results[,c(no.pairs, order)]), Rowv=NA, Colv=NA, scale='row', mar=c(10,5))
heatmap(t(gexp.results[,c(order)]), Rowv=NA, Colv=NA, scale='row', mar=c(10,5))

save.image('merfish.RData')

