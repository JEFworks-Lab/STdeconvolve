## Compare latent direchlet allocation
## with hierarchical direchlet allocation

## not sure if this is the best implementation
#devtools::install_github("nicolaroberts/hdp")
library(hdp)

######### Simulation
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
x3[5:10] = base[5:10] + 10

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
dim(train.data.sub)

heatmap(as.matrix(train.data.sub), scale='none', Rowv = NA, Colv = NA)

## HDA
set.seed(10)
quick_hdp <- hdp_quick_init(train.data.sub)
quick_hdp

## posterior sampling
#quick_chain <- hdp_posterior(quick_hdp, burnin=500, n=100, space=40, seed=1234)
quick_chain <- hdp_posterior(quick_hdp, burnin=500, n=1000, space=40, seed=1234)
quick_chain

plot_lik(quick_chain, bty="L")
plot_numcluster(quick_chain, bty="L")
plot_data_assigned(quick_chain, bty="L")

## merge
quick_chain <- hdp_extract_components(quick_chain)
quick_chain

## explore results
## row is number of components (cell-types identified?)
## column is number of categories (genes)
results <- comp_dp_distn(quick_chain)$mean
dim(results)
results <- results[-1,] ## bug it seems

par(mfrow=c(3,3))
plot(results[,1], train.data[, 'ctA'])
plot(results[,1], train.data[, 'ctB'])
plot(results[,1], train.data[, 'ctC'])

plot(results[,2], train.data[, 'ctA'])
plot(results[,2], train.data[, 'ctB'])
plot(results[,2], train.data[, 'ctC'])

plot(results[,3], train.data[, 'ctA'])
plot(results[,3], train.data[, 'ctB'])
plot(results[,3], train.data[, 'ctC'])

## Conclusion: wow near perfect results, much slower though
## Is the poor scalability due to implementation
## or just limitation inherent to method?


################# Real data
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

data <- as.matrix(t(cd.summary[od.genes,]))
quick_hdp <- hdp_quick_init(data)
quick_chain <- hdp_posterior(quick_hdp, burnin=500, n=1000, space=40, seed=1234)

plot_lik(quick_chain, bty="L")
plot_numcluster(quick_chain, bty="L")
plot_data_assigned(quick_chain, bty="L")

quick_chain <- hdp_extract_components(quick_chain)

results <- comp_dp_distn(quick_chain)$mean
dim(results)
dim(cd.summary)
head(results)
results <- results[-1,] ## bug it seems
rownames(results) <- colnames(cd.summary)

par(mfrow=c(3,3))
sapply(1:ncol(results), function(i){
  plotEmbedding(pos[rownames(results),], col=results[,i], cex=1, main=paste0('Proportion of Cell-Type:', i))
})
