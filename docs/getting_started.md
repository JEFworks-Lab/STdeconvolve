Load package

``` r
library(STdeconvolve)
```

Load built in data

``` r
data(mOB)
pos <- mOB$pos
cd <- mOB$counts
annot <- mOB$annot
```

Feature select

``` r
counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 100)
```

![](getting_started_files/figure-markdown_github/getting_started_feature-1.png)

``` r
corpus <- restrictCorpus(counts, alpha = 0.05, t = 1)
```

    ## [1] "Calculating variance fit ..."
    ## [1] "Using gam with k=5..."
    ## [1] "255 overdispersed genes ... "
    ## [1] "Removing genes present in more than 100% of datasets..."
    ## [1] "237 genes remaining..."

Choose optimal number of cell-types

``` r
lda <- fitLDA(corpus, Ks = seq(2, 9, by = 1), plot=TRUE)
```

    ## [1] "Time to train LDA models was 0.24mins"
    ## [1] "Computing perplexity for each fitted model..."
    ## final e step document 260
    ## final e step document 260
    ## final e step document 260
    ## final e step document 260
    ## final e step document 260
    ## final e step document 260
    ## final e step document 260
    ## final e step document 260

![](getting_started_files/figure-markdown_github/getting_started_opt-1.png)

Get best model results

``` r
results <- getBetaTheta(lda$models[["8"]])
```

Visualize deconvolved cell-type proportions

``` r
vizAllTopics(results$theta, pos, 
             groups = annot, 
             group_cols = rainbow(length(levels(annot))),
             r=0.4)
```

![](getting_started_files/figure-markdown_github/getting_started_proportions-1.png)

Visualize deconvolved cell-type gene expression

``` r
marker <- names(sort(results$beta[1, ], decreasing=TRUE))[1]

mat <- normalizeCounts(counts)
gexp <- MERINGUE:::winsorize(scale(mat[marker,])[,1])
MERINGUE::plotEmbedding(pos[colnames(mat),], col=gexp, cex=3)

## function needs to be simplified
df <- merge(as.data.frame(pos), as.data.frame(gexp), by = 0)
rownames(df) <- df$Row.names
df <- df[,-1]
colnames(df) <- c('x', 'y', marker)
vizGeneCounts(df = df,
              gene = marker,
              size = 5, stroke = 0.01,
              plotTitle = marker,
              showLegend = TRUE)
```
