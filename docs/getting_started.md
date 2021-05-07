In this tutorial, we will walk through some of the main functionalities
of `STdeconvolve`.

``` r
library(STdeconvolve)
```

Given a counts matrix from pixel-resolution spatial transcriptomics data
where each spatially resolved measurement may represents mixtures from
potentially multiple cell-types, STdeconvolve infers the putative
transcriptomic profiles of cell-types and their proportional
representation within each multi-cellular spatially resolved pixel. Such
a pixel-resolution spatial transcriptomics dataset of the mouse
olfactory bulb is built in and can be loaded.

``` r
data(mOB)
pos <- mOB$pos
cd <- mOB$counts
annot <- mOB$annot
```

`STdeconvolve` first feature selects for genes most likely to be
relevant for distinguishing between cell-types by looking for highly
overdispersed genes across ST pixels. Pixels with too few genes or genes
with too few reads can also be removed.

``` r
## remove pixels with too few genes
counts <- cleanCounts(cd)
```

    ## Converting to sparse matrix ...

    ## Filtering matrix with 262 cells and 15928 genes ...

    ## Resulting matrix has 260 cells and 14828 genes

![](getting_started_files/figure-markdown_github/getting_started_feature-1.png)

``` r
## feature select for genes
corpus <- restrictCorpus(counts, t1=0.05, t2 = 1)
```

    ## [1] "Restricting to overdispersed genes with alpha=0.05..."
    ## [1] "Calculating variance fit ..."
    ## [1] "Using gam with k=5..."
    ## [1] "345 overdispersed genes ... "
    ## [1] "Removing 81 genes present in less than 5% of datasets..."
    ## [1] "Removing 16 genes present in more than 100% of datasets..."
    ## [1] "248 genes remaining..."

`STdeconvolve` then applies Latent Dirichlet Allocation (LDA), a
generative statistical model commonly used in natural language
processing, to discover `K` latent cell-types. `STdeconvolve` fits a
range of LDA models to inform the choice of an optimal `K`.

``` r
lda <- fitLDA(corpus, Ks = seq(2, 9, by = 1), plot=TRUE, verbose=FALSE)
```

![](getting_started_files/figure-markdown_github/getting_started_opt-1.png)

In this example, we will use the model with the lowest model perplexity.
The resulting `theta` matrix can be interpreted as the proportion of
each deconvolved cell-type across each spatially resolved pixel. The
resulting `beta` matrix can be interpreted as the putative gene
expression profile for each deconvolved cell-type normalized to a
library size of 1. This `beta` matrix can be scaled by a depth factor
(ex. 1000) for interpretability.

``` r
results <- getBetaTheta(lda$models[[lda$kOpt2]])
deconProp <- results$theta
deconGexp <- results$beta*1000
```

We can now visualize the proportion of each deconvolved cell-type across
the original spatially resolved pixels.

``` r
vizAllTopics(deconProp, pos, 
             groups = annot, 
             group_cols = rainbow(length(levels(annot))),
             r=0.4)
```

![](getting_started_files/figure-markdown_github/getting_started_proportions-1.png)

We can also visualize the top marker genes for each deconvolved
cell-type. We will use deconvolved cell-types `7` and `2` as examples
here. We will define the top marker genes here as genes highly expressed
in the deconvolved cell-type (count \> 5) that also have the top 4
highest log2(fold change) when comparing the deconvolved cell-type’s
expression profile to the average of all other deconvolved cell-types’
expression profiles.

``` r
celltype <- 7
## highly expressed in cell-type of interest
highgexp <- names(which(deconGexp[celltype,] > 5))
## high log2(fold-change) compared to other deconvolved cell-types
log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
markers <- names(log2fc)[1:4]
## visualize
df <- merge(as.data.frame(pos), 
            as.data.frame(t(as.matrix(counts[markers,]))), 
            by = 0)
ps <- lapply(markers, function(marker) {
  vizGeneCounts(df = df,
              gene = marker,
              size = 2, stroke = 0.1,
              plotTitle = marker,
              winsorize = 0.05)
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2),
                        c(3, 4))
)
```

![](getting_started_files/figure-markdown_github/getting_started_expression-1.png)

``` r
celltype <- 2
## highly expressed in cell-type of interest
highgexp <- names(which(deconGexp[celltype,] > 5))
## high log2(fold-change) compared to other deconvolved cell-types
log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
markers <- names(log2fc)[1:4]
## visualize
df <- merge(as.data.frame(pos), 
            as.data.frame(t(as.matrix(counts[markers,]))), 
            by = 0)
ps <- lapply(markers, function(marker) {
  vizGeneCounts(df = df,
              gene = marker,
              size = 2, stroke = 0.1,
              plotTitle = marker,
              winsorize = 0.05)
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2),
                        c(3, 4))
)
```

![](getting_started_files/figure-markdown_github/getting_started_expression-2.png)
