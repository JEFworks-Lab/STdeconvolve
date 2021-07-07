In this tutorial, we will walk through some of the main functionalities
of `STdeconvolve`.

``` r
library(STdeconvolve)
```

Given a counts matrix from pixel-resolution spatial transcriptomics data
where each spatially resolved measurement may represent mixtures from
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
over dispersed genes across ST pixels. Pixels with too few genes or
genes with too few reads can also be removed.

``` r
## remove pixels with too few genes
counts <- cleanCounts(counts = cd,
                      min.lib.size = 100,
                      min.reads = 1,
                      min.detected = 1)
```

![](getting_started_files/figure-markdown_github/getting_started_feature-1.png)

``` r
## feature select for genes
corpus <- restrictCorpus(counts,
                         removeAbove=1.0,
                         removeBelow = 0.05,
                         alpha = 0.05,
                         plot = TRUE,
                         verbose = TRUE)
```

    ## [1] "Removing 124 genes present in 100% or more of pixels..."
    ## [1] "14704 genes remaining..."
    ## [1] "Removing 3009 genes present in 5% or less of pixels..."
    ## [1] "11695 genes remaining..."
    ## [1] "Restricting to overdispersed genes with alpha = 0.05..."
    ## [1] "Calculating variance fit ..."
    ## [1] "Using gam with k=5..."
    ## [1] "232 overdispersed genes ... "

![](getting_started_files/figure-markdown_github/getting_started_feature-2.png)

    ##  Using top 1000 overdispersed genes. 
    ##  number of top overdispersed genes available: 232

![](getting_started_files/figure-markdown_github/getting_started_feature-3.png)

`STdeconvolve` then applies Latent Dirichlet allocation (LDA), a
generative statistical model commonly used in natural language
processing, to discover `K` latent cell-types. `STdeconvolve` fits a
range of LDA models to inform the choice of an optimal `K`.

``` r
## Note: the input corpus needs to be an integer count matrix of pixels x genes
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2, 9, by = 1), plot=TRUE, verbose=TRUE)
```

    ## [1] "Time to fit LDA models was 0.32mins"
    ## [1] "Computing perplexity for each fitted model..."

![](getting_started_files/figure-markdown_github/getting_started_opt-1.png)

In this example, we will use the model with the lowest model perplexity.

The shaded region indicates where a fitted model for a given K had an
`alpha` \> 1. `alpha` is an LDA parameter that is solved for during
model fitting and corresponds to the shape parameter of a symmetric
Dirichlet distribution. In the model, this Dirichlet distribution
describes the cell-type proportions in the pixels. A symmetric Dirichlet
with alpha \> 1 would lead to more uniform cell-type distributions in
the pixels and difficulty identifying distinct cell-types. Instead, we
want models with alphas \< 1, resulting in sparse distributions where
only a few cell-types are represented in a given pixel.

The resulting `theta` matrix can be interpreted as the proportion of
each deconvolved cell-type across each spatially resolved pixel. The
resulting `beta` matrix can be interpreted as the putative gene
expression profile for each deconvolved cell-type normalized to a
library size of 1. This `beta` matrix can be scaled by a depth factor
(ex. 1000) for interpretability.

``` r
## select model with minimum perplexity
optLDA <- optimalModel(models = ldas, opt = "min")
## extract pixel cell-type proportions (theta) and cell-type gene expression profiles (beta) for the given dataset
results <- getBetaTheta(optLDA, corpus = t(as.matrix(corpus)))
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
cell-type. We will use deconvolved cell-types `5` and `1` as examples
here. We will define the top marker genes here as genes highly expressed
in the deconvolved cell-type (count \> 5) that also have the top 4
highest log2(fold change) when comparing the deconvolved cell-type’s
expression profile to the average of all other deconvolved cell-types’
expression profiles.

``` r
celltype <- 5
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
              winsorize = 0.05,
              showLegend = FALSE)
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2),
                        c(3, 4))
)
```

![](getting_started_files/figure-markdown_github/getting_started_expression-1.png)

``` r
celltype <- 1
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
              winsorize = 0.05,
              showLegend = FALSE)
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2),
                        c(3, 4))
)
```

![](getting_started_files/figure-markdown_github/getting_started_expression-2.png)
