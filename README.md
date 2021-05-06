<img src="https://github.com/JEFworks/STdeconvolve/blob/package/docs/img/STdeconvolve_logo.png?raw=true"/>

<!-- badges: start -->
[![R build status](https://github.com/JEFworks/STdeconvolve/workflows/R-CMD-check/badge.svg)](https://github.com/JEFworks/STdeconvolve/actions)
<!-- badges: end -->

`STdeconvolve` enables reference-free cell-type deconvolution of multi-cellular spatial transcriptomics data. The overall approach is detailed in the following publications: **COMING SOON**

## Overview

`STdeconvolve` is an unsupervised machine learning approach to deconvolve multi-cellular pixel resolution spatial transcriptomics datasets in order to recover the putative transcriptomic profiles of cell-types and their proportional representation within spatially resolved pixels without reliance on external single-cell transcriptomics references.

<img src="https://github.com/JEFworks/STdeconvolve/blob/package/docs/img/STdeconvolve_workflowforwebsite.png?raw=true"/>

## Installation

To install `STdeconvolve`, we recommend using `remotes`:

``` r
require(remotes)
remotes::install_github('JEFworks/STdeconvolve')
```

## Example

``` r
library(STdeconvolve)
## load built in data
data(mOB)
pos <- mOB$pos
cd <- mOB$counts
annot <- mOB$annot
## remove pixels with too few genes
counts <- cleanCounts(cd)
## feature select for genes
corpus <- restrictCorpus(counts, t1=0.05, t2 = 1)
## choose optimal number of cell-types
lda <- fitLDA(corpus, Ks = seq(2, 9, by = 1), plot=TRUE, verbose=FALSE)
## get best model results
results <- getBetaTheta(lda$models[[lda$kOpt2]])
deconProp <- results$theta
deconGexp <- results$beta*1000
## visualize deconvolved cell-type proportions
vizAllTopics(deconProp, pos, 
             groups = annot, 
             group_cols = rainbow(length(levels(annot))),
             r=0.4)	     
```

<img src="https://github.com/JEFworks/STdeconvolve/blob/package/docs/getting_started_files/figure-markdown_github/getting_started_proportions-1.png?raw=true"/>

## Tutorials
- [Getting Started with STdeconvolve](https://github.com/JEFworks/STdeconvolve/blob/package/docs/getting_started.md)