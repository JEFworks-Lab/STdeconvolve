<img src="https://github.com/JEFworks/STdeconvolve/blob/package/docs/img/STdeconvolve_logo.png?raw=true"/>

<!-- badges: start -->
[![R build status](https://github.com/JEFworks/STdeconvolve/workflows/R-CMD-check/badge.svg)](https://github.com/JEFworks/STdeconvolve/actions)
<!-- badges: end -->

`STdeconvolve` enables reference-free cell-type deconvolution of multi-cellular pixel-resolution spatial transcriptomics data. The overall approach is detailed on [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.06.15.448381v2)

## Overview

`STdeconvolve` is an unsupervised machine learning approach to deconvolve multi-cellular pixel-resolution spatial transcriptomics datasets in order to recover the putative transcriptomic profiles of cell-types and their proportional representation within spatially resolved pixels without reliance on external single-cell transcriptomics references.

<img src="https://github.com/JEFworks/STdeconvolve/blob/package/docs/img/STdeconvolve_workflowforwebsite_v2.png?raw=true"/>

## Installation

To install `STdeconvolve`, we recommend using `remotes`:

``` r
require(remotes)
remotes::install_github('JEFworks-Lab/STdeconvolve')
```

Installation should take a few minutes on a typical desktop computer.

## Example

``` r
library(STdeconvolve)
## load built in data
data(mOB)
pos <- mOB$pos
cd <- mOB$counts
annot <- mOB$annot
## remove pixels with too few genes
counts <- cleanCounts(cd, min.lib.size = 100)
## feature select for genes
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
## choose optimal number of cell-types
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2, 9, by = 1))
## get best model results
optLDA <- optimalModel(models = ldas, opt = "min")
## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta
## visualize deconvolved cell-type proportions
vizAllTopics(deconProp, pos,
             groups = annot, 
             group_cols = rainbow(length(levels(annot))),
             r=0.4)	  
```

<img src="https://github.com/JEFworks/STdeconvolve/blob/package/docs/getting_started_files/figure-markdown_github/getting_started_proportions-1.png?raw=true"/>

More details can be found in the tutorials.

## Tutorials
- [Getting started with `STdeconvolve`](https://github.com/JEFworks/STdeconvolve/blob/package/docs/getting_started.md)
- [Additional features with `STdeconvolve`](https://github.com/JEFworks/STdeconvolve/blob/package/docs/additional_features.md)
- [Annotating deconvolved cell-types](https://github.com/JEFworks/STdeconvolve/blob/package/docs/celltype_annotation.md)
- [Analysis of 10X Visium data](https://github.com/JEFworks/STdeconvolve/blob/package/docs/visium_10x.md)

## Preprocessing datasets

For commands to reproduce the preprocessing of certain datasets used in the manuscript, check out:

[https://jef.works/STdeconvolve/](https://jef.works/STdeconvolve/)