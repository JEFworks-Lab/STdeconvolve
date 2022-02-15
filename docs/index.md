# STdeconvolve

<!-- badges: start -->
[![R build status](https://github.com/JEFworks/STdeconvolve/workflows/R-CMD-check/badge.svg)](https://github.com/JEFworks/STdeconvolve/actions)
<!-- badges: end -->

`STdeconvolve` enables reference-free cell-type deconvolution of multi-cellular pixel-resolution spatially resolved transcriptomics data

<img src="https://github.com/JEFworks/STdeconvolve/blob/devel/docs/img/STdeconvolve_logo.png?raw=true"/>

The overall approach is detailed on [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.06.15.448381v2)

## Overview

`STdeconvolve` is an unsupervised machine learning approach to deconvolve multi-cellular pixel-resolution spatial transcriptomics datasets in order to recover the putative transcriptomic profiles of cell-types and their proportional representation within spatially resolved pixels without reliance on external single-cell transcriptomics references.

<img src="https://github.com/JEFworks/STdeconvolve/blob/devel/docs/img/STdeconvolve_workflowforwebsite_v2.png?raw=true"/>

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

<img src="https://github.com/JEFworks/STdeconvolve/blob/devel/docs/getting_started_files/figure-markdown_github/getting_started_proportions-1.png?raw=true"/>

More details can be found in the tutorials.

## Tutorials
- [Getting started with STdeconvolve](getting_started.md)
- [Additional features with STdeconvolve](additional_features.md)
- [Annotating deconvolved cell-types](celltype_annotation.md)
- [Analysis of 10X Visium data](visium_10x.md)

## Installation

To install `STdeconvolve`, we recommend using `remotes`:

``` r
require(remotes)
remotes::install_github('JEFworks-Lab/STdeconvolve')
```

## Contributing

We welcome any bug reports, enhancement requests, general questions, and other contributions. To submit a bug report or enhancement request, please use the `STdeconvolve` GitHub issues tracker. For more substantial contributions, please fork this repo, push your changes to your fork, and submit a pull request with a good commit message.

## Reproducing Analyses

Links below point to code outlining the preprocessing of datasets used for analyses in the manuscript.

- [MERFISH mouse medial preoptic area (Moffit et al. 2018)](process_mpoa_data.md)
- [Mouse Olfactory Bulb (Stahl et al. 2016)](process_mob_data.md)
- [Breast cancer sections (Yoosuf et al. 2020)](process_bcl_data.md)
- [DBiT-seq of E11 mouse lower embryo and tail (Liu et al. 2020)](process_dbitseq_data.md)
