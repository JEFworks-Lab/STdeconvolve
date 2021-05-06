<img src="https://github.com/JEFworks/STdeconvolve/blob/package/docs/img/STdeconvolve_logo.png?raw=true"/>

# STdeconvolve

<!-- badges: start -->
[![R build status](https://github.com/JEFworks/STdeconvolve/workflows/R-CMD-check/badge.svg)](https://github.com/JEFworks/STdeconvolve/actions)
<!-- badges: end -->

`STdeconvolve` enables reference-free cell-type deconvolution of multi-cellular spatial transcriptomics data. The overall approach is detailed in the following publications: *COMING SOON*

## Overview

`STdeconvolve` is an unsupervised machine learning approach to deconvolve multi-cellular pixel resolution spatial transcriptomics datasets in order to recover the putative transcriptomic profiles of cell-types and their proportional representation within spatially resolved pixels without reliance on external single-cell transcriptomics references.

<img src="https://github.com/JEFworks/STdeconvolve/blob/package/docs/img/STdeconvolve_workflowforwebsite.png?raw=true"/>

## Installation

To install `STdeconvolve`, we recommend using `remotes`:

``` r
require(remotes)
remotes::install_github('JEFworks/STdeconvolve')
```

