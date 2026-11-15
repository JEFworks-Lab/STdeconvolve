# Version 1.3.3 (2026-11-15)

+ `fitLDA()` fits the LDA models with a new C++ implementation of variational
  EM instead of `topicmodels::LDA()`. It fits the same model with the same
  initialization for a given `seed`, and returns objects of the same class
  ("LDA_VEM"), so results are equivalent (not bit-identical: the cell-type
  proportions of each pixel are optimized to a tighter tolerance), but it is
  30-50 times faster on one core.
+ `ncores` in `fitLDA()` is now the number of threads used to fit each model.
  The fitted models do not depend on the number of threads. `BiocParallel` is
  no longer used.
+ `fitLDA()` and `getBetaTheta()` accept sparse matrices, and
  `getBetaTheta()` with a new `corpus` uses the C++ implementation as well
  (new argument `ncores`).
+ `restrictCorpus()`, `preprocess()` and `getOverdispersedGenes()` keep the
  counts sparse (2-3 times faster, much less memory, same results), and
  `preprocess()` accepts sparse matrices.
+ `getBetaTheta()` and `perplexityPlot()` now stop with an error if a pixel in
  `corpus` has no counts of any gene in the model, e.g. when a genes x pixels
  matrix is given instead of pixels x genes. Previously this silently returned
  equal proportions of all cell-types for such pixels.
+ Fitting is faster when the number of cell-types is large: cell-types that
  are absent from a pixel are skipped while the proportions of that pixel are
  optimized.
+ `vizAllTopics()` draws the pie charts as polygons with `ggplot2` directly
  (about 10 times faster to render, same appearance). `scatterpie` is no
  longer needed. `groups` and `r` can be given for all pixels in `theta` when
  only some of them have a position.
+ Bug fixes:
    + `vizAllTopics()` with `overlay`: the limits of the image were dropped,
      so the pies were not aligned with the image.
    + `vizTopic()` did not accept a data.frame for `theta`.
    + `preprocess()` failed when `genes.to.remove` had more than one pattern,
      opened a graphics device with `plot = FALSE`, and continued after a
      path to a file that does not exist.
    + `fitLDA()` returned empty `kneedOptK` and `ctPropOptK` for fewer than 3
      values in `Ks`, so `optimalModel(opt = "kneed")` failed. The K with the
      lowest value is used in that case.
    + `annotateCellTypesGSEA()` did not order gene sets with equal adjusted
      p-values by their enrichment score.
    + `getCorrMtx()` with `thresh` printed its warnings instead of raising
      them, and `filterTheta()` (used by `getBetaTheta()`) returned a vector
      when a single cell-type remained.
    + `ggplot2` deprecation warnings in `fitLDA()`, `perplexityPlot()` and
      `correlationPlot()` (requires `ggplot2` >= 3.4.0).


# Version 1.3.2 (2024-11-07)

+ Remove `liger` as a dependency. Replace with `fgsea` Korotkevich G, Sukhov V, Sergushichev A (2019). “Fast gene set enrichment analysis.” bioRxiv. doi:10.1101/060012, http://biorxiv.org/content/early/2016/06/20/060012. 


# Version 1.3.1 (2022-12-10)

+ Based on `Pull requests` after being published, made some minor changes to help guide users:
+ Added new vignette: Discussion and examples of processing multiple datasets together or separately
+ Minor aesthetic updates to `VizAllTopics()`. Now lwd = 0.01 and "Topics" changes to "CellTypes"
+ Emphasis that functions `simulateBregmaSpots()` and `buildBregmaCorpus()` to generate MPOA multi-cellular spots are on the `devel` branch

# Version 0.99.12 (2022-04-22)

+ Accepted into Bioconductor 3.15 devel branch
+ Official paper out in Nat Comms 2022-04-29 !!

# Version 0.99.10 (2022-04-13)

+ R version for Bioconductor and `package` branch: (>= 4.1)
+ R version in `devel` branch in GitHub: (>= 3.6)

# Changes in version 0.99.0 (2022-02-15)

+ Submitted to Bioconductor
