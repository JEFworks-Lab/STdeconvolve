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
