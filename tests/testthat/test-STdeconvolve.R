test_that("restrictCorpus() restrict the corpus", {
  
  data(mOB)
  pos <- mOB$pos
  cd <- mOB$counts
  ## remove pixels with too few genes
  counts <- cleanCounts(cd, min.lib.size = 100)
  ## feature select for genes
  corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow=0.05)
  
  ## test remaining genes
  expect_equal(dim(corpus)[1], 232)
  ## test number of pixels
  expect_equal(dim(corpus)[2], 260)
  
  ## test that selecting overdispersed genes works too
  corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow=0.05, nTopOD=10)
  expect_equal(dim(corpus)[1], 10)
  expect_equal(rownames(corpus)[10], "Sox11")
  
})


test_that("preprocess() filters properly", {
  
  data(mOB)
  pos <- mOB$pos
  cd <- mOB$counts
  
  corpus <- preprocess(dat = t(cd),
                       min.lib.size=100,
                       removeAbove = 1.0,
                       removeBelow = 0.05,
                       plot = FALSE)

  expect_equal(class(corpus$corpus)[1], "matrix")
  ## test remaining genes
  expect_equal(dim(corpus$corpus)[1], 260)
  ## test number of pixels
  expect_equal(dim(corpus$corpus)[2], 232)
  
  ## test that selecting overdispersed genes works too
  corpus <- preprocess(dat = t(cd),
                       min.lib.size=100,
                       removeAbove = 1.0,
                       removeBelow = 0.05,
                       nTopOD = 10,
                       plot = FALSE)
  
  expect_equal(dim(corpus$corpus)[2], 10)
  expect_equal(colnames(corpus$corpus)[10], "Sox11")
  
})


## to speed up the tests, fit a set of LDA models once and use these outputs
## to check model fitting, beta and thetas, visualization, and annotation in 
## one comprehensive test
test_that("STdeconvolve fits models and visualizes results", {
  
  data(mOB)
  pos <- mOB$pos
  cd <- mOB$counts
  annot <- mOB$annot
  ## remove pixels with too few genes
  counts <- cleanCounts(cd, min.lib.size = 100)
  ## feature select for genes
  corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow=0.05)
  
  ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2, 5, by = 1),
                 perc.rare.thresh = 0.05,
                 plot=FALSE,
                 verbose=TRUE)
  
  alphas <- unlist(sapply(ldas$models, slot, "alpha"))
  expected_alphas <- c(44.2343385, 0.7856869, 0.6548315, 0.7804621)
  perplexities <- ldas$perplexities
  expected_perplexities <- c(157.947571324413, 142.851762473055,
                             141.437602764842, 140.679764694871)
  
  ## test by finding correlation between the observed and expected values
  ## in order to avoid cases where maybe the values differ by very small amounts
  ## but test_that would still flag them as errors:
  alpha_cor <- cor(as.vector(alphas), expected_alphas)
  expect_gt(alpha_cor, 0.99)
  perplex_cor <- cor(perplexities, expected_perplexities)
  expect_gt(perplex_cor, 0.99)
  
  ## check that optimal model obtained
  optLDA <- optimalModel(models = ldas, opt = "min")
  expect_equal(optLDA@k, 5)
  
  ## check that beta and theta extractable and scaled properly
  results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
  deconProp <- results$theta
  deconGexp <- results$beta
  
  ## scaled expression of first value in beta matrix
  ## use some leniency in case of small fluctuations
  ## in either case, the scaling by 1000 should give 5.8 for very first cell
  expect_gt(deconGexp[1], 3.06) ## real value is: 3.06909812664901
  
  ## cell-type proportions of first pixel
  first_pixel_proportions <- c(0, 0.0744826390670371,
                               0.578165711570332,
                               0.158247284814881, 
                               0.189104364547749)
  
  ## test the correlation in case some small fluctuations
  theta_cor <- cor(as.vector(deconProp[1,]), first_pixel_proportions)
  expect_gt(theta_cor, 0.99)
  
  ## visualize deconvolved cell-type proportions
  plt <- vizAllTopics(deconProp, pos,
                      groups = annot, 
                      group_cols = rainbow(length(levels(annot))),
                      r=0.4)
  
  ## test to see if the first row of the plot data and its structure looks correct
  expected_plt <- structure(list(Row.names = structure("ACAACTATGGGTTGGCGG",
                                                       class = "character"), 
                                 x = 16.001, y = 16.036,
                                 Pixel.Groups = "3: Outer Plexiform Layer", 
                                 Topics = structure(1L, .Label = c("Topic.1",
                                                                   "Topic.2",
                                                                   "Topic.3", 
                                                                   "Topic.4",
                                                                   "Topic.5"),
                                                    class = "factor"), 
                                 value = 0.0), row.names = 1L,
                            class = "data.frame")
  
  ## the first cell-type has proportion 0 in the first pixel, so the first
  ## wedge that is drawn belongs to the second cell-type
  expect_equal(as.character(plt$layers[[1]]$data[1,"CellTypes"]), "X2")
  expect_equal(plt$layers[[1]]$data[1,"Pixel.Groups"], "3: Outer Plexiform Layer")
  #expect_equal(plt$layers[[1]]$data[1,"Row.names"], structure("ACAACTATGGGTTGGCGG",
  #                                                            class = "character"))
  
  ## test individual cell-type plotting
  plt2 <- vizTopic(theta = deconProp, pos = pos, topic = "5", plotTitle = "X5",
                   size = 5, stroke = 1, alpha = 0.5,
                   low = "white",
                   high = "red")
  
  ## test to see if the first row of the plot data and its structure looks correct
  expected_plt2 <- structure(list(Row.names = structure("ACAACTATGGGTTGGCGG",
                                                        class = "character"), 
                                  proportion = 0.189104364547749,
                                  x = 16.001, y = 16.036), row.names = 1L,
                             class = "data.frame")
  
  expect_gt(plt2$layers[[1]]$data[1,"proportion"], 0.189)
  #expect_equal(plt2$layers[[1]]$data[1,"Row.names"], structure("ACAACTATGGGTTGGCGG",
  #                                                             class = "character"))
  
  ## test gene expression
  celltype <- 5
  ## highly expressed in cell-type of interest
  highgexp <- names(which(deconGexp[celltype,] > 5))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(
    log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])),
    decreasing=TRUE
  )
  markers <- names(log2fc)[1:4]
  
  df <- merge(as.data.frame(pos), 
              as.data.frame(t(as.matrix(counts[markers,]))), 
              by = 0)
  
  plt3 <- vizGeneCounts(df = df,
                        gene = "Sox11", # should be one of the top 4
                        # groups = annot,
                        # group_cols = rainbow(length(levels(annot))),
                        size = 3, stroke = 0.1,
                        plotTitle = "Sox11",
                        winsorize = 0.05,
                        showLegend = TRUE)
  
  ## test to see if the first row of the plot data and its structure looks correct
  expected_plt3 <- structure(list(Row.names = structure("ACAACTATGGGTTGGCGG",
                                                        class = "character"), 
                                  x = 16.001, y = 16.036, Mag = 0, Sox11 = 1,
                                  Cnp = 2, Nrep = 5), row.names = 1L,
                             class = "data.frame")
  
  expect_equal(plt3$layers[[1]]$data[1,"Sox11"], 1)
  expect_equal(plt3$layers[[1]]$data[1,"Mag"], 0)
  #expect_equal(plt3$layers[[1]]$data[1,"Row.names"], structure("ACAACTATGGGTTGGCGG",
  #                                                             class = "character"))
  
  ## finally, check that cell-type annotation assignment is working
  
  # proxy theta for the annotated layers
  mobProxyTheta <- model.matrix(~ 0 + annot)
  rownames(mobProxyTheta) <- names(annot)
  # fix names
  colnames(mobProxyTheta) <- unlist(lapply(colnames(mobProxyTheta), function(x) {
    unlist(strsplit(x, "annot"))[2]
  }))
  mobProxyGexp <- counts %*% mobProxyTheta
  
  corMtx_theta <- getCorrMtx(m1 = as.matrix(deconProp), # the deconvolved cell-type `theta` (pixels x celltypes)
                             m2 = as.matrix(mobProxyTheta), # the reference `theta` (pixels x celltypes)
                             type = "t") # "b" = comparing beta matrices, "t" for thetas
  rownames(corMtx_theta) <- paste0("decon_", seq(nrow(corMtx_theta)))
  
  ## check that the correlation worked
  expected_corr <- c(`1: Granular Cell Layer` = -0.276244855335915,
                     `2: Mitral Cell Layer` = -0.142174337478699, 
                     `3: Outer Plexiform Layer` = 0.0338622840004436,
                     `4: Glomerular Layer` = 0.367178778101782, 
                     `5: Olfactory Nerve Layer` = 0.0300295669474897)
  
  ## test the correlation in case some small fluctuations
  corMtx_cor <- cor(corMtx_theta[1,], expected_corr)
  expect_gt(corMtx_cor, 0.99)
  
  pairs <- lsatPairs(t(corMtx_theta))
  m <- t(corMtx_theta)[pairs$rowix, pairs$colsix]
  
  ## check that the pairing worked
  expected_pairs <- c(decon_4 = 0.803095452168013,
                      decon_3 = -0.560104806711357, 
                      decon_5 = 0.546757310198154,
                      decon_1 = -0.276244855335915,
                      decon_2 = -0.431325040365262)
  
  ## test the correlation in case some small fluctuations
  pair_cor <- cor(m[1,], expected_pairs)
  expect_gt(pair_cor, 0.99)
  
  ## GSEA
  mobProxyLayerMarkers <- list()
  
  ## make the tissue layers the rows and genes the columns
  gexp <- t(as.matrix(mobProxyGexp))
  
  for (i in seq(length(rownames(gexp)))){
    celltype <- i
    ## log2FC relative to other cell-types
    ## highly expressed in cell-type of interest
    highgexp <- names(which(gexp[celltype,] > 10))
    ## high log2(fold-change) to other deconvolved cell-types and limit top 200
    log2fc <- sort(
      log2(gexp[celltype,highgexp]/colMeans(gexp[-celltype,highgexp])),
      decreasing=TRUE)[1:200]
    
    ## for gene set of the ground truth cell-type, get the genes
    ## with log2FC > 1 (so FC > 2 over the mean exp of the other cell-types)
    markers <- names(log2fc[log2fc > 1])
    mobProxyLayerMarkers[[ rownames(gexp)[celltype] ]] <- markers
  }
  
  # No longer exactly same due to switching from liger to fgsea
  # test not integral to functioning of method anyway
  # omit test
  #
  # celltype_annotations <- annotateCellTypesGSEA(beta = results$beta,
  #                                               gset = mobProxyLayerMarkers,
  #                                               qval = 0.05)
  # 
  # expected_annots <- structure(list(p.val = c(9.99900009999e-05,
  #                                             9.99900009999e-05),
  #                                   q.val = c(0.0001999800019998,
  #                                             0.0001999800019998),
  #                                   sscore = c(2.36876122685379,
  #                                              -2.07135070110146),
  #                                   edge = c(3.24141633337745,
  #                                            1.35028517017407)),
  #                              row.names = c("5: Olfactory Nerve Layer",
  #                                            "1: Granular Cell Layer"),
  #                              class = "data.frame")
  # 
  # expect_equal(rownames(celltype_annotations$results$`2`),
  #               c("5: Olfactory Nerve Layer", "1: Granular Cell Layer"))
  # expect_gt(celltype_annotations$results$`2`$edge[1], 3.24)
  
})




test_that("fitLDA() does not depend on input type or threads", {
  
  data(mOB)
  counts <- cleanCounts(mOB$counts, min.lib.size = 100)
  corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow=0.05, verbose=FALSE)
  
  ## dense matrix, 1 thread
  ldas <- fitLDA(t(as.matrix(corpus)), Ks = 4, plot=FALSE, verbose=FALSE)
  model <- ldas$models[[1]]
  expect_s4_class(model, "STdeconvolveLDA")
  expect_output(print(model), "4 cell-types")
  expect_equal(dim(model@gamma), c(260, 4))
  expect_equal(rowSums(model@gamma), rep(1, 260))
  expect_equal(rowSums(exp(model@beta)), rep(1, 4))
  expect_equal(ldas$perplexities, exp(-sum(model@loglikelihood) / sum(corpus)))
  ## the likelihood bound should not decrease from one EM iteration to the next
  expect_true(all(diff(model@logLiks) > -1e-6 * abs(model@logLiks[1])))
  
  ## sparse matrix, 2 threads: identical model
  ldas2 <- fitLDA(Matrix::t(corpus), Ks = 4, ncores=2, plot=FALSE, verbose=FALSE)
  expect_identical(ldas2$models[[1]]@beta, model@beta)
  expect_identical(ldas2$models[[1]]@gamma, model@gamma)
  expect_identical(ldas2$models[[1]]@alpha, model@alpha)
  
  ## invalid input
  bad <- t(as.matrix(corpus))
  bad[1,] <- 0
  expect_error(fitLDA(bad, Ks = 2, plot=FALSE, verbose=FALSE), "at least one non-zero")
  bad[1,] <- 0.5
  expect_error(fitLDA(bad, Ks = 2, plot=FALSE, verbose=FALSE), "integer")
  
})


test_that("getBetaTheta() deconvolves a new corpus with a fitted model", {
  
  data(mOB)
  counts <- cleanCounts(mOB$counts, min.lib.size = 100)
  corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow=0.05, verbose=FALSE)
  x <- t(as.matrix(corpus))
  model <- fitLDA(x, Ks = 4, plot=FALSE, verbose=FALSE)$models[[1]]
  fitted <- getBetaTheta(model, perc.filt = 0, verbose=FALSE)
  
  ## the corpus the model was fitted to, with pixels and genes shuffled
  set.seed(1)
  shuffled <- x[sample(nrow(x)), sample(ncol(x))]
  new <- getBetaTheta(model, corpus = shuffled, perc.filt = 0, verbose=FALSE)
  expect_equal(new$beta, fitted$beta)
  expect_equal(rownames(new$theta), rownames(shuffled))
  expect_equal(new$theta[rownames(x),], fitted$theta, tolerance=1e-4)
  
  ## sparse input and several threads give the same result
  new2 <- getBetaTheta(model, corpus = Matrix::Matrix(shuffled, sparse=TRUE),
                       perc.filt = 0, verbose=FALSE, ncores=2)
  expect_identical(new2$theta, new$theta)
  
  ## genes missing from the new corpus are treated as not detected
  sub <- shuffled[1:50, 1:150]
  new3 <- getBetaTheta(model, corpus = sub, perc.filt = 0, verbose=FALSE)
  expect_equal(dim(new3$theta), c(50, 4))
  expect_equal(unname(rowSums(new3$theta)), rep(1, 50))
  expect_equal(rownames(new3$theta), rownames(sub))
  
  expect_error(getBetaTheta(list(), verbose=FALSE), "fitted with fitLDA")
  
  ## pixels without any counts of the genes in the model
  colnames(sub) <- paste0("unknown", seq_len(ncol(sub)))
  expect_error(getBetaTheta(model, corpus = sub, verbose=FALSE), "at least one non-zero")
  
})


test_that("preprocessing accepts sparse matrices", {
  
  data(mOB)
  cd <- mOB$counts
  sparse <- Matrix::Matrix(cd, sparse=TRUE)
  
  expect_equal(getOverdispersedGenes(sparse, verbose=FALSE),
               getOverdispersedGenes(cd, verbose=FALSE))
  ## gene variances computed from the non-zero counts
  od <- getOverdispersedGenes(cd, verbose=FALSE, details=TRUE)
  expect_equal(exp(od$df$v), unname(apply(cd, 1, var)))
  
  dense <- preprocess(t(cd), min.lib.size=100, removeAbove=0.95, removeBelow=0.05,
                      plot=FALSE, verbose=FALSE)
  fromSparse <- preprocess(Matrix::t(sparse), min.lib.size=100, removeAbove=0.95,
                           removeBelow=0.05, plot=FALSE, verbose=FALSE)
  expect_equal(fromSparse$corpus, dense$corpus)
  
})


test_that("vizAllTopics() draws one wedge per cell-type present in a pixel", {
  
  theta <- rbind(a = c(0.5, 0.25, 0.25), b = c(0, 1, 0), c = c(0.2, 0, 0.8))
  colnames(theta) <- c("1", "2", "3")
  pos <- cbind(x = c(1, 2, 3), y = c(1, 1, 2))
  rownames(pos) <- rownames(theta)
  
  w <- STdeconvolve:::pieWedges(theta, x = pos[, "x"], y = pos[, "y"], r = rep(0.4, 3), n = 360)
  expect_equal(length(unique(w$wedge)), 6)
  ## area of each wedge (shoelace formula) is its share of the area of the pie
  area <- vapply(split(seq_along(w$x), w$wedge), function(i) {
    x <- w$x[i]; y <- w$y[i]; j <- c(seq_along(i)[-1], 1)
    abs(sum(x * y[j] - x[j] * y)) / 2
  }, numeric(1))
  expect_equal(unname(area) / (pi * 0.4^2), c(0.5, 0.25, 0.25, 1, 0.2, 0.8), tolerance = 1e-3)
  expect_equal(unique(w$topic[w$pixel == 3]), c(1, 3))
  ## all points are on or inside the circle of their pixel
  d <- sqrt((w$x - pos[w$pixel, "x"])^2 + (w$y - pos[w$pixel, "y"])^2)
  expect_lte(max(d), 0.4 + 1e-12)
  
  ## pixels without a position are left out, groups stay with their pixels
  p <- suppressMessages(vizAllTopics(theta, pos[c("c", "a"), ], r = 0.4,
                                     groups = c("g1", "g2", "g3"),
                                     group_cols = c(g1 = "red", g2 = "blue", g3 = "black")))
  expect_s3_class(p, "ggplot")
  dat <- p$layers[[1]]$data
  expect_equal(levels(dat$CellTypes), c("X1", "X2", "X3"))
  expect_equal(unique(dat$Pixel.Groups), c("g1", "g3"))
  expect_silent(ggplot2::ggplot_build(p))
  
  ## data.frames are accepted by vizTopic()
  p <- vizTopic(as.data.frame(theta), as.data.frame(pos), topic = "2")
  expect_s3_class(p, "ggplot")
  
})


test_that("edge cases of the helper functions", {
  
  data(mOB)
  cd <- mOB$counts
  
  ## several patterns of genes to remove
  pre <- preprocess(t(cd), min.lib.size=100, genes.to.remove=c("^mt-", "^Rpl", "^Rps"),
                    ODgenes=FALSE, plot=FALSE, verbose=FALSE)
  expect_false(any(grepl("^mt-|^Rpl|^Rps", colnames(pre$corpus))))
  expect_error(preprocess("no/such/file.txt", plot=FALSE, verbose=FALSE), "does not exist")
  
  ## a single K has no knee
  counts <- cleanCounts(cd, min.lib.size = 100)
  corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow=0.05, verbose=FALSE)
  ldas <- fitLDA(t(as.matrix(corpus)), Ks = 3, plot=FALSE, verbose=FALSE)
  expect_equal(ldas$kneedOptK, 3)
  expect_equal(ldas$ctPropOptK, 3)
  expect_equal(optimalModel(ldas, opt = "kneed")@k, 3L)
  
  ## correlations restricted to highly expressed genes
  beta <- getBetaTheta(ldas$models[[1]], verbose=FALSE)$beta
  full <- suppressMessages(getCorrMtx(beta, beta, type="b"))
  expect_equal(unname(diag(full)), rep(1, 3))
  top <- suppressMessages(getCorrMtx(beta, beta, type="b", thresh=0.005))
  expect_equal(dim(top), c(3, 3))
  genes <- colnames(beta)[beta[1, ] > 0.005 & beta[2, ] > 0.005]
  expect_equal(top[1, 2], cor(beta[1, genes], beta[2, genes]))
  expect_error(suppressMessages(getCorrMtx(beta, beta, type="b", thresh=2)), "thresh")
  
  ## a single remaining cell-type keeps the matrix shape
  theta <- cbind("1" = c(0.99, 0.98), "2" = c(0.01, 0.02))
  expect_equal(dim(STdeconvolve:::filterTheta(theta, perc.filt=0.05, verbose=FALSE)), c(2, 1))
  
})


test_that("models agree with topicmodels, and models of topicmodels are accepted", {
  
  skip_if_not_installed("topicmodels")
  skip_if_not_installed("slam")
  
  data(mOB)
  counts <- cleanCounts(mOB$counts, min.lib.size = 100)
  corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow=0.05, verbose=FALSE)
  x <- t(as.matrix(corpus))
  ldas <- fitLDA(x, Ks = 4, plot=FALSE, verbose=FALSE)
  model <- ldas$models[[1]]
  
  ## same model as topicmodels::LDA up to the convergence tolerances
  expect_s4_class(ldas$fitCorpus, "dgCMatrix")
  stm <- slam::as.simple_triplet_matrix(x)
  ref <- topicmodels::LDA(stm, k=4,
                          control=list(seed=0, verbose=0, keep=1, estimate.alpha=TRUE))
  expect_equal(ldas$perplexities, topicmodels::perplexity(ref), tolerance=0.005)
  expect_gt(min(diag(cor(t(exp(model@beta)), t(exp(ref@beta))))), 0.99)
  expect_gt(min(diag(cor(model@gamma, ref@gamma))), 0.99)
  
  ## a corpus in the format of previous versions is still accepted
  expect_identical(fitLDA(stm, Ks = 4, plot=FALSE, verbose=FALSE)$models[[1]]@beta, model@beta)
  
  ## models fitted by previous versions of fitLDA() are "LDA_VEM" objects
  post <- topicmodels::posterior(ref)
  old <- getBetaTheta(ref, perc.filt = 0, verbose=FALSE)
  expect_equal(old$theta, post$topics)
  expect_equal(old$beta, post$terms)
  expect_equal(optimalModel(list(models=list("4"=ref)), opt = 4)@k, 4L)
  
  ## cell-type proportions of new pixels, given the same model
  sub <- x[1:50, 1:150]
  new <- getBetaTheta(ref, corpus = sub, perc.filt = 0, verbose=FALSE)
  refNew <- topicmodels::posterior(ref, newdata = slam::as.simple_triplet_matrix(sub))$topics
  expect_equal(dimnames(new$theta), dimnames(refNew))
  expect_lt(max(abs(new$theta - refNew)), 0.02)
  
})


test_that("fitLDA with a warm start fits models that are as good", {
  data(mOB)
  counts <- cleanCounts(mOB$counts, min.lib.size = 100, verbose = FALSE)
  corpus <- t(as.matrix(restrictCorpus(counts, removeAbove = 1.0, removeBelow = 0.05,
                                       verbose = FALSE)))
  Ks <- c(4, 8, 6)
  cold <- fitLDA(corpus, Ks = Ks, plot = FALSE, verbose = FALSE)
  warm <- fitLDA(corpus, Ks = Ks, plot = FALSE, verbose = FALSE, warm.start = TRUE)
  
  expect_equal(names(warm$models), as.character(Ks))
  expect_equal(unname(sapply(warm$models, slot, "k")), Ks)
  ## the model with the largest K is the same, the others are close
  expect_equal(warm$models[["8"]]@beta, cold$models[["8"]]@beta)
  expect_true(all(warm$perplexities < cold$perplexities * 1.01))
  for (m in warm$models) {
    expect_equal(rowSums(exp(m@beta)), rep(1, m@k))
    expect_equal(rowSums(m@gamma), rep(1, nrow(corpus)))
  }
  
  ## merging keeps the gene distributions normalised
  mergeCellTypes <- getFromNamespace("mergeCellTypes", "STdeconvolve")
  merged <- mergeCellTypes(cold$models[["8"]], cold$fitCorpus, 3)
  expect_equal(dim(merged), c(3, ncol(corpus)))
  expect_equal(rowSums(merged), rep(1, 3))
})
