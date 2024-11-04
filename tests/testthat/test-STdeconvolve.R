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
  
  expect_equal(plt$layers[[1]]$data[1,"value"], 0.0)
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
  
  celltype_annotations <- annotateCellTypesGSEA(beta = results$beta,
                                                gset = mobProxyLayerMarkers,
                                                qval = 0.05)
  
  expected_annots <- structure(list(p.val = c(9.99900009999e-05,
                                              9.99900009999e-05),
                                    q.val = c(0.0001999800019998,
                                              0.0001999800019998),
                                    sscore = c(2.36876122685379,
                                               -2.07135070110146),
                                    edge = c(3.24141633337745,
                                             1.35028517017407)),
                               row.names = c("5: Olfactory Nerve Layer",
                                             "1: Granular Cell Layer"),
                               class = "data.frame")
  
  expect_equal(rownames(celltype_annotations$results$`2`),
                c("5: Olfactory Nerve Layer", "1: Granular Cell Layer"))
  expect_gt(celltype_annotations$results$`2`$edge[1], 3.24)
  
})


