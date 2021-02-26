# Compute the cross correlation between top terms in each topic
# correlation based on the term counts across the corpus documents
#
# beta = topic-term distribution matrix
# corpus = document x term matrix. Counts of terms in each document
# thresh = the probability of the topic generating the term
#
# correlation = either "corpus" or "beta"
#        Choose which matrix to perform the correlation on.
#        Either wrt gene counts across the docs (corpus)
#        or the term probability for each topic (beta)
#
# returns:
# list of the resulting term cross correlation matrices for each topic
# also returns the correlation heat maps (optional)
#
topicTermCorrelationMats <- function(beta, corpus, thresh = 0.05, correlation = "corpus", plots = TRUE) {

  # get the top terms for each topic based on threshold
  top_topic_terms <- apply(beta, 1, function(x) {
    ordered_terms <- x[order(x, decreasing = TRUE)]
    top_terms <- names(which(ordered_terms > thresh))
    if (length(top_terms) < 2) {
      top_terms <- names(ordered_terms[1:2])
    }
    top_terms
  })

  termCorrelations <- lapply(1:length(top_topic_terms), function(topic) {
    topic_terms <- top_topic_terms[[topic]]
    if (correlation == "corpus") {
      term_corr <- cor(corpus[,topic_terms])
    } else if (correlation == "beta") {
      term_corr <- cor(beta[,topic_terms])
    } else {
      stop("correlation must equal either `corpus` or `beta`.")
    }


    if (plots) {
      par(mfrow=c(1,1), mar=c(8,8,2,2))
      h <- heatmap.2(term_corr,
                     density.info = "none",
                     trace = "none",
                     col = correlation_palette,
                     breaks = correlation_breaks,
                     cexRow=0.5,cexCol=0.5, margins=c(3,3),
                     main = topic,
                     lhei = c(1,5),
                     key.xlab = "Correlation",
                     key.title = NA)
    }

    term_corr
  })
  return(termCorrelations)
}


# Compute the "mean correlation" for each topic
# based on the topic term corss correlation matrices
#
# topicCorrList = list of the resulting term cross correlation matrices
#      for each topic. (topicTermCorrelationMats)
#
# returns:
# dataframe with the correlations for each topic (cor) and the number of
# topic terms used (numTerms)
#
topicTermCorrelation <- function(topicCorrList) {

  results <- lapply(topicCorrList, function(cor_mtx) {

    # sum the correlations, but subtract the auto-correlations (b/c all 1 and not interesting)
    # then divide by total correlations added
    num_terms <- nrow(cor_mtx)
    total_cor <- (sum(cor_mtx) - num_terms) / (num_terms**2 - num_terms)
    c(total_cor, num_terms)

  })

  df <- data.frame(matrix(unlist(results), nrow=length(results), byrow=T))
  colnames(df) <- c("cor", "numTerms")
  return(df)

}


# compute pvalue for topic correlation
# based on a null distribution of randomly sampled sets of terms
#
# null = ggplot transformed (long) dataframe
# ex:
# dat <- data.frame(draw = factor(rep(1:dim(termCorrNulls)[2], dim(termCorrNulls)[1])),
#            numTerms = factor(rep(2:(dim(termCorrNulls)[1]+1), each = dim(termCorrNulls)[2])),
#            corrs = as.vector(t(termCorrNulls)))
#
# was a matrix (numTerms x samples)
# so for 2:50 sizes of term sets, and 1000 random pairings to compute correlation on
# based on gene counts across documents in corpus
#
#     draw numTerms         corrs
# 1      1        2  0.0590920266
# 2      2        2 -0.4093004872
# 3      3        2 -0.0964689721
# 4      4        2  0.0803941294
# 5      5        2  0.0127734074
#
# predictions = `topicTermCorrelation()` output data frame
#
#            cor numTerms      zscore topic         pval
# 1   0.25002099       21  2.93820765     1 0.0074220276
# 2   0.29238589       27  4.83558037     2 0.0006501404
# 3   0.16600159       34  0.96881060     3 0.1725689210
# 4   0.27238757       30  4.71473737     4 0.0008689495
# 5   0.11743073       21 -0.66700411     5 0.7307249799
#
# procedure:
# gets the null distribution terms for a given term set size
# typically 1000 correlation values long.
# Computes a probability density distribution
# and integrates from the predicted correlation value
# from the predictions data frame.
#
# Thus returning a p-value for each predicted topic correlation
# based on the null distribution of randomly sampled gene set correlations.
#
computeCorrPval <- function(null, predictions, densityPlots = TRUE) {
  # null <- `dat`
  # predictions <- `LDA.k50.termCor.geneCounts`

  # zscores <- lapply((1:nrow(predictions)), function(i) {
  #   nterms <- predictions[i,"numTerms"]
  #   corr <- predictions[i,"cor"]
  #   nullsVals <- null[null$numTerms == nterms,]$corrs
  #   u <- mean(nullsVals)
  #   s <- sd(nullsVals)
  #   z <- (corr - u)/s
  #   z
  # })

  probs <- lapply((1:nrow(predictions)), function(i) {
    nterms <- predictions[i,"numTerms"]
    corr <- predictions[i,"cor"]
    nullsVals <- null[null$numTerms == nterms,]$corrs
    PDF <- approxfun(density(nullsVals), rule=2)

    if (densityPlots) {
      plot(PDF)
      points(corr, PDF(corr))
    }

    p <- integrate(PDF, lower=corr, upper=1)
    p$value
  })

  return(unlist(probs))
}


# get the correlation matrix between two topic-term betas
# The betas must have the same terms (columns) in the same order
#
# Typically the first beta is the topic-term beta for the predicted
# topics based on a fitted topic model. The second can be the same beta
# (for assessing cross-correlation of topics) or a ground truth beta to
# assess the concordance of predicted topics with known topics.
#
# thresh = can be used to only compare specific genes. Applied to select
# genes (i.e., terms) in the first beta.
# If a numeric fraction (0.01, let's say), then will be applied to only use
# genes with a probability in each given topic above that threshold.
#
correlationBetweenBetas <- function(beta1, beta2, thresh = NULL) {

  betaCorMtx <- do.call(rbind, lapply(1:nrow(beta1), function(i) {
    # if choosing top genes for a topic
    if (is.null(thresh)) {
      genes <- colnames(beta1)
    } else if (is.numeric(thresh)) {
      genes <- which(beta1[i,] > 0.01)
    } else {
      stop("thresh must be NULL or fraction > 0 and < 1")
    }
    # print(genes)
    sapply(1:nrow(beta2), function(j) {
      cor(beta1[i,genes], beta2[j,genes])
    })
  }))

  rownames(betaCorMtx) <- rownames(beta1)
  colnames(betaCorMtx) <- rownames(beta2)
  return(betaCorMtx)
}


# similar to `correlationBetweenBetas` but for the theta (doc-topic proportions)
# for the thetas the topics are the columns instead of the rows like in the betas
#
correlationBetweenThetas <- function(theta1, theta2, thresh = NULL) {

  thetaCorMtx <- do.call(rbind, lapply(1:ncol(theta1), function(i) {
    # # if choosing top genes for a topic
    # if (is.null(thresh)) {
    #   genes <- colnames(beta1)
    # } else if (is.numeric(thresh)) {
    #   genes <- which(beta1[i,] > 0.01)
    # } else {
    #   stop("thresh must be NULL or fraction > 0 and < 1")
    # }
    # print(genes)
    sapply(1:ncol(theta2), function(j) {
      cor(theta1[,i], theta2[,j])
    })
  }))

  rownames(thetaCorMtx) <- colnames(theta1)
  colnames(thetaCorMtx) <- colnames(theta2)
  return(thetaCorMtx)
}


# adjust values to 0-1 range
# x = can be vector or matrix
#
scale0_1 <- function(x) {
  xAdj <- (x - min(x)) / diff(range(x))
  return(xAdj)
}


# compute distance between 2 matrices
#
## frobenius norm of a matrix:
# sqrt(sum(M^2))
# or: norm(M, type = "F")
#
computeDistance <- function(mtx1, mtx2, normType = "F") {
  # d <- abs(norm(mtx1, type = normType) - norm(mtx2, type = normType))
  d <- Frobenius(mtx1, mtx2)
  return(d)
}




# preprocess ST count matrices
#
# input should be spot (row) x gene (columns) mtx with raw gene counts
#
# align = optional 3x3 alignment file to adjust spot coordinates to match
#         image pixels in order to properly align spots with image
# nTopGenes = number of top expressed genes to remove
#
# if the spot IDs are the positional information, extract:
#  ex: "31x20" extracted as 31, 20 for x and y
# extractPos = TRUE determines if this is done
#
# returns list with corpus (docs x genes) matrix,
# corpus in slm format, and spot positions.
#
# positions can be adjusted by the alignment matrices supplied if desired.
#
# remove = character vector of gene names to remove, based on grepl regular expression
#          ex: c("^mt-") or c("^MT", "^RPL", "^MRPL")
#
preprocess <- function(dat, alignFile = NA, extractPos = TRUE, nTopGenes = 5, remove = NA) {

  # if dat is a path to a file, read it, else assume matrix
  # definitely clean this up later..
  if (typeof(dat) == "character") {
    counts <- read.table(dat)
  } else {
    counts <- dat
  }


  # remove poor spots and genes
  countsClean <- MERINGUE::cleanCounts(counts = t(counts),
                                       min.reads = 100,
                                       min.lib.size = 100,
                                       plot=TRUE,
                                       verbose=TRUE)

  # get spot positions
  if (extractPos) {
    positions <- do.call(rbind, lapply(colnames(countsClean), function(spotID) {
      coords <- as.numeric(strsplit(spotID, "x")[[1]])
      coords
    }))
    colnames(positions) <- c("x", "y")
    rownames(positions) <- colnames(countsClean)
  } else {
    positions <- NULL
  }

  # adjust spot coordinates, optional
  if (is.na(alignFile) == FALSE) {
    align <- matrix(unlist(read.table(alignFile)), nrow = 3, ncol = 3)
    (positions[,"x"] * align[1,1]) - 290
    (positions[,"y"] * align[2,2]) - 290
  }

  # remove top expressed genes and mt genes
  top_expressed <- names(rowSums(countsClean)[order(rowSums(countsClean),
                                                    decreasing = TRUE)][1:nTopGenes])
  countsCleanFilt <- countsClean[rownames(countsClean) %in% top_expressed == FALSE,]
  print(paste("after removing top genes:", dim(countsCleanFilt)[1], dim(countsCleanFilt)[2]))

  # remove specific genes (if there are any)
  if (is.na(remove) == FALSE) {
    countsCleanFilt <- countsCleanFilt[!grepl(paste(remove, collapse="|"), rownames(countsCleanFilt)),]
    print(paste("after removing selected genes:", dim(countsCleanFilt)[1], dim(countsCleanFilt)[2]))
  }

  # remove genes that appear in every document
  countsCleanFilt_ <- countsCleanFilt
  countsCleanFilt_[which(countsCleanFilt_ > 0)] <- 1
  countsCleanFilt <- countsCleanFilt[which(rowSums(countsCleanFilt_) < 260),]
  print(paste("after removing genes present in every spot:", dim(countsCleanFilt)[1], dim(countsCleanFilt)[2]))

  # get variable genes
  par(mfrow=c(4,2), mar=c(1,1,1,1))
  countsCleanFiltOD <- getOverdispersedGenes(countsCleanFilt,
                                             plot = TRUE,
                                             details = TRUE)

  corpus <- t(as.matrix(countsCleanFilt[countsCleanFiltOD$ods,]))
  corpus_slm <- slam::as.simple_triplet_matrix(corpus)

  return(list(corpus = corpus,
              slm = corpus_slm,
              pos = positions))

}


# transform corpus raw gene counts into binned counts.
# Genes are binned based on their raw counts wrt the entire corpus
# and the bin is used as the new gene count value
#
binnedCorpus <- function(corpus, nbins = 50) {

  lin <- as.vector(corpus)

  b <- seq(from = min(lin), to = max(lin), by=max(lin)/nbins)
  b <- append(c(-Inf), b)

  l <- append(c(0), seq(nbins))

  binnedCorpus <- matrix(as.numeric(as.vector(cut(corpus, breaks = b, labels = l))),
                         nrow = dim(corpus)[1],
                         ncol = dim(corpus)[2])

  rownames(binnedCorpus) <- rownames(corpus)
  colnames(binnedCorpus) <- colnames(corpus)

  return(binnedCorpus)
}


