#' Restrict to informative words (genes) for topic modeling
#'
#' @export
#'
restrictCorpus <- function(counts,
                           t = 0.9,
                           alpha = 0.05,
                           plot = FALSE,
                           verbose = TRUE) {
  ## overdispersed genes only
  odGenes <- getOverdispersedGenes(counts,
                                   alpha = alpha,
                                   plot = plot,
                                   details = TRUE,
                                   verbose = verbose)
  countsFilt <- counts[odGenes$ods,]

  ## remove genes that are present in more than X% of spots
  vi <- rowSums(countsFilt > 0) >= ncol(countsFilt)*t
  if(verbose) {
    print(paste0('Removing genes present in more than ', t*100, '% of datasets...'))
  }
  countsFiltnoUnifGenes <- countsFilt[!vi,]

  if(verbose) {
    print(paste0(nrow(countsFiltnoUnifGenes), ' genes remaining...'))
  }
  if (plot) {
    par(mfrow=c(1,2), mar=rep(5,4))
    hist(log10(Matrix::colSums(countsFiltnoUnifGenes)+1), breaks=20, main='Genes Per Dataset')
    hist(log10(Matrix::rowSums(countsFiltnoUnifGenes)+1), breaks=20, main='Datasets Per Gene')
  }

  return(countsFiltnoUnifGenes)
}


#' Find the optimal number of topics K for topic modeling
#'
#' @param counts Gene expression counts with genes as rows
#' @param Ks vector of K parameters to search
#' @param seed Random seed
#' @param ncores Number of cores for parallelization
#' @plot Boolean for plotting
#'
#' @return A list that contains
#' \itemize{
#' \item models: each fitted LDA model
#' \item kOpt1: the optimal K based on Kneed algorithm
#' \item kOpt2: the optimal K based on minimum
#' \item perplexities: perplexity scores for each model
#' \item corpus: the corpus that was used to fit each model
#' }
#'
#' @export
#'
fitLDA <- function(counts, Ks = seq(2, 10, by = 2), seed = 0,
                   ncores = parallel::detectCores(logical = TRUE) - 1,
                   plot = TRUE) {

  # topic modeling
  corpus <- slam::as.simple_triplet_matrix(t(as.matrix(counts)))

  controls <- list(seed = seed,
                   verbose = 1, keep = 1,
                   alpha = 1, estimate.alpha = TRUE)

  fitted_models <- parallel::mclapply(Ks, function(k) {
    topicmodels::LDA(corpus, k=k, control = controls)
  },
  mc.cores = ncores
  )
  names(fitted_models) <- Ks

  pScores <- unlist(lapply(fitted_models, function(model){
    p <- topicmodels::perplexity(model, corpus)
    return(p)
  }))

  ## Kneed algorithm
  kOpt1 <- Ks[where.is.knee(pScores)]
  ## Min
  kOpt2 <- Ks[which(pScores == min(pScores))]

  if(plot) {
    plot(Ks, pScores)
    abline(v = kOpt1, col='blue')
    abline(v = kOpt2, col='red')
  }

  return(list(models = fitted_models,
              kOpt1 = kOpt1,
              kOpt2 = kOpt2,
              perplexities = pScores,
              fitCorpus = corpus))
}


#' Pull out topics and terms from fitted topic models from fitLDA
#' Jean: Brendan, the use of theta and beta vs topics and terms
#' across different functions is quite confusing
#' Please update for consistency
#'
#' @export
#'
getBetaTheta <- function(lda) {

  lda.tmResult <- topicmodels::posterior(lda)
  theta <- lda.tmResult$topics
  beta <- lda.tmResult$terms
  topicFreqsOverall <- colSums(theta) / length(lda@documents)

  return(list(beta = beta,
              theta = theta,
              topicFreq = topicFreqsOverall))

}


#' Cluster topics together using dynamic tree cutting.
#'
#' @param beta Beta matrix (topic-word distributions)
#' @param distance Distance measure to be used (default: euclidean)
#' @param clustering Clustering agglomeration method to be used (default: ward.D)
#' @param dynamic Dynamic tree cutting method to be used (default: hybrid)
#' @param deepSplit Dynamic tree cutting sensitivity parameter (default: 4)
#' @param plot Boolean for plotting
#'
#' @return A list that contains
#' \itemize{
#' \item order = vector of the dendrogram index order for the topics
#' \item clusters = factor of the topics (names) and their assigned cluster (levels)
#' \item dendro = dendrogram of the clusters
#' }
#'
#' @export
#'
clusterTopics <- function(beta,
                          #distance = "euclidean",
                          clustering = "ward.D",
                          dynamic = "hybrid",
                          deepSplit = 4,
                          plot = TRUE) {

  if (deepSplit == 4) {
    maxCoreScatter = 0.95
    minGap = (1 - maxCoreScatter) * 3/4
  } else if (deepSplit == 3) {
    maxCoreScatter = 0.91
    minGap = (1 - maxCoreScatter) * 3/4
  } else if (deepSplit == 2) {
    maxCoreScatter = 0.82
    minGap = (1 - maxCoreScatter) * 3/4
  } else if (deepSplit == 1) {
    maxCoreScatter = 0.73
    minGap = (1 - maxCoreScatter) * 3/4
  } else if (deepSplit == 0) {
    maxCoreScatter = 0.64
    minGap = (1 - maxCoreScatter) * 3/4
  }

  #d_ <- dist(beta, method = distance)
  ## Jean: use correlation instead
  d_ <- as.dist(1-cor(t(beta)))
  hc_ <- hclust(d_, method = clustering)

  groups <- cutreeDynamic(hc_,
                          method = dynamic,
                          distM = as.matrix(d_),
                          deepSplit = deepSplit,
                          minClusterSize=0,
                          maxCoreScatter = maxCoreScatter,
                          minGap = minGap,
                          maxAbsCoreScatter=NULL,
                          minAbsGap=NULL)

  names(groups) <- hc_$labels
  groups <- factor(groups)

  if (plot) {
    #plot(hc_)
    d2_ <- as.dist(1-cor(beta))
    rc_ <- hclust(d2_, method = clustering)
    heatmap(t(beta),
         Colv=as.dendrogram(hc_),
         Rowv=as.dendrogram(rc_),
         ColSideColors = fac2col(groups),
         col = correlation_palette)
  }

  return(list(clusters = groups,
              order = hc_$order,
              dendro = as.dendrogram(hc_)))

}


#' combine topics within same cluster.
#' the term beta distriubtions are added
#' or the theta document proportions are added
#'
#' mtx = needs to be topic (row) x feature (term or document) (column) format
#' clusters = factor of the topics (names) and their assigned cluster (levels)
#'
#' mtxType = either "b" (beta) or "t" (theta). Affects the adjustment to the
#'           combined topic vectors. "b" divides summed topic vectors by number
#'           of combined topics.
#'
combineTopics <- function(mtx, clusters, mtxType) {

  if (mtxType == "t") {
    mtx <- t(mtx)
  }

  combinedTopics <- do.call(rbind, lapply(levels(clusters), function(cluster) {
    # select the topics in the cluster in the same order of the dendrogram
    topics <- labels(clusters[which(clusters == cluster)])

    print(cluster)
    print(length(topics))

    # topic-term distribution reordered based on dendro and selected cluster topics
    mtx_slice <- mtx[topics,]

    if (length(topics) == 1) {
      topicVector <- mtx_slice
    } else if (length(topics) > 1) {
      mtx_slice <- as.data.frame(mtx_slice)

      # for the beta matrix, each row is a topic, each column is a gene.
      # The topic row is a distribution of terms that sums to 1.
      # So combining topic row vectors, these should be adjusted such that the
      # rowSum == 1. As in, take average of the terms after combining.

      # However, the theta matrix (after inversion) will have topic rows
      # and each column is a document. Because the topic can be represented
      # at various proportions in each doc and these will not necessarily add
      # to 1, should not take average. Just sum topic row vectors together.
      # This way, each document column still adds to 1 when considering the
      # proportion of each topic-cluster that makes up each the document.

      if (mtxType == "b") {
        topicVector <- colSums(mtx_slice) / length(topics)
      } else if (mtxType == "t") {
        topicVector <- colSums(mtx_slice)
      }

    }
    topicVector

  }))
  rownames(combinedTopics) <- levels(clusters)
  colnames(combinedTopics) <- colnames(mtx)
  print("topics combined.")

  if (mtxType == "t") {
    combinedTopics <- t(combinedTopics)
  }

  return(combinedTopics)
}


