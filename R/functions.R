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
#' @description The input for topicmodels::LDA needs to be a
#'     slam::as.simple_triplet_matrix (docs x words). Access a given model in
#'     the returned list via: lda$models[[k]][1]. The models are objects from
#'     the library `topicmodels`. The LDA models have slots with additional
#'     information.
#'
#' @param counts Gene expression counts with genes as rows or
#'     slam::simple_triplet_matrix where spots were rows and genes were columns
#' @param Ks vector of K parameters to search
#' @param seed Random seed
#' @param ncores Number of cores for parallelization
#' @plot Boolean for plotting
#'
#' @return A list that contains
#' \itemize{
#' \item models: each fitted LDA model for a given K
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

  if (slam::is.simple_triplet_matrix(counts) == TRUE){
    corpus <- counts
  } else {
    corpus <- slam::as.simple_triplet_matrix(t(as.matrix(counts)))
  }
  
  controls <- list(seed = seed,
                   verbose = 1, keep = 1, estimate.alpha = TRUE)

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


#' Pull out topic proportions across spots (theta) and
#' topic gene probabilities (beta) matrices from fitted topic models from fitLDA
#' 
#' @param lda an LDA model from `topicmodels`. From list of models returned by
#'     fitLDA
#'
#' @return A list that contains
#' \itemize{
#' \item beta: topic (rows) by gene (columns) distribution matrix.
#'     Each row is a probability distribution of a topic expressing each gene
#'     in the corpus
#' \item theta: spot (rows) by topics (columns) distribution matrix. Each row
#'     is the topic composition for a given spot
#' \item topicFreq: overall proportion of each topic in the entire corpus of
#'     spots
#' }
#'
#' @export
#'
getBetaTheta <- function(lda) {

  result <- topicmodels::posterior(lda)
  theta <- result$topics
  beta <- result$terms
  topicFreqsOverall <- colSums(theta) / length(lda@documents)

  return(list(beta = beta,
              theta = theta,
              topicFreq = topicFreqsOverall))
}


#' Cluster topics together using dynamic tree cutting.
#'
#' @param beta Beta matrix (topic-gene distribution matrix)
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


#' Collapse topics in the same cluster into a single topic
#' 
#' @description Note: for the beta matrix, each row is a topic, each column
#' is a gene. The topic row is a distribution of terms that sums to 1. So
#' combining topic row vectors, these should be adjusted such that the
#' rowSum == 1. As in, take average of the terms after combining. However, the
#' theta matrix (after inversion) will have topic rows and each column is a
#' document. Because the topic can be represented at various proportions in
#' each doc and these will not necessarily add to 1, should not take average.
#' Just sum topic row vectors together. This way, each document column still
#' adds to 1 when considering the proportion of each topic-cluster in the document.
#' 
#' @param mtx either a beta (topic-gene distribution matrix) or a
#'     theta (spot-topic distribution matrix)
#' @param clusters factor of the topics (names) and their assigned cluster (levels)
#' @param type either "t" or "b". Affects the adjustment to the combined
#'     topic vectors. "b" divides summed topic vectors by number of combined topics.
#' 
#' @return matrix where topics are now topic-clusters
#' 
#' @export
combineTopics <- function(mtx, clusters, type) {
  
  if (!type %in% c("t", "b")){
    stop("`type` must be either 't' or 'b'")
  }
  
  # if mtx is theta, transpose so topics are rows
  if (type == "t") {
    mtx <- t(mtx)
  }

  combinedTopics <- do.call(rbind, lapply(levels(clusters), function(cluster) {
    
    # get topics in a given cluster
    topics <- labels(clusters[which(clusters == cluster)])

    print(cluster)
    print(length(topics))

    # matrix slice for the topics in the cluster
    mtx_slice <- mtx[topics,]

    if (length(topics) == 1) {
      topicVector <- mtx_slice
    } else if (length(topics) > 1) {
      mtx_slice <- as.data.frame(mtx_slice)
      if (type == "b") {
        topicVector <- colSums(mtx_slice) / length(topics)
      } else if (type == "t") {
        topicVector <- colSums(mtx_slice)
      }
    }
    topicVector

  }))
  rownames(combinedTopics) <- levels(clusters)
  colnames(combinedTopics) <- colnames(mtx)
  print("topics combined.")
  
  # if theta, make topics the columns again
  if (type == "t") {
    combinedTopics <- t(combinedTopics)
  }
  return(combinedTopics)
}


#' Get the optimal LDA model
#' 
#' @param models list returned from fitLDA
#' @param opt either "kneed" (kOpt1) or "min" (kOpt2)
#' 
#' @return optimal LDA model fitted to the K based on `opt`
#' 
#' @export
optimalModel <- function(models, opt) {
  
  if (!opt %in% c("kneed", "min")){
    stop("`opt` must be either 'kneed' or 'min'")
  }
  if (opt == "kneed"){
    m <- models$models[[which(sapply(models$models, slot, "k") == models$kOpt1)]]
  }
  if (opt == "min"){
    m <- models$models[[which(sapply(models$models, slot, "k") == models$kOpt2)]]
  }
  return(m)
}


#' Wrapper to extract beta (topic-gene distribution matrix),
#' theta (spot-topic distribution) for individual topic and combined topic-clusters
#' for an LDA model in `fitLDA` output list.
#' 
#' @param LDAmodel LDA model from fitLDA
#' @param deepSplit parameter for `clusterTopics` for dynamic tree splitting
#'     when clustering topics (default: 4)
#' @param colorScheme color scheme for generating colors assigned to topic
#'     clusters for visualizing. Either "rainbow" or "ggplot" (default: "rainbow")
#' 
#' @return A list that contains
#' \itemize{
#' \item beta: topic (rows) by gene (columns) distribution matrix.
#'     Each row is a probability distribution of a topic expressing each gene
#'     in the corpus
#' \item theta: spot (rows) by topics (columns) distribution matrix. Each row
#'     is the topic composition for a given spot
#' \item topicFreq: overall proportion of each topic in the entire corpus of
#'     spots
#' \item clusters = factor of the topics (names) and their assigned cluster (levels)
#' \item dendro = dendrogram of the clusters
#' \item cols = factor of colors for each topic where colors correspond to their assigned cluster
#' \item betaCombn = topic (rows) by gene (columns) distribution matrix for combined topic-clusters
#' \item thetaCombn = spot (rows) by topic (columns) distribution matrix for combined topic-clusters
#' }
#'
#' @export
buildLDAobject <- function(LDAmodel,
                           deepSplit = 4,
                           colorScheme = "rainbow"){
  
  # get beta and theta list object from the LDA model
  m <- getBetaTheta(LDAmodel)
  
  # cluster topics
  clust <- clusterTopics(beta = m$beta,
                         deepSplit = deepSplit)
  
  # add cluster information to the list
  m$clusters <- clust$clusters
  m$dendro <- clust$dendro
  
  # colors for the topics. Essentially colored by the cluster they are in
  cols <- m$clusters
  if (colorScheme == "rainbow"){
    levels(cols) <- rainbow(length(levels(cols)))
  }
  if (colorScheme == "ggplot"){
    levels(cols) <- gg_color_hue(length(levels(cols)))
  }
  m$cols <- cols
  
  # construct beta and thetas for the topic clusters
  m$betaCombn <- combineTopics(m$beta, clusters = m$clusters, type = "b")
  m$thetaCombn <- combineTopics(m$theta, clusters = m$clusters, type = "t")
  
  # colors for the topic clusters
  # separate factor for ease of use with vizTopicClusters and others
  # note that these color assignments are different than the
  # cluster color assignments in the levels of `cols`
  clusterCols <- as.factor(colnames(m$thetaCombn))
  names(clusterCols) <- colnames(m$thetaCombn)
  levels(clusterCols) <- levels(m$cols)
  m$clustCols <- clusterCols
  
  m$k <- LDAmodel@k
  
  return(m)
  
}




