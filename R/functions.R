#' Restrict to informative words (genes) for topic modeling
#'
#' @description identifies over dispersed genes across pixels to use as informative words
#'     (genes) in topic modeling. Also allows ability to restrict over dispersed genes
#'     to those that occur in more than and/or less than selected fractions of pixels in corpus.
#'
#' @param counts genes x pixels gene count matrix
#' @param removeAbove remove over dispersed genes that are present in more than this fraction of pixels (default: 1.0)
#' @param removeBelow remove over dispersed genes that are present in less than this fraction of pixels (default: 0.05)
#' @param alpha alpha parameter for `getOverdispersedGenes`.
#'     Higher = less stringent and more overdispersed genes returned (default: 0.05)
#' @param plot return histogram plots of genes per pixel and pixels per genes
#'     for over dispersed genes and after corpus restriction. (default: FALSE)
#' @param verbose (default: TRUE)
#'
#' @export
restrictCorpus <- function(counts,
                           removeAbove = 1.0,
                           removeBelow = 0.05,
                           alpha = 0.05,
                           plot = FALSE,
                           verbose = TRUE) {
  
  ## overdispersed genes only
  if(verbose) {
    print(paste0('Restricting to overdispersed genes with alpha=', alpha, '...'))
  }
  odGenes <- getOverdispersedGenes(counts,
                                   alpha = alpha,
                                   plot = plot,
                                   details = TRUE,
                                   verbose = verbose)
  countsFilt <- counts[odGenes$ods,]
  
  ## remove genes that are present in more than X% of pixels
  vi <- rowSums(as.matrix(countsFilt) > 0) >= ncol(countsFilt)*removeAbove
  if(verbose) {
    print(paste0('Removing ', sum(vi), ' genes present in ', removeAbove*100, '% or more of pixels...'))
  }
  countsFilt_ <- countsFilt[!vi,]
  if(verbose) {
    print(paste0(nrow(countsFilt_), ' genes remaining...'))
  }
  
  ## remove genes that are present in less than X% of pixels
  vi <- rowSums(as.matrix(countsFilt_) > 0) <= ncol(countsFilt_)*removeBelow
  if(verbose) {
    print(paste0('Removing ', sum(vi), ' genes present in ', removeBelow*100, '% or less of pixels...'))
  }
  countsFiltRestricted <- countsFilt_[!vi,]
  if(verbose) {
    print(paste0(nrow(countsFiltRestricted), ' genes remaining...'))
  }
  if (plot) {
    par(mfrow=c(1,2), mar=rep(5,4))
    hist(log10(Matrix::colSums(countsFiltRestricted)+1), breaks=20, main='Genes Per Pixel')
    hist(log10(Matrix::rowSums(countsFiltRestricted)+1), breaks=20, main='Pixels Per Gene')
  }
  
  return(countsFiltRestricted)
}


#' Find the optimal number of cell-types K for the LDA model
#'
#' @description The input for topicmodels::LDA needs to be a
#'     slam::as.simple_triplet_matrix (docs x words). Access a given model in
#'     the returned list via: lda$models[[k]][1]. The models are objects from
#'     the library `topicmodels`. The LDA models have slots with additional
#'     information.
#'
#' @param counts Gene expression counts with pixels as rows and genes as columns
#' @param Ks vector of K parameters to search
#' @param seed Random seed
#' @param testSize fraction of pixels to set aside for test corpus when computing perplexity (default: NULL)
#'    Either NULL or decimal between 0 and 1.
#' @param perc.rare.thresh the number of deconvolved cell-types with mean pixel proportion below this fraction used to assess
#'     performance of fitted models for each K. Recorded for each K. (default: 0.05)
#' @param ncores Number of cores for parallelization
#' @plot Boolean for plotting
#'
#' @return A list that contains
#' \itemize{
#' \item models: each fitted LDA model for a given K
#' \item kneedOptK: the optimal K based on Kneed algorithm
#' \item minOptK: the optimal K based on minimum
#' \item ctPropOptK: Suggested upper bound on K. K in which number of returned cell-types
#'     with mean proportion < 5% starts to increases steadily.
#' \item numRare: number of cell-types with mean pixel proportion < 5% for each K
#' \item perplexities: perplexity scores for each model
#' \item corpus: the corpus that was used to fit each model
#' }
#'
#' @export
fitLDA <- function(counts, Ks = seq(2, 10, by = 2),
                   seed = 0, testSize = NULL, perc.rare.thresh = 0.05,
                   ncores = parallel::detectCores(logical = TRUE) - 1,
                   plot = TRUE, verbose = TRUE) {
  
  counts <- as.matrix(counts)
  
  if (is.null(testSize)){
    set.seed(seed)
    testingPixels <- seq(nrow(counts))
    fittingPixels <- seq(nrow(counts))
  } else if ((0 < testSize) & (testSize < 1.0)){
    cat("Splitting pixels into", testSize*100, "% and", 100-testSize*100, "% testing and fitting corpuses", "\n")
    set.seed(seed)
    testingPixels <- sample(nrow(counts), round(nrow(counts)*testSize))
    fittingPixels <- seq(nrow(counts))[-testingPixels]
  } else {
    stop("`testSize` must be NULL or decimal between 0 and 1")
  }
  
  ## counts must be pixels (rows) x genes (cols) matrix
  corpus <- slam::as.simple_triplet_matrix((as.matrix(counts)))
  corpusFit <- slam::as.simple_triplet_matrix((as.matrix(counts[fittingPixels,])))
  corpusTest <- slam::as.simple_triplet_matrix((as.matrix(counts[testingPixels,])))
  
  # if (slam::is.simple_triplet_matrix(counts) == TRUE){
  #   corpus <- counts
  # } else {
  #   corpus <- slam::as.simple_triplet_matrix(t(as.matrix(counts)))
  # }
  
  if (verbose == TRUE){
    verbose <- 1
  } else {
    verbose <- 0
  }
  
  controls <- list(seed = seed,
                   verbose = verbose, keep = 1, estimate.alpha = TRUE)
  
  start_time <- Sys.time()
  fitted_models <- parallel::mclapply(Ks, function(k) {
    topicmodels::LDA(corpusFit, k=k, control = controls)
  },
  mc.cores = ncores
  )
  names(fitted_models) <- Ks
  
  if(verbose) {
    total_t <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
    print(sprintf("Time to fit LDA models was %smins", total_t))
  }
  
  if(verbose) {
    print("Computing perplexity for each fitted model...")
  }
  
  pScores <- unlist(lapply(fitted_models, function(model){
    p <- topicmodels::perplexity(model, newdata = corpusTest)
    return(p)
  }))
  
  ## Kneed algorithm
  kOpt1 <- Ks[where.is.knee(pScores)]
  ## Min
  kOpt2 <- Ks[which(pScores == min(pScores))]
  
  ## check number of predicted cell-types at low proportions
  out <- lapply(1:length(Ks), function(i) {
    apply(getBetaTheta(fitted_models[[i]], corpus = corpus)$theta, 2, mean)
  })
  ## number of cell-types present at fewer than `perc.rare.thresh` on average across pixels
  numrare <- unlist(lapply(out, function(x) sum(x < perc.rare.thresh)))
  kOpt3 <- Ks[where.is.knee(numrare)]
  
  if(plot) {
    # plot(Ks, pScores,
    #      ylab = "perplexity",
    #      xlab = "K",
    #      main = "K vs perplexity")
    # abline(v = kOpt1, col='blue')
    # abline(v = kOpt2, col='red')
    # legend(x = "top",
    #        legend = c("kneed", "min"), col = c("blue", "red"), lty = 1, lwd = 1)
    #
    # plot(Ks, numrare,
    #      ylab = paste0("number of cell-types", "\n","with < 5% mean proportion"),
    #      xlab = "K",
    #      main = "K vs number of rare predicted cell-types",
    #      ylim = c(min=0, max=round((max(numrare)+1)*1.5)) )
    # abline(v = kOpt3, col='red')
    # legend(x = "top",
    #        legend = c("kneed"), col = c("red"), lty = 1, lwd = 1)
    
    dat <- data.frame(K = as.double(Ks),
                      rareCts = numrare,
                      perplexity = pScores,
                      rareCtsAdj = scale0_1(numrare),
                      perplexAdj = scale0_1(pScores))
    
    prim_ax_labs <- seq(min(dat$rareCts), max(dat$rareCts))
    prim_ax_breaks <- scale0_1(prim_ax_labs)
    ## if number rareCts stays constant, then only one break. scale0_1(prim_ax_labs) would be NaN so change to 0
    if(length(prim_ax_labs) == 1){
      prim_ax_breaks <- 0
      ## also the rareCtsAdj <- scale0_1(rareCts) would be NaN, so set to 0, so at same position as the tick,
      ## and its label will still be set to the constant value of rareCts
      dat$rareCtsAdj <- 0
    }
    if(max(dat$rareCts) < 1){
      sec_ax_labs <- seq(min(dat$perplexity), max(dat$perplexity), (max(dat$perplexity)-min(dat$perplexity))/1)
    } else {
      sec_ax_labs <- seq(min(dat$perplexity), max(dat$perplexity), (max(dat$perplexity)-min(dat$perplexity))/max(dat$rareCts))
    }
    sec_ax_breaks <- scale0_1(sec_ax_labs)
    
    plt <- ggplot2::ggplot(dat, ggplot2::aes(x=K)) +
      ggplot2::geom_line(ggplot2::aes(y=rareCtsAdj), col="blue", lwd = 2) +
      ggplot2::geom_line(ggplot2::aes(y=perplexAdj), col="red", lwd = 2) +
      ggplot2::scale_y_continuous(name=paste0("# cell-types with mean proportion < ", round(perc.rare.thresh*100, 2), "%"), breaks = prim_ax_breaks, labels = prim_ax_labs,
                                  sec.axis=ggplot2::sec_axis(~ ., name="perplexity", breaks = sec_ax_breaks, labels = round(sec_ax_labs, 2))) +
      ggplot2::scale_x_continuous(breaks = min(dat$K):max(dat$K)) +
      ggplot2::ggtitle("Fitted model K's vs deconvolved cell-types and perplexity") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size=15, face=NULL),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(color = "black", size = 0.1),
        axis.title.y.left = ggplot2::element_text(color="blue", size = 13),
        axis.text.y.left = ggplot2::element_text(color="blue", size = 13),
        axis.title.y.right = ggplot2::element_text(color="red", size = 15, vjust = 1.5),
        axis.text.y.right = ggplot2::element_text(color="red", size = 13),
        axis.text.x = ggplot2::element_text(angle = 0, size = 13),
        axis.title.x = ggplot2::element_text(size=13)
      )
    print(plt)
  }
  
  return(list(models = fitted_models,
              kneedOptK = kOpt1,
              minOptK = kOpt2,
              ctPropOptK = kOpt3,
              numRare = numrare,
              perplexities = pScores,
              fitCorpus = corpusFit,
              testCorpus = corpusTest))
}


#' Pull out cell-type proportions across pixels (theta) and
#' cell-type gene probabilities (beta) matrices from fitted LDA models from fitLDA
#'
#' @param lda an LDA model from `topicmodels`. From list of models returned by
#'     fitLDA
#' @param corpus Gene expression counts with pixels as rows and genes as columns
#'     for which to get predicted cell-type proportions and cell-type gene expression profiles.
#'     Typically same one that was used to fit the LDA model.
#'
#' @return A list that contains
#' \itemize{
#' \item beta: cell-type (rows) by gene (columns) distribution matrix.
#'     Each row is a probability distribution of a cell-type expressing each gene
#'     in the corpus
#' \item theta: pixel (rows) by cell-types (columns) distribution matrix. Each row
#'     is the cell-type composition for a given pixel
#' }
#'
#' @export
getBetaTheta <- function(lda, corpus) {
  
  result <- topicmodels::posterior(lda, newdata = corpus)
  theta <- result$topics
  beta <- result$terms
  # topicFreqsOverall <- colSums(theta) / length(lda@documents)
  
  return(list(beta = beta,
              theta = theta))
}


#' Aggregate cell-types together using dynamic tree cutting.
#'
#' @param beta Beta matrix (cell-type gene distribution matrix)
#' @param distance Distance measure to be used (default: euclidean)
#' @param clustering Clustering agglomeration method to be used (default: ward.D)
#' @param dynamic Dynamic tree cutting method to be used (default: hybrid)
#' @param deepSplit Dynamic tree cutting sensitivity parameter (default: 4)
#' @param plot Boolean for plotting
#'
#' @return A list that contains
#' \itemize{
#' \item order: vector of the dendrogram index order for the cell-types
#' \item clusters: factor of the cell-types (names) and their assigned cluster (levels)
#' \item dendro: dendrogram of the clusters
#' }
#'
#' @export
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
  hc_ <- stats::hclust(d_, method = clustering)
  
  groups <- dynamicTreeCut::cutreeDynamic(hc_,
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
            col = correlation_palette,
            xlab = "Cell-types",
            ylab = "Genes",
            main = "Predicted cell-types and clusters")
  }
  
  return(list(clusters = groups,
              order = hc_$order,
              dendro = as.dendrogram(hc_)))
  
}


#' Aggregate cell-types in the same cluster into a single cell-type
#'
#' @description Note: for the beta matrix, each row is a cell-type, each column
#'     is a gene. The cell-type row is a distribution of terms that sums to 1. So
#'     combining cell-type row vectors, these should be adjusted such that the
#'     rowSum == 1. As in, take average of the terms after combining. However, the
#'     theta matrix (after inversion) will have cell-type rows and each column is a
#'     document. Because the cell-type can be represented at various proportions in
#'     each doc and these will not necessarily add to 1, should not take average.
#'     Just sum cell-type row vectors together. This way, each document column still
#'     adds to 1 when considering the proportion of each cell-type-cluster in the document.
#'
#' @param mtx either a beta (cell-type gene distribution matrix) or a
#'     theta (pixel-cell type distribution matrix)
#' @param clusters factor of the cell-types (names) and their assigned cluster (levels)
#' @param type either "t" or "b". Affects the adjustment to the combined
#'     cell-type vectors. "b" divides summed cell-type vectors by number of aggregated cell-types
#'
#' @return matrix where cell-types are now cell-type-clusters
#'
#' @export
combineTopics <- function(mtx, clusters, type) {
  
  if (!type %in% c("t", "b")){
    stop("`type` must be either 't' or 'b'")
  }
  
  # if mtx is theta, transpose so cell-types are rows
  if (type == "t") {
    mtx <- t(mtx)
  }
  
  combinedTopics <- do.call(rbind, lapply(levels(clusters), function(cluster) {
    
    # get cell-types in a given cluster
    topics <- labels(clusters[which(clusters == cluster)])
    
    # print(cluster)
    # print(length(topics))
    
    # matrix slice for the cell-types in the cluster
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
  print("cell-types combined.")
  
  # if theta, make cell-types the columns again
  if (type == "t") {
    combinedTopics <- t(combinedTopics)
  }
  return(combinedTopics)
}


#' Get the optimal LDA model
#'
#' @param models list returned from fitLDA
#' @param opt either "kneed" (kOpt1) or "min" (kOpt2), or designate a specific K.
#'     "kneed" = K vs perplexity inflection point.
#'     "min" = K corresponding to minimum perplexity
#'     "proportion" = K vs number of cell-type with mean proportion < 5% inflection point
#'
#' @return optimal LDA model fitted to the K based on `opt`
#'
#' @export
optimalModel <- function(models, opt) {
  
  if (opt == "kneed"){
    m <- models$models[[which(sapply(models$models, slot, "k") == models$kneedOptK)]]
  } else if (opt == "min"){
    m <- models$models[[which(sapply(models$models, slot, "k") == models$minOptK)]]
  } else if (opt == "proportion"){
    m <- models$models[[which(sapply(models$models, slot, "k") == models$ctPropOptK)]]
  } else if (opt %in% sapply(models$models, slot, "k")){
    m <- models$models[[which(sapply(models$models, slot, "k") == opt)]]
  } else {
    stop("`opt` must be either 'kneed',  'min', 'proportion', or int of a fitted K")
  }
  return(m)
}


#' Wrapper to extract beta (cell type-gene distribution matrix),
#' theta (pixel cell-type distribution) for individual cell-types and aggregated cell-type-clusters
#' for an LDA model in `fitLDA` output list.
#'
#' @description Wrapper that combines the functions `getBetaTheta`, `clusterTopics`,
#'     `combineTopics` and slots of the topicmodels::LDA object to return a list
#'     that contains the most relevant components of a given LDA model for ease
#'     of analysis and visualization.
#'
#' @param LDAmodel LDA model from fitLDA
#' @param corpus Gene expression counts with pixels as rows and genes as columns
#'     for which to get predicted cell-type proportions and cell-type gene expression profiles.
#'     Typically same one that was used to fit the LDA model.
#' @param clustering Clustering agglomeration method to be used (default: ward.D)
#' @param dynamic Dynamic tree cutting method to be used (default: hybrid)
#' @param deepSplit parameter for `clusterTopics` for dynamic tree splitting
#'     when aggregating cell-types (default: 4)
#' @param colorScheme color scheme for generating colors assigned to cell-type
#'     clusters for visualizing. Either "rainbow" or "ggplot" (default: "rainbow")
#' @param plot Boolean for plotting
#'
#' @return A list that contains
#' \itemize{
#' \item beta: cell-type (rows) by gene (columns) distribution matrix.
#'     Each row is a probability distribution of a cell-type expressing each gene
#'     in the corpus
#' \item theta: pixel (rows) by cell-type (columns) distribution matrix. Each row
#'     is the cell-type composition for a given pixel
#' \item clusters: factor of the cell-types (names) and their assigned cluster (levels)
#' \item dendro: dendrogram of the clusters. Returned from `stats::hclust()` in `clusterTopics`
#' \item cols: factor of colors for each cell-type where colors correspond to their assigned cluster
#' \item betaCombn: cell-type (rows) by gene (columns) distribution matrix for combined cell-type-clusters
#' \item thetaCombn: pixel (rows) by cell-type (columns) distribution matrix for aggregated cell-type-clusters
#' \item clustCols: factor of colors for each cell-type-cluster
#' \item k: number of cell-types K of the model
#' }
#'
#' @export
buildLDAobject <- function(LDAmodel,
                           corpus,
                           clustering = "ward.D",
                           dynamic = "hybrid",
                           deepSplit = 4,
                           colorScheme = "rainbow",
                           plot = TRUE){
  
  # get beta and theta list object from the LDA model
  m <- getBetaTheta(LDAmodel, corpus = corpus)
  
  # cluster cell-types
  clust <- clusterTopics(beta = m$beta,
                         clustering = clustering,
                         dynamic = dynamic,
                         deepSplit = deepSplit,
                         plot = plot)
  
  # add cluster information to the list
  m$clusters <- clust$clusters
  m$dendro <- clust$dendro
  
  # colors for the cell-types. Essentially colored by the cluster they are in
  cols <- m$clusters
  if (colorScheme == "rainbow"){
    levels(cols) <- rainbow(length(levels(cols)))
  }
  if (colorScheme == "ggplot"){
    levels(cols) <- gg_color_hue(length(levels(cols)))
  }
  m$cols <- cols
  
  # construct beta and thetas for the cell-type-clusters
  m$betaCombn <- combineTopics(m$beta, clusters = m$clusters, type = "b")
  m$thetaCombn <- combineTopics(m$theta, clusters = m$clusters, type = "t")
  
  # colors for the cell-type-clusters
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


#' Function to get Hungarian sort pairs via clue::lsat
#'
#' @description Finds best matches between cell-types that correlate between
#'     beta or theta matrices that have been compared via `getCorrMtx`.
#'     Each row is paired with a column in the output matrix from `getCorrMtx`.
#'     If there are less rows than columns, then some columns will not be
#'     matched and not part of the output.
#'
#' @param mtx output correlation matrix from `getCorrMtx`. Must not have more rows
#'     than columns
#'
#' @return A list that contains
#' \itemize{
#' \item pairs: output of clue::solve_LSAP. A vectorized object where for each
#'     position the first element is a row and the second is the paired column.
#' \item rowix: the indices of the rows. Essentially seq_along(pairing)
#' \item colsix: the indices of each column paired to each row
#' }
#'
#' @export
lsatPairs <- function(mtx){
  # must have equal or more rows than columns
  # values in matrix converted to 0-1 scale relative to all values in mtx
  pairing <- clue::solve_LSAP(scale0_1(mtx), maximum = TRUE)
  # clue::lsat returns vector where for each position the first element is a row
  # and the second is the paired column
  rowsix <- seq_along(pairing)
  colsix <- as.numeric(pairing)
  
  return(list(pairs = pairing,
              rowix = rowsix,
              colsix = colsix))
}

