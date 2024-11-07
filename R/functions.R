#' Restrict to informative words (genes) for topic modeling
#'
#' @description identifies over dispersed genes across pixels to use as informative words
#'     (genes) in topic modeling. Also allows ability to restrict over dispersed genes
#'     to those that occur in more than and/or less than selected fractions of pixels in corpus.
#'     Limits to the top 1000 overdispersed genes in order to keep the corpus to a reasonable size.
#'
#' @param counts genes x pixels gene count matrix
#' @param removeAbove remove over dispersed genes that are present in more than this fraction of pixels (default: 1.0)
#' @param removeBelow remove over dispersed genes that are present in less than this fraction of pixels (default: 0.05)
#' @param alpha alpha parameter for getOverdispersedGenes().
#'     Higher = less stringent and more overdispersed genes returned (default: 0.05)
#' @param nTopOD number of top over dispersed genes to use. int (default: 1000).
#'     If the number of overdispersed genes is less then this number will use all of them,
#'     or set to NA to use all overdispersed genes.
#' @param plot return histogram plots of genes per pixel and pixels per genes
#'     for over dispersed genes and after corpus restriction. (default: FALSE)
#' @param verbose (default: TRUE)
#' 
#' @return a gene by pixel matrix where the remaining genes have been filtered
#' 
#' @examples 
#' data(mOB)
#' corpus <- restrictCorpus(counts = mOB$counts)
#' corpus
#'
#' @export
restrictCorpus <- function(counts,
                           removeAbove=1.0,
                           removeBelow=0.05,
                           alpha=0.05,
                           nTopOD=1000,
                           plot=FALSE,
                           verbose=TRUE) {
  
  ## remove genes that are present in more than X% of pixels
  vi <- rowSums(as.matrix(counts) > 0) >= ncol(counts)*removeAbove
  if(verbose) {
    message(paste0('Removing ', sum(vi), ' genes present in ',
                   removeAbove*100, '% or more of pixels...'))
  }
  counts <- counts[!vi,]
  if(verbose) {
    message(paste0(nrow(counts), ' genes remaining...'))
  }
  
  ## remove genes that are present in less than X% of pixels
  vi <- rowSums(as.matrix(counts) > 0) <= ncol(counts)*removeBelow
  if(verbose) {
    message(paste0('Removing ', sum(vi), ' genes present in ',
                   removeBelow*100, '% or less of pixels...'))
  }
  counts <- counts[!vi,]
  if(verbose) {
    message(paste0(nrow(counts), ' genes remaining...'))
  }
  
  ## overdispersed genes
  if(verbose) {
    message(paste0('Restricting to overdispersed genes with alpha = ', alpha, '...'))
  }
  OD <- getOverdispersedGenes(counts,
                              alpha=alpha,
                              plot=plot,
                              details=TRUE,
                              verbose=verbose)

  # option to select just the top n overdispersed genes based on
  # log p-val adjusted
  if (!is.na(nTopOD)){
    if(verbose){
      message(" Using top ", nTopOD, " overdispersed genes.", "\n")
    }
    OD_filt <- OD$df[OD$ods,]
    # check if actual number of OD genes less than `nTopOD`
    if (dim(OD_filt)[1] < nTopOD){
      if(verbose){
        message(" number of top overdispersed genes available: ", dim(OD_filt)[1], "\n")
      }
      od_genes <- rownames(OD_filt)
    } else {
      od_genes <- rownames(OD_filt[order(OD_filt$lpa),][1:nTopOD,])
    }
  } else {
    od_genes <- OD$ods
  }
  countsFiltRestricted <- counts[od_genes,]
  
  if (plot) {
    par(mfrow=c(1,2), mar=rep(5,4))
    hist(log10(Matrix::colSums(countsFiltRestricted)+1), breaks=20, main='Genes Per Pixel')
    hist(log10(Matrix::rowSums(countsFiltRestricted)+1), breaks=20, main='Pixels Per Gene')
  }
  
  if (dim(countsFiltRestricted)[1] > 1000){
    message("Genes in corpus > 1000 (", dim(countsFiltRestricted)[1],
        "). This may cause model fitting to take a while. Consider reducing the number of genes.", "\n")
  }
  return(countsFiltRestricted)
}


#' Pre-process ST pixel gene count matrices to construct corpus for input into LDA
#'
#' @description Takes pixel (row) x gene (columns) matrix and filters out poor genes
#'              and pixels. Then selects for genes to be included in final corpus for input into LDA.
#'              If the pixel IDs are made up of their positions in "XxY" these
#'              can be extracted as the pixel position coordinates (a characteristic of Stahl datasets).
#'
#'              Order of filtering options:
#'              1. Selection to use specific genes only
#'              2. `cleanCounts` to remove poor pixels and genes
#'              3. Remove top expressed genes in matrix
#'              4. Remove specific genes based on grepl pattern matching
#'              5. Remove genes that appear in more/less than a percentage of pixels
#'              6. Use the over dispersed genes computed from the remaining genes
#'                 after filtering steps 1-5 (if selected)
#'              7. Choice to use the top over dispersed genes based on -log10(p.adj)
#'
#' @param dat pixel (row) x gene (columns) mtx with gene counts OR path to it
#' @param extractPos Boolean to extract pixel positional coordinates from pixel name
#'     names
#'    (default: FALSE)
#' @param selected.genes vector of gene names to use specifically for the corpus
#'    (default: NA)
#' @param nTopGenes integer for number of top expressed genes to remove
#'     (default: NA)
#' @param genes.to.remove vector of gene names or patterns for matching to genes
#'     to remove (default: NA). ex: c("^mt-") or c("^MT", "^RPL", "^MRPL")
#' @param removeAbove non-negative numeric <=1 to use as a percentage.
#'    Removes genes present in this fraction or more of pixels (default: NA)
#' @param removeBelow non-negative numeric <=1 to use as a percentage.
#'    Removes genes present in this fraction or less of pixels (default: NA)
#' @param min.reads cleanCounts() param; minimum number of reads to keep
#'     a gene (default: 1)
#' @param min.lib.size cleanCounts() param; minimum number of counts a
#'     pixel needs to keep (default: 1)
#' @param min.detected cleanCounts() param; minimum number of pixels a gene
#'     needs to have been detected in to keep (default: 1)
#' @param ODgenes Boolean to use getOverdispersedGenes() for the corpus
#'    genes (default: TRUE)
#' @param nTopOD number of top over dispersed genes to use. int (default: 1000).
#'     If the number of overdispersed genes is less then this number will use all of them,
#'     or set to NA to use all overdispersed genes.
#' @param od.genes.alpha alpha parameter for getOverdispersedGenes().
#'     Higher = less stringent and more over dispersed genes returned (default: 0.05)
#' @param gam.k gam.k parameter for getOverdispersedGenes(). Dimension
#'     of the "basis" functions in the GAM used to fit, higher = "smoother"
#'     (default: 5)
#' @param verbose control verbosity (default: TRUE)
#' @param plot control if plots are returned (default: TRUE)
#'
#' @return A list that contains
#' \itemize{
#' \item corpus: (pixels x genes) matrix of the counts of the selected genes
#' \item slm: slam::as.simple_triplet_matrix(corpus); required format for topicmodels::LDA input
#' \item positions: matrix of x and y coordinates of pixels rownames = pixels, colnames = "x", "y"
#' }
#' 
#' @examples 
#' data(mOB)
#' cd <- mOB$counts
#' corpus <- preprocess(t(cd), removeAbove = 0.95, removeBelow = 0.05)
#' 
#' @importFrom utils read.table
#'
#' @export
preprocess <- function(dat,
                       extractPos=FALSE,
                       selected.genes=NA,
                       nTopGenes=NA,
                       genes.to.remove=NA,
                       removeAbove=NA,
                       removeBelow=NA,
                       min.reads=1,
                       min.lib.size=1,
                       min.detected=1,
                       ODgenes=TRUE,
                       nTopOD=1000,
                       od.genes.alpha=0.05,
                       gam.k=5,
                       verbose=TRUE,
                       plot=TRUE) {
  
  if (typeof(dat) == "character") {
    if (file.exists(dat)){
      counts <- read.table(dat)
    } else {
      message("path to file does not exist", "\n")
    }
  } else if (is.matrix(dat)){
    counts <- dat
  } else {
    stop("`dat` is not a viable path or matrix")
  }
  
  if(verbose) {
    message("Initial genes: ",
            dim(t(counts))[1],
            " Initial pixels: ",
            dim(t(counts))[2], "\n")
  }
  
  # use specific genes in the corpus
  if (!is.na(selected.genes[1])) {
    if(verbose){
      message("- Using genes in `selected.genes` for corpus.", "\n")
    }
    counts <- counts[,selected.genes]
    if(verbose){
      message(" ", dim(counts)[2], " genes are present in dataset.", "\n")
    }
  }
  # remove poor pixels and genes
  if(verbose){
    message("- Removing poor pixels with <= ",
            min.lib.size, " reads", "\n")
    message("- Removing genes with <= ",
            min.reads,
            " reads across pixels and detected in <= ",
            min.detected, " pixels", "\n")
  }
  countsClean <- cleanCounts(counts=t(counts),
                             min.reads=min.reads,
                             min.lib.size=min.lib.size,
                             min.detected=min.detected,
                             plot=plot,
                             verbose=FALSE)
  countsClean <- as.matrix(countsClean)
  if(verbose){
    message("  Remaining genes: ", dim(countsClean)[1],
            " and remaining pixels: ", dim(countsClean)[2], "\n")
  }
  
  # remove top expressed genes (nTopGenes needs to be integer or NA)
  if (!is.na(nTopGenes) & is.numeric(nTopGenes)) {
    nTopGenes <- round(nTopGenes)
    if(verbose){
      message("- Removing the top ", nTopGenes, " expressed genes.", "\n")
    }
    top_expressed <- names(Matrix::rowSums(countsClean)[order(Matrix::rowSums(countsClean),
                                                              decreasing = TRUE)][1:nTopGenes])
    countsClean <- countsClean[rownames(countsClean) %in% top_expressed == FALSE,]
  }
  
  # remove specific genes (if there are any). Use grepl to search for gene name pattern matches
  if (!is.na(genes.to.remove)) {
    countsClean <- countsClean[!grepl(paste(genes.to.remove, collapse="|"), rownames(countsClean)),]
    if(verbose){
      message("- After filtering for `genes.to.remove`:", "\n",
              " Remaining genes: ", dim(countsClean)[1], "\n")
    }
  }
  
  # remove genes that appear in more than certain percentage of the pixels
  if (!is.na(removeAbove) & is.numeric(removeAbove)){
    if (removeAbove >= 0 & removeAbove <= 1){
      # matrix where all positive gene counts are set to 1
      countsClean_ <- countsClean
      countsClean_[which(countsClean_ > 0)] <- 1
      # convert the percentage to number of pixels
      numberSpots <- (removeAbove * ncol(countsClean_))
      # rowSums of countsClean_ equate to number of pixels where each gene is expressed
      # remove genes expressed in "numberSpots" or more pixels
      countsClean <- countsClean[which(Matrix::rowSums(countsClean_) < numberSpots),]
      if(verbose){
        message("- Removed genes present in ",
                as.character(removeAbove*100), "% or more of pixels", "\n",
                " Remaining genes: ", dim(countsClean)[1], "\n")
      }
    } else {
      warning("Warning: `removeAbove` must be a numeric from 0 to 1. Skipping `removeAbove` gene filter.", "\n")
    }
  }
  # remove genes that appear in less than certain percentage of the pixels
  if (!is.na(removeBelow) & is.numeric(removeBelow)){
    if (removeBelow >= 0 & removeBelow <= 1){
      # matrix where all positive gene counts are set to 1
      countsClean_ <- countsClean
      countsClean_[which(countsClean_ > 0)] <- 1
      # convert the percentage to number of pixels
      numberSpots <- (removeBelow * ncol(countsClean_))
      # rowSums of countsClean_ equate to number of pixels where each gene is expressed
      # remove genes expressed in "numberSpots" or less pixels
      countsClean <- countsClean[which(Matrix::rowSums(countsClean_) > numberSpots),]
      if(verbose){
        message("- Removed genes present in ",
                as.character(removeBelow*100), "% or less of pixels", "\n",
                " Remaining genes: ", dim(countsClean)[1], "\n")
      }
    } else {
      warning("Warning: `removeBelow` must be a numeric from 0 to 1. Skipping `removeBelow` gene filter.", "\n")
    }
  }
  
  # use overdispersed variable genes for corpus
  if (ODgenes) {
    if(verbose){
      message("- Capturing only the overdispersed genes...", "\n")
    }
    par(mfrow=c(4,2), mar=c(1,1,1,1))
    OD <- getOverdispersedGenes(countsClean,
                                alpha=od.genes.alpha,
                                gam.k=gam.k,
                                plot=plot,
                                details=TRUE)
    
    # option to select just the top n overdispersed genes based on
    # log p-val adjusted
    if (!is.na(nTopOD)){
      if(verbose){
        message("- Using top ", nTopOD, " overdispersed genes.", "\n")
      }
      OD_filt <- OD$df[OD$ods,]
      # check if actual number of OD genes less than `nTopOD`
      if (dim(OD_filt)[1] < nTopOD){
        if(verbose){
          message(" number of top overdispersed genes available: ",
                  dim(OD_filt)[1], "\n")
        }
        od_genes <- rownames(OD_filt)
      } else {
        od_genes <- rownames(OD_filt[order(OD_filt$lpa),][1:nTopOD,])
      }
    } else {
      od_genes <- OD$ods
    }
    
    countsClean <- countsClean[od_genes,]
  }
  
  if (dim(countsClean)[1] > 1000){
    message("Genes in corpus > 1000 (", dim(countsClean)[1],
            "). This may cause model fitting to take a while. Consider reducing the number of genes.", "\n")
  }
  
  corpus <- t(as.matrix(countsClean))
  
  # last filter: each row must have at least 1 non-zero entry
  # to be compatible with `topicmodels`.
  if(verbose){
    message("- Check that each pixel has at least 1 non-zero gene count entry..", "\n")
  }
  corpus <- corpus[which(!Matrix::rowSums(corpus) == 0),]
  corpus_slm <- slam::as.simple_triplet_matrix(corpus)
  message("Final corpus:", "\n")
  print(corpus_slm)
  
  # get pixel positions if pixel colnames contain the positions
  # this is the case for some ST datasets like Stahl 2016 sets
  # ex: "20x30" -> 20 x, 30 y positions
  if (extractPos) {
    if(verbose){
      message("Extracting positions from pixel names.", "\n")
    }
    positions <- do.call(rbind, lapply(rownames(corpus), function(spotID) {
      coords <- as.numeric(strsplit(spotID, "x")[[1]])
      coords
    }))
    colnames(positions) <- c("x", "y")
    rownames(positions) <- rownames(corpus)
  } else {
    positions <- NULL
  }
  
  message("Preprocess complete.", "\n")
  return(list(corpus = corpus,
              slm = corpus_slm,
              pos = positions))
}


#' Find the optimal number of cell-types K for the LDA model
#'
#' @description The input for topicmodels::LDA needs to be a
#'     slam::as.simple_triplet_matrix (docs x words). Access a given model in
#'     the returned list via: lda$models$k. The models are objects from
#'     the R package "topicmodels". The LDA models have slots with additional
#'     information.
#'
#' @param counts Gene expression counts with pixels as rows and genes as columns
#' @param Ks vector of K parameters, or number of cell-types, to fit models with
#' @param seed Random seed
#' @param perc.rare.thresh the number of deconvolved cell-types with mean pixel proportion below this fraction used to assess
#'     performance of fitted models for each K. Recorded for each K. (default: 0.05)
#' @param ncores Number of cores for parallelization (default: 1). Suggest: parallel::detectCores()
#' @param plot Boolean for plotting (default: TRUE)
#' @param verbose Boolean for verbosity (default: TRUE)
#'
#' @return A list that contains
#' \itemize{
#' \item models: each fitted LDA model for a given K
#' \item kneedOptK: the optimal K based on Kneed algorithm
#' \item minOptK: the optimal K based on minimum
#' \item ctPropOptK: Suggested upper bound on K. K in which number of returned cell-types
#'     with mean proportion < perc.rare.thresh starts to increase steadily.
#' \item numRare: number of cell-types with mean pixel proportion < perc.rare.thresh for each K
#' \item perplexities: perplexity scores for each model
#' \item fitCorpus: the corpus that was used to fit each model
#' \item testCorpus: the corpus used to compute model perplexity.
#' }
#' 
#' @examples 
#' data(mOB)
#' pos <- mOB$pos
#' cd <- mOB$counts
#' counts <- cleanCounts(cd, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = 3, ncores=2)
#' 
#' @importFrom methods slot
#' 
#' @export
fitLDA <- function(counts,
                   Ks=seq(2, 10, by=2),
                   seed=0,
                   perc.rare.thresh=0.05,
                   ncores=1,
                   plot=TRUE,
                   verbose=TRUE) {
  
  if ( min(Ks) < 2 | !isTRUE(all(Ks == floor(Ks))) ){
    stop("`Ks` must be a vector of integers greater than 2.")
  }
  
  counts <- as.matrix(counts)
  
  if (length(which(rowSums(counts) == 0))){
    stop("Each row (pixel) of `counts` needs to contain at least one non-zero entry")
  }
  
  if ( !isTRUE(all(counts == floor(counts))) ){
    stop("`counts` must contain integer gene counts")
  }
  
  testingPixels <- seq(nrow(counts))
  fittingPixels <- seq(nrow(counts))
  
  ## counts must be pixels (rows) x genes (cols) matrix
  corpus <- slam::as.simple_triplet_matrix((as.matrix(counts)))
  corpusFit <- slam::as.simple_triplet_matrix((as.matrix(counts[fittingPixels,])))
  corpusTest <- slam::as.simple_triplet_matrix((as.matrix(counts[testingPixels,])))
  
  verbose <- as.integer(verbose)
  
  controls <- list(seed=seed,
                   verbose=0, keep=1, estimate.alpha=TRUE)
  
  start_time <- Sys.time()
  fitted_models <- BiocParallel::bplapply(Ks, function(k){
    if(verbose) system(paste("echo 'now fitting LDA model with K =", k, "'"))
    # if(verbose) print(paste("now fitting LDA model with K =", k))
    topicmodels::LDA(corpusFit, k=k, control=controls)
  }, BPPARAM=BiocParallel::SnowParam(workers=ncores))
  names(fitted_models) <- Ks
  
  if(verbose) {
    total_t <- round(difftime(Sys.time(), start_time, units="mins"), 2)
    message(sprintf("Time to fit LDA models was %s mins", total_t))
  }
  
  if(verbose) {
    message("Computing perplexity for each fitted model...")
  }
  
  start_time <- Sys.time()
  pScores <- BiocParallel::bplapply(Ks, function(k){
    if(verbose) system(paste("echo 'computing perplexity for LDA model with K =", k, "'"))
    # if(verbose) print(paste("computing perplexity for LDA model with K =", k))
    model <- fitted_models[[as.character(k)]]
    ## if no splitting, then same corpus and thus no need to refit in order to save time
    topicmodels::perplexity(model)
    
  }, BPPARAM=BiocParallel::SnowParam(workers=ncores))
  pScores <- unlist(pScores)
  
  if(verbose) {
    total_t <- round(difftime(Sys.time(), start_time, units="mins"), 2)
    message(sprintf("Time to compute perplexities was %s mins", total_t))
  }
  
  ## Kneed algorithm
  kOpt1 <- Ks[where.is.knee(pScores)]
  ## Min
  kOpt2 <- Ks[which(pScores == min(pScores))]
  
  ## check number of predicted cell-types at low proportions
  if(verbose) {
    message("Getting predicted cell-types at low proportions...")
  }
  
  start_time <- Sys.time()
  out <- lapply(1:length(Ks), function(i) {
    ## if no splitting, then already fit to entire corpus and no need to refit
    
    ## note that we are including all cell-types when computing theta here because we are trying to
    ## assess the best model by the number of rare cell-types predicted.
    apply(getBetaTheta(fitted_models[[i]], corpus=NULL, perc.filt=0,
                       verbose=FALSE)$theta, 2, mean)
    
  })
  ## number of cell-types present at fewer than `perc.rare.thresh` on average across pixels
  numrare <- unlist(lapply(out, function(x) sum(x < perc.rare.thresh)))
  kOpt3 <- Ks[where.is.knee(numrare)]
  
  if(verbose) {
    total_t <- round(difftime(Sys.time(), start_time, units="mins"), 2)
    message(sprintf("Time to compute cell-types at low proportions was %s mins",
                    total_t))
  }
  
  if(plot) {
    
    if(verbose) {
      message("Plotting...")
    }
    
    dat <- data.frame(K=as.double(Ks),
                      rareCts=numrare,
                      perplexity=pScores,
                      rareCtsAdj=scale0_1(numrare),
                      perplexAdj=scale0_1(pScores),
                      alphas=unlist(sapply(fitted_models, slot, "alpha")))
    dat[["alpha < 1"]] <- ifelse(dat$alphas < 1, 'gray90', 'gray50')
    dat$alphaBool <- ifelse(dat$alphas < 1, 0, 1)

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
      sec_ax_labs <- seq(min(dat$perplexity), max(dat$perplexity),
                         (max(dat$perplexity)-min(dat$perplexity))/1)
    } else {
      sec_ax_labs <- seq(min(dat$perplexity), max(dat$perplexity),
                         (max(dat$perplexity)-min(dat$perplexity))/max(dat$rareCts))
    }
    sec_ax_breaks <- scale0_1(sec_ax_labs)

    if(length(sec_ax_labs) == 1){
      sec_ax_breaks <- 0
      dat$perplexAdj <- 0
    }

    plt <- ggplot2::ggplot(dat) +
      ggplot2::geom_point(ggplot2::aes(y=rareCtsAdj, x=K), col="blue", lwd=2) +
      ggplot2::geom_point(ggplot2::aes(y=perplexAdj, x=K), col="red", lwd=2) +
      ggplot2::geom_line(ggplot2::aes(y=rareCtsAdj, x=K), col="blue", lwd=2) +
      ggplot2::geom_line(ggplot2::aes(y=perplexAdj, x=K), col="red", lwd=2) +
      ggplot2::geom_bar(ggplot2::aes(x=K, y=alphaBool), fill=dat$`alpha < 1`,
                        stat="identity", width=1, alpha=0.5) +
      ggplot2::scale_y_continuous(name=paste0("# cell-types with mean proportion < ",
                                              round(perc.rare.thresh*100, 2), "%"),
                                  breaks=prim_ax_breaks,
                                  labels=prim_ax_labs,
                                  sec.axis=ggplot2::sec_axis(~ ., name="perplexity",
                                                             breaks=sec_ax_breaks,
                                                             labels=round(sec_ax_labs, 2))) +
      ggplot2::scale_x_continuous(breaks=min(dat$K):max(dat$K)) +
      ggplot2::labs(title="Fitted model K's vs deconvolved cell-types and perplexity",
                    subtitle="LDA models with alpha > 1 shaded") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        panel.background=ggplot2::element_blank(),
        plot.title=ggplot2::element_text(size=15, face=NULL),
        # plot.subtitle = ggplot2::element_text(size=13, face=NULL),
        panel.grid.minor=ggplot2::element_blank(),
        panel.grid.major=ggplot2::element_line(color="black", size=0.1),
        panel.ontop=TRUE,
        axis.title.y.left=ggplot2::element_text(color="blue", size=13),
        axis.text.y.left=ggplot2::element_text(color="blue", size=13),
        axis.title.y.right=ggplot2::element_text(color="red", size=15, vjust=1.5),
        axis.text.y.right=ggplot2::element_text(color="red", size=13),
        axis.text.x=ggplot2::element_text(angle=0, size=13),
        axis.title.x=ggplot2::element_text(size=13)
      )
    print(plt)
    
  }
  
  return(list(models=fitted_models,
              kneedOptK=kOpt1,
              minOptK=kOpt2,
              ctPropOptK=kOpt3,
              numRare=numrare,
              perplexities=pScores,
              fitCorpus=corpusFit,
              testCorpus=corpusTest))
}


#' Pull out cell-type proportions across pixels (theta) and
#' cell-type gene probabilities (beta) matrices from fitted LDA models from fitLDA
#'
#' @param lda an LDA model from "topicmodels" R package. From list of models returned by
#'     fitLDA
#' @param corpus If corpus is NULL, then it will use the original corpus that
#'     the model was fitted to. Otherwise, compute deconvolved topics from this
#'     new corpus. Needs to be pixels x genes and nonnegative integer counts. 
#'     Each row needs at least 1 nonzero entry (default: NULL)
#' @param perc.filt proportion threshold to remove cell-types in pixels (default: 0.05)
#' @param betaScale factor to scale the predicted cell-type gene expression profiles (default: 1)
#' @param verbose Boolean for verbosity (default: TRUE)
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
#' @examples 
#' data(mOB)
#' pos <- mOB$pos
#' cd <- mOB$counts
#' counts <- cleanCounts(cd, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = 3, ncores=2)
#' optLDA <- optimalModel(models = ldas, opt = 3)
#' results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
#' head(results$theta)
#' head(results$beta)
#' 
#' @export
getBetaTheta <- function(lda, corpus=NULL, perc.filt=0.05, betaScale=1, verbose=TRUE) {
  
  ## If corpus is NULL, then it will use the original fitted data
  ## otherwise, use new corpus to predict topic proportions based on the 
  ## model. pixels x genes, with nonzero integer counts. Also converted
  ## to slam matrix format to be compatible with topicmodels package.
  if(is.null(corpus)){
    result <- topicmodels::posterior(lda)
  } else {
    
    corpus <- as.matrix(corpus)
    
    if (length(which(rowSums(corpus) == 0))){
      stop("Each row (pixel) of `corpus` needs to contain at least one non-zero entry")
    }
    if ( !isTRUE(all(corpus == floor(corpus))) ){
      stop("`corpus` must contain integer gene counts")
    }
    
    corpus <- slam::as.simple_triplet_matrix((as.matrix(corpus)))
    result <- topicmodels::posterior(lda, newdata = corpus)
  }
  
  theta <- result$topics
  beta <- result$terms
  
  ## filter out cell-types with low proportions in pixels
  if(verbose){
    message("Filtering out cell-types in pixels that contribute less than ",
            perc.filt, " of the pixel proportion.", "\n")
  }
  theta <- filterTheta(theta, perc.filt=perc.filt, verbose=verbose)
  
  ## scale the beta
  beta <- beta * betaScale
    
  return(list(beta=beta,
              theta=theta))
}


#' Find Pearson's correlations between topics (cell-types) with respect to their
#' proportions across documents (pixels), i.e. thetas, or gene probabilities,
#' i.e. betas.
#'
#' @param m1 first matrix
#' @param m2 second matrix
#' @param type must be either "t" (theta; cell-type proportions across pixels) or "b" (beta; cell-type gene expression profiles)
#' @param thresh if comparing betas, use to compare genes above this probability (e.g., expression level).
#'               NULL or 0 < numeric < 1.0 (default: NULL)
#' @param verbose control the verbosity (default: TRUE)
#'
#' @return matrix of Pearson's correlations; m1 (rows) by m2 (cols)
#' 
#' @examples 
#' data(mOB)
#' pos <- mOB$pos
#' cd <- mOB$counts
#' counts <- cleanCounts(cd, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = 3)
#' optLDA <- optimalModel(models = ldas, opt = 3)
#' results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
#' deconProp <- results$theta
#' corMtx <- getCorrMtx(m1 = as.matrix(deconProp), m2 = as.matrix(deconProp), type = "t")
#' rownames(corMtx) <- paste0("X", seq(nrow(corMtx)))
#' colnames(corMtx) <- paste0("X", seq(ncol(corMtx)))
#' head(corMtx)
#' 
#' @importFrom stats cor
#' 
#' @export
getCorrMtx <- function(m1, m2, type=c("t", "b"), thresh=NULL, verbose=TRUE) {
  
  type <- match.arg(type)
  
  if (!is.matrix(m1) | !is.matrix(m2)){
    stop("`m1` and `m2` must be matrices")
  }
  
  if(type == "t"){
    message("NOTE: using type='t' and comparing thetas where the cell-types are
            the columns (pixels x cell-types)")
  }
  if(type == "b"){
    message("NOTE: using type='b' and comparing betas where the cell-types are
            the rows (cell-types x genes)")
  }
  
  
  # if comparing thetas the cell-types are the columns (pixels x cell-types)
  if (type == "t"){
    
    # make sure the same pixels are being compared
    keep_spots <- intersect(rownames(m1), rownames(m2))
    if(verbose){
      message("cell-type correlations based on ",
              length(keep_spots),
              " shared pixels between m1 and m2.", "\n")
    }
    corMtx <- do.call(rbind, lapply(1:ncol(m1), function(i) {
      sapply(1:ncol(m2), function(j) {
        cor(m1[keep_spots,i], m2[keep_spots,j])
      })
    }))
    rownames(corMtx) <- colnames(m1)
    colnames(corMtx) <- colnames(m2)
    return(corMtx)
  }
  
  # if comparing betas the cell-types are the rows (cell-types x genes)
  if (type == "b"){
    
    # make sure the same genes are being compared
    keep_genes <- intersect(colnames(m1), colnames(m2))
    if(verbose){
      message("cell-type correlations based on ",
              length(keep_genes),
              " shared genes between m1 and m2.", "\n")
    }
    
    if (is.numeric(thresh)){
      if(verbose){
        message("comparing genes with cell-type probability > ", thresh, "\n")
      }
    }
    
    corMtx <- do.call(rbind, lapply(1:nrow(m1), function(i) {
      
      # if choosing top genes for a cell-type using thresh
      if (is.null(thresh)){
        genes <- keep_genes
      } else if (is.numeric(thresh) & thresh > 0 & thresh < 1){
        m1genes <- keep_genes[which(m1[i,keep_genes] > thresh)]
      } else {
        stop("thresh must be NULL or numeric > 0 and < 1")
      }
      
      sapply(1:nrow(m2), function(j) {
        if (is.null(thresh)){
          genes <- keep_genes
        } else {
          # at this point thresh should be between 0 and 1 and m1genes selected
          m2genes <- keep_genes[which(m2[j,keep_genes] > thresh)]
          genes <- intersect(m1genes, m2genes)
          if (length(genes) < 2){
            warning(cat("WARNING: 0 or 1 shared genes >", thresh,
                        "for both m1 cell-type", i, "and m2 cell-type", j,
                        "; corr will be NA.", "\n"))
          }
        }
        cor(m1[i,genes], m2[j,genes])
      })
    }))
    rownames(corMtx) <- rownames(m1)
    colnames(corMtx) <- rownames(m2)
    return(corMtx)
  }
}


#' Match deconvolved cell-types to ground truth cell-types based on transcriptional profiles
#'
#' @description Match deconvolved cell-types to ground truth cell-types by testing for
#'     enrichment of ground truth marker gene sets in the deconvolved transcriptional profiles.
#'     Uses fgsea Korotkevich G, Sukhov V, Sergushichev A (2019). “Fast gene set enrichment analysis.” bioRxiv. doi:10.1101/060012, http://biorxiv.org/content/early/2016/06/20/060012.
#'
#' @param beta cell-type (rows) x gene (columns) matrix of deconvolved cell-type transcriptional profiles
#' @param gset named list where each entry is a vector of marker genes for a given ground truth cell-type.
#' @param qval adjusted p-value threshold (default: 0.05)
#' @param ... additional options for fgsea
#'
#' @return A list that contains
#' \itemize{
#' \item results: A named list that contains sorted matrices for each deconvolved cell-type.
#'     The matrix rows are the ground truth cell-types ordered by significance, edge-score, and enrichment score
#'     of their gene sets in the deconvolved transcriptional profile of a given deconvolved cell-type.
#' \item predictions: a named vector where the names are the deconvolved cell-types and the values
#'     are the best matched ground truth cell-type that is also positively enriched.
#' }
#'
#' @export
annotateCellTypesGSEA <- function(beta, gset, qval=0.05, ...) {
  
  results <- list()
  top.pos.enrich <- c()
  
  for (i in seq(nrow(beta))){
    celltype <- i
    vals <- sort(beta[celltype,], decreasing=TRUE)
    
    gsea.results <- fgsea::fgsea(stats=vals, pathways=gset, scoreType = "pos", ...)
    
    # filter for top hits
    gsea.sig <- gsea.results[gsea.results$padj < qval,]
    
    ## order of selection:
    ## 1. adjusted p-value
    ## 2. leadingEdge (the leading edge subset of a gene set is the subset of genes that contribute most to the Enrichment Score)
    ## 3. ES (Enrichment Score)
    if(nrow(gsea.sig) > 1) {
      gsea.sig <- gsea.sig[order(gsea.sig$padj, rev(gsea.sig$leadingEdge), rev(gsea.sig$ES)), ]
    }
    
    results[[ rownames(beta)[celltype] ]] <- gsea.sig
    
    ## the top entry that is also positiviely enriched in the txn profile is predicted to be the best matching
    gsea.sig.pos <- gsea.sig[which(gsea.sig$ES > 0), ]
    top.pos.enrich <- append(top.pos.enrich, gsea.sig.pos[1]$pathway)
    
  }
  
  #names(top.pos.enrich) <- rownames(beta)
  
  return(list(results=results,
              predictions=top.pos.enrich))
  
}


#' Get the optimal LDA model
#'
#' @param models list returned from fitLDA
#' @param opt either "kneed" (kOpt1) or "min" (kOpt2), or designate a specific K.
#'     "kneed" = K vs perplexity inflection point.
#'     "min" = K corresponding to minimum perplexity
#'     "proportion" = K vs number of cell-type with mean proportion < 5% inflection point
#'
#' @return optimal LDA model fitted to the K based on opt parameter
#' 
#' @examples 
#' data(mOB)
#' pos <- mOB$pos
#' cd <- mOB$counts
#' counts <- cleanCounts(cd, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2,4), ncores=2)
#' optLDA <- optimalModel(models = ldas, opt = "min")
#' 
#' @importFrom methods slot
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
    stop("`opt` must be either 'kneed',  'min', 'proportion', or integer of a fitted K")
  }
  return(m)
}


#' Function to get Hungarian sort pairs via clue::lsat
#'
#' @description Finds best matches between cell-types that correlate between
#'     beta or theta matrices that have been compared via getCorrMtx().
#'     Each row is paired with a column in the output matrix from getCorrMtx().
#'     If there are less rows than columns, then some columns will not be
#'     matched and not part of the output.
#'
#' @param mtx output correlation matrix from getCorrMtx(). Must not have more rows
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
#' @examples 
#' data(mOB)
#' pos <- mOB$pos
#' cd <- mOB$counts
#' counts <- cleanCounts(cd, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = 3)
#' optLDA <- optimalModel(models = ldas, opt = 3)
#' results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
#' deconProp <- results$theta
#' corMtx <- getCorrMtx(m1 = as.matrix(deconProp), m2 = as.matrix(deconProp), type = "t")
#' pairs <- lsatPairs(corMtx)
#' pairs
#'
#' @export
lsatPairs <- function(mtx){
  # must have equal or more rows than columns
  # values in matrix converted to 0-1 scale relative to all values in mtx
  pairing <- clue::solve_LSAP(scale0_1(mtx), maximum=TRUE)
  # clue::lsat returns vector where for each position the first element is a row
  # and the second is the paired column
  rowsix <- seq_along(pairing)
  colsix <- as.numeric(pairing)
  
  return(list(pairs=pairing,
              rowix=rowsix,
              colsix=colsix))
}


#' Function to filter out cell-types in pixels below a certain proportion
#' 
#' @description Sets cell-types in each pixel to 0 that are below a given proportion.
#'     Then renormalizes the pixel proportions to sum to 1.
#'     Cell-types that result in 0 in all pixels after this filtering step are removed.
#'
#' @param theta pixel (rows) by cell-types (columns) distribution matrix. Each row
#'     is the cell-type composition for a given pixel
#' @param perc.filt proportion threshold to remove cell-types in pixels (default: 0.05)
#' @param verbose Boolean for verbosity (default: TRUE)
#' 
#' @return A filtered pixel (rows) by cell-types (columns) distribution matrix.
#' 
#' @noRd
filterTheta <- function(theta, perc.filt=0.05, verbose=TRUE){
  ## remove rare cell-types in pixels
  theta[theta < perc.filt] <- 0
  ## re-adjust pixel proportions to 1
  theta <- theta/rowSums(theta)
  ## if NAs because all cts are 0 in a spot, replace with 0
  theta[is.na(theta)] <- 0
  ## drop any cts that are 0 for all pixels
  dropped_cts <- names(which(colSums(theta) == 0))
  if(length(dropped_cts) > 0){
    if(verbose){
      message("Cell-types with no proportions in pixels after filtering dropped:\n",
              dropped_cts, "\n")
    }
  }
  theta <- theta[,which(colSums(theta) > 0)]
  
  empty_pixels <- names(which(rowSums(theta) == 0))
  if(length(empty_pixels) > 0){
    if(verbose){
      message(length(empty_pixels), " pixels with no cell-types after filtering.", "\n")
    }
  }
  
  return(theta)
}


#' Returns top n genes of each deconvolved cell-type for a given beta matrix
#'
#' @description For a given beta matrix (cell-type gene distribution matrix),
#'     returns the top n genes based on their probability.
#' 
#' @param beta beta matrix (cell-type gene distribution matrix)
#' @param n number of top genes for each deconvolved cell-type to return (default: 10)
#' 
#' @return a list where each item is a vector of the top genes and their associated probabilities for
#'     a given deconvolved cell-type
#'     
#' @examples 
#' data(mOB)
#' pos <- mOB$pos
#' cd <- mOB$counts
#' counts <- cleanCounts(cd, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = 3, ncores=2)
#' optLDA <- optimalModel(models = ldas, opt = 3)
#' results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
#' deconGexp <- results$beta
#' genes <- topGenes(deconGexp)
#' 
#' @export
topGenes <- function(beta, n=10){
  topgenes <- lapply(seq(nrow(beta)), function(ct){
    sort(beta[ct,], decreasing=TRUE)[1:n]
  })
  names(topgenes) <- rownames(beta)
  return(topgenes)
}


#' Helper function to scale values to 0-1 range relative to each other. For use
#' with lsatPairs()
#'
#' @param x vector or matrix
#' @return vector or matrix with all values adjusted 0-1 scale relative to each other.
#' 
#' @noRd
scale0_1 <- function(x) {
  xAdj <- (x - min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE))
  return(xAdj)
}

