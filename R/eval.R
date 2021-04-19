#' Find Pearson correlations between topics (cell types) with respect to their
#' proportions across documents (spots), i.e. thetas, or gene probabilities,
#' i.e. betas.
#' 
#' @param m1 first matrix
#' @param m2 second matrix
#' @param type must be either "t" (theta) or "b" (beta)
#' @param thresh if comparing betas, use to compare genes above this probability.
#'               NULL or 0 < numeric < 1.0 (default: NULL)
#' 
#' @return matrix of pearson correlations between topics m1 (rows) by topics m2 (cols)
#' @export
getCorrMtx <- function(m1, m2, type, thresh = NULL) {
  
  if (is.matrix(m1) == FALSE | is.matrix(m2) == FALSE){
    stop("`m1` and `m2` must be matrices")
  }
  if (!type %in% c("t", "b")){
    stop("`type` must be either 't' or 'b'")
  }
  
  # if comparing thetas the topics are the columns (spots x topics)
  if (type == "t"){
    
    # make sure the same spots are being compared
    keep_spots <- intersect(rownames(m1), rownames(m2))
    cat("topic correlations based on", length(keep_spots), "shared spots between m1 and m2.", "\n")
    
    corMtx <- do.call(rbind, lapply(1:ncol(m1), function(i) {
      sapply(1:ncol(m2), function(j) {
      cor(m1[keep_spots,i], m2[keep_spots,j])
      })
    }))
    rownames(corMtx) <- colnames(m1)
    colnames(corMtx) <- colnames(m2)
    return(corMtx)
    }
  
  # if comparing betas the topics are the rows (topics x genes)
  if (type == "b"){
    
    # make sure the same genes are being compared
    keep_genes <- intersect(colnames(m1), colnames(m2))
    cat("topic correlations based on", length(keep_genes), "shared genes between m1 and m2.", "\n")
    
    if (is.numeric(thresh)){
      cat("comparing genes with topic probability >", thresh, "\n")
    }
    
    corMtx <- do.call(rbind, lapply(1:nrow(m1), function(i) {
      
      # if choosing top genes for a topic using thresh
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
          # cat(length(genes), "genes compared for m1 topic", i, "and m2 topic", j, "\n")
          if (length(genes) < 2){
            warning(cat("WARNING: 0 or 1 shared genes >", thresh, "for both m1 topic", i, "and m2 topic", j, "; corr will be NA.", "\n"))
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


#' Helper function to scale values to 0-1 range relative to each other. For use
#' with `lsatPairs`
#' 
#' @param x vector or matrix
#' @return vector or matrix with values adjusted to 0-1 scale
#' @export
scale0_1 <- function(x) {
  xAdj <- (x - min(x, na.rm = T)) / diff(range(x, na.rm = T))
  return(xAdj)
}


#' Pre-process ST count matrices to construct corpus for input into LDA
#' 
#' @description Takes spot (row) x gene (columns) mtx and filters out poor genes
#'              and spots. Then selects for genes to be included in corpus.
#'              If the spot IDs are made up of their positions in "XxY" these
#'              can be extracted as the spot position coordinates.
#'              
#'              Order of filtering options:
#'              1. MERINGUE::cleanCounts to remove poor spots and genes
#'              2. Selection to use specific genes only
#'              3. Remove top expressed genes in matrix
#'              4. Remove specific genes based on grepl pattern matching
#'              5. Remove genes that appear in more than a percentage of spots
#'              6. Use the overdispersed genes computed from the remaining genes
#'                 after filtering steps 1-5 (if selected)
#'              
#' @param dat spot (row) x gene (columns) mtx with gene counts OR path to it
#' @param alignFile path to 3x3 alignment file to adjust spot coordinates
#'     (optional).
#'    Stahl 2016 datasets (default: NA)
#' @param extractPos Boolean to extract spot positional coordinates from spot
#'     names
#'    (default: FALSE)
#' @param selected.genes vector of gene names to use specifically for the corpus
#'    (default: NA)
#' @param nTopGenes integer for number of top expressed genes to remove
#'     (default: NA)
#' @param genes.to.remove vector of gene names or patterns for matching to genes
#'     to remove (default: NA). ex: c("^mt-") or c("^MT", "^RPL", "^MRPL")
#' @param perc.spots non-negative numeric <=1 to use as a percentage.
#'    Removes genes present in this fraction of spots (default: NA)
#' @param min.reads MERINGUE::cleanCounts param; minimum number of reads to keep
#'     a gene (default: 100)
#' @param min.lib.size MERINGUE::cleanCounts param; minimum number of counts a
#'     spot needs to keep (default: 100)
#' @param min.detected MERIGNUE::cleanCounts param; minimum number of spots a gene
#'     needs to have been detected in to keep (default: 1)
#' @param ODgenes Boolean to use MERINGUE::getOverdispersedGenes for the corpus
#'    genes (default: TRUE)
#' @param nTopOD number of top OD genes to use. int (default: NA)
#' @param od.genes.alpha alpha parameter for MERINGUE::getOverdispersedGenes.
#'     Higher = less stringent and more OD genes returned (default: 0.05)
#' @param gam.k gam.k parameter for MERINGUE::getOverdispersedGenes. Dimension
#'     of the "basis" functions in the GAM used to fit, higher = "smoother"
#'     (default: 5)
#' 
#' @return A list that contains
#' \itemize{
#' \item corpus: (spots x genes) matrix of the counts of the selected genes 
#' \item slm: slam::as.simple_triplet_matrix(corpus); required format for topicmodels::LDA input
#' \item positions: matrix of x and y coordinates of spots. rownames = spots, colnames = "x", "y"
#' }
#' 
#' @export
preprocess <- function(dat,
                       alignFile = NA,
                       extractPos = FALSE,
                       selected.genes = NA,
                       nTopGenes = NA,
                       genes.to.remove = NA,
                       perc.spots = NA,
                       min.reads = 100,
                       min.lib.size = 100,
                       min.detected = 1,
                       ODgenes = TRUE,
                       nTopOD = NA,
                       od.genes.alpha = 0.05,
                       gam.k = 5) {
  
  if (typeof(dat) == "character") {
    if (file.exists(dat) == TRUE){
      counts <- read.table(dat)
    } else {
      cat("path to file does not exist", "\n")
    }
  } else if (is.matrix(dat) == TRUE){
    counts <- dat
  } else {
    stop("`dat` is not a viable path or matrix")
  }
  
  cat("Initial genes:", dim(t(counts))[1], "Initial spots:", dim(t(counts))[2], "\n")
  
  # use specific genes in the corpus
  if (is.na(selected.genes[1]) == FALSE) {
    cat("- Using genes in `selected.genes` for corpus.", "\n")
    counts <- counts[,selected.genes]
    cat(" ", dim(counts)[2], "genes are present in dataset.", "\n")
  }
  
  # remove poor spots and genes
  cat("- Removing poor spots with <=", min.lib.size, "reads", "\n")
  cat("- Removing genes with <=", min.reads, "reads across spots and detected in <=", min.detected, "spots.", "\n")
  countsClean <- cleanCounts(counts = t(counts), 
                             min.reads = min.reads, 
                             min.lib.size = min.lib.size, 
                             min.detected = min.detected,
                             plot=TRUE,
                             verbose=FALSE)
  
  cat("  Remaining genes:", dim(countsClean)[1], "and remaining spots:", dim(countsClean)[2], "\n")
  
  # adjust spot coordinates, optional.
  # based on Stahl 2016 ST data with alignment matrices
  if (is.na(alignFile) == FALSE) {
    if (file.exists(alignFile) == TRUE){
      cat("- Adjusting spot positions based on alignment file.")
      align <- matrix(unlist(read.table(alignFile)), nrow = 3, ncol = 3)
      (positions[,"x"] * align[1,1]) - 290 # note that I found they were off by one spot distance in pixels
      (positions[,"y"] * align[2,2]) - 290
    } else {
      cat("Warning: `alignFile` path does not exists. Skipping position adjustments.", "\n")
    }
  }
  
  # remove top expressed genes (nTopGenes needs to be integer or NA)
  if (is.na(nTopGenes) == FALSE & is.numeric(nTopGenes) == TRUE) {
    nTopGenes <- round(nTopGenes)
    cat("- Removing the top", nTopGenes, "expressed genes.", "\n")
    top_expressed <- names(rowSums(countsClean)[order(rowSums(countsClean),
                                                      decreasing = TRUE)][1:nTopGenes])
    countsClean <- countsClean[rownames(countsClean) %in% top_expressed == FALSE,]
    
    # print(paste("after removing top ",
    #             as.character(nTopGenes),
    #             " genes:", dim(countsClean)[1], "genes remain."))
  }
  
  # remove specific genes (if there are any). Use grepl to search for gene name pattern matches
  if (is.na(genes.to.remove) == FALSE) {
    countsClean <- countsClean[!grepl(paste(genes.to.remove, collapse="|"), rownames(countsClean)),]
    cat("- After filtering for `genes.to.remove`:", "\n",
        " Remaining genes:", dim(countsClean)[1], "\n")
  }
  
  # remove genes that appear in certain percentage of the spots
  if (is.na(perc.spots) == FALSE & is.numeric(perc.spots) == TRUE){
    if (perc.spots >= 0 & perc.spots <= 1){
      # matrix where all positive gene counts are set to 1
      countsClean_ <- countsClean
      countsClean_[which(countsClean_ > 0)] <- 1
      # convert the percentage to number of spots
      numberSpots <- (perc.spots * ncol(countsClean_))
      # rowSums of countsClean_ equate to number of spots where each gene is expressed
      # remove genes expressed in "numberSpots" or more spots
      countsClean <- countsClean[which(rowSums(countsClean_) < numberSpots),]
      cat("- Removed genes present in", 
                  as.character(perc.spots*100), "% or more of spots", "\n",
                  " Remaining genes:", dim(countsClean)[1], "\n")
    } else {
      cat("Warning: `perc.spots` must be a numeric from 0 to 1. Skipping perc.spots gene filter.", "\n")
    }
  }
  
  # use overdispersed variable genes for corpus
  if (ODgenes == TRUE) {
    cat("- Capturing only the overdispersed genes...", "\n")
    par(mfrow=c(4,2), mar=c(1,1,1,1))
    OD <- getOverdispersedGenes(countsClean,
                                alpha = od.genes.alpha,
                                gam.k = gam.k,
                                plot = TRUE,
                                details = TRUE)
    
    # option to select just the top n OD genes based on
    # log p-val adjusted
    if (is.na(nTopOD) == FALSE){
      cat("- Using top", nTopOD, "overdispersed genes.", "\n")
      OD_filt <- OD$df[OD$ods,]
      # check if actual number of OD genes less than `nTopOD`
      if (dim(OD_filt)[1] < nTopOD){
        cat(" number of top OD genes available:", dim(OD_filt)[1], "\n")
        od_genes <- rownames(OD_filt)
      } else {
        od_genes <- rownames(OD_filt[order(OD_filt$lpa),][1:nTopOD,])
      }
    } else {
      od_genes <- OD$ods
    }
    
    countsClean <- countsClean[od_genes,]
  }
  
  corpus <- t(as.matrix(countsClean))
  
  # last filter: each row must have at least 1 non-zero entry
  # to be compatible with `topicmodels`.
  cat("- Check that each spot has at least 1 non-zero gene count entry..", "\n")
  corpus <- corpus[which(!rowSums(corpus) == 0),]
  corpus_slm <- slam::as.simple_triplet_matrix(corpus)
  cat("Final corpus:", "\n")
  print(corpus_slm)
  
  # get spot positions if spot colnames contain the positions
  # this is the case for some ST datasets like Stahl 2016 sets
  # ex: "20x30" -> 20 x, 30 y positions
  if (extractPos) {
    cat("Extracting positions from spot names.", "\n")
    positions <- do.call(rbind, lapply(rownames(corpus), function(spotID) {
      coords <- as.numeric(strsplit(spotID, "x")[[1]])
      coords
    }))
    colnames(positions) <- c("x", "y")
    rownames(positions) <- rownames(corpus)
  } else {
    positions <- NULL
  }
  
  cat("Preprocess complete.", "\n")
  return(list(corpus = corpus,
              slm = corpus_slm,
              pos = positions))
}



nrmse_func <-  function(obs, pred, type = "sd") {
  
  squared_sums <- sum((obs - pred)^2)
  mse <- squared_sums/length(obs)
  rmse <- sqrt(mse)
  if (type == "sd") nrmse <- rmse/sd(obs)
  if (type == "mean") nrmse <- rmse/mean(obs)
  if (type == "maxmin") nrmse <- rmse/ (max(obs) - min(obs))
  if (type == "iq") nrmse <- rmse/ (quantile(obs, 0.75) - quantile(obs, 0.25))
  if (!type %in% c("mean", "sd", "maxmin", "iq")) message("Wrong type!")
  nrmse <- round(nrmse, 3)
  return(nrmse)
  
}


