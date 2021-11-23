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
#' @export
getCorrMtx <- function(m1, m2, type, thresh = NULL, verbose = TRUE) {
  
  if (is.matrix(m1) == FALSE | is.matrix(m2) == FALSE){
    stop("`m1` and `m2` must be matrices")
  }
  if (!type %in% c("t", "b")){
    stop("`type` must be either 't' or 'b'")
  }
  
  # if comparing thetas the cell-types are the columns (pixels x cell-types)
  if (type == "t"){
    
    # make sure the same pixels are being compared
    keep_spots <- intersect(rownames(m1), rownames(m2))
    if(verbose){
      cat("cell-type correlations based on", length(keep_spots), "shared pixels between m1 and m2.", "\n")
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
      cat("cell-type correlations based on", length(keep_genes), "shared genes between m1 and m2.", "\n")
    }
    
    if (is.numeric(thresh)){
      if(verbose){
        cat("comparing genes with cell-type probability >", thresh, "\n")
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
          # cat(length(genes), "genes compared for m1 topic", i, "and m2 topic", j, "\n")
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


#' Helper function to scale values to 0-1 range relative to each other. For use
#' with `lsatPairs`
#'
#' @param x vector or matrix
#' @return vector or matrix with all values adjusted 0-1 scale relative to each other.
#' @export
scale0_1 <- function(x) {
  xAdj <- (x - min(x, na.rm = T)) / diff(range(x, na.rm = T))
  return(xAdj)
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
#' @param alignFile path to 3x3 alignment file to adjust pixel coordinates
#'     (optional).
#'    Stahl 2016 datasets (default: NA)
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
#' @param min.reads `cleanCounts` param; minimum number of reads to keep
#'     a gene (default: 1)
#' @param min.lib.size `cleanCounts` param; minimum number of counts a
#'     pixel needs to keep (default: 1)
#' @param min.detected `cleanCounts` param; minimum number of pixels a gene
#'     needs to have been detected in to keep (default: 1)
#' @param ODgenes Boolean to use ``getOverdispersedGenes`` for the corpus
#'    genes (default: TRUE)
#' @param nTopOD number of top over dispersed genes to use. int (default: 1000).
#'     If the number of overdispersed genes is less then this number will use all of them,
#'     or set to NA to use all overdispersed genes.
#' @param od.genes.alpha alpha parameter for `getOverdispersedGenes`.
#'     Higher = less stringent and more over dispersed genes returned (default: 0.05)
#' @param gam.k gam.k parameter for `getOverdispersedGenes`. Dimension
#'     of the "basis" functions in the GAM used to fit, higher = "smoother"
#'     (default: 5)
#' @param verbose control verbosity (default: TRUE)
#'
#' @return A list that contains
#' \itemize{
#' \item corpus: (pixels x genes) matrix of the counts of the selected genes
#' \item slm: slam::as.simple_triplet_matrix(corpus); required format for topicmodels::LDA input
#' \item positions: matrix of x and y coordinates of pixels rownames = pixels, colnames = "x", "y"
#' }
#'
#' @export
preprocess <- function(dat,
                       alignFile = NA,
                       extractPos = FALSE,
                       selected.genes = NA,
                       nTopGenes = NA,
                       genes.to.remove = NA,
                       removeAbove = NA,
                       removeBelow = NA,
                       min.reads = 1,
                       min.lib.size = 1,
                       min.detected = 1,
                       ODgenes = TRUE,
                       nTopOD = 1000,
                       od.genes.alpha = 0.05,
                       gam.k = 5,
                       verbose = TRUE) {
  
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
  
  if(verbose) {
    cat("Initial genes:", dim(t(counts))[1], "Initial pixels:", dim(t(counts))[2], "\n")
  }
  
  # use specific genes in the corpus
  if (is.na(selected.genes[1]) == FALSE) {
    if(verbose){
      cat("- Using genes in `selected.genes` for corpus.", "\n")
    }
    counts <- counts[,selected.genes]
    if(verbose){
      cat(" ", dim(counts)[2], "genes are present in dataset.", "\n")
    }
  }
  # remove poor pixels and genes
  if(verbose){
    cat("- Removing poor pixels with <=", min.lib.size, "reads", "\n")
    cat("- Removing genes with <=", min.reads, "reads across pixels and detected in <=", min.detected, "pixels", "\n")
  }
  countsClean <- cleanCounts(counts = t(counts),
                             min.reads = min.reads,
                             min.lib.size = min.lib.size,
                             min.detected = min.detected,
                             plot=TRUE,
                             verbose=FALSE)
  countsClean <- as.matrix(countsClean)
  if(verbose){
    cat("  Remaining genes:", dim(countsClean)[1], "and remaining pixels:", dim(countsClean)[2], "\n")
  }
  
  # adjust pixel coordinates, optional.
  # based on Stahl 2016 ST data with alignment matrices
  if (is.na(alignFile) == FALSE) {
    if (file.exists(alignFile) == TRUE){
      if(verbose){
        cat("- Adjusting pixel positions based on alignment file.")
      }
      align <- matrix(unlist(read.table(alignFile)), nrow = 3, ncol = 3)
      (positions[,"x"] * align[1,1]) - 290 # note that I found they were off by one pixel distance in pixels
      (positions[,"y"] * align[2,2]) - 290
    } else {
      cat("Warning: `alignFile` path does not exists. Skipping position adjustments.", "\n")
    }
  }
  
  # remove top expressed genes (nTopGenes needs to be integer or NA)
  if (is.na(nTopGenes) == FALSE & is.numeric(nTopGenes) == TRUE) {
    nTopGenes <- round(nTopGenes)
    if(verbose){
      cat("- Removing the top", nTopGenes, "expressed genes.", "\n")
    }
    top_expressed <- names(Matrix::rowSums(countsClean)[order(Matrix::rowSums(countsClean),
                                                      decreasing = TRUE)][1:nTopGenes])
    countsClean <- countsClean[rownames(countsClean) %in% top_expressed == FALSE,]
    
    # print(paste("after removing top ",
    #             as.character(nTopGenes),
    #             " genes:", dim(countsClean)[1], "genes remain."))
  }
  
  # remove specific genes (if there are any). Use grepl to search for gene name pattern matches
  if (is.na(genes.to.remove) == FALSE) {
    countsClean <- countsClean[!grepl(paste(genes.to.remove, collapse="|"), rownames(countsClean)),]
    if(verbose){
      cat("- After filtering for `genes.to.remove`:", "\n",
          " Remaining genes:", dim(countsClean)[1], "\n")
    }
  }
  
  # remove genes that appear in more than certain percentage of the pixels
  if (is.na(removeAbove) == FALSE & is.numeric(removeAbove) == TRUE){
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
        cat("- Removed genes present in",
            as.character(removeAbove*100), "% or more of pixels", "\n",
            " Remaining genes:", dim(countsClean)[1], "\n")
      }
    } else {
      cat("Warning: `removeAbove` must be a numeric from 0 to 1. Skipping `removeAbove` gene filter.", "\n")
    }
  }
  # remove genes that appear in less than certain percentage of the pixels
  if (is.na(removeBelow) == FALSE & is.numeric(removeBelow) == TRUE){
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
        cat("- Removed genes present in",
            as.character(removeBelow*100), "% or less of pixels", "\n",
            " Remaining genes:", dim(countsClean)[1], "\n")
      }
    } else {
      cat("Warning: `removeBelow` must be a numeric from 0 to 1. Skipping `removeBelow` gene filter.", "\n")
    }
  }
  
  # use overdispersed variable genes for corpus
  if (ODgenes == TRUE) {
    if(verbose){
      cat("- Capturing only the overdispersed genes...", "\n")
    }
    par(mfrow=c(4,2), mar=c(1,1,1,1))
    OD <- getOverdispersedGenes(countsClean,
                                alpha = od.genes.alpha,
                                gam.k = gam.k,
                                plot = TRUE,
                                details = TRUE)
    
    # option to select just the top n overdispersed genes based on
    # log p-val adjusted
    if (is.na(nTopOD) == FALSE){
      if(verbose){
        cat("- Using top", nTopOD, "overdispersed genes.", "\n")
      }
      OD_filt <- OD$df[OD$ods,]
      # check if actual number of OD genes less than `nTopOD`
      if (dim(OD_filt)[1] < nTopOD){
        if(verbose){
          cat(" number of top overdispersed genes available:", dim(OD_filt)[1], "\n")
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
    cat("Genes in corpus > 1000 (", dim(countsClean)[1],
        "). This may cause model fitting to take a while. Consider reducing the number of genes.", "\n")
  }
  
  corpus <- t(as.matrix(countsClean))
  
  # last filter: each row must have at least 1 non-zero entry
  # to be compatible with `topicmodels`.
  if(verbose){
    cat("- Check that each pixel has at least 1 non-zero gene count entry..", "\n")
  }
  corpus <- corpus[which(!Matrix::rowSums(corpus) == 0),]
  corpus_slm <- slam::as.simple_triplet_matrix(corpus)
  cat("Final corpus:", "\n")
  print(corpus_slm)
  
  # get pixel positions if pixel colnames contain the positions
  # this is the case for some ST datasets like Stahl 2016 sets
  # ex: "20x30" -> 20 x, 30 y positions
  if (extractPos) {
    if(verbose){
      cat("Extracting positions from pixel names.", "\n")
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
  
  cat("Preprocess complete.", "\n")
  return(list(corpus = corpus,
              slm = corpus_slm,
              pos = positions))
}


#' Match deconvolved cell-types to ground truth cell-types based on transcriptional profiles
#'
#' @description Match deconvolved cell-types to ground truth cell-types by testing for
#'     enrichment of ground truth marker gene sets in the deconvolved transcriptional profiles.
#'     Uses `liger::iterative.bulk.gsea`.
#'
#' @param beta cell-type (rows) x gene (columns) matrix of deconvolved cell-type transcriptional profiles
#' @param gset named list where each entry is a vector of marker genes for a given ground truth cell-type.
#' @param qval adjusted p-value threshold (default: 0.05)
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
annotateCellTypesGSEA <- function(beta, gset, qval = 0.05) {
  
  results <- list()
  top.pos.enrich <- c()
  
  for (i in seq(nrow(beta))){
    celltype <- i
    vals <- sort(beta[celltype,], decreasing = TRUE)
    
    gsea.results <- liger::iterative.bulk.gsea(values=vals, set.list=gset, rank=TRUE)
    
    # filter for top hits
    gsea.sig <- gsea.results[gsea.results$q.val < qval,]
    
    ## order of selection:
    ## 1. q-val
    ## 2. edge (the leading edge subset of a gene set is the subset of genes that contribute most to the Expression Score)
    ## 3. sscore (Expression Score)
    gsea.sig <- gsea.sig[order(gsea.sig$q.val, rev(gsea.sig$edge), rev(gsea.sig$sscore)), ]
    
    results[[ rownames(beta)[celltype] ]] <- gsea.sig
    
    ## the top entry that is also positiviely enriched in the txn profile is predicted to be the best matching
    gsea.sig.pos <- gsea.sig[which(gsea.sig$sscore > 0), ]
    top.pos.enrich <- append(top.pos.enrich, rownames(gsea.sig.pos)[1])
    
  }
  
  names(top.pos.enrich) <- rownames(beta)
  
  return(list(results = results,
              predictions = top.pos.enrich))
  
}

