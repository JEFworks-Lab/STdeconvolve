# Additional functions specific for generating figures but not necessary for 
# the actual STdeconvolve pipeline


#' Build hash table of simulated spots for each bregma layer of a given
#' MERFISH experiment.
#' 
#' @description uses `hash` to create a hash table of each bregma for a given
#'     input MERFISH experiment. Goal is to generate simulated spots for each bregma
#'     by partitioning cells into spots based on their Centroid coordinates. The
#'     size of the simulated spots is chosen as well but default is 100um x 100um.
#'     Note that for a given experiment/dataset there are multiple bregma layers. 
#'     Important: "Cell_class" column is where the cell labels will be pulled
#'     from. These will be used to construct the: "cellTypeTable", where downstream
#'     in `buildBregmaCorpus`, will be used to make "gtSpotTopics" and "gtCtGenes".
#' 
#' @param cellCentroidsAndClass an "annot.table" data.frame from
#'     "mpoa_merfish_clean.RData" that has been filtered for cells (the rows)
#'     of a given merfish experiment and has the following columns:
#'     c('Centroid_X', 'Centroid_Y', 'Bregma', "Cell_class", "Neuron_cluster_ID")
#' @param counts cell x gene count matrix for the individual cells in all the
#'     bregmas for which spots will be simulated
#' @param patch_size dimension in um to make the simulated spots for a MERFISH
#'     experiment. The Centroid coords are already in um. (default: 100)
#' 
#' @return a hash table where each key is a bregma ID (ie [["-0.04]])
#'     and the returned object is a list that contains
#' \itemize{
#' \item bregmaFullDf: the "annot.table" data.frame for the specific bregma with
#'     a new column "patch_id" indicating which patch a cell is assigned to. Written
#'     as "xcoord_ycoord". Simulated spots on the edges are dropped so some cells
#'     are not assigned to a patch and their patch ID is ""
#' \item cellTypeTable: table of counts of different cell types in each simulated patch
#' \item patchTotalCells: vector of total cells in each simulated patch
#' \item cellTypeCount: vector of counts of unique cell types in each simulated patch
#' \item cellGexp: individual cell x gene counts matrix for the specific cells
#'     in the bregma
#' \item patchGexp: patch x gene count matrix; ie the simulated gene counts for
#'     each simulated patch
#' }
#' 
#' @export
#' 
simulateBregmaSpots <- function (cellCentroidsAndClass, counts, patch_size = 100) {
  
  # dictionary hash table
  h <- hash::hash()
  
  data <- cellCentroidsAndClass
  
  for (bregma in unique(data$Bregma)) {
    
    bregma_key <- as.character(bregma)
    print(bregma_key)
    
    selected_bregma <- data[which(data$Bregma == bregma),]
    
    # 1. Get patch edge coordinates:
    
    # Sequence of X-coord positions for left edge of each patch:
    x_edges <- seq(min(selected_bregma$Centroid_X), max(selected_bregma$Centroid_X), patch_size)
    # drop first and last to avoid any issues with the edges of the whole region
    inner_x_edges <- x_edges[2:(length(x_edges)-1)]
    # Sequence of Y-coord positions for bottom edge of each patch:
    y_edges <- seq(min(selected_bregma$Centroid_Y), max(selected_bregma$Centroid_Y), patch_size)
    inner_y_edges <- y_edges[2:(length(y_edges)-1)]
    
    selected_bregma$patch_id <- character(length(rownames(selected_bregma)))
    
    # 2. add patch IDs to bregma cells, for the patch they belong to:
    
    for (x in inner_x_edges) {
      for (y in inner_y_edges) {
        patch_id <- paste0(as.character(x), "_", as.character(y))
        patch <- selected_bregma[which( (selected_bregma$Centroid_X > x) &
                                          (selected_bregma$Centroid_X < x+patch_size) &
                                          (selected_bregma$Centroid_Y > y) &
                                          (selected_bregma$Centroid_Y < y+patch_size) ),]
        
        if (length(rownames(patch)) > 0) {
          selected_bregma[rownames(patch),]$patch_id <- patch_id
        }
      }
    }
    
    # 3. get table of counts of cell types for each patch in bregma
    patch_cell_types <- selected_bregma[which(selected_bregma$patch_id != ""),
                                               c("patch_id", "Cell_class")]
    cellTypeTable <- table(patch_cell_types[])
    
    # 4. total cell counts in each patch
    patchTotalCells <- rowSums(cellTypeTable)
    
    # 5. counts of unique cell types in each patch
    cellTypeCount <- c()
    for (i in seq_len(length(rownames(cellTypeTable)))) {
      patch_num_cell_types <- length(which(cellTypeTable[i,] != 0))
      cellTypeCount <- append(cellTypeCount, patch_num_cell_types)
    }
    
    # 6. collapse gene counts for cells in same spot to make simulation
    patches <- unique(selected_bregma$patch_id[which(!selected_bregma$patch_id == "")])
    patchGexp <- do.call(rbind, lapply(patches, function(patch){
      cells <- rownames(selected_bregma[which(selected_bregma$patch_id == patch),])
      mat <- as.matrix(counts[cells,])
      if (length(cells) == 1){
        patch_counts <- as.vector(mat)
      } else if (length(cells) > 1){
        patch_counts <- colSums(mat)
      } else if (length(cells) == 0){
        cat("WARNING:", bregma, "patch", patch, "had no cells in `counts` to use for simulated gene expression", "\n")
      }
      patch_counts
    }))
    rownames(patchGexp) <- patches
    
    # 7. gene count matrix for individual cells in the bregma
    # this also includes cells not assigned to patches
    bregma_cells <- rownames(selected_bregma)
    cellGexp <- counts[bregma_cells,]
    
    # 8. combine data objects and append to hash table
    h[[bregma_key]] <- list(bregmaFullDf = selected_bregma,
                            cellTypeTable = cellTypeTable,
                            patchTotalCells = patchTotalCells,
                            cellTypeCount = cellTypeCount,
                            cellGexp = cellGexp,
                            patchGexp = patchGexp)
  }
  return(h)
}



#' Generate simulated corpus as input into `topicmodels` and the matched ground
#' truth spot-topic proportion and topic-gene proportion matrices for a given
#' simulated bregma spot dataset (built with `simulateBregmaSpots`)
#' 
#' @description Important: "Cell_class" column in "bregmaFullDf" as input for `simulateBregmaSpots`
#'     is where the cell labels will be pulled from. These will be used to construct
#'     the: "cellTypeTable", where downstream in `buildBregmaCorpus`, will be used
#'     to make "gtSpotTopics" and "gtCtGenes".
#' 
#' @param hashTable output hash table of simulated bregma spots from
#'     `simulateBregmaSpots`
#' @param bregmaID ID to reference a given bregma in the hashTable. ex: "-0.04"
#' 
##' @return a list that contains
#' \itemize{
#' \item sim: slam::as.simple_triplet_matrix(corpus); required format for topicmodels::LDA input
#' \item gtSpotTopics: the ground truth of major cell class proportions in each simulated patch.
#'     Use as ground truth "theta" (document x topix proportions)
#' \item gtCtGenes: the normalized gene count proportions of each of the major cell classes.
#'     Use as ground truth "beta" (topics x word frequencies)
#' \item cellCounts: data.frame that has patch names, "x", and "y" centroid coordinates,
#'     and counts of total cells in each patch
#' \item classColors: vector of colors for each cell class
#' \item annotDf: data.frame of the individual cells of the bregma,
#'     their position, cell type label, and assigned patch
#' }
#' 
#' @export
#'
buildBregmaCorpus <- function (hashTable, bregmaID) {
  
  bregmaID <- as.character(bregmaID)
  
  # corpus in slam format
  sim <- hashTable[[bregmaID]][["patchGexp"]]
  sim <- sim[order(rownames(sim)),]
  # print(sim)
  # remove "Blanks" from data:
  sim <- sim[,!grepl("Blank", colnames(sim))]
  sim <- slam::as.simple_triplet_matrix(sim)
  
  # ground truth spot - cell type proportions
  gtSpotTopics <- hashTable[[bregmaID]][["cellTypeTable"]]/rowSums(hashTable[[bregmaID]][["cellTypeTable"]])
  gtSpotTopics <- as.data.frame.matrix(gtSpotTopics[order(rownames(sim)),])
  
  # reformat `gtDocTopic` proportions into data frame with spot coordinates
  # tmp_positions <- do.call(rbind, lapply(rownames(gtSpotTopics), function(x){
  #   coords <- strsplit(x, "_")[[1]]
  #   as.numeric(coords)
  # }))
  # colnames(tmp_positions) <- c("x", "y")
  # rownames(tmp_positions) <- rownames(gtSpotTopics)
  # tmp_proportions <- lapply(colnames(gtSpotTopics), function(i) {
  #   gtSpotTopics[,i]
  # })
  # names(tmp_proportions) <- colnames(gtSpotTopics)
  # gtSpotTopics <- merge(tmp_positions, as.data.frame(tmp_proportions), by="row.names")
  
  # the individual cell annotation df but only for cells assigned to spots
  df <- hashTable[[bregmaID]][["bregmaFullDf"]]
  df <- df[which(df$patch_id != ""),]
  
  # cells and their assigned class (only cells in simulated patches)
  cellTypes <- df[,c("Cell_class")]
  cells <- rownames(df)
  
  # ground truth cell type - gene expression frequencies 
  mat <- hashTable[[bregmaID]][["cellGexp"]][cells,]
  mm <- model.matrix(~ 0 + factor(cellTypes))
  colnames(mm) <- levels(factor(cellTypes))
  gtCtGenes <- t(t(as.matrix(mat)) %*% mm)
  # remove "Blanks" from data:
  gtCtGenes <- gtCtGenes[,!grepl("Blank", colnames(gtCtGenes))]
  gtCtGenes <- gtCtGenes/rowSums(gtCtGenes)
  
  # number of total cells in each spot
  cell_counts <- hashTable[[bregmaID]]$patchTotalCells
  count_df <- do.call(rbind, lapply(names(cell_counts), function(x){
    coords <- strsplit(x, "_")[[1]]
    as.numeric(coords)
  }))
  colnames(count_df) <- c("x", "y")
  rownames(count_df) <- names(cell_counts)
  count_df <- as.data.frame(count_df)
  count_df$counts <- cell_counts
  # same order as simulated spots
  count_df <- count_df[order(rownames(sim)),]
  
  # colors for each cell class
  classColors <- gg_color_hue(length(unique(df$Cell_class)))
  names(classColors) <- names(unique(df$Cell_class))
  
  bregma <- list(sim = sim,
                 gtSpotTopics = gtSpotTopics,
                 gtCtGenes = gtCtGenes,
                 cellCounts = count_df,
                 classColors = classColors,
                 annotDf = df)
  
  return(bregma)
}


#' Wrapper around functions in `SPOTlight::spotlight_deconvolution()`
#' that take place subsequently after training the NMF model.
#' 
#' @description Allows for applying a trained NMF model to additional
#'     ST datasets without retraining a new NMF in spotlight, which can take
#'     a very long time (hours depending on number of genes and clusters in scRNAseq ref).
#'     
#'     Additionally, other useful matrices and information can also be collected
#'     that are not returned by `SPOTlight::spotlight_deconvolution()`.
#'     
#'     This does come with a caveat that one must be careful about:
#'     WARNING: SPOTlight functions crash the R session if you use a corpus with a gene set
#'     that is not in the NMF W. I think this happens in `SPOTlight::predict_spatial_mixtures()`
#'     during the nnls step at the end. Because it generates the full W matrix after selecting
#'     genes in the ST mtrx. When the NMF was trained, there was a step that made it only
#'     train on the intersecting genes with the ST data set. So under normal circumstances this would not be noticed
#'     because all genes should match between W and ST mtrx. But if one were to use a different
#'     ST mtrx on the trained NMF, then I bet it would crash because the matrices end up not being comparable shapes.
#' 
#' 
#' @param nmfRef the list returned from `SPOTlight::train_nmf()`. In `SPOTlight::spotlight_deconvolution()`,
#'     this returned list actually has this nmf input list as its first element
#' @param stCounts ST count matrix to deconvolve, genes x spots
#' 
#' @return a list that contains
#' \itemize{
#' \item betaTopics: t(w), # topic x gene weights that sum to 1 (all genes and topics)
#' \item betaCt: ct_beta, # CellType x gene weights that sum to 1 (all genes and topics)
#' \item ctTopicProps: ct_topics, # topic proportion for each cell type
#' \item thetaTopics: t(topics_in_spots_norm), # proportion of each topic in each spot (all topics)
#' \item thetaCt: ct_in_spots_clean)) # spot x cell type proportions (all spots and cell types)
#' }
#'
#' @export
SPOTlightPredict <- function(nmfRef, stCounts) {
  
  # get basis matrix W [genes x topics]
  w <- basis(nmfRef[[1]])
  colnames(w) <- paste("Topic", 1:ncol(w), sep = "_")
  # get coefficient matrix H [topics x cells]
  h <- coef(nmfRef[[1]])
  
  
  # reference for which cell type(s) a topic represents
  # ct_topic_profiles = [topics x CellTypes] (seurat clusters)
  # values are coefficients
  # uses the H matrix and the cluster labels to get a new mtx where topic x cluster,
  # Cell types can be made up of multiple topics...
  ct_topic_profiles <- topic_profile_per_cluster_nmf(h = h,
                                                     train_cell_clust = nmfRef[[2]])
  
  # convert to topic proportions for each cell type
  ct_topic_profiles_r <- round(ct_topic_profiles, 4)
  ct_topics <- do.call(cbind, lapply(seq(ncol(ct_topic_profiles_r)), function(i){
    ct <- ct_topic_profiles_r[,i]
    ct/sum(ct)
  }))
  rownames(ct_topics) <- paste("Topic", 1:nrow(ct_topic_profiles_r), sep = "_")
  colnames(ct_topics) <- colnames(ct_topic_profiles_r)
  
  # now that we have the proportion of each topic for a cell type,
  # use this information to get gene weights (sum to 1) for each cell type
  # (w is gene weights for each topic)
  # do this by getting the topic-gene vectors from w for topics that make up a given cell type and
  # merge the gene values of these vectors based on the proportion that the topic contributes to the cell type
  ct_beta <- do.call(rbind, lapply(seq(ncol(ct_topics)), function(ct){
    # topics for the given cell type in ct_topics
    ts <- ct_topics[,ct]
    # topics that contribute to the cell type and their proportions
    ts <- ts[which(ts > 0)]
    # w topics-gene vectors (rows = genes, cols = topics for given ct)
    g <- w[,names(ts)]
    if (length(ts) > 1){
      # if more than 1 topic contributes to ct, g will be mtx with multiple columns (one for each topic)
      # mtx multiply to multiply first column (topic) gene values by first (topic proportion),
      # second col by second proportion, etc
      g <- g %*% diag(ts)
    }
    # after adjusting gene values for each topic column by the proportion the given topic contributes to a cell type,
    # sum the adjusted topic gene vectors together
    if (is.null(dim(g)) == FALSE){
      g <- rowSums(g)
    }
    
    # final adjustment to make topic gene proportions sum to 1
    g <- g/sum(g)
    g
  }))
  rownames(ct_beta) <- colnames(ct_topics)
  
  
  # [topics x spots]; topic coefficients for each spot
  # This function appears to account for that only using genes in ST that are also in W
  # NNLS to get topic coefficients for each spot based on genes associated with each topic (W)
  topics_in_spots <- predict_spatial_mixtures_nmf(nmf_mod = nmfRef[[1]],
                                                  mixture_transcriptome = stCounts,
                                                  transf = "uv")
  # normalize coefficients to get topic proportions in each spot
  topics_in_spots_norm <- do.call(cbind, lapply(seq(ncol(topics_in_spots)), function(i){
    ct <- topics_in_spots[,i]
    ct/sum(ct)
  }))
  rownames(topics_in_spots_norm) <- paste("Topic", 1:nrow(topics_in_spots), sep = "_")
  colnames(topics_in_spots_norm) <- colnames(topics_in_spots)
  
  # returns [spot x CellType]; proportion of each cell type in each spot
  # topics below `min_cont` are not counted in a spot
  # within this function, `topics_in_spots` is also made but not returned
  ct_in_spots <- mixture_deconvolution_nmf(nmf_mod = nmfRef[[1]],
                                           mixture_transcriptome = stCounts,
                                           transf = "uv", 
                                           reference_profiles = ct_topics, 
                                           min_cont = 0.0) # only keep topics if 9% or more in a spot
  # note that last column is an additional columns for the residual error
  # cleanup to get actual spot-celltype predictions:
  rownames(ct_in_spots) <- colnames(stCounts)
  ct_in_spots_clean <- ct_in_spots[,1:(ncol(ct_in_spots)-1)] # last column is residuals
  # note: all spots and cell types.
  # Possible some cell type will not be detected at all and be 0 for all spots
  
  return(list(betaTopics = t(w), # topic x gene weights that sum to 1 (all genes and topics)
              betaCt = ct_beta, # CellType x gene weights that sum to 1 (all genes and topics)
              ctTopicProps = ct_topics, # topic proportion for each cell type
              thetaTopics = t(topics_in_spots_norm), # proportion of each topic in each spot (all topics)
              thetaCt = ct_in_spots_clean)) # spot x celltype proportions (all spots and cell types)
  
  # if made with original raw `stCounts` then will have thousands of genes and all spots.
  # so when comparing, will need to make sure the betas or thetas have same spots and genes as
  # the `stCounts` they might actually be compared to (corpus for a given topic model)
  # however, a raw `stCounts` can be used to make the SPOTlight predictions at first, and
  # after the betas and theta can be filtered. OR, a processed corpus could be used for
  # SPOTlight, which will probably affect the predicted topics because the gene set used
  # will be smaller. Worth comparing both cases.
  
  # WARNING: SPOTlight functions crash the session if use a corpus with a gene set
  # that is not in the NMF. I think this happens in `SPOTlight::predict_spatial_mixtures()`
  # during the nnls step at the end. Because it generates the full W matrix after selecting
  # genes in the ST mtrx. When the NMF was trained, there was a step that made it only
  # train on the intersecting genes. So under normal circumstances this would not be noticed
  # because all genes should match between W and ST mtrx.
  # But if one were to use a different ST mtrx on the trained NMF, then I bet it would crash
  # because the matrices end up not being comparable shapes.
  
  # looking at the code, it seems the cd is bigger than w. So actually cd is being trimmed
  # such that it's genes equal W. If other way around it all crashes
  
  # remember that during training of the NMF, it got cluster genes, and then only kept
  # those that were in cd.
  # with the corpuses, genes could come up as variable genes that were not picked
  # as cluster genes for the NMF, and thus would not be in the NMF at all
  
}



