# Additional functions specific for generating figures but not necessary for 
# the actual STdeconvolve pipeline


#' Aggregate cell-types together using dynamic tree cutting.
#'
#' @param beta Beta matrix (cell-type gene distribution matrix)
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
  message("cell-types combined.")
  
  # if theta, make cell-types the columns again
  if (type == "t") {
    combinedTopics <- t(combinedTopics)
  }
  return(combinedTopics)
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
#' @param corpus If corpus is NULL, then it will use the original corpus that
#'     the model was fitted to. Otherwise, compute deconvolved topics from this
#'     new corpus. Needs to be pixels x genes and nonnegative integer counts. 
#'     Each row needs at least 1 nonzero entry (default: NULL)
#' @param perc.filt proportion threshold to remove cell-types in pixels (default: 0.05)
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
buildLDAobject <- function(LDAmodel,
                           corpus = NULL,
                           perc.filt = 0.05,
                           clustering = "ward.D",
                           dynamic = "hybrid",
                           deepSplit = 4,
                           colorScheme = "rainbow",
                           plot = TRUE){
  
  # get beta and theta list object from the LDA model
  m <- getBetaTheta(LDAmodel, corpus = corpus, perc.filt = perc.filt)
  
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
    levels(cols) <- grDevices::rainbow(length(levels(cols)))
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


#' Generate heatmap of correlations
#' 
#' @description Visualize the correlations between topics stored in a matrix, typically one
#'     returned via `getCorrMtx()`. This function uses gplots::heatmap.2.
#'
#' @param mat matrix with correlation values from -1 to 1
#' @param rowLabs y-axis label for plot. These are the rows of the matrix, or specifically m1 from getCorrMtx. (default: NULL)
#' @param colLabs x-axis label for plot. These are the columns of the matrix, or specifically m2 from getCorrMtx. (default: NULL)
#' @param rowv cluster order for the rows to make dendrogram (RowV in heatmap.2()). (default: NA)
#' @param colv cluster order for the columns to make dendrogram (Colv in heatmap.2()). (default: NA)
#' @param margins set margins of the plot. (default: c(6,8))
#' @param textSize set the text size for both x-axis and y-axis labels. (default: 0.9)
#' 
correlationPlot_2 <- function(mat, rowLabs = NA, colLabs = NA, rowv = NA, colv = NA, margins = c(6,8), textSize = 0.9) {
  
  correlation_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(n = 209)
  correlation_breaks <- c(seq(-1,-0.01,length=100),
                          seq(-0.009,0.009,length=10),
                          seq(0.01,1,length=100))
  
  plt <- gplots::heatmap.2(x = mat,
                           density.info = "none",
                           trace = "none",
                           Rowv = rowv,
                           Colv = colv,
                           col = correlation_palette,
                           breaks = correlation_breaks,
                           key.xlab = "Correlation",
                           xlab = colLabs,
                           ylab = rowLabs,
                           key.title = NA,
                           cexRow = textSize,
                           cexCol = textSize,
                           # lhei = c(1,3),
                           margins = margins)
  
  return(plt)
}


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
#' @return a hash table where each key is a bregma ID
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
  mm <- stats::model.matrix(~ 0 + factor(cellTypes))
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
#'     WARNING: SPOTlight functions crash the R session if you use a corpus with gene(s)
#'     that is not in the NMF W. I think this happens in `SPOTlight::predict_spatial_mixtures()`
#'     during the nnls step at the end. Because it generates the full W matrix after selecting
#'     genes in the ST mtrx. When the NMF was trained, there was a step that made it only
#'     train on the intersecting genes with the ST data set. So under normal circumstances this would not be noticed
#'     because all genes should match between W and ST mtrx. But if one were to use a different
#'     ST mtrx on the trained NMF, then I bet it would crash because the matrices end up not being comparable shapes.
#' 
#' 
#' @param nmfRef the list returned from `SPOTlight::train_nmf()`. For `SPOTlight::spotlight_deconvolution()`,
#'     this returned list actually has this nmfRef input list as its first element
#' @param stCounts ST count matrix to deconvolve, genes x spots
#' @param min_cont param of `mixture_deconvolution_nmf()`. remove topics less than this percent in a spot
#' @param normCtTopicProfiles if TRUE, uses the normalized cell-type topic-proportions instead of the 
#'     unnormalized cell-type topic coefficients that is typically used by SPOTlight. (default: FALSE)
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
SPOTlightPredict <- function(nmfRef, stCounts, min_cont = 0.0, normCtTopicProfiles = FALSE) {
  
  # get basis matrix W [genes x topics]
  w <- NMF::basis(nmfRef[[1]])
  colnames(w) <- paste("Topic", 1:ncol(w), sep = "_")
  # get coefficient matrix H [topics x cells]
  h <- NMF::coef(nmfRef[[1]])
  
  
  # reference for which cell type(s) a topic represents
  # ct_topic_profiles = [topics x CellTypes] (seurat clusters)
  # values are coefficients
  # uses the H matrix and the cluster labels to get a new mtx where topic x cluster,
  # Cell types can be made up of multiple topics...
  ct_topic_profiles <- SPOTlight::topic_profile_per_cluster_nmf(h = h,
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
  # This function appears only using genes in ST that are also in W
  # NNLS to get topic coefficients for each spot based on genes associated with each topic (W)
  topics_in_spots <- SPOTlight::predict_spatial_mixtures_nmf(nmf_mod = nmfRef[[1]],
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
  # celltypes below `min_cont` are not counted in a spot
  # within this function, `topics_in_spots` is also made but not returned
  # uses the raw unnormalized celltype topic coefficients (ct_topic_profiles)
  # could use the normalized and adjusted cell-type topic proportions (ct_topics) but 
  # this will cause a slight change in the reported pixel cell-type proportions
  # SPOTlight code uses the unnormalized cell-type topic coefficients
  if(normCtTopicProfiles){
    topic_profiles <- ct_topics
  } else {
    topic_profiles <- ct_topic_profiles
  }
  ct_in_spots <- SPOTlight::mixture_deconvolution_nmf(nmf_mod = nmfRef[[1]],
                                           mixture_transcriptome = stCounts,
                                           transf = "uv", 
                                           reference_profiles = topic_profiles, 
                                           min_cont = min_cont)
  # note that last column is an additional column for the residual error
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

