#' Collection of functions for STDeconvolve with Topic Modeling
#' 
#' 
#' 

#' convert a matrix to a long form data frame for easy plotting via ggplot2
#' 
#' mtx = the matrix to be converted
#' colLab = df column label for the mtx columns, whatever they may represent
#' rowLab = df row label for the row labels, whatever they may represent
#' cellLab = label for the cell values, whatever they may represent
#' 
#' example:
#'        [,1] [,2] [,3]
#' [1,]    6    9    1
#' [2,]    7   10    2
#' [3,]    8    3   12
#' [4,]    4   11    5
#' becomes:
#'    columns rows cells
#'  1        1    1     6
#'  2        2    1     9
#'  3        3    1     1
#'  4        1    2     7
#'  5        2    2    10
#'  6        3    2     2
#'  7        1    3     8
#'  8        2    3     3
#'  9        3    3    12
#'  10       1    4     4
#'  11       2    4    11
#'  12       3    4     5
#' 
mtx2ggplotDf <- function(mtx,
                         colLab = "columns", rowLab = "rows", cellLab = "cells") {
  
  if (is.null(colnames(mtx))) {
    colvals <- 1:ncol(mtx)
  } else { 
    colvals <- colnames(mtx)
  }
  
  if (is.null(rownames(mtx))) {
    rowvals <- 1:nrow(mtx)
  } else { 
    rowvals <- rownames(mtx)
  }
  
  dat <- data.frame(columns = factor(rep(colvals, dim(mtx)[1])), # columns
                    rows = factor(rep(rowvals, each = dim(mtx)[2])), # rows
                    cells = as.vector(t(mtx))) # cell values
  
  colnames(dat) <- c(as.character(colLab),
                     as.character(rowLab),
                     as.character(cellLab))
  return(dat)
}


#' cluster topics together using dynamic tree cutting.
#' Topics clustered by beta matrix (topic-word distributions)
#'
#' returns:
#' a list that contains:
#' order = vector of the dendrogram index order for the topics
#' clusters = factor of the topics (names) and their assigned cluster (levels)
#' dendro = dendrogram of the clusters
#'
clusterTopics <- function(beta,
                          distance = "euclidean", clustering = "ward.D", dynamic = "hybrid",
                          deepSplit = 4,
                          plotDendro = TRUE) {
  
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
  
  d_ <- dist(beta, method = distance)
  hc_ <- hclust(d_, method = clustering)
  
  if (plotDendro) {
    plot(hc_)
  }
  
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
  
  return(list(clusters = groups,
              order = hc_$order,
              dendro = as.dendrogram(hc_)))

}


#' Train vanilla LDA model from `topicmodels`
#' 
#` corpus - slam::as.simple_triplet_matrix; docs x words
#` Ks - vector of K parameters to search
#'
#' returns list of each fitted LDA model, the optimal K
#' perplexity scores for each model, and the corpus
#' that was used to fit them.
#' 
#' Note:
#' access given model via: lda$models[[k]][1]
#' models are objects from the library `topicmodels`
#' LDA models have slots with additional information.
#' 
fitLDA <- function(corpus, Ks, seed = 0) {

  controls <- list(seed = seed,
                   verbose = 1, keep = 1, estimate.alpha = TRUE)
  
  fitted_models <- mclapply(Ks, function(k) {
    topicmodels::LDA(corpus, k=k, control = controls)
  },
  mc.cores = detectCores(logical = TRUE) - 1
  )
  
  pScores <- unlist(lapply(seq(length(Ks)), function(k){
    p <- topicmodels::perplexity(fitted_models[[k]], corpus)
    p
  }))
  
  # pScores <- c()
  # for (k in seq(length(Ks))) {
  #   p <- topicmodels::perplexity(fitted_models[[k]], corpus)
  #   pScores <- append(pScores, p)
  # }
  
  names(pScores) <- Ks
  pScoresNorm <- (pScores-min(pScores))/(max(pScores)-min(pScores))
  
  par(mfrow=c(2,1), mar=c(3,5,1,1))
  barplot(pScores)
  barplot(pScoresNorm)
  
  kOpt <- Ks[which(pScores == min(pScores))]
  # lda_model <- topicmodels::LDA(corpus, k=kOpt, control = controls)
  
  return(list(models = fitted_models,
              kOpt = kOpt,
              perplexities = pScores,
              fitCorpus = corpus))
}


#' get the optimal LDA model
#' 
#' models = list returned from `fitLDA`
#' 
optimalModel <- function(models) {
  models$models[[which(sapply(models$models, slot, "k") == models$kOpt)]]
}


#' custom correlation color range for heatmap.2 correlation plots
correlation_palette <- colorRampPalette(c("blue", "white", "red"))(n = 209)
correlation_breaks = c(seq(-1,-0.01,length=100),
                       seq(-0.009,0.009,length=10),
                       seq(0.01,1,length=100))

#' lighten and darken a color
lighten <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- col2rgb(color)
  col <- col + (255 - col)*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}


#' color palette to replicate ggplot2
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#' Visualize topic proportions across spots with `scatterpie`
#' 
#` theta = document x topic proportion matrix
#` pos = position of documents, x and y columns
#` topicOrder = order of topics based on dendrogram; a numeric vector
#'             from: (clusterTopics$order)
#'
#` cluster_cols = vector of colors for each of the topics
#'     can use factor of clusters for each topic (via clusterTopics$clusters)
#'     as long as the cluster levels/values have been converted to colors.
#'     This gets reordered wrt the topicOrder fyi. 
#'
#' groups = color spot piecharts group labels or cell layer they belong to.
#'          Needs to be a character vector. Ex: c("0", "1", "0", ...).
#' group_cols = color labels for the groups. Ex: c("0" = "gray", "1" = "red")
#'
#' r = radius of the circles. Adjust based on size of spots.
#'     40 = merfish; 0.4 = mOB
#'     
#' lwd = width of lines of the pie charts
#'
vizAllTopics <- function(theta, pos, topicOrder, cluster_cols,
                         groups = NA,
                         group_cols = NA,
                         r = 40,
                         lwd = 0.5,
                         showLegend = TRUE,
                         plotTitle = NA) {
  
  # reorder colors of topics wrt the dendrogram
  colors_ordered <- as.vector(cluster_cols[topicOrder])
  
  # doc-topic distribution reordered based on topicOrder
  theta_ordered <- theta[, topicOrder]
  theta_ordered <- as.data.frame(theta_ordered)
  # colnames(theta_ordered) <- paste0("Topic.", topicOrder)
  colnames(theta_ordered) <- paste0("Topic.", colnames(theta_ordered))
  
  # add columns with document positions
  # colnames(pos) <- c("x", "y")
  theta_ordered_pos <- merge(data.frame(theta_ordered),
                                 data.frame(pos), by=0)
  
  # column names of topics in DF, in order of topicOrder
  # topicColumns <- paste0("Topic.", topicOrder)
  
  # first column after merge is "Row.names", last two are "x" and "y"
  # problem is that data frame will replace "-" and " " with "."
  topicColumns <- colnames(theta_ordered_pos)[2:(dim(theta_ordered_pos)[2]-2)]
  
  # color of piechart groups (lines of piechart):
  if (is.na(groups[1]) == TRUE) {
    groups <- rep("0", dim(theta_ordered_pos)[1])
    theta_ordered_pos$groups <- groups
  } else {
    theta_ordered_pos$groups <- as.character(groups)
  }
  if (is.na(group_cols[1]) == TRUE) {
    group_cols <- c("0" = "gray")
  }
  
  p <- ggplot() +
    theme(panel.background = element_rect(fill = "black"),
          panel.grid = element_blank()) +
    # theme_classic() +
    geom_scatterpie(aes(x=x, y=y, group=Row.names, r=r, color = groups),
                    lwd = lwd,
                    data = theta_ordered_pos,
                    cols = topicColumns,
                    legend_name = "Topics") +
    scale_fill_manual(values=colors_ordered) +
    scale_color_manual(values = group_cols)
    coord_equal()
  
  if (showLegend == FALSE) {
    p <- p + guides(fill=FALSE)
  }
    
  if (is.na(plotTitle) == FALSE) {
    p <- p + ggtitle(plotTitle)
  }
  
  print(p)
}


#' Visualize proportions of topic clusters separately
#' 
#` theta - document topic proportion matrix
#` pos - position of documents, x and y columns
#` topicOrder - order of topics based on dendrogram; a numeric vector NOT NEEDED
#` clusters - factor of the color (topic cluster) each cluster is assigned to
#'            In this case, the levels should be colors. In `vizAllTopics`,
#'            clusters is "cluster_cols" and can just be a vector of colors.
#'            
#' groups = color spot piecharts based on a group or cell layer they belong to.
#'          Needs to be a character vector. Ex: c("0", "1", "0", ...).
#' group_cols = color labels for the groups. Ex: c("0" = "gray", "1" = "red")
#' 
#' r = radius of the circles. Adjust based on size of spots.
#'     40 = merfish; 0.4 = mOB
#'     
#' lwd = width of lines of the pie charts
#' 
vizTopicClusters <- function(theta, pos, clusters,
                             groups = NA,
                             group_cols = NA,
                             sharedCol = FALSE,
                             r = 40,
                             lwd = 0.5,
                             plotTitle = NA) {
  
  print("Topic cluster members:")
  # produce a plot for each topic cluster:
  for (cluster in levels(clusters)) {
    
    # select the topics in the cluster in the same order of the dendrogram
    topics <- labels(clusters[which(clusters == cluster)])
    
    cat(cluster, ":", topics, "\n")
    
    # doc-topic distribution reordered based on topicOrder and selected cluster topics
    theta_ <- theta[, topics]
    
    # get percentage of other topics not in cluster
    if (is.null(dim(theta_))) {
      other <- 1 - theta_
    } else {
      other <- 1 - rowSums(theta_)
    }
    
    theta_ <- as.data.frame(theta_)
    colnames(theta_) <- paste0("Topic.", topics)
    theta_$other <- other
    
    # if any topics not represented at all, drop them
    # Apparently if a topic is 0 for all pie charts, it is not plotted
    # and doesn't appear in the legend. So it messes with the colors.
    # "other" takes one of the colors of the topics and is not gray
    theta_ <- theta_[,which(!colSums(theta_) == 0)]
    # theta_ <- as.data.frame(theta_)
    # print(theta_)
    # add columns with document positions
    # rownames(theta_) <- rownames(pos)
    theta_pos <- merge(data.frame(theta_),
                                   data.frame(pos), by=0)
    
    # first column after merge is "Row.names", last two are "x" and "y"
    # problem is that data frame will replace "-" and " " with "."
    topicColumns <- colnames(theta_pos)[2:(dim(theta_pos)[2]-2)]
    
    # get a hue of colors representing the cluster color
    if (sharedCol){
      color_ramp <- colorRampPalette(c(cluster, cluster))
    } else {
      color_ramp <- colorRampPalette(c(lighten(cluster, factor = 0.8), darken(cluster, factor = 1)))
    }
    
    # topic_colors <- viridis_pal(option = "C")(length(blue_cluster))
    topic_colors <- color_ramp(ncol(theta_) - 1) # don't count "other"
    topic_colors <- append(topic_colors, c("gray"))
    
    # color of piechart groups (lines of piechart):
    if (is.na(groups[1]) == TRUE) {
      groups <- rep("0", dim(theta_pos)[1])
      theta_pos$groups <- groups
    } else {
      theta_pos$groups <- as.character(groups)
    }
    if (is.na(group_cols[1]) == TRUE) {
      group_cols <- c("0" = "gray")
    }
    
    p <- ggplot() +
      theme(panel.background = element_rect(fill = "black"),
            panel.grid = element_blank()) +
      # theme_classic() +
      geom_scatterpie(aes(x=x, y=y, group=Row.names, r = r, color = groups), # r=40 for MERFISH 0.4 mOB
                      lwd = lwd,
                      data=theta_pos,
                      cols = topicColumns,
                      legend_name = "Topics") +
      coord_equal() +
      scale_fill_manual(values=topic_colors) +
      scale_color_manual(values = group_cols)
    
    if (is.na(plotTitle) == FALSE) {
      p <- p + ggtitle(plotTitle)
    }
    
    print(p)
  }
}


#' Compute the cross correlation between top terms in each topic
#' correlation based on the term counts across the corpus documents
#'
#' beta = topic-term distribution matrix
#' corpus = document x term matrix. Counts of terms in each document
#' thresh = the probability of the topic generating the term
#' 
#' correlation = either "corpus" or "beta"
#'        Choose which matrix to perform the correlation on.
#'        Either wrt gene counts across the docs (corpus)
#'        or the term probability for each topic (beta)
#' 
#' returns:
#' list of the resulting term cross correlation matrices for each topic
#' also returns the correlation heat maps (optional)
#' 
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


#' Compute the "mean correlation" for each topic
#' based on the topic term correlation matrices
#'
#' topicCorrList = list of the resulting term cross correlation matrices
#'      for each topic. (topicTermCorrelationMats)
#'
#' returns:
#' dataframe with the correlations for each topic (cor) and the number of
#' topic terms used (numTerms)
#'
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


#' build hash table of different bregma experiments to use as input for
#' `extractBregmaCorpus`
#'
#'
build_bregma_hash_table <- function (cellCentroidsAndClass, patch_size) {
  
  # cellCentroidsAndClass -> `annot.table` from `mpoa_merfish_clean.RData` that has been filtered for cells of a given
  # merfish experiment and has: c('Centroid_X', 'Centroid_Y', 'Bregma', "Cell_class", "Neuron_cluster_ID")
  # patch_size -> size of patch in um. Centroid coords are already in um
  
  # dictionary hash table
  h <- hash()
  
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
    selected_bregma_patches <- selected_bregma[which(selected_bregma$patch_id != ""),
                                               c("patch_id", "Cell_class")]
    selected_bregma_patch_cells <- table(selected_bregma_patches[])
    
    # 4. total cell counts in each patch
    bregma_cell_counts <- rowSums(selected_bregma_patch_cells)
    
    # 5. counts of unique cell types in each patch
    unique_types_per_patch <- c()
    for (i in seq_len(length(rownames(selected_bregma_patch_cells)))) {
      patch_num_cell_types <- length(which(selected_bregma_patch_cells[i,] != 0))
      unique_types_per_patch <- append(unique_types_per_patch, patch_num_cell_types)
    }
    
    # 6. combine data objects and append to hash table
    h[[bregma_key]] <- list(bregmaFullDf = selected_bregma,
                            cellTypeTable = selected_bregma_patch_cells,
                            totalCells = bregma_cell_counts,
                            cellTypeCount = unique_types_per_patch)
  }
  
  return(h)
  
}


#' Extract corpus and ground truth doc-topic/topic-term matrices from a 
#' MERFISH bregma in the hash table
#'
#' For FN7_2_M22_M26_hash[["-0.04"]] specifically:
#' 
#' bregmas: [1] "-0.04" "-0.09" "-0.14" "-0.19" "-0.24" "-0.29"
#' 
#' what is in the hash table:
#' h[[bregma_key]] <- list(bregmaFullDf = selected_bregma, # use with cellGeneCounts to get topic-word proportions
#'                             cellTypeTable = selected_bregma_patch_cells, # gt document-topic proportions
# '                            totalCells = bregma_cell_counts,
# '                            cellTypeCount = unique_types_per_patch,
#'                             cellGeneCounts = cellGeneCounts, # use to get gt topic-word proportions
# '                            patchGeneCounts = patchGeneCounts) # the simulation
#'
extractBregmaCorpus <- function (hashTable, bregmaID) {
  
  bregmaID <- as.character(bregmaID)
  
  # corpus in slam format
  sim <- hashTable[[bregmaID]][["patchGeneCounts"]]
  sim <- sim[order(rownames(sim)),]
  # print(sim)
  # remove "Blanks" from data:
  sim <- sim[,!grepl("Blank", colnames(sim))]
  sim <- slam::as.simple_triplet_matrix(sim)
  
  # ground truth spot cell type proportions
  gtDocTopics <- hashTable[[bregmaID]][["cellTypeTable"]]/rowSums(hashTable[[bregmaID]][["cellTypeTable"]])
  gtDocTopics <- gtDocTopics[order(rownames(sim)),]
  
  # reformat `gtDocTopic` proportions into data frame with spot coordinates
  tmp_positions <- do.call(rbind, lapply(rownames(gtDocTopics), function(x){
    coords <- strsplit(x, "_")[[1]]
    as.numeric(coords)
  }))
  colnames(tmp_positions) <- c("x", "y")
  rownames(tmp_positions) <- rownames(gtDocTopics)
  tmp_proportions <- lapply(colnames(gtDocTopics), function(i) {
    gtDocTopics[,i]
  })
  names(tmp_proportions) <- colnames(gtDocTopics)
  gtDocTopics <- merge(tmp_positions, as.data.frame(tmp_proportions), by="row.names")
  
  
  df <- hashTable[[bregmaID]][["bregmaFullDf"]]
  df <- df[which(df$patch_id != ""),]
  cellTypes <- df[,c("Cell_class")]
  cells <- rownames(df)
  
  # ground truth topic words
  mat <- hashTable[[bregmaID]][["cellGeneCounts"]][cells,]
  mm <- model.matrix(~ 0 + factor(cellTypes))
  colnames(mm) <- levels(factor(cellTypes))
  gtTopicWords <- t(t(as.matrix(mat)) %*% mm)
  gtTopicWords <- gtTopicWords/rowSums(gtTopicWords)
  
  # number of total cells in each spot
  cell_counts <- hashTable[[bregmaID]]$totalCells
  count_df <- do.call(rbind, lapply(names(cell_counts), function(x){
    coords <- strsplit(x, "_")[[1]]
    as.numeric(coords)
  }))
  colnames(count_df) <- c("x", "y")
  rownames(count_df) <- names(cell_counts)
  count_df <- as.data.frame(count_df)
  count_df$counts <- cell_counts
  
  # colors for each cell class
  classColors <- gg_color_hue(length(unique(df$Cell_class)))
  names(classColors) <- names(tmp_proportions)
  
  bregma <- list(sim = sim,
                 gtDocTopics = gtDocTopics,
                 gtTopicWords = gtTopicWords,
                 cellCounts = count_df,
                 cols = classColors,
                 df = df)
  
  return(bregma)
  
}


#' compute pvalue for topic correlation
#' based on a null distribution of randomly sampled sets of terms
#' 
#' null = ggplot transformed (long) dataframe
#' ex:
#' dat <- data.frame(draw = factor(rep(1:dim(termCorrNulls)[2], dim(termCorrNulls)[1])),
#'            numTerms = factor(rep(2:(dim(termCorrNulls)[1]+1), each = dim(termCorrNulls)[2])),
#'            corrs = as.vector(t(termCorrNulls)))
#' 
#' was a matrix (numTerms x samples)
#' so for 2:50 sizes of term sets, and 1000 random pairings to compute correlation on
#' based on gene counts across documents in corpus
#'            
#'     draw numTerms         corrs
#' 1      1        2  0.0590920266
#' 2      2        2 -0.4093004872
#' 3      3        2 -0.0964689721
#' 4      4        2  0.0803941294
#' 5      5        2  0.0127734074
#'
#' predictions = `topicTermCorrelation()` output data frame
#' 
#'            cor numTerms      zscore topic         pval
#' 1   0.25002099       21  2.93820765     1 0.0074220276
#' 2   0.29238589       27  4.83558037     2 0.0006501404
#' 3   0.16600159       34  0.96881060     3 0.1725689210
#' 4   0.27238757       30  4.71473737     4 0.0008689495
#' 5   0.11743073       21 -0.66700411     5 0.7307249799
#' 
#' procedure:
#' gets the null distribution terms for a given term set size
#' typically 1000 correlation values long.
#' Computes a probability density distribution
#' and integrates from the predicted correlation value
#' from the predictions data frame.
#' 
#' Thus returning a p-value for each predicted topic correlation
#' based on the null distribution of randomly sampled gene set correlations.
#'
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


#' get the correlation matrix between two topic-term betas
#' The betas must have the same terms (columns) in the same order
#'
#' Typically the first beta is the topic-term beta for the predicted
#' topics based on a fitted topic model. The second can be the same beta
#' (for assessing cross-correlation of topics) or a ground truth beta to
#' assess the concordance of predicted topics with known topics.
#'
#' thresh = can be used to only compare specific genes. Applied to select
#' genes (i.e., terms) in the first beta.
#' If a numeric fraction (0.01, let's say), then will be applied to only use
#' genes with a probability in each given topic above that threshold.
#' 
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


#' similar to `correlationBetweenBetas` but for the theta (doc-topic proportions)
#' for the thetas the topics are the columns instead of the rows like in the betas
#' 
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


#' adjust values to 0-1 range
#' x = can be vector or matrix
#' 
#' note that this will adjust based on all of the values in the object
#' so if it is a matrix, all the values are adjusted based on the highest
#' and lowest values in the entire matrix.
#' 
scale0_1 <- function(x) {
  xAdj <- (x - min(x)) / diff(range(x))
  return(xAdj)
}


#' compute distance between 2 matrices
#'
#'# frobenius norm of a matrix:
# sqrt(sum(M^2))
# or: norm(M, type = "F")
#'
computeDistance <- function(mtx1, mtx2, normType = "F") {
  # d <- abs(norm(mtx1, type = normType) - norm(mtx2, type = normType))
  d <- Frobenius(mtx1, mtx2)
  return(d)
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
    
    # print(cluster)
    # print(length(topics))
    
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
      # proportion of each topic-cluster that makes up each document.
      
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


#' visualize gene counts in spots in space. Also see group assignment of
#' spots.
#'
#' df = data.frame where rows are spots and columns must be at least:
#'      "x", "y" for spot positions in space
#'      "gene" counts of a gene for each spot
#'
#' gene = column name of the gene counts in df to be visualized
#' 
#' groups = a character vector of labels assigning spots to different groups.
#'          Ex: c("0", "1", "0", ...). 
#' group_cols = a vector designating the spot border color to be used for
#'              group assignment. Ex: c("0" = "gray", "1" = "red").
#'
vizGeneCounts <- function(df, gene,
                           groups = NA, group_cols = NA,
                           size = 7, stroke = 2,
                           plotTitle = NA,
                           showLegend = TRUE) {
  
  counts <- df[,gene]
  
  # color spots by group:
  if (is.na(groups[1]) == TRUE) {
    groups <- " "
    stroke <- 0.5
  } else {
    groups <- as.character(groups)
  }
  if (is.na(group_cols[1]) == TRUE) {
    group_cols <- c(" " = "gray")
  }
  
  p <- ggplot() +
    geom_point(data = df, aes(x=x, y=y, fill=counts, color = groups),
               shape = 21,
               stroke = stroke, size = size) +
    scale_fill_viridis(option = "A", direction = -1) +
    scale_color_manual(values = group_cols)
  
  if (showLegend == FALSE) {
    p <- p + guides(fill=FALSE)
  }
  if (is.na(plotTitle) == FALSE) {
    p <- p + ggtitle(plotTitle)
  }
  
  p <- p + theme_classic()
  print(p)
}


#' preprocess ST count matrices
#' 
#' input should be spot (row) x gene (columns) mtx with raw gene counts
#' 
#' align = optional 3x3 alignment file to adjust spot coordinates to match
#'         image pixels in order to properly align spots with image
#' nTopGenes = number of top expressed genes to remove
#' 
#' if the spot IDs are the positional information, extract:
#'  ex: "31x20" extracted as 31, 20 for x and y
#' extractPos = TRUE determines if this is done
#' 
#' perc.spots = genes with counts in this percent of spots are removed
#' 
#' min.reads = MERINGUE::cleanCounts param; minimum number of reads a gene needs
#' min.lib.size = MERINGUE::cleanCounts param; minimum number of counts a spot needs to have
#' 
#' od.genes.alpha MERINGUE::getOverdispersedGenes param
#' 
#' selected.genes = can supply a character vector of gene names to use specifically for the 
#'                  for the corpus.
#'                  
#' ODgenes = boolean whether to use OD genes or not
#' 
#' returns list with corpus (docs x genes) matrix,
#' corpus in slm format, and spot positions.
#' 
#' positions can be adjusted by the alignment matrices supplied if desired.
#'
#' genes.to.remove = character vector of gene names to remove, based on grepl regular expression
#'          ex: c("^mt-") or c("^MT", "^RPL", "^MRPL")
#'
preprocess <- function(dat, alignFile = NA, extractPos = TRUE,
                       selected.genes = NA, # can supply a vector of genes to specifically use for corpus.
                       nTopGenes = 5, # set to NA if do not want to remove top expressed genes
                       genes.to.remove = NA,
                       perc.spots = 1.0, # set to NA if do not want to remove genes based on this threshold
                       min.reads = 100,
                       min.lib.size = 100,
                       ODgenes = TRUE,
                       od.genes.alpha = 0.05,
                       gam.k = 5) {
  
  # if dat is a path to a file, read it, else assume matrix
  # definitely clean this up later..
  if (typeof(dat) == "character") {
    counts <- read.table(dat)
  } else {
    counts <- dat
  }
  
  # remove poor spots and genes
  countsClean <- MERINGUE::cleanCounts(counts = t(counts), 
                                       min.reads = min.reads, 
                                       min.lib.size = min.lib.size, 
                                       plot=TRUE,
                                       verbose=TRUE)
  
  if (is.na(selected.genes) == FALSE) {
    countsClean <- countsClean[selected.genes,]
  }
  
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
  
  # adjust spot coordinates, optional.
  # (based on Stahl ST data with alignment matrices)
  if (is.na(alignFile) == FALSE) {
    align <- matrix(unlist(read.table(alignFile)), nrow = 3, ncol = 3)
    (positions[,"x"] * align[1,1]) - 290
    (positions[,"y"] * align[2,2]) - 290
  }
  
  # remove top expressed genes (nTopGenes needs to be integer or NA)
  if (is.na(nTopGenes) == FALSE) {
    top_expressed <- names(rowSums(countsClean)[order(rowSums(countsClean),
                                                      decreasing = TRUE)][1:nTopGenes])
    countsClean <- countsClean[rownames(countsClean) %in% top_expressed == FALSE,]
    
    print(paste("after removing top ",
                as.character(nTopGenes),
                " genes:", dim(countsClean)[1], dim(countsClean)[2]))
  }
  
  # remove specific genes (if there are any)
  if (is.na(genes.to.remove) == FALSE) {
    countsClean <- countsClean[!grepl(paste(genes.to.remove, collapse="|"), rownames(countsClean)),]
    print(paste("after removing selected genes:", dim(countsClean)[1], dim(countsClean)[2]))
  }
  
  # remove genes that appear in every document
  if (is.na(perc.spots) == FALSE){
    countsClean_ <- countsClean
    countsClean_[which(countsClean_ > 0)] <- 1
    numberSpots <- (perc.spots * ncol(countsClean_))
    countsClean <- countsClean[which(rowSums(countsClean_) < numberSpots),]
    print(paste("after removing genes present in more than ", 
                as.character(perc.spots*100), "% of spots: ",
                dim(countsClean)[1], dim(countsClean)[2]))
  }
  
  # get variable genes
  if (ODgenes == TRUE) {
    par(mfrow=c(4,2), mar=c(1,1,1,1))
    OD <- getOverdispersedGenes(countsClean,
                                alpha = od.genes.alpha,
                                gam.k = gam.k,
                                plot = TRUE,
                                details = TRUE)
    
    countsClean <- countsClean[OD$ods,]
  }
  
  corpus <- t(as.matrix(countsClean))
  corpus_slm <- slam::as.simple_triplet_matrix(corpus)
  
  return(list(corpus = corpus,
              slm = corpus_slm,
              pos = positions))
  
}


#' transform corpus raw gene counts into binned counts.
#' Genes are binned based on their raw counts wrt the entire corpus
#' and the bin is used as the new gene count value
#'
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

#' return the beta and theta matrices for a given lda model
#'
getBetaTheta <- function(lda) {
  
  lda.tmResult <- posterior(lda)
  theta <- lda.tmResult$topics
  beta <- lda.tmResult$terms
  topicFreqsOverall <- colSums(theta) / length(lda@documents)
  
  return(list(beta = beta,
              theta = theta,
              topicFreq = topicFreqsOverall))
  
}


#' wrapper for functions to extract beta and theta and cluster topics
#' for a LDA model
#' 
#' LDAmodel is an LDA model from `fitLDA` list.
#' can use: optimalModel(fitLDAoutput) to input the optimal model
#' 
#' color schemes for coloring topics:
#' "rainbow", "ggplot"
#'
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
  m$betaCombn <- combineTopics(m$beta, clusters = m$clusters, mtxType = "b")
  m$thetaCombn <- combineTopics(m$theta, clusters = m$clusters, mtxType = "t")
  
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


#' function to properly construct input ST count matrix
#' from Moncada et al GSE111672
#' 
#' returns spot x gene matrix where row spot names are the x by y spatial coordinates
#' such that extractPos option in `preprocess` can be set to TRUE
#' 
#' path is path to file
#'
loadPDACfile <- function(path) {
  
  dat <- read.table(path, header = TRUE, sep = "\t")
  
  spotids <- unlist(lapply(colnames(dat)[2:ncol(dat)], function(i) {
    ix <- strsplit(i, "X")[[1]][2]
    ix
  }))
  
  t <- t(dat)
  colnames(t) <- dat$Genes
  t <- t[2:nrow(t),]
  dat <- apply(t, 2,FUN = as.numeric)
  rownames(dat) <- spotids
  
  return(dat)
  
}


#' function to get Hungarian sort pairs via clue::lsat
#' 
#' pairs each row with a column
#' So must have equal or less rows than number of columns.
#' If less rows, then some columns will not be matched and not part of
#' the output.
#' 
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



#' attempt to optimize search time for optimal K
#' 
#' based on optimal search algorithm to split data into intervals
#' and choosing the interval of k's that produce the lower perplexity.
#' Idea is to not have to fit a model to every k
#'
#' interval = c(i, j) where i and j indicate the lower and upper bounds of a
#' sequence of integer K's to search
#'
fitLDA2 <- function(interval, corpus = corpus,
                    K = c(), Plx = c(), models = c(),
                    lower = min(interval),
                    mid = round((min(interval) + max(interval))/2),
                    upper = max(interval),
                    tol = .Machine$double.eps^0.25,
                    seed = 0) {
  
  # params for fitting model
  params <- list(seed = seed,
                 verbose = 0, keep = 1, estimate.alpha = TRUE)
  
  print("begin")
  cat("interval:", lower, mid, upper, "\n")
  
  # if optimal found:
  if(lower == mid | mid == upper) {
    # final adjustment in case the optimum was actually at end of the interval
    # because when computing mid, it rounds. 
    # ex: for a final interval 6,7, it rounds the mid to 6 and ends search, labeling 6 as optimal
    # although 7 may have actually been the lower value
    kopt <- K[which.min(Plx)]
    
    cat("optimal K estimated to be", kopt, "\n")
    
    # plot
    oix <- order(K) # order by K's
    kopt_ix <- which(K[oix] %in% kopt) # index of reordered k's that is optimal
    
    plot.new()
    plot(K[oix], Plx[oix], type="l")
    points(kopt, Plx[oix][kopt_ix], col = "red")
    
    return(list(models = models,
                kOpt = kopt,
                perplexities = Plx,
                Ks = K,
                fitCorpus = corpus))
  }
  
  # 1. check if perplexities have already been computed for k's.
  
  # if not, compute them.
  # Use mclapply to get values at least for first iteration when there will be multiple k's
  ks <- c(lower, mid, upper)
  notComputed <- ks[which(!ks %in% K)]
  cat("computing perplexities for:", notComputed, "\n")
  
  # returns list of lists, where each list is: [[1]]perplexity, [[2]]k, [[3]]LDAmodel
  fitmodels <- mclapply(notComputed, function(k){
    model <- topicmodels::LDA(corpus, k=k, control = params)
    p <- topicmodels::perplexity(model, corpus)
    return(c(p, k, model))
  }, mc.cores = detectCores(logical = TRUE) - 1)
  
  # append newly obtained perplexities, corresponding k's, and models. Need to be in same order
  for (i in seq(length(fitmodels))) {
    m <- fitmodels[[i]]
    Plx <- c(Plx, m[[1]])
    K <- c(K, m[[2]])
    models <- c(models, m[[3]])
  }
  # can't get these global lists to update
  # lapply(fitmodels, function(m){
  #   perp_vals <- c(perp_vals, m[1])
  #   k_coords <- c(k_coords, m[2])
  # })
  
  # 2. determine new interval
  
  # get perplexity values for the current k interval values (lower, mid, upper)
  plxs <- sapply(ks, function(k) {
    # get index of k and perplexity value for each element of current interval `ks`
    ix <- which(K %in% c(k))
    p <- Plx[ix]
    p
  })
  
  # perplexity values for the current interval `ks`
  fl <- plxs[1]
  fm <- plxs[2]
  fu <- plxs[3]
  
  cat("lower k", lower, "\n")
  cat("perplexity:", fl, "\n")
  cat("mid k", mid, "\n")
  cat("perplexity", fm, "\n")
  cat("upper k", upper, "\n")
  cat("perplexity", fu, "\n")
  
  if (fl < fu) {
    print('searching lower half')
    fitLDA2(interval = c(lower, mid),
            corpus = corpus,
            K = K,
            Plx = Plx,
            models = models)
  } else {
    print('searching upper half')
    fitLDA2(interval = c(mid, upper),
            corpus = corpus,
            K = K,
            Plx = Plx,
            models = models)
  }
}



#' Another alternative to optimizing search for K
#' 
#' fit a small number of K's evenly spread out along a chosen interval of
#' positive integer K's and fit models to them.
#' 
#' Then fit a model to the K's and resulting perplexities and try to
#' predict which K gives the lower perplexity.
#' 
#' With mclapply parallelization, the initial starting K's have all be fit at the same time.
#' Then, the predicted optimal K can be fitted if it is not one of the initial Ks.
#' Otherwise you're done.
#' 
#' This may be useful because the limiting step for every fitting is that you must wait for the
#' model with the largest K to be fit, which takes the longest.
#' 
#' So do this first and in parallel with the smaller initial K's but few enough K's such that
#' the largest K will be assigned to a thread from the get-go
#' (instead of waiting for a thread to free up in the original `fitLDA`)
#' 
#' Any extra time will be from fitting one more K at the end if it is different from one of the initial K's
#' but will still be less time than what was needed for the largest K in the interval range chosen.
#' 
#' This way it has the potential to be faster than `fitLDA`, but also possible to be a little bit slower
#' But `fitLDA2`, moving to each segment cannot be threaded, so you are iterating through each segment
#' and this cannot be parallelized. So if you keep moving in direction of larger K's, this will be ver slow.
#' 
#' So this version may be the best compromise. 
#' 
#' p = number of perplexities to solve for (k's to fit) initially
#' 
#' 
fitLDA3 <- function(interval, p = 5, corpus, seed = 0){
  
  # # sequence across interval
  # s <- seq(min(interval), max(interval))
  # # split sequence into p segments
  # segs <- split(s, cut(s, p))
  # 
  # # for each segments, take first element,
  # # or if last one, take last elements
  # # selected elements are the Ks to fit
  # Ks <- sapply(seq(p), function(seg){
  #   if (seg == tail(seq(p), 1)) {
  #     k <- max(segs[[seg]])
  #   } else {
  #     k <- min(segs[[seg]])
  #   }
  #   k
  # })
  
  Ks <- round(seq(min(interval), max(interval), length = p))
  
  # params for fitting model
  params <- list(seed = seed,
                 verbose = 0, keep = 1, estimate.alpha = TRUE)
  
  # fit model and compute perplexity for each k
  # returns list of lists, where each list is: [[1]]perplexity, [[2]]LDAmodel
  cat("Training models for following K's:", Ks, "\n")
  fitmodels <- mclapply(Ks, function(k){
    model <- topicmodels::LDA(corpus, k=k, control = params)
    p <- topicmodels::perplexity(model, corpus)
    return(c(p, model))
  }, mc.cores = detectCores(logical = TRUE) - 1)
  
  plxs <- c() # y coords for model fitting
  models <- c()
  for (i in seq(length(fitmodels))) {
    m <- fitmodels[[i]]
    plxs <- c(plxs, m[[1]])
    models <- c(models, m[[2]])
  }
  
  # fit models and pick best one
  y <- plxs
  x <- Ks
  # first degree 
  fit  <- lm(y~x)
  # second degree
  fit2 <- lm(y~poly(x,2,raw=TRUE))
  # third degree
  fit3 <- lm(y~poly(x,3,raw=TRUE))
  # fourth degree
  fit4 <- lm(y~poly(x,4,raw=TRUE))
  
  fits <- list(fit, fit2, fit3, fit4)
  
  xx <- seq(min(Ks), max(Ks), length=50)
  plot(x, y, pch=19)
  lines(xx, predict(fit, data.frame(x=xx)), col="red")
  lines(xx, predict(fit2, data.frame(x=xx)), col="green")
  lines(xx, predict(fit3, data.frame(x=xx)), col="blue")
  lines(xx, predict(fit4, data.frame(x=xx)), col="purple")
  
  # pick best fit by highest r.squared value:
  rsqds <- sapply(fits, function(f){
    s <- summary(f)
    r2 <- s[["r.squared"]]
    r2
    #s
  })
  best <- fits[[which.max(rsqds)]]
  
  # find lowest predicted K
  xx <- seq(min(Ks), max(Ks))
  predictions <- predict(best, data.frame(x=xx))
  optK <- xx[which.min(as.vector(predictions))]
  
  # if optK not one of the initial sampled Ks, then a model has not been fitted yet
  # so fit a model:
  if (!optK %in% Ks){
    cat("Training putative optimal K", optK, "\n")
    model <- topicmodels::LDA(corpus, k=optK, control = params)
    p <- topicmodels::perplexity(model, corpus)
    
    Ks <- c(Ks, optK)
    plxs <- c(plxs, p)
    models <- c(models, model)
  }
  
  # final adjustment in case one of the sampled Ks was actually better
  # than the predicted
  optK <- Ks[which.min(plxs)]
  cat("Determined optimal K:", optK, "\n")
  
  # plot
  oix <- order(Ks) # order by K's
  kopt_ix <- which(Ks[oix] %in% optK) # index of reordered k's that is optimal
  
  plot(Ks[oix], plxs[oix], type="l")
  points(optK, plxs[oix][kopt_ix], col = "red")
  
  return(list(models = models,
              kOpt = optK,
              perplexities = plxs,
              Ks = Ks,
              fitCorpus = corpus))
  
}



#' nmfRef = the list from SPOTlight::train_nmf()
#' stCounts = ST count matrix, genes x spots
#'
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
  spotlight_ct_topics <- do.call(cbind, lapply(seq(ncol(ct_topic_profiles_r)), function(i){
    ct <- ct_topic_profiles_r[,i]
    ct/sum(ct)
  }))
  rownames(spotlight_ct_topics) <- paste("Topic", 1:nrow(ct_topic_profiles_r), sep = "_")
  colnames(spotlight_ct_topics) <- colnames(ct_topic_profiles_r)
  
  # now that we have the proportion of each topic for a cell type,
  # use this information to get gene weights (sum to 1) for each cell type
  # (w is gene weights for each topic)
  # do this by getting the topic-gene vectors from w for each cell type and
  # merging them based on the topic proportions for the cell type
  spotlight_beta <- do.call(rbind, lapply(seq(ncol(spotlight_ct_topics)), function(ix){
    # topics for the given cell type in spotlight_ct_topics
    ts <- spotlight_ct_topics[,ix]
    # topics that make up the cell type and their proportions
    ts <- ts[which(ts > 0)]
    # w topics and the gene associations, multiplied by the topic proportions of the cell type
    g <- w[,names(ts)]
    if (length(ts) > 1){
      g <- g %*% diag(ts) # mtx multi. to multiply first column by first element, second col by second element, etc
    }
    
    # if multiple topics, add the adjusted topic gene vectors together
    if (is.null(dim(g)) == FALSE){
      g <- rowSums(g)
    }
    
    # final adjustment to make topic gene proportions sum to 1
    g <- g/sum(g)
    g
  }))
  rownames(spotlight_beta) <- colnames(spotlight_ct_topics)
  
  
  # [topics x spots]; topic coefficients for each spot
  # This function appears to account for that only using genes in ST that are also in W
  # NNLS to get topic coefficients for each spot based on genes associated with each topic (W)
  profile_mtrx <- predict_spatial_mixtures_nmf(nmf_mod = nmfRef[[1]],
                                               mixture_transcriptome = stCounts,
                                               transf = "uv")
  # normalize to get topic proportions in each spot
  spotlight_spot_topics <- do.call(cbind, lapply(seq(ncol(profile_mtrx)), function(i){
    ct <- profile_mtrx[,i]
    ct/sum(ct)
  }))
  rownames(spotlight_spot_topics) <- paste("Topic", 1:nrow(profile_mtrx), sep = "_")
  colnames(spotlight_spot_topics) <- colnames(profile_mtrx)
  
  
  # returns [spot x CellType]; proportion of each cell type in each spot
  # topics below `min_cont` are not counted in a spot
  # within this function, `profile_mtrx` is also made but not returned
  decon_mtrx <- mixture_deconvolution_nmf(nmf_mod = nmfRef[[1]],
                                          mixture_transcriptome = stCounts,
                                          transf = "uv", 
                                          reference_profiles = ct_topic_profiles, 
                                          min_cont = 0.0) # only keep topics if 9% or more in a spot
  # note that last column is an additional columns for the residual error
  # cleanup to get actual spot-celltype predictions:
  rownames(decon_mtrx) <- colnames(stCounts)
  spotlight <- decon_mtrx[,1:(ncol(decon_mtrx)-1)] # last column is residuals
  # note: all spots and cell types.
  # Possible some cell type will not be detected at all and be 0 for all spots
  
  return(list(topicBeta = t(w), # topic x gene weights that sum to 1 (all genes and topics)
              ctBeta = spotlight_beta, # CellType x gene weights that sum to 1 (all genes and topics)
              ctTopicProps = spotlight_ct_topics, # topic proportion for each cell type
              spotTopicTheta = t(spotlight_spot_topics), # proportion of each topic in each spot (all topics)
              spotCtTheta = spotlight)) # spot x celltype proportions (all spots and cell types)
  
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


#' Get filtered spotCtTheta, spotTopicTheta, ctBeta, topicBeta sl matrices
#' wrt a given lda model such that they can be compared equally in terms
#' of shared genes, shared spots, and cell types that were detected
#' to be present in the sl predictions. Also generated colors for the
#' cell type clusters and topics of the sl object.
#'
#'
filterSPOTlightMtxs <- function(sl, lda_model){
  
  # Thetas:
  # make sure the spots are the same
  spotCtThetaFilt <- sl$spotCtTheta[rownames(lda_model$theta),]
  spotTopicThetaFilt <- sl$spotTopicTheta[rownames(lda_model$theta),]
  
  # drop any clusters that were not detected at all in the ST data
  spotCtThetaFilt <- spotCtThetaFilt[,which(!colSums(spotCtThetaFilt) == 0)]
  spotTopicThetaFilt <- spotTopicThetaFilt[,which(!colSums(spotTopicThetaFilt) == 0)]
  
  # Betas:
  # make sure genes are consistent between the topic model lda beta
  keep_genes <- colnames(sl$ctBeta)[colnames(sl$ctBeta) %in% colnames(lda_model$beta)]
  
  # also filter for cell types that were detected in the spots 
  ctBetaFilt <- sl$ctBeta[colnames(spotCtThetaFilt),keep_genes]
  
  # also filter for topics that were detected in the spots 
  topicBetaFilt <- sl$topicBeta[colnames(spotTopicThetaFilt),keep_genes]
  
  # Before subetting genes, topicBeta was just W. ctBeta combined W topic genes vectors
  # proportional to topics in each cell type). So gene values would add to 1.
  # But after subsetting to have same genes as corpus to compare, the genes no longer add to 1.
  # So adjust such that gene values sum to 1
  ctBetaFilt <- ctBetaFilt/rowSums(ctBetaFilt)
  topicBetaFilt <- topicBetaFilt/rowSums(topicBetaFilt)
  
  # colors for the predicted and filtered cell types
  ctcols <- as.factor(colnames(spotCtThetaFilt))
  names(ctcols) <- colnames(spotCtThetaFilt)
  levels(ctcols) <- gg_color_hue(length(colnames(spotCtThetaFilt)))
  
  # colors for the topics
  topiccols <- as.factor(colnames(spotTopicThetaFilt))
  names(topiccols) <- colnames(spotTopicThetaFilt)
  levels(topiccols) <- gg_color_hue(length(colnames(spotTopicThetaFilt)))
  
  return(list(ctThetaFilt = spotCtThetaFilt,
              topicThetaFilt = spotTopicThetaFilt,
              ctBetaFilt = ctBetaFilt,
              topicBetaFilt = topicBetaFilt,
              ctcols = ctcols,
              topiccols = topiccols,
              shared_genes = keep_genes))
  
}






