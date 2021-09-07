#' Visualize all topic proportions across pixels with `scatterpie`
#'
#' @description Note: visualizes all cell-types in theta at once (could be
#'     individual cell-types or cell-type-clusters) so for accuracy of the proportions
#'     of each cell-type in a pixel, the row (pixel) should sum to 1.
#'
#' @param theta document (pixel) x cell-type proportion matrix
#' @param pos position of pixels, as data.frame with `x` and `y` columns
#' @param topicOrder order of topics in theta to visualize as a numeric vector
#'     and same length as topicCols (default: seq(ncol(theta)))
#' @param topicCols Vector of colors for each of the cell-types to be visualized.
#'     Same length and order as topicOrder (default: rainbow(ncol(theta)))
#' @param groups Indicates color of the scatterpie strokes with the goal of coloring them
#'     by their assigned group. This can be a vector or factor indicating the group of each
#'     scatterpie. Needs to be in the same order as the pixel rows in "theta" (default: NA)
#' @param group_cols Color labels for the groups. Can be a vector or factor. (default: NA)
#' @param r Radius of the scatterpie circles. Adjust based on positions of pixels (default: max(0.4, max(pos)/nrow(pos)*4))
#' @param lwd Width of lines of the pie charts. Increasing helps visualize
#'     group_cols if being used.
#' @param showLegend Boolean to show the legend indicating cell-types and their color
#' @param plotTitle add title to the resulting plot (default: NA)
#' @param overlay plot the scatterpies on top of a raster image of the H&E tissue
#'     (default: NA)
#'
#' @export
vizAllTopics <- function(theta, pos,
                         topicOrder=seq(ncol(theta)),
                         topicCols=rainbow(ncol(theta)),
                         groups = NA,
                         group_cols = NA,
                         r = max(0.4, max(pos)/nrow(pos)*4),
                         lwd = 0.5,
                         showLegend = TRUE,
                         plotTitle = NA,
                         overlay = NA) {
  
  # pixel cell-type distribution reordered based on topicOrder
  theta_ordered <- theta[, topicOrder]
  theta_ordered <- as.data.frame(theta_ordered)
  # colnames(theta_ordered) <- paste0("Topic.", topicOrder)
  colnames(theta_ordered) <- paste0("Topic.", colnames(theta_ordered))
  
  # add columns "x", "y" with document positions from `pos`
  theta_ordered_pos <- merge(data.frame(theta_ordered),
                             data.frame(pos), by=0)
  
  # column names of cell-types, in order of topicOrder
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
  
  if (is.na(overlay[1]) == FALSE){
    p <- ggplot2::ggplot(mapping = ggplot2::aes(x = 0:dim(overlay)[2], y = 0:dim(overlay)[1])) +
      ggplot2::coord_equal(xlim = c(0,dim(overlay)[2]), ylim = c(0, dim(overlay)[1]), expand = FALSE) +
      ggplot2::theme(
        #panel.background = element_rect(fill = "white"),
        panel.grid = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.background= ggplot2::element_blank()) +
      # geom_point(aes(x = c(0,dim(overlay)[2]), y = c(0, dim(overlay)[1]))) +
      ggplot2::annotation_raster(overlay, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      # theme_classic() +
      scatterpie::geom_scatterpie(ggplot2::aes(x=x, y=y, group=Row.names, r=r, color = groups),
                                  lwd = lwd,
                                  data = theta_ordered_pos,
                                  cols = topicColumns,
                                  legend_name = "Topics") +
      ggplot2::scale_fill_manual(values = as.vector(topicCols)) +
      ggplot2::scale_color_manual(values = group_cols)
  } else {
    p <- ggplot2::ggplot() +
      ggplot2::theme(
        #panel.background = element_rect(fill = "white"),
        panel.grid = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank()) +
      # theme_classic() +
      scatterpie::geom_scatterpie(ggplot2::aes(x=x, y=y, group=Row.names, r=r, color = groups),
                                  lwd = lwd,
                                  data = theta_ordered_pos,
                                  cols = topicColumns,
                                  legend_name = "Topics") +
      ggplot2::scale_fill_manual(values = as.vector(topicCols)) +
      ggplot2::scale_color_manual(values = group_cols)
  }
  
  if (showLegend == FALSE) {
    # p <- p + ggplot2::guides(fill=FALSE)
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  if (is.na(plotTitle) == FALSE) {
    p <- p + ggplot2::ggtitle(plotTitle)
  }
  
  return(p)
}


#' Visualize proportions of cell-types or aggregated cell-type-clusters individually
#' 
#' @description Similar to `vizAllTopics` but will generate a separate plot for
#'     each cell-type or cell-type-cluster where the other cell-types or clusters will be
#'     colored gray. In this way, the actual proportions of each cell-type in a pixel
#'     will be maintained such that pixel cell-type proportions still sum to 1.
#'
#' @param theta document (pixel) x cell-type proportion matrix
#' @param pos position of pixels, as data.frame with `x` and `y` columns
#' @param clusters factor of colors that each cluster (i.e., cell-type-cluster) is
#'     assigned to. In this case, the levels should be colors. In `vizAllTopics`,
#'     clusters is "topicCols" and can just be a vector of colors.
#' @param sharedCol Boolean indicating if the cell-types in a cluster will be plotted with
#'     the same color or if each celll-type will be colored by its own shade to also
#'     show how the cell-types in a cluster are distributed in space wrt each other.
#' @param groups colors the pixel scatterpie lines based on a group or cell layer
#'     they belong to. Needs to be a character vector in the order of the spot
#'     rows in theta. Ex: c("0", "1", "0", ...)
#' @param group_cols color labels for the groups. Ex: c("0" = "gray", "1" = "red")
#' @param r = radius of the scatterpies. Adjust based on position of pixels (default: 1)
#' @param lwd = width of lines of the scatterpies. Increasing helps visualize
#'     group_cols if being used.
#' @param showLegend Boolean to show the legend indicating topics and their color
#' @param plotTitle add title to the resulting plot (default: NA)
#' @param overlay plot the scatterpies on top of a raster image of the H&E tissue
#'     (default: NA)
#' @param fig_path path so save output figures for each plotted cluster (not in use)
#' @param fig_prefix prefix to name each output figure for each plotted cluster (not in use)
#'     
#' @export
vizTopicClusters <- function(theta, pos, clusters,
                             sharedCol = FALSE,
                             groups = NA,
                             group_cols = NA,
                             r = 1,
                             lwd = 0.5,
                             showLegend = TRUE,
                             plotTitle = NA,
                             overlay = NA,
                             fig_path = "./",
                             fig_prefix = NA) {
  
  print("Topic cluster members:")
  # produce a plot for each cell-type-cluster:
  for (cluster in levels(clusters)) {
    
    # select the cell-types in the cluster
    topics <- labels(clusters[which(clusters == cluster)])
    
    cat(cluster, ":", topics, "\n")
    
    # pixel cell-type distribution reordered based on topicOrder and selected cluster cell-types
    theta_ordered <- theta[, topics]
    
    # get percentage of other cell-types not in cluster to maintain actual
    # pixel cell-type proportions that sum to 1
    if (is.null(dim(theta_ordered))) {
      other <- 1 - theta_ordered
    } else {
      other <- 1 - rowSums(theta_ordered)
    }
    
    theta_ordered <- as.data.frame(theta_ordered)
    colnames(theta_ordered) <- paste0("Topic.", topics)
    theta_ordered$other <- other
    
    # if any cell-types not represented at all, drop them
    # Apparently if proportion of a  cell-type is 0 for all pixels, it is not plotted
    # and doesn't appear in the legend and messes with the colors such that
    # "other" takes one of the colors of the cell-types and is not gray
    if ( length(which(colSums(theta_ordered) == 0)) > 0 ) {
      missing_topics <- colnames(theta_ordered)[which(colSums(theta_ordered) == 0)]
      cat("NOTE:", missing_topics, "not present in any pixels and will be dropped.", "\n")
      theta_ordered <- theta_ordered[,which(!colSums(theta_ordered) == 0)]
      # if the entire cluster/topic is represented in no spots, skip
      if (is.null(dim(theta_ordered))){
        cat("No pixels contain this cell-type. Skipping", "\n")
        next
      }
    }
    
    # add columns with pixel positions
    rownames(theta_ordered) <- rownames(pos)
    theta_ordered_pos <- merge(data.frame(theta_ordered),
                               data.frame(pos), by=0)
    
    # first column after merge is "Row.names", last two are "x" and "y"
    # problem is that data frame will replace "-" and " " with "."
    topicColumns <- colnames(theta_ordered_pos)[2:(dim(theta_ordered_pos)[2]-2)]
    
    # get a hue of colors for each cell-type in cell-type-cluster
    if (sharedCol){
      color_ramp <- grDevices::colorRampPalette(c(cluster, cluster))
    } else {
      color_ramp <- grDevices::colorRampPalette(c(lighten(cluster, factor = 0.5), darken(cluster, factor = 2)))
    }
    
    topic_colors <- color_ramp(ncol(theta_ordered) - 1) # don't count "other" here
    # topic_colors <- append(topic_colors, c("gray")) # add gray to other here
    topic_colors <- append(topic_colors, c(transparentCol("white", percent = 60)))
    
    # color of scatterpie groups (lines of scatterpies):
    if (is.na(groups[1]) == TRUE) {
      groups <- rep("0", dim(theta_ordered_pos)[1])
      theta_ordered_pos$groups <- groups
    } else {
      theta_ordered_pos$groups <- as.character(groups)
    }
    if (is.na(group_cols[1]) == TRUE) {
      group_cols <- c("0" = "gray")
      # group_cols <- c("0" = transparentCol("white", percent = 90))
    }
    
    if (is.na(overlay[1]) == FALSE){
      p <- ggplot2::ggplot(mapping = ggplot2::aes(x = 0:dim(overlay)[2], y = 0:dim(overlay)[1])) +
        ggplot2::coord_equal(xlim = c(0,dim(overlay)[2]), ylim = c(0, dim(overlay)[1]), expand = FALSE) +
        ggplot2::theme(
          #panel.background = element_rect(fill = "white"),
          panel.grid = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank()) +
        # geom_point(aes(x = c(0,dim(overlay)[2]), y = c(0, dim(overlay)[1]))) +
        ggplot2::annotation_raster(overlay, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
        # theme_classic() +
        scatterpie::geom_scatterpie(ggplot2::aes(x=x, y=y, group=Row.names, r = r, color = groups),
                                    lwd = lwd,
                                    data=theta_ordered_pos,
                                    cols = topicColumns,
                                    legend_name = "Topics") +
        # coord_equal() +
        ggplot2::scale_fill_manual(values=topic_colors) +
        ggplot2::scale_color_manual(values = group_cols)
    } else {
      p <- ggplot2::ggplot() +
        ggplot2::theme(
          #panel.background = element_rect(fill = "white"),
          panel.grid = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank()) +
        # theme_classic() +
        scatterpie::geom_scatterpie(ggplot2::aes(x=x, y=y, group=Row.names, r = r, color = groups),
                                    lwd = lwd,
                                    data=theta_ordered_pos,
                                    cols = topicColumns,
                                    legend_name = "Topics") +
        # coord_equal() +
        ggplot2::scale_fill_manual(values=topic_colors) +
        ggplot2::scale_color_manual(values = group_cols)
    }
    
    if (showLegend == FALSE) {
      # p <- p + ggplot2::guides(fill=FALSE)
      p <- p + ggplot2::theme(legend.position = "none")
    }
    
    if (is.na(plotTitle) == FALSE) {
      p <- p + ggplot2::ggtitle(plotTitle)
    }
    
    if (is.na(fig_prefix) == FALSE) {
      fig_name <- paste0(fig_prefix, "_", topics, ".pdf")
    } else {
      fig_name <- paste0(topics, ".pdf")
    }
    
    print(p)
    # ggsave(filename = fig_name,
    #        device = "pdf",
    #        path = fig_path,
    #        scale = 1.5,
    #        width = 5,
    #        height = 4,
    #        units = c("in"),
    #        dpi = 600)
  }
}


#' Visualize gene counts in pixels in space. Can also see group assignment of
#' spots.
#' 
#' @description Note: visualized one gene at a time. Can set up a loop to plot
#'    a different gene column in df individually.
#'
#' @param df data.frame where rows are spots and columns must be at least:
#'      "x", "y" for spot positions in space and "gene" column that is counts
#'      of a gene for each spot.
#' @param gene column name of the gene counts in df to be visualized
#' @param groups colors the spots lines based on a group or cell layer
#'     they belong to. Needs to be a character vector in the order of the spot
#'     rows in df. Ex: c("0", "1", "0", ...)
#' @param group_cols color labels for the groups. Ex: c("0" = "gray", "1" = "red")
#' @param winsorize Winsorization quantile
#' @param size size of the geom_points to plot (default: 7)
#' @param stroke thickness of the geom_point lines to help in emphasizing groups
#'     (default: 2)
#' @param alpha alpha value of colored pixels (default: 1)
#' @param plotTitle option to add a title to the plot
#' @param showLegend Boolean to show the plot legend
#' 
#' @export
vizGeneCounts <- function(df, gene,
                          groups = NA,
                          group_cols = NA,
                          winsorize = 0,
                          size = 7, stroke = 0.5,
                          alpha = 1,
                          plotTitle = NA,
                          showLegend = TRUE) {
  
  counts <- df[,gene]
  
  ## winsorize
  counts <- winsorize(counts, qt=winsorize)
  
  # color spots by group:
  if (is.na(groups[1]) == TRUE) {
    groups <- " "
    # stroke <- 0.5
  } else {
    groups <- as.character(groups)
  }
  if (is.na(group_cols[1]) == TRUE) {
    group_cols <- c(" " = "white")
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = df, ggplot2::aes(x=x, y=y, fill=counts, color = groups),
               shape = 21,
               stroke = stroke, size = size, 
               alpha = alpha) +
    viridis::scale_fill_viridis(option = "A", direction = -1) +
    ggplot2::scale_color_manual(values = group_cols)
  
  p <- p +
    ggplot2::theme(
      #panel.background = element_rect(fill = "white"),
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank())
  # theme_classic()
  
  if (showLegend == FALSE) {
    # p <- p + ggplot2::guides(fill=FALSE)
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  if (is.na(plotTitle) == FALSE) {
    p <- p + ggplot2::ggtitle(plotTitle)
  }
  
  return(p)
}


# custom correlation color range for heatmap.2 correlation plots
correlation_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(n = 209)
correlation_breaks <- c(seq(-1,-0.01,length=100),
                       seq(-0.009,0.009,length=10),
                       seq(0.01,1,length=100))


#' lighten color
#' 
#' @description lighten a color
#'
#' @param color color
#' @param factor how much to lighten (default: 0.5)
#' 
#' @export
lighten <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- col2rgb(color)
  col <- col + (255 - col)*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}


#' darken color
#' 
#' @description darken a color
#'
#' @param color color
#' @param factor how much to darken (default: 1.4)
#' 
#' @export
darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}


#' color palette to replicate ggplot2
#' 
#' @description get a ggplot2 color palette of length n
#'
#' @param n length of color palette
#' 
#' @export
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



#' transparent version of color.
#' 
#' @description adjust transparency of a color
#'
#' @param color color name
#' @param percent % transparency (default: 50)
#' @param name an optional name for the color (default: NULL)
#' 
#' @export
transparentCol <- function(color, percent = 50, name = NULL) {
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}


#' Generate heatmap of correlations
#' 
#' @description Visualize the correlations between topics stored in a matrix, typically one
#'     returned via `getCorrMtx()`. This function uses ggplot2::geom_tile.
#'
#' @param mat matrix with correlation values from -1 to 1
#' @param rowLabs y-axis label for plot. These are the rows of the matrix, or specifically m1 from getCorrMtx. (default: NULL)
#' @param colLabs x-axis label for plot. These are the columns of the matrix, or specifically m2 from getCorrMtx. (default: NULL)
#' @param title title of the plot. (default: NULL)
#' 
#' @export
correlationPlot <- function(mat, colLabs = NA, rowLabs = NA, title = NA){
  
  correlation_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(n = 209)
  correlation_breaks <- c(seq(-1,-0.01,length=100),
                          seq(-0.009,0.009,length=10),
                          seq(0.01,1,length=100))
  
  dat <- reshape2::melt(mat)
  plt <- ggplot2::ggplot(data = dat) +
    ggplot2::geom_tile(ggplot2::aes(x = Var1, y = Var2, fill=value)) +
    ggplot2::scale_fill_gradientn(colors = correlation_palette, breaks = correlation_breaks, limits = c(-1,1),
                                  guide = ggplot2::guide_colorbar(title = "correlation", ticks = FALSE, label = FALSE)) +
    ggplot2::scale_y_discrete(breaks = dat$Var2, labels = dat$Var2) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -90, size = 10),
                   axis.text.y = ggplot2::element_text(angle = 0, size = 10),
                   axis.title.y = ggplot2::element_text(size = 12),
                   axis.title.x = ggplot2::element_text(size = 12),
                   plot.title = ggplot2::element_text(size = 15),
                   panel.grid = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank()
                   # legend.position = "bottom"
    ) +
    # ggplot2::guides(color = guide_colourbar(nbin = 10, raster = F, direction = "horizontal", barwidth = 20, barheight = 0.8)) +
    ggplot2::coord_fixed()
  
  if (is.na(colLabs) == FALSE){
    plt <- plt + ggplot2::xlab(colLabs)
  }
  if (is.na(rowLabs) == FALSE){
    plt <- plt + ggplot2::ylab(rowLabs)
  }
  if (is.na(title) == FALSE){
    plt <- plt + ggplot2::ggtitle(title)
  }
  
  return(plt)
  
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
#' @export
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


#' Plot the perplexity and rare cell-types versus fitted Ks
#' 
#' @description the same plot returned by fitLDA() but now callable as a 
#'     separate function. 
#'     
#' @param models models output from fitLDA()
#' @param corpus If corpus is NULL, then it will use the original corpus that
#'     the model was fitted to. Otherwise, compute deconvolved topics from this
#'     new corpus. Needs to be pixels x genes and nonnegative integer counts. 
#'     Each row needs at least 1 nonzero entry (default: NULL)
#' @param perc.rare.thresh the number of deconvolved cell-types with mean pixel proportion below this fraction used to assess
#'     performance of fitted models for each K. Recorded for each K. (default: 0.05)
#' 
#' @export
perplexityPlot <- function(models, corpus = NULL, perc.rare.thresh = 0.05){
  
  fitted_models <- models$models
  Ks <- as.vector(unlist(sapply(fitted_models, slot, "k")))
  pScores <- models$perplexities
  
  out <- lapply(1:length(Ks), function(i) {
    apply(getBetaTheta(fitted_models[[i]], corpus = corpus)$theta, 2, mean)
  })
  ## number of cell-types present at fewer than `perc.rare.thresh` on average across pixels
  numrare <- unlist(lapply(out, function(x) sum(x < perc.rare.thresh)))
  
  dat <- data.frame(K = as.double(Ks),
                    rareCts = numrare,
                    perplexity = pScores,
                    rareCtsAdj = scale0_1(numrare),
                    perplexAdj = scale0_1(pScores),
                    alphas = unlist(sapply(fitted_models, slot, "alpha")))
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
    sec_ax_labs <- seq(min(dat$perplexity), max(dat$perplexity), (max(dat$perplexity)-min(dat$perplexity))/1)
  } else {
    sec_ax_labs <- seq(min(dat$perplexity), max(dat$perplexity), (max(dat$perplexity)-min(dat$perplexity))/max(dat$rareCts))
  }
  sec_ax_breaks <- scale0_1(sec_ax_labs)
  
  if(length(sec_ax_labs) == 1){
    sec_ax_breaks <- 0
    dat$perplexAdj <- 0
  }
  
  print(dat)
  print(prim_ax_labs)
  print(prim_ax_breaks)
  print(sec_ax_labs)
  print(sec_ax_breaks)
  
  plt <- ggplot2::ggplot(dat) +
    ggplot2::geom_point(ggplot2::aes(y=rareCtsAdj, x=K), col="blue", lwd = 2) +
    ggplot2::geom_point(ggplot2::aes(y=perplexAdj, x=K), col="red", lwd = 2) +
    ggplot2::geom_line(ggplot2::aes(y=rareCtsAdj, x=K), col="blue", lwd = 2) +
    ggplot2::geom_line(ggplot2::aes(y=perplexAdj, x=K), col="red", lwd = 2) +
    ggplot2::geom_bar(ggplot2::aes(x = K, y = alphaBool), fill = dat$`alpha < 1`, stat = "identity", width = 1, alpha = 0.5) +
    ggplot2::scale_y_continuous(name=paste0("# cell-types with mean proportion < ", round(perc.rare.thresh*100, 2), "%"), breaks = prim_ax_breaks, labels = prim_ax_labs,
                                sec.axis= ggplot2::sec_axis(~ ., name="perplexity", breaks = sec_ax_breaks, labels = round(sec_ax_labs, 2))) +
    ggplot2::scale_x_continuous(breaks = min(dat$K):max(dat$K)) +
    ggplot2::labs(title = "Fitted model K's vs deconvolved cell-types and perplexity",
                  subtitle = "") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size=15, face=NULL),
      # plot.subtitle = ggplot2::element_text(size=13, face=NULL),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "black", size = 0.1),
      panel.ontop = TRUE,
      axis.title.y.left = ggplot2::element_text(color="blue", size = 13),
      axis.text.y.left = ggplot2::element_text(color="blue", size = 13),
      axis.title.y.right = ggplot2::element_text(color="red", size = 15, vjust = 1.5),
      axis.text.y.right = ggplot2::element_text(color="red", size = 13),
      axis.text.x = ggplot2::element_text(angle = 0, size = 13),
      axis.title.x = ggplot2::element_text(size=13)
    )
  print(plt)
}

