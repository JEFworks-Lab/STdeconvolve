#' Visualize all topic proportions across pixels with `scatterpie`
#'
#' @description Note: visualizes all cell-types in theta at once (could be
#'     individual cell-types or cell-type-clusters) so for accuracy of the proportions
#'     of each cell-type in a pixel, the row (pixel) should sum to 1.
#'
#' @param theta document (pixel) x cell-type proportion matrix
#' @param pos position of pixels, as data.frame with x and y columns
#' @param topicOrder order of topics in theta to visualize as a numeric vector
#'     and same length as topicCols (default: seq(ncol(theta)))
#' @param topicCols Vector of colors for each of the cell-types to be visualized.
#'     Same length and order as topicOrder (default: rainbow(ncol(theta)))
#' @param groups Indicates color of the scatterpie strokes (borders) with the goal of coloring them
#'     by their assigned group. This can be a vector or factor indicating the group of each
#'     pixel. Needs to be in the same order as the pixel rows in "theta" (default: NA)
#' @param group_cols Color labels for the groups. Can be a vector or factor. (default: NA)
#' @param r Radius of the scatterpie circles. Adjust based on positions of pixels (default: max(0.4, max(pos)/nrow(pos)*4))
#' @param lwd Width of lines of the pie charts. Increasing helps visualize
#'     group_cols if being used.
#' @param showLegend Boolean to show the legend indicating cell-types and their color
#' @param plotTitle add title to the resulting plot (default: NA)
#' @param overlay raster image of an H&E tissue (for example) to plot the scatterpies on top of
#'     (default: NA)
#' 
#' @return a plot of scatterpies, where each scatterpie represents
#'     a pixel in space based on the x,y coordinates and the components
#'     represent the proportion of each cell-type at that pixel.
#' 
#' @examples 
#' data(mOB)
#' pos <- mOB$pos
#' cd <- mOB$counts
#' annot <- mOB$annot
#' counts <- cleanCounts(cd, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = 8)
#' optLDA <- optimalModel(models = ldas, opt = 8)
#' results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
#' deconProp <- results$theta
#' vizAllTopics(deconProp,pos, groups = annot, group_cols = rainbow(length(levels(annot))), r=0.4)
#' 
#' @importFrom grDevices rainbow
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
  
  ## check that theta and pos are either data.frames or matrices
  if( !is.matrix(theta) & !is.data.frame(theta) ){
    stop("`theta` must be a matrix or data.frame.")
  }
  if( !is.matrix(pos) & !is.data.frame(pos) ){
    stop("`pos` must be a matrix or data.frame with exactly 2 columns named `x` and `y`.")
  }
  
  if( (any(!colnames(pos) %in% c("x", "y")) == TRUE) | (dim(pos)[2] != 2) ){
    stop("`pos` must have exactly 2 columns named `x` and `y`.")
  }
  
  # pixel cell-type distribution reordered based on topicOrder
  theta_ordered <- theta[, topicOrder]
  theta_ordered <- as.data.frame(theta_ordered)
  colnames(theta_ordered) <- paste0("Topic.", colnames(theta_ordered))
  
  # ensure that `theta` and `pos` pixel rownames maintain same order
  # after the merge so as to not mess up the order of `groups`
  # if provided
  # make sure only using the shared pixels
  pixels <- intersect(rownames(theta_ordered), rownames(pos))
  pixels <- rownames(theta_ordered)[which(rownames(theta_ordered) %in% pixels)]
  
  # add columns "x", "y" with document positions from `pos`
  theta_ordered_pos <- merge(data.frame(theta_ordered),
                             data.frame(pos), by=0)
  rownames(theta_ordered_pos) <- theta_ordered_pos[,"Row.names"]
  ## make sure pixels in the original order before the merge
  theta_ordered_pos <- theta_ordered_pos[pixels,]
  
  # first column after merge is "Row.names", last two are "x" and "y"
  # problem is that data frame will replace "-" and " " with "."
  topicColumns <- colnames(theta_ordered_pos)[2:(dim(theta_ordered_pos)[2]-2)]
  
  # color of piechart groups (lines of piechart):
  if (is.na(groups[1])) {
    groups <- rep("0", dim(theta_ordered_pos)[1])
    theta_ordered_pos$Pixel.Groups <- groups
  } else {
    theta_ordered_pos$Pixel.Groups <- as.character(groups)
  }
  if (is.na(group_cols[1])) {
    group_cols <- c("0" = "gray")
  }
  
  message("Plotting scatterpies for ", dim(theta_ordered_pos)[1], " pixels with ", length(topicColumns),
      " cell-types...this could take a while if the dataset is large.", "\n")
  
  if (!is.na(overlay[1])){
    p <- ggplot2::ggplot(mapping = ggplot2::aes(x = 0:dim(overlay)[2], y = 0:dim(overlay)[1])) +
      ggplot2::coord_equal(xlim = c(0,dim(overlay)[2]), ylim = c(0, dim(overlay)[1]), expand = FALSE) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 12, colour = "black"),
        legend.title = ggplot2::element_text(size = 12, colour = "black")
      ) +
      ggplot2::annotation_raster(overlay, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      scatterpie::geom_scatterpie(ggplot2::aes(x=x, y=y, group=Row.names, r=r, color = Pixel.Groups),
                                  lwd = lwd,
                                  data = theta_ordered_pos,
                                  cols = topicColumns,
                                  legend_name = "Topics") +
      ggplot2::scale_fill_manual(values = as.vector(topicCols)) +
      ggplot2::scale_color_manual(values = group_cols)
  } else {
    p <- ggplot2::ggplot() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 12, colour = "black"),
        legend.title = ggplot2::element_text(size = 12, colour = "black")
        ) +
      scatterpie::geom_scatterpie(ggplot2::aes(x=x, y=y, group=Row.names, r=r, color = Pixel.Groups),
                                  lwd = lwd,
                                  data = theta_ordered_pos,
                                  cols = topicColumns,
                                  legend_name = "Topics") +
      ggplot2::scale_fill_manual(values = as.vector(topicCols)) +
      ggplot2::scale_color_manual(values = group_cols)
  }
  
  if (!showLegend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  if (!is.na(plotTitle)) {
    p <- p + ggplot2::ggtitle(plotTitle)
  }
  
  p <- p + ggplot2::coord_equal()
  
  return(p)
}


#' Visualize pixel proportions of a single cell-type.
#' 
#' @description Visualize the pixel proportions of a single topic.
#'
#' @param theta document (pixel) x cell-type proportion matrix
#' @param pos position of pixels, as data.frame with `x` and `y` columns
#' @param topic the index of the topic
#' @param groups colors the pixel border lines based on a group or cell layer
#'     they belong to. Needs to be a character or named vector of assigned groups for each pixel
#'     Ex: c("0", "1", "0", ...)
#' @param group_cols color labels for the groups. Ex: c("0" = "gray", "1" = "red")
#' @param size size of the geom_points to plot (default: 2)
#' @param stroke thickness of the geom_point lines to help in emphasizing groups
#'     (default: 0.5)
#' @param alpha alpha value of colored pixels (default: 1)
#' @param low sets the color for the low end of the topic proportion color scale (default: "white")
#' @param high sets the color the the high end of the topic proportion color scale (default: "red")
#' @param plotTitle option to add a title to the plot (character)
#' @param showLegend Boolean to show the plot legend
#' 
#' @return a plot where each point is a pixel colored by the proportion of the selected cell-type
#' 
#' @examples 
#' data(mOB)
#' pos <- mOB$pos
#' cd <- mOB$counts
#' counts <- cleanCounts(cd, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = 8)
#' optLDA <- optimalModel(models = ldas, opt = 8)
#' results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
#' deconProp <- results$theta
#' vizTopic(theta = deconProp, pos = pos, topic = "5", plotTitle = "X5",
#'     size = 5, stroke = 1, alpha = 0.5, low = "white", high = "red")
#' 
#' @export
vizTopic <- function(theta, pos, topic,
                     groups = NA,
                     group_cols = NA,
                     size = 2,
                     stroke = 0.3,
                     alpha = 1,
                     low = "white",
                     high = "red",
                     plotTitle = NA,
                     showLegend = TRUE) {
  
  ## check that theta and pos are either data.frames or matrices
  if( !is.matrix(theta) & !is.data.frame(theta) == FALSE ){
    stop("`theta` must be a matrix or data.frame.")
  }
  if( !is.matrix(pos) == FALSE & !is.data.frame(pos) == FALSE ){
    stop("`pos` must be a matrix or data.frame with exactly 2 columns named `x` and `y`.")
  }
  
  if( (any(!colnames(pos) %in% c("x", "y")) == TRUE) | (dim(pos)[2] != 2) ){
    stop("`pos` must have exactly 2 columns named `x` and `y`.")
  }
  
  # ensure that `theta` and `pos` pixel rownames maintain same order
  # after the merge so as to not mess up the order of `groups`
  # if provided
  # make sure only using the shared pixels
  pixels <- intersect(rownames(theta), rownames(pos))
  pixels <- rownames(theta)[which(rownames(theta) %in% pixels)]
  
  proportion <- theta[,topic]
  dat <- merge(data.frame(proportion),
               data.frame(pos), by=0)
  
  rownames(dat) <- dat[,"Row.names"]
  ## make sure pixels in the original order before the merge
  dat <- dat[pixels,]
  
  # color spots by group:
  if (is.na(groups[1])) {
    Groups <- " "
  } else {
    Groups <- as.character(groups)
  }
  if (is.na(group_cols[1])) {
    group_cols <- c(" " = "black")
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = dat, ggplot2::aes(x=x, y=y, fill=proportion, color = Groups),
                        shape = 21,
                        stroke = stroke, size = size, 
                        alpha = alpha) +
    ggplot2::scale_color_manual(values = group_cols)
  
  p <- p +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 12, colour = "black"),
      legend.title = ggplot2::element_text(size = 12, colour = "black")
    ) +
    
    ggplot2::scale_fill_gradientn(limits = c(0, 1.0),
                                  breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                                  colors=(grDevices::colorRampPalette(c(low, high)))(n = 209)
    ) +
    
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Proportion",
                                                   title.position = "left",
                                                   title.hjust = 0.5,
                                                   ticks.colour = "black",
                                                   ticks.linewidth = 2,
                                                   frame.colour= "black",
                                                   frame.linewidth = 2,
                                                   label.hjust = 0,
                                                   title.theme = ggplot2::element_text(angle = 90)
    ))
  
  if (!showLegend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  if (!is.na(plotTitle)) {
    p <- p + ggplot2::ggtitle(plotTitle)
  }
  
  p <- p + ggplot2::coord_equal()
  
  return(p)
  
}


#' Visualize gene counts for a given gene in the pixels. Can also see group assignment of
#' spots.
#' 
#' @description Visualize one gene at a time.
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
#'     (default: 0.5)
#' @param alpha alpha value of colored pixels (default: 1)
#' @param plotTitle option to add a title to the plot
#' @param showLegend Boolean to show the plot legend
#' 
#' @return a plot where each point is a pixel colored by the expression level of the selected gene
#' 
#' @examples 
#' data(mOB)
#' pos <- mOB$pos
#' cd <- mOB$counts
#' counts <- cleanCounts(cd, min.lib.size = 100) 
#' df <- merge(as.data.frame(pos), as.data.frame(t(as.matrix(counts))), by = 0)
#' vizGeneCounts(df = df, gene = "Sox11",
#'     size = 3, stroke = 0.1, plotTitle = "Sox11",
#'     winsorize = 0.05, showLegend = TRUE)
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
  if (is.na(groups[1])) {
    groups <- " "
  } else {
    groups <- as.character(groups)
  }
  if (is.na(group_cols[1])) {
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
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()) + 
    
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Counts",
                                                   title.position = "left",
                                                   title.hjust = 0.5,
                                                   ticks.colour = "black",
                                                   ticks.linewidth = 2,
                                                   frame.colour= "black",
                                                   frame.linewidth = 2,
                                                   label.hjust = 0
    ))
  
  if (!showLegend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  if (!is.na(plotTitle)) {
    p <- p + ggplot2::ggtitle(plotTitle)
  }
  
  p <- p + ggplot2::coord_equal()
  
  return(p)
}


# custom correlation color range for heatmap.2 correlation plots
correlation_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(n = 209)
correlation_breaks <- c(seq(-1,-0.01,length=100),
                       seq(-0.009,0.009,length=10),
                       seq(0.01,1,length=100))


#' Generate heatmap of correlations
#' 
#' @description Visualize the correlations between topics stored in a matrix, typically one
#'     returned via `getCorrMtx()`. This function uses ggplot2::geom_tile.
#'
#' @param mat matrix with correlation values from -1 to 1
#' @param rowLabs y-axis label for plot. These are the rows of the matrix, or specifically m1 from getCorrMtx. (default: NULL)
#' @param colLabs x-axis label for plot. These are the columns of the matrix, or specifically m2 from getCorrMtx. (default: NULL)
#' @param title title of the plot. (default: NULL)
#' @param annotation Boolean to show the correlation values in the squares of the heatmap (default; FALSE)
#' 
#' @return a heatmap of the values in the input mat
#' 
#' @examples 
#' data(mOB)
#' pos <- mOB$pos
#' cd <- mOB$counts
#' counts <- cleanCounts(cd, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = 8)
#' optLDA <- optimalModel(models = ldas, opt = 8)
#' results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
#' deconProp <- results$theta
#' corMtx <- getCorrMtx(m1 = as.matrix(deconProp), m2 = as.matrix(deconProp), type = "t")
#' rownames(corMtx) <- paste0("X", seq(nrow(corMtx)))
#' colnames(corMtx) <- paste0("X", seq(ncol(corMtx)))
#' correlationPlot(mat = corMtx, title = "Proportional correlation", annotation = TRUE) +
#'     ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))
#' 
#' @export
correlationPlot <- function(mat, colLabs = NA, rowLabs = NA, title = NA, annotation = FALSE){
  
  correlation_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(n = 209)
  correlation_breaks <- c(seq(-1,-0.01,length=100),
                          seq(-0.009,0.009,length=10),
                          seq(0.01,1,length=100))
  
  dat <- reshape2::melt(mat)
  plt <- ggplot2::ggplot(data = dat) +
    ggplot2::geom_tile(ggplot2::aes(x = Var1, y = Var2, fill=value)) +
    
    ggplot2::scale_y_discrete(breaks = as.character(dat$Var2), labels = as.character(dat$Var2))
    
    ## if correlation values are to be plotted in squares:
    if(annotation){
      plt <- plt + ggplot2::geom_text(ggplot2::aes(x = as.character(Var1), y = as.character(Var2), label = format(round(value, 2), nsmall = 2) ))
    }
    
    plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_text(size=12, color = "black", hjust = 0, vjust = 0.5),
                   axis.text.y = ggplot2::element_text(size=12, color = "black"),
                   axis.title.y = ggplot2::element_text(size=13),
                   axis.title.x = ggplot2::element_text(size=13),
                   plot.title = ggplot2::element_text(size=15),
                   legend.text = ggplot2::element_text(size = 15, colour = "black"),
                   legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 90),
                   panel.background = ggplot2::element_blank(),
                   ## border around plot
                   panel.border = ggplot2::element_rect(fill = NA, color = "black", size = 2),
                   plot.background = ggplot2::element_blank()
                   # legend.position="none"
      ) +
      ## fix up colorbar legend
      ggplot2::scale_fill_gradientn(limits = c(-1,1),
                                    breaks = c(-1,0,1),
                                    colors=(grDevices::colorRampPalette(c("blue","white","red")))(n = 209)
      ) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Correlation",
                                                     title.position = "left",
                                                     title.hjust = 0.5,
                                                     ticks.colour = "black",
                                                     ticks.linewidth = 2,
                                                     frame.colour= "black",
                                                     frame.linewidth = 2,
                                                     label.hjust = 0
      )) +
      
      ggplot2::coord_fixed()
  
  if (!is.na(colLabs)){
    plt <- plt + ggplot2::xlab(colLabs)
  }
  if (!is.na(rowLabs)){
    plt <- plt + ggplot2::ylab(rowLabs)
  }
  if (!is.na(title)){
    plt <- plt + ggplot2::ggtitle(title)
  }
  
  return(plt)
  
}


#' Plot the perplexity and rare cell-types versus fitted Ks
#' 
#' @description the same plot returned by fitLDA() but now callable as a 
#'     separate function. 
#'     
#' @param models list returned from fitLDA
#' @param corpus If corpus is NULL, then it will use the original corpus that
#'     the model was fitted to. Otherwise, compute deconvolved topics from this
#'     new corpus. Needs to be pixels x genes and nonnegative integer counts. 
#'     Each row needs at least 1 nonzero entry (default: NULL)
#' @param perc.rare.thresh the number of deconvolved cell-types with mean pixel proportion below this fraction used to assess
#'     performance of fitted models for each K. Recorded for each K. (default: 0.05)
#'     
#' @return a plot indicating the perplexity and number of rare cell-types of a list of fitted LDA models
#' 
#' @examples 
#' data(mOB)
#' pos <- mOB$pos
#' cd <- mOB$counts
#' counts <- cleanCounts(cd, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2,6))
#' perplexityPlot(models = ldas, corpus = corpus)
#' 
#' @importFrom methods slot
#' 
#' @export
perplexityPlot <- function(models, corpus = NULL, perc.rare.thresh = 0.05){
  
  fitted_models <- models$models
  Ks <- as.vector(unlist(sapply(fitted_models, slot, "k")))
  pScores <- models$perplexities
  
  out <- lapply(1:length(Ks), function(i) {
    apply(getBetaTheta(fitted_models[[i]], corpus = corpus, verbose = FALSE)$theta, 2, mean)
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
                  subtitle = "LDA models with \u03b1 > 1 shaded") +
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

