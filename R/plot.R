#' Convert a matrix to a long form data frame for easy plotting via ggplot2
#' Jean: Brendan, please check if this function is necessary
#' or could be replaced with reshape2::melt
#'
#' @param mtx the matrix to be converted
#' @param colLab df column label for the mtx columns
#' @param rowLab df row label for the row labels
#' @param cellLab label for the cell values
#'
mtx2ggplotDf <- function(mtx, colLab = "columns", rowLab = "rows", cellLab = "cells") {

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


#' Visualize topic proportions across spots with `scatterpie`
#'
#' @param theta document x topic proportion matrix
#' @param pos position of documents, x and y columns
#' @param topicOrder order of topics based on dendrogram; a numeric vector
#'             from: (clusterTopics$order)
#' @param cluster_cols vector of colors for each of the topics
#'     can use factor of clusters for each topic (via clusterTopics$clusters)
#'     as long as the cluster levels/values have been converted to colors.
#'     This gets reordered wrt the topicOrder fyi.
#' @param groups color spot piecharts based on a group or cell layer they belong to.
#'          Needs to be a character vector. Ex: c("0", "1", "0", ...).
#' @param group_cols color labels for the groups. Ex: c("0" = "gray", "1" = "red")
#' @param r = radius of the circles. Adjust based on size of spots.
#' @param lwd = width of lines of the pie charts
#'
#' @export
#'
vizAllTopics <- function(theta, pos,
                         topicOrder=seq(ncol(theta)),
                         cluster_cols=rainbow(ncol(theta)),
                         groups = NA,
                         group_cols = NA,
                         r = 1,
                         lwd = 0.5,
                         showLegend = TRUE,
                         plotTitle = NA) {

  # reorder colors of topics wrt the dendrogram
  #colors_ordered <- as.vector(cluster_cols[topicOrder])
  colors_ordered <- cluster_cols

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
  if (is.na(groups) == TRUE) {
    groups <- rep("0", dim(theta_ordered_pos)[1])
    theta_ordered_pos$groups <- groups
  } else {
    theta_ordered_pos$groups <- as.character(groups)
  }
  if (is.na(group_cols) == TRUE) {
    group_cols <- c("0" = "gray")
  }

  p <- ggplot() +
    theme(
      #panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      panel.background=element_blank()) +
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


#' Visualize proportions of topic clusters
#'
#' @param theta document topic proportion matrix
#' @param pos position of documents, x and y columns
#' @param topicOrder order of topics based on dendrogram; a numeric vector
#' @param clusters factor of the color (topic cluster) each cluster is assigned to
#'            In this case, the levels should be colors. In `vizAllTopics`,
#'            clusters is "cluster_cols" and can just be a vector of colors.
#' @param groups color spot piecharts based on a group or cell layer they belong to.
#'          Needs to be a character vector. Ex: c("0", "1", "0", ...).
#' @param group_cols color labels for the groups. Ex: c("0" = "gray", "1" = "red")
#' @param r = radius of the circles. Adjust based on size of spots.
#'     40 = merfish; 0.4 = mOB
#' @param lwd = width of lines of the pie charts
#'
#' @export
#'
vizTopicClusters <- function(theta, pos, topicOrder, clusters,
                             groups = NA,
                             group_cols = NA,
                             sharedCol = FALSE,
                             r = 40,
                             lwd = 0.5,
                             plotTitle = NA) {

  # reorder factor wrt the dendrogram
  clusters_ordered <- clusters[topicOrder]

  print("Topic cluster members:")
  # produce a plot for each topic cluster:
  for (cluster in levels(clusters)) {

    # select the topics in the cluster in the same order of the dendrogram
    topics <- labels(clusters_ordered[which(clusters_ordered == cluster)])

    cat(cluster, ":", topics, "\n")

    # doc-topic distribution reordered based on topicOrder and selected cluster topics
    theta_ordered <- theta[, topics]

    # get percentage of other topics not in cluster
    if (is.null(dim(theta_ordered))) {
      other <- 1 - theta_ordered
    } else {
      other <- 1 - rowSums(theta_ordered)
    }

    theta_ordered <- as.data.frame(theta_ordered)
    colnames(theta_ordered) <- paste0("Topic.", topics)
    theta_ordered$other <- other

    # if any topics not represented at all, drop them
    # Apparently if a topic is 0 for all pie charts, it is not plotted
    # and doesn't appear in the legend. So it messes with the colors.
    # "other" takes one of the colors of the topics and is not gray
    theta_ordered <- theta_ordered[,which(!colSums(theta_ordered) == 0)]

    # add columns with document positions
    rownames(theta_ordered) <- rownames(pos)
    theta_ordered_pos <- merge(data.frame(theta_ordered),
                               data.frame(pos), by=0)

    # first column after merge is "Row.names", last two are "x" and "y"
    # problem is that data frame will replace "-" and " " with "."
    topicColumns <- colnames(theta_ordered_pos)[2:(dim(theta_ordered_pos)[2]-2)]

    # get a hue of colors representing the cluster color
    if (sharedCol){
      color_ramp <- colorRampPalette(c(cluster, cluster))
    } else {
      color_ramp <- colorRampPalette(c(lighten(cluster, factor = 0.8), darken(cluster, factor = 1)))
    }

    # topic_colors <- viridis_pal(option = "C")(length(blue_cluster))
    topic_colors <- color_ramp(ncol(theta_ordered) - 1) # don't count "other"
    topic_colors <- append(topic_colors, c("gray"))

    # color of piechart groups (lines of piechart):
    if (is.na(groups) == TRUE) {
      groups <- rep("0", dim(theta_ordered_pos)[1])
      theta_ordered_pos$groups <- groups
    } else {
      theta_ordered_pos$groups <- as.character(groups)
    }
    if (is.na(group_cols) == TRUE) {
      group_cols <- c("0" = "gray")
    }

    p <- ggplot() +
      theme(panel.background = element_rect(fill = "black"),
            panel.grid = element_blank()) +
      # theme_classic() +
      geom_scatterpie(aes(x=x, y=y, group=Row.names, r = r, color = groups), # r=40 for MERFISH 0.4 mOB
                      lwd = lwd,
                      data=theta_ordered_pos,
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


# visualize gene counts in spots in space. Also see group assignment of
# spots.
#
# df = data.frame where rows are spots and columns must be at least:
#      "x", "y" for spot positions in space
#      "gene" counts of a gene for each spot
#
# gene = column name of the gene counts in df to be visualized
#
# groups = a character vector of labels assigning spots to different groups.
#          Ex: c("0", "1", "0", ...).
# group_cols = a vector designating the spot border color to be used for
#              group assignment. Ex: c("0" = "gray", "1" = "red").
#
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
    scale_color_manual(values = group_cols) +
    ggtitle("Granule txn C1")

  if (showLegend == FALSE) {
    p <- p + guides(fill=FALSE)
  }
  if (is.na(plotTitle) == FALSE) {
    p <- p + ggtitle(plotTitle)
  }

  p <- p + theme_classic()
  print(p)
}



# custom correlation color range for heatmap.2 correlation plots
correlation_palette <- colorRampPalette(c("blue", "white", "red"))(n = 209)
correlation_breaks = c(seq(-1,-0.01,length=100),
                       seq(-0.009,0.009,length=10),
                       seq(0.01,1,length=100))

# lighten and darken a color
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


# color palette to replicate ggplot2
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


