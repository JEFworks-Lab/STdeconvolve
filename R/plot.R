#' Visualize all topic proportions across spots with `scatterpie`
#'
#' @description Note: visualizes all topics in theta at once (could be
#'     individual topics or topic clusters) so for accuracy of the proportions
#'     of each topic in a spot, the row (spot) should sum to 1.
#'
#' @param theta document (spot) x topic proportion matrix
#' @param pos position of documents, x and y columns
#' @param topicOrder order of topics in theta to visualize as a numeric vector
#'     and same length as topicCols ( default: seq(ncol(theta)) )
#' @param topicCols vector of colors for each of the topics to be visualized.
#'     Same length and order as topicOrder ( default: rainbow(ncol(theta)) )
#' @param groups colors the spot piechart lines based on a group or cell layer
#'     they belong to. Needs to be a character vector in the order of the spot
#'     rows in theta. Ex: c("0", "1", "0", ...)
#' @param group_cols color labels for the groups. Ex: c("0" = "gray", "1" = "red")
#' @param r radius of the circles. Adjust based on size of spots. (default: 1)
#' @param lwd width of lines of the pie charts. Increasing helps visualize
#'     group_cols if being used.
#' @param showLegend Boolean to show the legend indicating topics and their color
#' @plotTitle add title to the resulting plot (default: NA)
#'
#' @export
#'
vizAllTopics <- function(theta, pos,
                         topicOrder=seq(ncol(theta)),
                         topicCols=rainbow(ncol(theta)),
                         groups = NA,
                         group_cols = NA,
                         r = 1,
                         lwd = 0.5,
                         showLegend = TRUE,
                         plotTitle = NA) {

  # doc-topic distribution reordered based on topicOrder
  theta_ordered <- theta[, topicOrder]
  theta_ordered <- as.data.frame(theta_ordered)
  # colnames(theta_ordered) <- paste0("Topic.", topicOrder)
  colnames(theta_ordered) <- paste0("Topic.", colnames(theta_ordered))

  # add columns "x", "y" with document positions from `pos`
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
    scatterpie::geom_scatterpie(aes(x=x, y=y, group=Row.names, r=r, color = groups),
                    lwd = lwd,
                    data = theta_ordered_pos,
                    cols = topicColumns,
                    legend_name = "Topics") +
    scale_fill_manual(values = topicCols) +
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


#' Visualize proportions of topic clusters individually
#' 
#' @description Similar to `vizAllTopics` but will generate a separate plot for
#'     each topic or topic-cluster where the other topics or clusters will be
#'     colored gray. In this way, the actual proportions of each topic in a spot
#'     will be maintained such that spot topic proportions still sum to 1.
#'
#' @param theta document (spot) x topic proportion matrix
#' @param pos position of documents, x and y columns
#' @param clusters factor of the color (i.e., topic cluster) each cluster is
#'     assigned to. In this case, the levels should be colors. In `vizAllTopics`,
#'     clusters is "topicCols" and can just be a vector of colors.
#' @sharedCol Boolean indicating if the topics in a cluster will be plotted with
#'     the same color or if each topic will be colored by its own shade to also
#'     show how the topics in a cluster are distributed in space wrt each other.
#' @param groups colors the spot piechart lines based on a group or cell layer
#'     they belong to. Needs to be a character vector in the order of the spot
#'     rows in theta. Ex: c("0", "1", "0", ...)
#' @param group_cols color labels for the groups. Ex: c("0" = "gray", "1" = "red")
#' @param r = radius of the circles. Adjust based on size of spots. (default: 1)
#' @param lwd = width of lines of the pie charts. Increasing helps visualize
#'     group_cols if being used.
#' @param showLegend Boolean to show the legend indicating topics and their color
#' @plotTitle add title to the resulting plot (default: NA)
#'
#' @export
#'
vizTopicClusters <- function(theta, pos, clusters,
                             sharedCol = FALSE,
                             groups = NA,
                             group_cols = NA,
                             r = 1,
                             lwd = 0.5,
                             showLegend = TRUE,
                             plotTitle = NA) {

  print("Topic cluster members:")
  # produce a plot for each topic cluster:
  for (cluster in levels(clusters)) {

    # select the topics in the cluster
    topics <- labels(clusters[which(clusters == cluster)])

    cat(cluster, ":", topics, "\n")

    # doc-topic distribution reordered based on topicOrder and selected cluster topics
    theta_ordered <- theta[, topics]

    # get percentage of other topics not in cluster to maintain actual
    # spot proportions that sum to 1
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
    # and doesn't appear in the legend and messes with the colors such that
    # "other" takes one of the colors of the topics and is not gray
    if ( length(which(colSums(theta_ordered) == 0)) > 0 ) {
      missing_topics <- colnames(theta_ordered)[which(colSums(theta_ordered) == 0)]
      cat("NOTE:", missing_topics, "not present in any spots and will be dropped.", "\n")
      theta_ordered <- theta_ordered[,which(!colSums(theta_ordered) == 0)]
    }
    
    # add columns with document positions
    rownames(theta_ordered) <- rownames(pos)
    theta_ordered_pos <- merge(data.frame(theta_ordered),
                               data.frame(pos), by=0)

    # first column after merge is "Row.names", last two are "x" and "y"
    # problem is that data frame will replace "-" and " " with "."
    topicColumns <- colnames(theta_ordered_pos)[2:(dim(theta_ordered_pos)[2]-2)]

    # get a hue of colors for each topic in topic-cluster
    if (sharedCol){
      color_ramp <- colorRampPalette(c(cluster, cluster))
    } else {
      color_ramp <- colorRampPalette(c(lighten(cluster, factor = 0.5), darken(cluster, factor = 2)))
    }

    topic_colors <- color_ramp(ncol(theta_ordered) - 1) # don't count "other" here
    topic_colors <- append(topic_colors, c("gray")) # add gray to other here

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
      # theme(panel.background = element_rect(fill = "black"),
      #       panel.grid = element_blank()) +
      # theme_classic() +
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
      scatterpie::geom_scatterpie(aes(x=x, y=y, group=Row.names, r = r, color = groups),
                      lwd = lwd,
                      data=theta_ordered_pos,
                      cols = topicColumns,
                      legend_name = "Topics") +
      coord_equal() +
      scale_fill_manual(values=topic_colors) +
      scale_color_manual(values = group_cols)

    if (showLegend == FALSE) {
      p <- p + guides(fill=FALSE)
    }
    
    if (is.na(plotTitle) == FALSE) {
      p <- p + ggtitle(plotTitle)
    }

    print(p)
  }
}


#' Visualize gene counts in spots in space. Can also see group assignment of
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
#' @param size size of the geom_points to plot (default: 7)
#' @param stroke thickness of the geom_point lines to help in emphasizing groups
#'     (default: 2)
#' @param plotTitle option to add a title to the plot
#' @param showLegend Boolean to show the plot legend
#' 
#' @export
#'
vizGeneCounts <- function(df, gene,
                          groups = NA,
                          group_cols = NA,
                          size = 7, stroke = 2,
                          plotTitle = NA,
                          showLegend = TRUE) {

  counts <- df[,gene]

  # color spots by group:
  if (is.na(groups[1]) == TRUE) {
    groups <- " "
    # stroke <- 0.5
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
  
  p <- p +
    theme(
      #panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      panel.background=element_blank())
    # theme_classic()
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


