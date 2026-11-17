## columns of the data that are referred to in ggplot2::aes()
utils::globalVariables(c("panel", "celltype", "gene", "expression", "relative",
                         "label", "row", "column"))


## colors of the cell-types: 8 hues that are distinct for the common types of
## color blindness, followed by a light and a dark version of each of them
celltype_palette <- c(
  "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948",
  "#9fc2ec", "#f6bba3", "#98dbc3", "#f6d48c", "#f4c3d6", "#8cc78c", "#b8b0ec", "#f2adac",
  "#174275", "#81391c", "#0e6043", "#825800", "#7f435a", "#004800", "#352c75", "#7c2827")

## one hue from light to dark, for proportions and expression levels
sequential_palette <- c("#f0efec", "#b7d3f6", "#5598e7", "#256abf", "#0d366b")

## two hues with a neutral midpoint, for correlations
diverging_palette <- c("#256abf", "#f0efec", "#d03b3b")


#' Colors of the cell-types
#'
#' @description The color of a cell-type depends only on its position among
#'     the cell-types of the model, so each cell-type has the same color in all
#'     figures of this package. The first 8 colors are distinct for people
#'     with the common types of color blindness. Colors 9 to 16 are light
#'     versions and colors 17 to 24 are dark versions of the first 8, so with
#'     more than 8 cell-types some pairs of colors are hard to tell apart.
#'     With more than 24 cell-types, evenly spaced hues are used instead.
#'
#' @param celltypes the number of cell-types, their names, a theta matrix
#'     (pixels x cell-types) or the list returned by `getBetaTheta()`
#'
#' @return character vector of colors, named by cell-type
#'
#' @examples
#' cellTypeColors(5)
#' cellTypeColors(c("neurons", "glia"))
#'
#' @export
cellTypeColors <- function(celltypes) {
  if (is.list(celltypes) || is.matrix(celltypes)) {
    celltypes <- colnames(thetaOf(celltypes))
  } else if (is.numeric(celltypes) && length(celltypes) == 1) {
    celltypes <- as.character(seq_len(celltypes))
  }
  celltypes <- as.character(celltypes)
  n <- length(celltypes)
  if (n <= length(celltype_palette)) {
    cols <- celltype_palette[seq_len(n)]
  } else {
    cols <- grDevices::hcl(h = seq(15, 375, length.out = n + 1)[-1],
                           c = 80, l = rep_len(c(45, 65, 80), n))
  }
  stats::setNames(cols, celltypes)
}


#' The pixel x cell-type proportions of the result of getBetaTheta()
#'
#' @param x list with an element `theta`, or a theta matrix or data.frame
#' @return matrix
#' @noRd
thetaOf <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    if (is.null(x$theta)) {
      stop("`theta` must be a matrix or data.frame, or the list returned by getBetaTheta().")
    }
    x <- x$theta
  }
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("`theta` must be a matrix or data.frame, or the list returned by getBetaTheta().")
  }
  x <- as.matrix(x)
  if (is.null(colnames(x))) {
    colnames(x) <- as.character(seq_len(ncol(x)))
  }
  x
}


#' The cell-type x gene matrix of the result of getBetaTheta()
#'
#' @param x list with an element `beta`, or a beta matrix or data.frame
#' @param required if FALSE, return NULL if `x` is not a list with a beta
#' @return matrix with rows that sum to 1
#' @noRd
betaOf <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    if (is.null(x$beta)) {
      stop("`beta` must be a matrix, or the list returned by getBetaTheta().")
    }
    x <- x$beta
  }
  x <- as.matrix(x)
  if (is.null(rownames(x))) {
    rownames(x) <- as.character(seq_len(nrow(x)))
  }
  if (is.null(colnames(x))) {
    stop("The columns of `beta` must be named by gene.")
  }
  x / rowSums(x)
}


#' Log2 fold change of each gene in each cell-type versus the other cell-types
#'
#' @param beta cell-type x gene matrix
#' @return matrix of the same size
#' @noRd
geneSpecificity <- function(beta) {
  beta <- beta / rowSums(beta)
  k <- nrow(beta)
  ## a tenth of the probability of a gene if all genes were expressed equally
  eps <- 0.1 / ncol(beta)
  if (k == 1) {
    others <- beta
  } else {
    others <- (rep(colSums(beta), each = k) - beta) / (k - 1)
  }
  log2((beta + eps) / (others + eps))
}


#' Values of pixels at their positions, as points or as cells of an image
#'
#' @description With many pixels, the pixels are averaged in the cells of a
#'     regular grid, which is drawn as an image.
#'
#' @param values matrix with a row for each pixel and a column for each panel
#' @param pos positions of the pixels with columns `x` and `y`
#' @param raster whether to average in grid cells. NULL: with more than 20000 pixels.
#' @param resolution largest number of grid cells along an axis
#' @param weights weight of each pixel in the average of a grid cell
#'
#' @return list with the positions `x`, `y`, the matrix `values` with a row
#'     for each point or grid cell, and `raster`
#' @noRd
mapData <- function(values, pos, raster = NULL, resolution = 400, weights = NULL) {
  if (!is.matrix(pos) && !is.data.frame(pos)) {
    stop("`pos` must be a matrix or data.frame with columns named `x` and `y`.")
  }
  if (!all(c("x", "y") %in% colnames(pos))) {
    stop("`pos` must have columns named `x` and `y`.")
  }
  if (!is.null(rownames(values)) && !is.null(rownames(pos)) &&
      !identical(rownames(values), rownames(pos))) {
    shared <- rownames(values) %in% rownames(pos)
    if (!any(shared)) {
      stop("None of the pixels (rows) of `theta` have a position in `pos`.")
    }
    values <- values[shared, , drop = FALSE]
    weights <- weights[shared]
    pos <- pos[rownames(values), , drop = FALSE]
  } else if (nrow(values) != nrow(pos)) {
    stop("`pos` must have a row for each pixel.")
  }
  x <- as.numeric(pos[, "x"])
  y <- as.numeric(pos[, "y"])
  n <- length(x)
  if (is.null(raster)) {
    raster <- n > 20000
  }
  if (!raster) {
    return(list(x = x, y = y, values = values, raster = FALSE))
  }

  ## size of the grid cells: a multiple of the spacing of the pixels if they
  ## are on a regular grid, otherwise large enough for about 4 pixels per cell
  step <- function(v) {
    d <- diff(sort(unique(v)))
    if (length(d) == 0) return(NA)
    s <- stats::median(d)
    if (s > 0 && all(abs(d / s - round(d / s)) < 0.05)) s else NA
  }
  sx <- step(x)
  sy <- step(y)
  xr <- range(x)
  yr <- range(y)
  if (!is.na(sx) && !is.na(sy)) {
    f <- max(1, ceiling(max(diff(xr) / sx, diff(yr) / sy) / resolution))
    w <- f * sx
    h <- f * sy
    x0 <- xr[1] - sx / 2
    y0 <- yr[1] - sy / 2
  } else {
    res <- max(1, min(resolution, floor(sqrt(n) / 2)))
    w <- h <- max(diff(xr), diff(yr), 1e-12) / res * (1 + 1e-9)
    x0 <- xr[1]
    y0 <- yr[1]
  }
  ix <- floor((x - x0) / w)
  iy <- floor((y - y0) / h)
  nx <- max(ix) + 1
  cell <- ix + iy * nx
  if (is.null(weights)) {
    weights <- rep(1, n)
  }
  total <- rowsum(weights, cell)[, 1]
  values <- rowsum(values * weights, cell) / pmax(total, 1e-300)
  id <- as.numeric(rownames(values))
  ix <- id %% nx
  iy <- id %/% nx
  rownames(values) <- NULL
  ## an image needs all columns and rows of the grid: empty ones get a cell
  ## without a value, which is not drawn
  gapx <- setdiff(seq(0, max(ix)), ix)
  gapy <- setdiff(seq(0, max(iy)), iy)
  if (length(gapx) + length(gapy) > 0) {
    ix <- c(ix, gapx, rep(ix[1], length(gapy)))
    iy <- c(iy, rep(iy[1], length(gapx)), gapy)
    values <- rbind(values, matrix(NA_real_, length(gapx) + length(gapy), ncol(values)))
  }
  list(x = x0 + (ix + 0.5) * w, y = y0 + (iy + 0.5) * h, values = values, raster = TRUE)
}


#' Layer of a map, as an image or as points
#'
#' @param df data.frame with the columns `x`, `y` and `value`
#' @param raster whether `df` has the cells of a regular grid
#' @param size size of the points. NULL: depending on the number of points.
#' @param npanel number of columns of panels the plot is divided in
#' @noRd
mapLayer <- function(df, raster, size = NULL, npanel = 1) {
  if (raster) {
    ggplot2::geom_raster(ggplot2::aes(x = x, y = y, fill = value), data = df, na.rm = TRUE)
  } else {
    if (is.null(size)) {
      n <- nrow(df) / max(1, length(unique(df$panel)))
      size <- max(0.1, min(6, 90 / sqrt(n))) / sqrt(npanel)
    }
    ggplot2::geom_point(ggplot2::aes(x = x, y = y, colour = value), data = df,
                        shape = 16, size = size)
  }
}


## the theme of the maps: no axes
mapTheme <- function() {
  ggplot2::theme_void() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 11, colour = "black",
                                         margin = ggplot2::margin(2, 2, 4, 2)),
      legend.text = ggplot2::element_text(size = 10, colour = "black"),
      legend.title = ggplot2::element_text(size = 11, colour = "black"),
      plot.title = ggplot2::element_text(size = 13, colour = "black",
                                         margin = ggplot2::margin(2, 2, 6, 2)))
}


## pixel x panel matrix as a data.frame with one row per pixel and panel
longMap <- function(m) {
  panels <- colnames(m$values)
  data.frame(x = rep(m$x, length(panels)),
             y = rep(m$y, length(panels)),
             panel = factor(rep(panels, each = length(m$x)), levels = panels),
             value = as.vector(m$values))
}


#' Maps of the proportions of all cell-types, one panel per cell-type
#'
#' @description Unlike the pie charts of `vizAllTopics()`, this also works
#'     for data with very many pixels: with more than 20000 pixels the
#'     proportions are averaged in the cells of a grid and drawn as an image.
#'
#' @param theta pixel x cell-type proportion matrix, or the list returned by
#'     `getBetaTheta()`
#' @param pos positions of the pixels, as a matrix or data.frame with columns
#'     `x` and `y` and the pixels as row names
#' @param topics names or indices of the cell-types to show (default: all)
#' @param scale "shared": the same color scale from 0 to 1 in all panels.
#'     "free": each panel is scaled to the largest proportion of its
#'     cell-type, which is stated in its title. This shows the pattern of
#'     cell-types with low proportions. (default: "shared")
#' @param ncol number of columns of panels (default: chosen by ggplot2)
#' @param size size of the points (default: depends on the number of pixels)
#' @param raster whether to draw an image of averages in grid cells instead of
#'     a point for each pixel (default: with more than 20000 pixels)
#' @param resolution largest number of grid cells along an axis of the image
#'     (default: 400)
#' @param colors colors of the proportions, from low to high
#' @param showLegend Boolean to show the legend
#' @param plotTitle title of the plot (default: NA)
#'
#' @return a ggplot object
#'
#' @examples
#' data(mOB)
#' counts <- cleanCounts(mOB$counts, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove = 1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = 3)
#' results <- getBetaTheta(optimalModel(ldas, opt = 3))
#' vizTopicMaps(results, mOB$pos)
#'
#' @export
vizTopicMaps <- function(theta, pos, topics = NULL,
                         scale = c("shared", "free"),
                         ncol = NULL, size = NULL,
                         raster = NULL, resolution = 400,
                         colors = sequential_palette,
                         showLegend = TRUE, plotTitle = NA) {
  scale <- match.arg(scale)
  theta <- thetaOf(theta)
  if (!is.null(topics)) {
    theta <- theta[, topics, drop = FALSE]
  }
  m <- mapData(theta, pos, raster = raster, resolution = resolution)
  if (scale == "free") {
    top <- apply(m$values, 2, max, na.rm = TRUE)
    m$values <- t(t(m$values) / pmax(top, 1e-300))
    colnames(m$values) <- paste0("Cell-type ", colnames(theta),
                                 " (max ", format(round(top, 2), nsmall = 2), ")")
    legend <- "Proportion\nrelative to\nthe maximum\nof the panel"
  } else {
    colnames(m$values) <- paste("Cell-type", colnames(theta))
    legend <- "Proportion"
  }
  df <- longMap(m)
  if (is.null(ncol)) {
    ncol <- ggplot2::wrap_dims(ncol(theta))[2]
  }
  p <- ggplot2::ggplot() +
    mapLayer(df, m$raster, size = size, npanel = ncol) +
    ggplot2::facet_wrap(~ panel, ncol = ncol) +
    ggplot2::scale_fill_gradientn(name = legend, colors = colors, limits = c(0, 1),
                                  na.value = NA, aesthetics = c("fill", "colour")) +
    ggplot2::coord_equal() +
    mapTheme()
  if (!showLegend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  if (!is.na(plotTitle)) {
    p <- p + ggplot2::ggtitle(plotTitle)
  }
  p
}


#' Map of the cell-type with the highest proportion in each pixel
#'
#' @description A summary of all cell-types in one map, which is fast to
#'     draw and readable with any number of pixels. It hides the mixture of
#'     cell-types within a pixel: use `vizTopicMaps()` or, for few pixels,
#'     `vizAllTopics()` to see them.
#'
#' @inheritParams vizTopicMaps
#' @param topicCols colors of the cell-types, in the order of the columns of
#'     `theta` (default: `cellTypeColors(theta)`)
#'
#' @return a ggplot object
#'
#' @examples
#' data(mOB)
#' counts <- cleanCounts(mOB$counts, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove = 1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = 3)
#' results <- getBetaTheta(optimalModel(ldas, opt = 3))
#' vizDominantTopic(results, mOB$pos)
#'
#' @export
vizDominantTopic <- function(theta, pos,
                             topicCols = cellTypeColors(theta),
                             size = NULL, raster = NULL, resolution = 400,
                             showLegend = TRUE, plotTitle = NA) {
  theta <- thetaOf(theta)
  topicCols <- stats::setNames(as.vector(topicCols), colnames(theta))
  m <- mapData(theta, pos, raster = raster, resolution = resolution)
  ## pixels without any cell-type, and empty grid cells, have no color
  present <- !is.na(m$values[, 1]) & rowSums(m$values) > 0
  dominant <- rep(NA_character_, length(present))
  dominant[present] <- colnames(theta)[max.col(m$values[present, , drop = FALSE],
                                               ties.method = "first")]
  df <- data.frame(x = m$x, y = m$y, value = factor(dominant, levels = colnames(theta)))
  if (!m$raster) {
    df <- df[present, , drop = FALSE]
  }
  p <- ggplot2::ggplot() +
    mapLayer(df, m$raster, size = size) +
    ggplot2::scale_fill_manual(name = "Cell-type", values = topicCols, drop = FALSE,
                               breaks = colnames(theta), na.value = NA,
                               aesthetics = c("fill", "colour")) +
    ggplot2::coord_equal() +
    mapTheme()
  if (!m$raster) {
    ## legend keys that are large enough to see the color
    p <- p + ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 4)))
  }
  if (!showLegend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  if (!is.na(plotTitle)) {
    p <- p + ggplot2::ggtitle(plotTitle)
  }
  p
}


#' Bar chart of the abundance of each cell-type
#'
#' @description The length of a bar is the mean proportion of the cell-type
#'     across all pixels. The label also states the percentage of the pixels
#'     that contain the cell-type.
#'
#' @inheritParams vizDominantTopic
#' @param sort whether to order the cell-types by their mean proportion
#'     (default: TRUE)
#'
#' @return a ggplot object
#'
#' @examples
#' data(mOB)
#' counts <- cleanCounts(mOB$counts, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove = 1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = 3)
#' results <- getBetaTheta(optimalModel(ldas, opt = 3))
#' vizTopicAbundance(results)
#'
#' @export
vizTopicAbundance <- function(theta, topicCols = cellTypeColors(theta),
                              sort = TRUE, plotTitle = NA) {
  theta <- thetaOf(theta)
  topicCols <- stats::setNames(as.vector(topicCols), colnames(theta))
  df <- data.frame(celltype = colnames(theta),
                   value = colMeans(theta),
                   present = colMeans(theta > 0))
  df$label <- paste0(format(round(100 * df$value, 1), nsmall = 1, trim = TRUE), "%, in ",
                     round(100 * df$present), "% of pixels")
  ## the first cell-type is at the top
  levels <- if (sort) df$celltype[order(df$value)] else rev(df$celltype)
  df$celltype <- factor(df$celltype, levels = levels)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = value, y = celltype)) +
    ggplot2::geom_col(ggplot2::aes(fill = celltype), width = 0.7, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = label), hjust = 0, nudge_x = max(df$value) * 0.01,
                       size = 3.5, colour = "black") +
    ggplot2::scale_fill_manual(values = topicCols) +
    ggplot2::scale_x_continuous(labels = function(v) paste0(100 * v, "%"),
                                expand = ggplot2::expansion(mult = c(0, 0.35))) +
    ggplot2::labs(x = "Mean proportion across pixels", y = "Cell-type") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 11, colour = "black"),
                   axis.title = ggplot2::element_text(size = 12),
                   axis.ticks.y = ggplot2::element_blank())
  if (!is.na(plotTitle)) {
    p <- p + ggplot2::ggtitle(plotTitle)
  }
  p
}


#' Dot plot of the most specific genes of each cell-type
#'
#' @description Shows what each cell-type is. The genes are the `n` genes of
#'     each cell-type with the highest fold change of their expression in the
#'     cell-type versus their mean expression in the other cell-types (see
#'     `topGenes(beta, specific = TRUE)`). The area of a dot is the expression
#'     of the gene in the cell-type relative to the cell-type with the highest
#'     expression of the gene, so each gene has one dot of the largest size.
#'     The color of a dot is the expression of the gene in the cell-type, as
#'     counts per 1000 counts of the cell-type, on a log scale.
#'
#' @param beta cell-type x gene expression matrix, or the list returned by
#'     `getBetaTheta()`
#' @param n number of genes of each cell-type (default: 5)
#' @param topics names or indices of the cell-types whose genes are shown
#'     (default: all)
#' @param genes genes to show instead of the most specific genes, e.g. known
#'     marker genes (default: NULL)
#' @param colors colors of the relative expression, from low to high
#' @param plotTitle title of the plot (default: NA)
#'
#' @return a ggplot object
#'
#' @examples
#' data(mOB)
#' counts <- cleanCounts(mOB$counts, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove = 1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = 3)
#' results <- getBetaTheta(optimalModel(ldas, opt = 3))
#' vizTopGenes(results, n = 5)
#'
#' @export
vizTopGenes <- function(beta, n = 5, topics = NULL, genes = NULL,
                        colors = sequential_palette, plotTitle = NA) {
  beta <- betaOf(beta)
  if (is.null(genes)) {
    top <- topGenes(beta, n = n, specific = TRUE)
    if (!is.null(topics)) {
      top <- top[topics]
    }
    genes <- unique(unlist(lapply(top, names)))
  } else {
    missing <- setdiff(genes, colnames(beta))
    if (length(missing) > 0) {
      warning("Genes that are not in `beta` are left out: ", paste(missing, collapse = ", "))
    }
    genes <- unique(genes[genes %in% colnames(beta)])
    if (length(genes) == 0) {
      stop("None of the `genes` are in `beta`.")
    }
  }
  sub <- beta[, genes, drop = FALSE]
  df <- data.frame(celltype = factor(rep(rownames(sub), length(genes)), levels = rownames(sub)),
                   ## the first gene is at the top
                   gene = factor(rep(genes, each = nrow(sub)), levels = rev(genes)),
                   ## the color scale starts at 0.1 counts per 1000
                   expression = log10(pmax(as.vector(sub) * 1000, 0.1)),
                   relative = as.vector(t(t(sub) / apply(sub, 2, max))))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = celltype, y = gene)) +
    ggplot2::geom_point(ggplot2::aes(size = relative, fill = expression),
                        shape = 21, colour = "#52514e", stroke = 0.3) +
    ggplot2::scale_size_area(name = "Expression\nrelative to the\nhighest cell-type",
                             max_size = 5.5, limits = c(0, 1), breaks = c(0.25, 0.5, 1)) +
    ggplot2::scale_fill_gradientn(name = "Expression\n(per 1000 counts\nof the cell-type)",
                                  colors = colors, limits = c(-1, max(2, df$expression)),
                                  breaks = -1:2, labels = c("0.1 or less", "1", "10", "100")) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(x = "Cell-type", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "#e1e0d9", linewidth = 0.3),
                   axis.text = ggplot2::element_text(size = 10, colour = "black"),
                   axis.text.y = ggplot2::element_text(face = "italic"),
                   axis.title = ggplot2::element_text(size = 12))
  if (!is.na(plotTitle)) {
    p <- p + ggplot2::ggtitle(plotTitle)
  }
  p
}


#' Map of a cell-type next to maps of its most specific genes
#'
#' @description A check that a cell-type is real: the pixels with a high
#'     proportion of the cell-type should be the pixels where its genes are
#'     expressed. Each panel has its own scale, from 0 to the largest value
#'     in the panel, which is stated in its title.
#'
#' @inheritParams vizTopicMaps
#' @param theta the list returned by `getBetaTheta()`, or a pixel x cell-type
#'     proportion matrix (then `genes` are required)
#' @param counts gene x pixel count matrix, e.g. the output of `cleanCounts()`
#' @param topic name or index of the cell-type
#' @param genes genes to show (default: the `n` most specific genes of the
#'     cell-type that are in `counts`)
#' @param n number of genes (default: 3)
#' @param normalize whether to show the counts of a gene as a fraction of all
#'     counts of the pixel in `counts`, instead of the counts (default: TRUE)
#' @param winsorize the largest values of a gene, this fraction of the
#'     pixels, are set to the next largest value, so that a few pixels do not
#'     determine the color scale (default: 0.005)
#'
#' @return a ggplot object
#'
#' @examples
#' data(mOB)
#' counts <- cleanCounts(mOB$counts, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove = 1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = 3)
#' results <- getBetaTheta(optimalModel(ldas, opt = 3))
#' vizTopicMarkers(results, mOB$pos, counts, topic = 1)
#'
#' @export
vizTopicMarkers <- function(theta, pos, counts, topic, genes = NULL, n = 3,
                            normalize = TRUE, winsorize = 0.005,
                            ncol = NULL, size = NULL,
                            raster = NULL, resolution = 400,
                            colors = sequential_palette,
                            showLegend = TRUE, plotTitle = NA) {
  results <- theta
  theta <- thetaOf(results)
  if (length(topic) != 1) {
    stop("`topic` must be the name or index of one cell-type.")
  }
  if (is.numeric(topic)) {
    topic <- colnames(theta)[topic]
  }
  topic <- as.character(topic)
  if (is.na(topic) || !topic %in% colnames(theta)) {
    stop("`topic` is not a cell-type (column) of `theta`.")
  }
  if (!inherits(counts, "Matrix") && !is.matrix(counts)) {
    counts <- as.matrix(counts)
  }
  if (is.null(genes)) {
    if (!is.list(results) || is.data.frame(results) || is.null(results$beta)) {
      stop("Provide `genes`, or the list returned by getBetaTheta() as `theta`.")
    }
    beta <- betaOf(results)
    beta <- beta[, colnames(beta) %in% rownames(counts), drop = FALSE]
    if (ncol(beta) == 0) {
      stop("`counts` must be a gene x pixel matrix with the genes of the model as row names.")
    }
    genes <- names(topGenes(beta, n = min(n, ncol(beta)), specific = TRUE)[[topic]])
  } else if (!all(genes %in% rownames(counts))) {
    stop("Genes that are not in the rows of `counts`: ",
         paste(setdiff(genes, rownames(counts)), collapse = ", "))
  }

  pixels <- rownames(theta)
  if (is.null(pixels) || is.null(colnames(counts))) {
    if (ncol(counts) != nrow(theta)) {
      stop("`counts` must have a column for each pixel (row) of `theta`.")
    }
    cols <- seq_len(nrow(theta))
  } else {
    pixels <- pixels[pixels %in% colnames(counts)]
    if (length(pixels) == 0) {
      stop("None of the pixels (rows) of `theta` are columns of `counts`.")
    }
    theta <- theta[pixels, , drop = FALSE]
    cols <- pixels
  }
  expr <- t(as.matrix(counts[genes, cols, drop = FALSE]))
  total <- if (normalize) Matrix::colSums(counts)[cols] else rep(1, nrow(expr))
  rownames(expr) <- rownames(theta)

  ## the genes are averaged with the pixels weighted by their counts, so that
  ## a grid cell shows the counts of the gene over all counts in the cell
  mt <- mapData(theta[, topic, drop = FALSE], pos, raster = raster, resolution = resolution)
  mg <- mapData(expr / pmax(total, 1), pos, raster = raster, resolution = resolution,
                weights = if (normalize) total)
  if (!identical(mt$x, mg$x) || !identical(mt$y, mg$y)) {
    stop("Internal error: the maps of the cell-type and of the genes differ.")
  }
  values <- cbind(mt$values, mg$values)
  top <- vapply(seq_len(ncol(values)), function(j) {
    v <- values[!is.na(values[, j]), j]
    q <- if (j == 1) max(v) else stats::quantile(v, 1 - winsorize, names = FALSE)
    ## genes that are expressed in fewer pixels than are winsorized
    if (q <= 0) max(v) else q
  }, numeric(1))
  values <- pmin(t(t(values) / pmax(top, 1e-300)), 1)
  unit <- if (normalize) paste0(signif(100 * top[-1], 2), "% of counts") else
    paste(signif(top[-1], 2), "counts")
  colnames(values) <- c(paste0("Cell-type ", topic, " (max ", format(round(top[1], 2), nsmall = 2), ")"),
                        paste0(genes, " (max ", unit, ")"))
  mt$values <- values
  df <- longMap(mt)
  if (is.null(ncol)) {
    ncol <- ggplot2::wrap_dims(ncol(values))[2]
  }
  p <- ggplot2::ggplot() +
    mapLayer(df, mt$raster, size = size, npanel = ncol) +
    ggplot2::facet_wrap(~ panel, ncol = ncol) +
    ggplot2::scale_fill_gradientn(name = "Relative to\nthe maximum\nof the panel",
                                  colors = colors, limits = c(0, 1),
                                  na.value = NA, aesthetics = c("fill", "colour")) +
    ggplot2::coord_equal() +
    mapTheme()
  if (!showLegend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  if (!is.na(plotTitle)) {
    p <- p + ggplot2::ggtitle(plotTitle)
  }
  p
}


#' Heatmap of the similarity of the cell-types of a model
#'
#' @description Cell-types that are very similar may be one cell-type that
#'     was split because the number of cell-types K is too high.
#'
#' @param results the list returned by `getBetaTheta()`. A beta matrix (for
#'     type "beta") or a theta matrix (for type "theta") is accepted as well.
#' @param type "beta": Pearson's correlation of the gene expression profiles
#'     of the cell-types, on a log scale. "theta": Pearson's correlation of
#'     the proportions of the cell-types across pixels, which shows the
#'     cell-types that occur in the same pixels. (default: "beta")
#' @param cluster whether to order the cell-types by hierarchical clustering,
#'     so that similar cell-types are next to each other (default: TRUE)
#' @param annotation whether to print the correlations in the heatmap
#'     (default: with at most 12 cell-types)
#' @param plotTitle title of the plot (default: depends on `type`)
#'
#' @return a ggplot object
#'
#' @examples
#' data(mOB)
#' counts <- cleanCounts(mOB$counts, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove = 1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = 4)
#' results <- getBetaTheta(optimalModel(ldas, opt = 4))
#' vizTopicCorrelation(results, type = "beta")
#' vizTopicCorrelation(results, type = "theta")
#'
#' @export
vizTopicCorrelation <- function(results, type = c("beta", "theta"), cluster = TRUE,
                                annotation = NULL, plotTitle = NULL) {
  type <- match.arg(type)
  if (type == "beta") {
    beta <- betaOf(results)
    profiles <- log(beta + 0.1 / ncol(beta))
    title <- "Similarity of the gene expression of the cell-types"
  } else {
    profiles <- t(thetaOf(results))
    title <- "Co-occurrence of the cell-types in pixels"
  }
  k <- nrow(profiles)
  mat <- suppressWarnings(stats::cor(t(profiles)))
  ## cell-types that are absent from all pixels have no correlation
  mat[is.na(mat)] <- 0
  diag(mat) <- 1
  ord <- seq_len(k)
  if (cluster && k > 2) {
    ord <- stats::hclust(stats::as.dist(1 - mat), method = "average")$order
  }
  names <- rownames(profiles)[ord]
  mat <- mat[ord, ord, drop = FALSE]
  df <- data.frame(row = factor(rep(names, k), levels = rev(names)),
                   column = factor(rep(names, each = k), levels = names),
                   value = as.vector(mat))
  if (is.null(annotation)) {
    annotation <- k <= 12
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = column, y = row)) +
    ggplot2::geom_tile(ggplot2::aes(fill = value), colour = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradientn(name = "Correlation", colors = diverging_palette,
                                  limits = c(-1, 1), breaks = c(-1, 0, 1)) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(x = "Cell-type", y = "Cell-type") +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = 11, colour = "black"),
                   axis.title = ggplot2::element_text(size = 12))
  if (annotation) {
    df$label <- format(round(df$value, 2), nsmall = 2)
    p <- p + ggplot2::geom_text(ggplot2::aes(label = label), data = df, size = 3,
                                colour = ifelse(abs(df$value) > 0.6, "white", "black"))
  }
  if (is.null(plotTitle)) {
    plotTitle <- title
  }
  if (!is.na(plotTitle)) {
    p <- p + ggplot2::ggtitle(plotTitle)
  }
  p
}


#' The standard set of figures of a deconvolution
#'
#' @description Makes the figures that show a result from different
#'     perspectives, and optionally writes them to a PDF file with one figure
#'     per page:
#' \itemize{
#' \item abundance: the mean proportion of each cell-type (`vizTopicAbundance()`)
#' \item topGenes: the most specific genes of each cell-type (`vizTopGenes()`)
#' \item maps: the proportions of each cell-type in the tissue (`vizTopicMaps()`)
#' \item dominant: the cell-type with the highest proportion in each pixel
#'     (`vizDominantTopic()`)
#' \item pies: the proportions in each pixel as pie charts (`vizAllTopics()`),
#'     only with at most `max.pies` pixels
#' \item similarity: the correlation of the gene expression of the cell-types
#'     (`vizTopicCorrelation()`)
#' \item cooccurrence: the correlation of the proportions of the cell-types
#'     (`vizTopicCorrelation()`)
#' \item markers1, markers2, ...: each cell-type next to its most specific
#'     genes (`vizTopicMarkers()`), only if `counts` are given
#' }
#'
#' @param results the list returned by `getBetaTheta()`
#' @param pos positions of the pixels, as a matrix or data.frame with columns
#'     `x` and `y` and the pixels as row names
#' @param counts gene x pixel count matrix, e.g. the output of
#'     `cleanCounts()`, for the maps of the genes (default: NULL)
#' @param file path of the PDF file to write. NULL: no file is written.
#'     (default: NULL)
#' @param width,height size of the pages of the PDF file in inches
#'     (default: 10 and 8)
#' @param n number of genes of each cell-type in the dot plot (default: 5)
#' @param max.pies largest number of pixels for which the pie charts are
#'     drawn (default: 3000)
#' @param topicCols colors of the cell-types (default: `cellTypeColors(results)`)
#' @param verbose Boolean for verbosity (default: TRUE)
#'
#' @return a named list of ggplot objects, invisibly if a file is written.
#'     Each of them can be changed or saved separately, e.g. with
#'     `ggplot2::ggsave()`.
#'
#' @examples
#' data(mOB)
#' counts <- cleanCounts(mOB$counts, min.lib.size = 100)
#' corpus <- restrictCorpus(counts, removeAbove = 1.0, removeBelow = 0.05)
#' ldas <- fitLDA(t(as.matrix(corpus)), Ks = 3)
#' results <- getBetaTheta(optimalModel(ldas, opt = 3))
#' plots <- vizReport(results, mOB$pos, counts, verbose = FALSE)
#' names(plots)
#' plots$topGenes
#' \dontrun{
#' vizReport(results, mOB$pos, counts, file = "report.pdf")
#' }
#'
#' @export
vizReport <- function(results, pos, counts = NULL, file = NULL,
                      width = 10, height = 8, n = 5, max.pies = 3000,
                      topicCols = cellTypeColors(results), verbose = TRUE) {
  if (!is.list(results) || is.data.frame(results) ||
      is.null(results$theta) || is.null(results$beta)) {
    stop("`results` must be the list with `beta` and `theta` returned by getBetaTheta().")
  }
  theta <- thetaOf(results)
  topicCols <- stats::setNames(as.vector(topicCols), colnames(theta))
  pos <- pos[, c("x", "y"), drop = FALSE]

  plots <- list(
    abundance = vizTopicAbundance(theta, topicCols = topicCols,
                                  plotTitle = "Abundance of the cell-types"),
    topGenes = vizTopGenes(results, n = n,
                           plotTitle = "Most specific genes of each cell-type"),
    maps = vizTopicMaps(theta, pos, scale = "free",
                        plotTitle = "Proportion of each cell-type"),
    dominant = vizDominantTopic(theta, pos, topicCols = topicCols,
                                plotTitle = "Cell-type with the highest proportion in each pixel"))
  if (sum(rownames(theta) %in% rownames(pos)) <= max.pies) {
    plots$pies <- suppressMessages(
      vizAllTopics(theta, pos, topicCols = topicCols, lwd = 0,
                   plotTitle = "Proportions of the cell-types in each pixel") +
        ggplot2::guides(colour = "none"))
  }
  plots$similarity <- vizTopicCorrelation(results, type = "beta")
  plots$cooccurrence <- vizTopicCorrelation(results, type = "theta")
  if (!is.null(counts)) {
    for (topic in colnames(theta)) {
      plots[[paste0("markers", topic)]] <- vizTopicMarkers(
        results, pos, counts, topic = topic,
        plotTitle = paste("Cell-type", topic, "and its most specific genes"))
    }
  }

  if (is.null(file)) {
    return(plots)
  }
  grDevices::pdf(file, width = width, height = height)
  ## the file is closed even if a figure fails
  on.exit(grDevices::dev.off())
  for (p in plots) {
    print(p)
  }
  if (verbose) {
    message("Wrote ", length(plots), " figures to ", file)
  }
  invisible(plots)
}
