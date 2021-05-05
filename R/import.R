## from https://rdrr.io/cran/gama/src/R/bestk.R
# Calculates an aproximation of the second derivative of a set of points
# the maximum second derivative will be a good choice for the inflexion point (the elbow or knee)
# https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
# https://raghavan.usc.edu/papers/kneedle-simplex11.pdf (Finding a "Kneedle" in a Haystack:
# Detecting Knee Points in System Behavior)
where.is.knee <- function(points = NULL) {

  lower.limit <- 2
  upper.limit <- length(points) -1

  second.derivative <- sapply(lower.limit:upper.limit, function(i) { points[i+1] + points[i-1] - (2 * points[i]) })

  w.max <- which.max(second.derivative)
  best.k <- w.max +1

  return(best.k)

}


#' Filter a counts matrix
#'
#' @description Filter a counts matrix based on gene (row) and cell (column)
#'      requirements.
#'
#' @param counts A sparse read count matrix. The rows correspond to genes,
#'      columns correspond to individual cells
#' @param min.lib.size Minimum number of genes detected in a cell. Cells with
#'      fewer genes will be removed (default: 1)
#' @param max.lib.size Maximum number of genes detected in a cell. Cells with
#'      more genes will be removed (default: Inf)
#' @param min.reads Minimum number of reads per gene. Genes with fewer reads
#'      will be removed (default: 1)
#' @param min.detected Minimum number of cells a gene must be seen in. Genes
#'      not seen in a sufficient number of cells will be removed (default: 1)
#' @param verbose Verbosity (default: FALSE)
#' @param plot Whether to plot (default: TRUE)
#'
#' @return a filtered read count matrix
#'
#' @export
#'
#' @importFrom Matrix Matrix colSums rowSums
#'
cleanCounts <- function (counts, min.lib.size = 1, max.lib.size = Inf, min.reads = 1, min.detected = 1, verbose = FALSE, plot=TRUE) {
  if (!any(class(counts) %in% c("dgCMatrix", "dgTMatrix"))) {
    if (verbose) {
      message("Converting to sparse matrix ...")
    }
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }
  if (verbose) {
    message("Filtering matrix with ", ncol(counts), " cells and ",
            nrow(counts), " genes ...")
  }
  ix_col <- Matrix::colSums(counts)
  ix_col <- ix_col > min.lib.size & ix_col < max.lib.size
  counts <- counts[, ix_col]
  counts <- counts[Matrix::rowSums(counts) > min.reads, ]
  counts <- counts[Matrix::rowSums(counts > 0) > min.detected, ]
  if (verbose) {
    message("Resulting matrix has ", ncol(counts), " cells and ", nrow(counts), " genes")
  }
  if (plot) {
    par(mfrow=c(1,2), mar=rep(5,4))
    hist(log10(Matrix::colSums(counts)+1), breaks=20, main='Genes Per Dataset')
    hist(log10(Matrix::rowSums(counts)+1), breaks=20, main='Datasets Per Gene')
  }
  return(counts)
}

#' Normalize gene expression variance relative to transcriptome-wide expectations
#' (Modified from SCDE/PAGODA2 code)
#'
#' Normalizes gene expression magnitudes to with respect to its ratio to the
#' transcriptome-wide expectation as determined by local regression on expression magnitude
#'
#' @param counts Read count matrix. The rows correspond to genes, columns correspond to individual cells
#' @param gam.k Generalized additive model parameter; the dimension of the basis used to represent the smooth term (default: 5)
#' @param alpha Significance threshold (default: 0.05)
#' @param plot Whether to plot the results (default: FALSE)
#' @param use.unadjusted.pvals If true, will apply BH correction (default: FALSE)
#' @param do.par Whether to adjust par for plotting if plotting (default: TRUE)
#' @param max.adjusted.variance Ceiling on maximum variance after normalization to prevent infinites (default: 1e3)
#' @param min.adjusted.variance Floor on minimum variance after normalization (default: 1e-3)
#' @param verbose Verbosity (default: TRUE)
#' @param details If true, will return data frame of normalization parameters. Else will return list of overdispersed genes. (default: FALSE)
#'
#' @return If details is true, will return data frame of normalization parameters. Else will return list of overdispersed genes.
#'
#' @importFrom mgcv s
#'
#' @export
#'
getOverdispersedGenes <- function(counts, gam.k=5, alpha=0.05, plot=FALSE, use.unadjusted.pvals=FALSE, do.par=TRUE, max.adjusted.variance=1e3, min.adjusted.variance=1e-3, verbose=TRUE, details=FALSE) {

  if (!any(class(counts) %in% c("dgCMatrix", "dgTMatrix"))) {
    if(verbose) {
      message('Converting to sparse matrix ...')
    }
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }

  mat <- Matrix::t(counts) ## make rows as cells, cols as genes

  if(verbose) {
    print("Calculating variance fit ...")
  }
  dfm <- log(Matrix::colMeans(mat))
  dfv <- log(apply(mat, 2, var))
  names(dfm) <- names(dfv) <- colnames(mat)
  df <- data.frame(m=dfm, v=dfv)

  vi <- which(is.finite(dfv))

  if(length(vi)<gam.k*1.5) { gam.k=1 } ## too few genes

  if(gam.k<2) {
    if(verbose) {
      print("Using lm ...")
    }
    m <- lm(v ~ m, data = df[vi,])
  } else {
    if(verbose) {
      print(paste0("Using gam with k=", gam.k, "..."))
    }
    fm <- as.formula(sprintf("v ~ s(m, k = %s)", gam.k))
    m <- mgcv::gam(fm, data = df[vi,])
  }
  df$res <- -Inf;  df$res[vi] <- resid(m,type='response')
  n.cells <- ncol(mat)
  n.obs <- nrow(mat)
  df$lp <- as.numeric(pf(exp(df$res),n.obs,n.obs,lower.tail=F,log.p=T))
  df$lpa <- bh.adjust(df$lp,log=TRUE)
  df$qv <- as.numeric(qchisq(df$lp, n.cells-1, lower.tail = FALSE,log.p=TRUE)/n.cells)

  if(use.unadjusted.pvals) {
    ods <- which(df$lp<log(alpha))
  } else {
    ods <- which(df$lpa<log(alpha))
  }
  if(verbose) {
    print(paste0(length(ods), ' overdispersed genes ... ' ))
  }

  df$gsf <- geneScaleFactors <- sqrt(pmax(min.adjusted.variance,pmin(max.adjusted.variance,df$qv))/exp(df$v));
  df$gsf[!is.finite(df$gsf)] <- 0;

  if(plot) {
    if(do.par) {
      par(mfrow=c(1,2), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
    }
    smoothScatter(df$m,df$v,main='',xlab='log10[ magnitude ]',ylab='log10[ variance ]')
    grid <- seq(min(df$m[vi]),max(df$m[vi]),length.out=1000)
    lines(grid,predict(m,newdata=data.frame(m=grid)),col="blue")
    if(length(ods)>0) {
      points(df$m[ods],df$v[ods],pch='.',col=2,cex=1)
    }
    smoothScatter(df$m[vi],df$qv[vi],xlab='log10[ magnitude ]',ylab='',main='adjusted')
    abline(h=1,lty=2,col=8)
    if(is.finite(max.adjusted.variance)) { abline(h=max.adjusted.variance,lty=2,col=1) }
    points(df$m[ods],df$qv[ods],col=2,pch='.')
  }

  ## variance normalize
  norm.mat <- counts*df$gsf
  if(!details) {
    return(rownames(mat)[ods])
  } else {
    ## return normalization factor
    return(list(mat=norm.mat, ods=colnames(mat)[ods], df=df))
  }
}

# Multiple testing correction
bh.adjust <- function(x, log = FALSE) {
  nai <- which(!is.na(x))
  ox <- x
  x <- x[nai]
  id <- order(x, decreasing = FALSE)
  if(log) {
    q <- x[id] + log(length(x)/seq_along(x))
  } else {
    q <- x[id]*length(x)/seq_along(x)
  }
  a <- rev(cummin(rev(q)))[order(id)]
  ox[nai] <- a
  ox
}


#' Normalizes counts to CPM
#'
#' @description Normalizes raw counts to log10 counts per million with pseudocount
#'
#' @param counts Read count matrix. The rows correspond to genes, columns
#'      correspond to individual cells
#' @param normFactor Normalization factor such as cell size. If not provided
#'      column sum as proxy for library size will be used
#' @param depthScale Depth scaling. Using a million for CPM (default: 1e6)
#' @param pseudo Pseudocount for log transform (default: 1)
#' @param log Whether to apply log transform
#' @param verbose Verbosity (default: TRUE)
#'
#' @return a normalized matrix
#'
#' @export
#'
#' @importFrom Matrix Matrix colSums t
#'
normalizeCounts <- function (counts, normFactor = NULL, depthScale = 1e+06, pseudo = 1,
                             log = TRUE, verbose = TRUE) {
  if (!class(counts) %in% c("dgCMatrix", "dgTMatrix")) {
    if (verbose) {
      message("Converting to sparse matrix ...")
    }
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }
  if (verbose) {
    message("Normalizing matrix with ", ncol(counts), " cells and ",
            nrow(counts), " genes.")
  }
  if (is.null(normFactor)) {
    if (verbose) {
      message("normFactor not provided. Normalizing by library size.")
    }
    normFactor <- Matrix::colSums(counts)
  }
  if (verbose) {
    message(paste0("Using depthScale ", depthScale))
  }
  counts <- Matrix::t(Matrix::t(counts)/normFactor)
  counts <- counts * depthScale
  if (log) {
    if (verbose) {
      message("Log10 transforming with pseudocount ",
              pseudo, ".")
    }
    counts <- log10(counts + pseudo)
  }
  return(counts)
}

# from factors to colors
fac2col <- function(x,s=1,v=1,shuffle=FALSE,min.group.size=1,return.details=F,unclassified.cell.color='lightgrey',level.colors=NULL) {
  x <- as.factor(x);
  if(min.group.size>1) {
    x <- factor(x,exclude=levels(x)[unlist(tapply(rep(1,length(x)),x,length))<min.group.size])
    x <- droplevels(x)
  }
  if(is.null(level.colors)) {
    col <- rainbow(length(levels(x)),s=s,v=v);
  } else {
    col <- level.colors[1:length(levels(x))];
  }
  names(col) <- levels(x);

  if(shuffle) col <- sample(col);

  y <- col[as.integer(x)]; names(y) <- names(x);
  y[is.na(y)] <- unclassified.cell.color;
  if(return.details) {
    return(list(colors=y,palette=col))
  } else {
    return(y);
  }
}
