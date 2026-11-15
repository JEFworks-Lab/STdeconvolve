#' @useDynLib STdeconvolve, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL


#' Coerce a pixels x genes count matrix to a sparse corpus and validate it
#'
#' @param counts pixels (rows) x genes (columns) counts. A base matrix, a
#'     data.frame, a sparse matrix from the Matrix package or a
#'     slam::simple_triplet_matrix.
#' @param what name of the argument, used in error messages
#'
#' @return a dgCMatrix without explicit zeros
#'
#' @noRd
asCorpus <- function(counts, what="counts") {
  if (inherits(counts, "simple_triplet_matrix")) {
    counts <- Matrix::sparseMatrix(i=counts$i, j=counts$j, x=as.numeric(counts$v),
                                   dims=c(counts$nrow, counts$ncol),
                                   dimnames=if (is.null(counts$dimnames)) list(NULL, NULL) else counts$dimnames)
  } else if (methods::is(counts, "Matrix")) {
    counts <- methods::as(methods::as(methods::as(counts, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  } else {
    counts <- as.matrix(counts)
    storage.mode(counts) <- "double"
    counts <- methods::as(counts, "CsparseMatrix")
  }
  counts <- Matrix::drop0(counts)
  
  if (any(counts@x != floor(counts@x)) || any(counts@x < 0)) {
    stop("`", what, "` must contain non-negative integer gene counts")
  }
  if (any(Matrix::rowSums(counts) == 0)) {
    stop("Each row (pixel) of `", what, "` needs to contain at least one non-zero entry")
  }
  counts
}


#' Fitted LDA model
#'
#' @description Latent Dirichlet allocation model fitted by variational EM, as
#'     returned by [fitLDA()]. Use [getBetaTheta()] to get the cell-type gene
#'     expression profiles and the cell-type proportions of the pixels.
#'     The slots have the same names and content as those of the "LDA_VEM"
#'     class of the `topicmodels` package.
#'
#' @slot call the call that fitted the model
#' @slot Dim number of pixels and genes of the corpus
#' @slot control list of the settings used to fit the model
#' @slot k number of cell-types
#' @slot terms names of the genes
#' @slot documents names of the pixels
#' @slot beta cell-types x genes matrix, logarithm of the probability of each
#'     gene in each cell-type
#' @slot gamma pixels x cell-types matrix of cell-type proportions
#' @slot loglikelihood lower bound of the log-likelihood of each pixel
#' @slot iter number of EM iterations
#' @slot logLiks lower bound of the log-likelihood of the corpus at each EM iteration
#' @slot n total number of counts in the corpus
#' @slot alpha parameter of the Dirichlet prior of the cell-type proportions
#'
#' @exportClass STdeconvolveLDA
methods::setClass("STdeconvolveLDA",
                  representation(call="call", Dim="integer", control="list",
                                 k="integer", terms="character", documents="character",
                                 beta="matrix", gamma="matrix",
                                 loglikelihood="numeric", iter="integer",
                                 logLiks="numeric", n="integer", alpha="numeric"))

methods::setMethod("show", "STdeconvolveLDA", function(object) {
  cat("LDA model with ", object@k, " cell-types fitted to ", object@Dim[1],
      " pixels and ", object@Dim[2], " genes (alpha = ", signif(object@alpha, 4),
      ", ", object@iter, " EM iterations)\n", sep="")
})


#' Check that an object is a fitted LDA model
#'
#' @description Models are "STdeconvolveLDA" objects. "LDA_VEM" objects of
#'     the `topicmodels` package, as returned by previous versions of
#'     fitLDA(), have the same slots and are accepted as well.
#'
#' @param lda object to check
#'
#' @return `lda`, invisibly
#'
#' @noRd
checkModel <- function(lda) {
  slots <- c("k", "alpha", "beta", "gamma", "terms", "documents", "Dim")
  if (!isS4(lda) || !all(slots %in% names(attributes(lda)))) {
    stop("`lda` must be a model fitted with fitLDA()")
  }
  invisible(lda)
}


#' Fit an LDA model, or infer cell-type proportions with a fitted model
#'
#' @description Variational EM for latent Dirichlet allocation, implemented in
#'     C++. It fits the same model, with the same initialization for a given
#'     seed and the same stopping rule for the EM iterations, as
#'     `topicmodels::LDA(method = "VEM")`, but is much faster and can use
#'     several threads.
#'
#' @param corpus pixels x genes dgCMatrix, see asCorpus()
#' @param k number of cell-types
#' @param model a fitted model. If given, its cell-type gene
#'     distributions (and alpha) are kept fixed and only the cell-type
#'     proportions of the pixels in `corpus` are inferred.
#' @param seed seed for the random initialization
#' @param alpha initial value of the Dirichlet parameter of the cell-type
#'     proportions (default: 50/k)
#' @param estimate.alpha estimate alpha (default: TRUE) or keep it fixed
#' @param em.max.iter,em.tol maximum number of EM iterations and tolerance on
#'     the relative change of the likelihood bound
#' @param var.max.iter,var.tol maximum number of iterations to optimize the
#'     cell-type proportions of a single pixel and tolerance on the change
#'     of these proportions (sum of absolute differences) between iterations
#' @param ncores number of threads
#'
#' @return an "STdeconvolveLDA" object
#'
#' @noRd
ldaVEM <- function(corpus, k, model=NULL, seed=0, alpha=50/k, estimate.alpha=TRUE,
                   em.max.iter=1000L, em.tol=1e-4,
                   var.max.iter=500L, var.tol=1e-6,
                   ncores=1) {
  
  if (!is.null(model)) {
    k <- model@k
    alpha <- model@alpha
    estimate.alpha <- FALSE
    em.max.iter <- 0L
    init <- exp(model@beta)
  } else {
    init <- NULL
  }
  
  fit <- lda_vem_cpp(corpus@p, corpus@i, corpus@x,
                     nrow(corpus), ncol(corpus), as.integer(k),
                     alpha, estimate.alpha, init, as.integer(seed),
                     as.integer(em.max.iter), em.tol,
                     as.integer(var.max.iter), var.tol,
                     as.integer(max(1, ncores)))
  
  control <- list(seed=as.integer(seed), alpha=alpha, estimate.alpha=estimate.alpha,
                  em=list(iter.max=as.integer(em.max.iter), tol=em.tol),
                  var=list(iter.max=as.integer(var.max.iter), tol=var.tol))
  
  methods::new("STdeconvolveLDA",
               call=match.call(),
               Dim=dim(corpus),
               control=control,
               k=as.integer(k),
               terms=as.character(colnames(corpus)),
               documents=as.character(rownames(corpus)),
               beta=fit$beta,
               gamma=fit$gamma / rowSums(fit$gamma),
               loglikelihood=fit$loglikelihood,
               iter=fit$iter,
               logLiks=fit$logLiks,
               n=as.integer(sum(corpus@x)),
               alpha=fit$alpha)
}


#' Arrange the genes of a new corpus to match the genes of a fitted model
#'
#' @param corpus pixels x genes dgCMatrix
#' @param model fitted LDA model
#'
#' @return dgCMatrix with one column for each gene in the model. Genes that
#'     are not in the model are dropped, genes missing in the corpus are 0.
#'
#' @noRd
matchGenes <- function(corpus, model) {
  genes <- model@terms
  if (length(genes) == 0 || is.null(colnames(corpus))) {
    if (ncol(corpus) != model@Dim[2]) {
      stop("The number of genes in `corpus` and in the fitted model need to match")
    }
    return(corpus)
  }
  if (identical(colnames(corpus), genes)) {
    return(corpus)
  }
  corpus <- corpus[, colnames(corpus) %in% genes, drop=FALSE]
  tm <- methods::as(corpus, "TsparseMatrix")
  Matrix::sparseMatrix(i=tm@i + 1L, j=match(colnames(corpus), genes)[tm@j + 1L], x=tm@x,
                       dims=c(nrow(corpus), length(genes)),
                       dimnames=list(rownames(corpus), genes))
}
