## Corpora (pixels x genes, sparse) used by the benchmark scripts.
##  mob:  mouse olfactory bulb data of the package, 260 pixels x 232 genes
##  sim:  simulated mixtures of 10 cell-types, 3000 pixels x 1000 genes,
##        about 5000 counts per pixel
##  simlow: simulated mixtures of 8 cell-types, 5000 pixels x 500 genes,
##        about 300 counts per pixel
bench_corpus <- function(name = c("mob", "sim", "simlow")) {
  name <- match.arg(name)
  if (name == "mob") {
    e <- new.env(); utils::data("mOB", package = "STdeconvolve", envir = e)
    counts <- cleanCounts(e$mOB$counts, min.lib.size = 100)
    corpus <- suppressMessages(restrictCorpus(counts, removeAbove = 1, removeBelow = 0.05,
                                              verbose = FALSE))
    return(Matrix::t(corpus))
  }
  par <- list(sim = list(N = 3000, G = 1000, K = 10, depth = 5000),
              simlow = list(N = 5000, G = 500, K = 8, depth = 300))[[name]]
  set.seed(1)
  beta <- matrix(stats::rgamma(par$K * par$G, 0.1), par$K); beta <- beta / rowSums(beta)
  theta <- matrix(stats::rgamma(par$N * par$K, 0.3), par$N); theta <- theta / rowSums(theta)
  x <- matrix(stats::rpois(par$N * par$G, par$depth * theta %*% beta), par$N, par$G,
              dimnames = list(paste0("p", seq_len(par$N)), paste0("g", seq_len(par$G))))
  keep <- rowSums(x) > 0
  x <- methods::as(x[keep, ], "CsparseMatrix")
  attr(x, "truth") <- list(beta = beta, theta = theta[keep, ])
  x
}
