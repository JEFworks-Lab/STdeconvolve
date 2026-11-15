## The accelerated E-step must find the same models as the plain fixed-point
## iteration. Compares both over datasets, numbers of cell-types and seeds.
##   Rscript bench/check_acceleration.R [ncores]
suppressMessages(library(STdeconvolve, lib.loc = Sys.getenv("STLIB", "/tmp/claude-stlib")))
source("bench/datasets.R")
ncores <- as.integer(c(commandArgs(TRUE), 1)[1])
cpp <- getFromNamespace("lda_vem_cpp", "STdeconvolve")
cases <- list(mob = c(3, 5, 8, 12, 20), sim = c(6, 10, 15), simlow = c(4, 8, 12))
worst <- 0
for (nm in names(cases)) {
  x <- bench_corpus(nm); n <- sum(x@x)
  for (k in cases[[nm]]) for (seed in c(0L, 2L, 4L)) {
    fit <- function(accel) {
      el <- system.time(r <- cpp(x@p, x@i, x@x, nrow(x), ncol(x), k, 50 / k, TRUE, NULL, seed,
                                 1000L, 1e-4, 500L, 1e-6, ncores, accel))[["elapsed"]]
      c(s = el, updates = r$var_iters / nrow(x) / (r$iter + 1), em = r$iter,
        perp = exp(-sum(r$loglikelihood) / n), alpha = r$alpha)
    }
    a <- fit(TRUE); p <- fit(FALSE)
    d <- a[["perp"]] / p[["perp"]] - 1; worst <- max(worst, d)
    cat(sprintf("%-6s K=%2d seed %d | plain %6.2fs %5.0f upd EM %3.0f perp %9.4f alpha %.4f | accel %6.2fs %5.0f upd EM %3.0f perp %9.4f alpha %.4f | perp diff %+.4f%%\n",
                nm, k, seed, p[["s"]], p[["updates"]], p[["em"]], p[["perp"]], p[["alpha"]],
                a[["s"]], a[["updates"]], a[["em"]], a[["perp"]], a[["alpha"]], 100 * d))
  }
}
cat(sprintf("largest increase of perplexity with acceleration: %.4f%%\n", 100 * worst))
