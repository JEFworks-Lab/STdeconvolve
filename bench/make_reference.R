## Fit reference models with topicmodels::LDA (the pre-C++ implementation) so the
## new fitter can be compared against them. Run once from the package root:
##   Rscript bench/make_reference.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("bench/datasets.R")
ctl <- list(seed = 0, verbose = 0, keep = 1, estimate.alpha = TRUE)
ref <- list()
for (nm in c("mob", "sim")) {
  x <- bench_corpus(nm)
  stm <- slam::as.simple_triplet_matrix(as.matrix(x))
  Ks <- if (nm == "mob") c(2, 3, 4, 5, 8, 12) else c(6, 10)
  for (k in Ks) {
    el <- system.time(m <- topicmodels::LDA(stm, k = k, control = ctl))[["elapsed"]]
    message(sprintf("%s K=%d  %.1fs  perplexity %.3f alpha %.4f iter %d",
                    nm, k, el, topicmodels::perplexity(m), m@alpha, m@iter))
    ref[[paste(nm, k, sep = "_")]] <- list(
      dataset = nm, k = k, elapsed = el, alpha = m@alpha, iter = m@iter,
      perplexity = topicmodels::perplexity(m), logLik = sum(m@loglikelihood),
      beta = exp(m@beta), gamma = m@gamma)
  }
}
saveRDS(ref, "bench/reference_topicmodels.rds")
