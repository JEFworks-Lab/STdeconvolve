## Compare the C++ fitter with the topicmodels reference fits.
##   R CMD INSTALL --library=/tmp/claude-stlib . && Rscript bench/compare.R [ncores]
suppressMessages(library(STdeconvolve, lib.loc = Sys.getenv("STLIB", "/tmp/claude-stlib")))
source("bench/datasets.R")
ncores <- as.integer(c(commandArgs(TRUE), 1)[1])
ref <- readRDS("bench/reference_topicmodels.rds")
ldaVEM <- getFromNamespace("ldaVEM", "STdeconvolve")
corpora <- list()
cat(sprintf("%-8s %7s %7s %6s | %9s %9s | %8s %8s | %4s %4s | %9s %9s\n", "", "ref_s", "new_s", "x",
            "ref_perp", "new_perp", "ref_alph", "new_alph", "it", "it", "beta_cor", "theta_cor"))
for (r in ref) {
  if (is.null(corpora[[r$dataset]])) corpora[[r$dataset]] <- bench_corpus(r$dataset)
  x <- corpora[[r$dataset]]
  el <- system.time(m <- ldaVEM(x, k = r$k, seed = 0, ncores = ncores))[["elapsed"]]
  ## match topics by correlation of their gene distributions
  cm <- cor(t(exp(m@beta)), t(r$beta))
  pr <- as.integer(clue::solve_LSAP(pmax(cm, 0), maximum = TRUE))
  bc <- min(cm[cbind(seq_len(r$k), pr)])
  tc <- min(diag(cor(m@gamma, r$gamma[, pr])))
  cat(sprintf("%-8s %7.2f %7.3f %6.0f | %9.3f %9.3f | %8.4f %8.4f | %4d %4d | %9.5f %9.5f %s\n",
              paste0(r$dataset, " K=", r$k), r$elapsed, el, r$elapsed / el,
              r$perplexity, exp(-sum(m@loglikelihood) / m@n), r$alpha, m@alpha, r$iter, m@iter, bc, tc,
              if (all(pr == seq_len(r$k))) "" else "(permuted)"))
}
