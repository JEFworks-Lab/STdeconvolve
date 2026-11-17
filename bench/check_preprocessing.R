## Check that the preprocessing functions give the same results as a reference
## installation of the package (e.g. the previous release), and time both.
##   Rscript bench/check_preprocessing.R <lib with new version> <lib with old version>
libs <- commandArgs(TRUE)
run <- function(lib) {
  code <- '
    suppressMessages(library(STdeconvolve, lib.loc = commandArgs(TRUE)[1]))
    data(mOB); set.seed(1); N <- 3000; G <- 20000; K <- 10
    beta <- matrix(rgamma(K * G, 0.1), K); beta <- beta / rowSums(beta)
    theta <- matrix(rgamma(N * K, 0.3), N); theta <- theta / rowSums(theta)
    big <- matrix(rpois(N * G, 5000 * theta %*% beta), G, N, byrow = TRUE,
                  dimnames = list(paste0("g", 1:G), paste0("p", 1:N)))
    out <- list(); tm <- c()
    for (nm in c("mob", "big")) {
      cd <- if (nm == "mob") mOB$counts else big
      t0 <- function(lab, e) { el <- system.time(v <- e)[["elapsed"]]; tm[[paste(nm, lab)]] <<- el; v }
      cc <- t0("cleanCounts", cleanCounts(cd, min.lib.size = 100))
      out[[nm]] <- list(
        clean = as.matrix(cc),
        od = t0("getOverdispersedGenes", getOverdispersedGenes(cc, verbose = FALSE, details = TRUE)$df),
        restrict = as.matrix(t0("restrictCorpus", suppressMessages(restrictCorpus(cc, removeAbove = 1, removeBelow = 0.05, verbose = FALSE)))),
        restrict10 = rownames(suppressMessages(restrictCorpus(cc, removeAbove = 0.9, removeBelow = 0.1, nTopOD = 10, verbose = FALSE))),
        pre = t0("preprocess", suppressMessages(preprocess(t(cd), min.lib.size = 100, removeAbove = 0.95, removeBelow = 0.05,
                                           nTopGenes = 5, plot = FALSE, verbose = FALSE))[c("corpus", "slm")]))
    }
    th <- matrix(runif(3000 * 12), 3000, dimnames = list(paste0("p", 1:3000), 1:12))
    be <- matrix(runif(12 * 1000), 12, dimnames = list(1:12, paste0("g", 1:1000)))
    out$cor_t <- suppressMessages(getCorrMtx(th, th[, 1:8], type = "t"))
    out$cor_b <- suppressMessages(getCorrMtx(be, be[1:8, ], type = "b"))
    out$cor_bt <- suppressMessages(getCorrMtx(be, be[1:8, ], type = "b", thresh = 0.3))
    saveRDS(list(out = out, tm = unlist(tm)), commandArgs(TRUE)[2])'
  f <- tempfile(fileext = ".rds"); s <- tempfile(fileext = ".R"); writeLines(code, s)
  system2(file.path(R.home("bin"), "Rscript"), c(s, lib, f), stdout = FALSE, stderr = FALSE)
  readRDS(f)
}
new <- run(libs[1]); old <- run(libs[2])
print(round(cbind(old_s = old$tm, new_s = new$tm, speedup = old$tm / new$tm), 2))
cat("results equal:", isTRUE(all.equal(old$out, new$out, tolerance = 1e-10)), "\n")
if (!isTRUE(all.equal(old$out, new$out, tolerance = 1e-10))) print(all.equal(old$out, new$out, tolerance = 1e-10))
