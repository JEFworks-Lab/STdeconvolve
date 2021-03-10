# Load individual Rds files into list

# absolute path but maybe way to just use relative path when working from repo?
p <- "/Users/brendan/Desktop/PostDoc/work/STDeconvolve/repos/STdeconvolve/data/"

rds_names <- list.files(path = p, pattern = "*.Rds", full.names = FALSE)
data <- lapply(rds_names, function(f){
  dat <- readRDS(paste0(p,f))
  dat
})
names(data) <- unlist(strsplit(rds_names, ".Rds"))


# alternatively, just load in Rds files as the variables they were originally assigned
# rds_names <- list.files(path = p, pattern = "*.Rds", full.names = TRUE)
# lapply(rds_names, load, .GlobalEnv)