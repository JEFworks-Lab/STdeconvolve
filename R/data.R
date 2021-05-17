#' Spatial transcriptomics transcriptomics of the mouse olfactory bulb
#'
#' @format List where 'counts' is a sparse matrix with columns as voxels and rows as genes and
#'                    'pos' is a data frame of x and y position values per voxel
#'
#' @source \url{https://science.sciencemag.org/content/353/6294/78}
"mOB"

#' MERFISH data of the mouse pre-optic region for a female naive animal (FN7)
#'
#' @format List where 'mat' is a sparse matrix with columns as cells and rows as genes
#'                          where expression values have already been normalized by volume
#'                    'pos' is a data frame of x, y, z position values per cell
#'                          and brain position as 6 slice indices from anterior to posterior
#'
#' @source \url{https://science.sciencemag.org/content/362/6416/eaau5324/}
"mPOA"

#' Spatial transcriptomics transcriptomics of 4 breast cancer biopsy sections
#'
#' @format List where 'counts' is a sparse matrix with columns as voxels and rows as genes and
#'                    'pos' is a data frame of x and y position values per voxel
#'                          and slice index for 4 consecutive slices
#'
#' @source \url{https://science.sciencemag.org/content/353/6294/78}
"BCL"


# Load individual Rds files into list

# absolute path but maybe way to just use relative path when working from repo?
# p <- "/Users/brendan/Desktop/PostDoc/work/STDeconvolve/repos/STdeconvolve/data/"

# rds_names <- list.files(path = p, pattern = "*.Rds", full.names = FALSE)
# data <- lapply(rds_names, function(f){
#   dat <- readRDS(paste0(p,f))
#   dat
# })
# names(data) <- unlist(strsplit(rds_names, ".Rds"))


# alternatively, just load in Rds files as the variables they were originally assigned
# rds_names <- list.files(path = p, pattern = "*.Rds", full.names = TRUE)
# lapply(rds_names, load, .GlobalEnv)
