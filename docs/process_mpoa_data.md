---
editor_options: 
  markdown: 
    wrap: 80
---

The following code shows how to generate the gene counts and metadata tables of
individual cells profiled by MERFISH from the Moffit et al. 2018 raw data.

After, the gene counts and metadata tables can be used as input in functions
included in `STdeconvolve` functions in `R/misc.R` (on the `devel` branch) to
generate simulated ST pixels.

## Process raw data

``` {.r}
## Moffit et al. 2018 raw data downloaded from: https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248/
moffit <- read.csv2(file = "Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv",
                    sep = ",")
## split into metadata
annot.table_ <- moffit[,c(1:9)]
## and gene counts
counts_ <- moffit[,-c(1:9)]

## convert gene counts to numerics
counts_ <- as.matrix(data.frame(apply(counts_, 2, function(x) as.numeric(as.character(x)))))

rownames(annot.table_) <- annot.table_[,"Cell_ID"]
rownames(counts_) <- annot.table_[,"Cell_ID"]
```

``` {.r}
## Moffit et al. 2018 Supplemental Table S6 contains the genes profiled
merfish_genes <- readxl::read_excel(path = "aau5324_Moffitt_Table-S6.xlsx",
                                    range = "A2:D163")
## find the smFISH genes
smFISH_genes <- dplyr::filter(merfish_genes, Barcode == "Sequential stain")
```

``` {.r}
## filter for MERFISH profiled genes
counts_ <- counts_[, which(!colnames(counts_) %in% smFISH_genes$`Gene name`)]
## remove Blank measurements
counts_ <- counts_[, !grepl("Blank", colnames(counts_))]
```

``` {.r}
## transform the expression values to discrete counts
## divide by 1000 and multiply by each cell's volume
```

``` {.r}
## "171023_FN7_1_M22_M26" (posterior) and "171021_FN7_2_M22_M26" (anterior) are datasets for Female Naive "Animal_ID" = 2.
## Together they account for all 12 bregma tissue sections.

## collect the pertinent metadata for the cells in the tissue sections
annot.table_ <- annot.table_[annot.table_$Animal_ID == 2,
                             c('Centroid_X', 'Centroid_Y', 'Bregma', "Cell_class", "Neuron_cluster_ID")]
```

``` {.r}
## collapse the Oligodendrocytes and Endothelial cell-types
annot.table_[grep(pattern = "OD Mature",
                                x = annot.table_$Cell_class),]$Cell_class <- "OD Mature"
annot.table_[grep(pattern = "OD Immature",
                                x = annot.table_$Cell_class),]$Cell_class <- "OD Immature"
annot.table_[grep(pattern = "Endothelial",
                                x = annot.table_$Cell_class),]$Cell_class <- "Endothelial"

## remove "Ambiguous" cell-types
annot.table_ <- annot.table_[annot.table_$Cell_class != "Ambiguous",]
dim(annot.table_)
```

    ## [1] 59651     5

``` {.r}
## remove rows with NA
annot.table_ <- na.omit(annot.table_) 
```

``` {.r}
## change the type of some of the metadata columns to numerics
annot.table_$Centroid_X <- as.numeric(as.character(annot.table_$Centroid_X))
annot.table_$Centroid_Y <- as.numeric(as.character(annot.table_$Centroid_Y))
```

``` {.r}
## filter the gene counts matrix to only include the filtered cells
counts_ <- counts_[rownames(annot.table_),]
dim(counts_)
```

    ## [1] 59651   135

## Generate simulated pixels

``` {.r}
library(STdeconvolve)
```

``` {.r}
## get list of major cell-types of each cell
majorCellTypes <- annot.table_$Cell_class
length(majorCellTypes)
```

    ## [1] 59651

``` {.r}
## get a similar list but with neuronal cell-types expanded to subtypes
neuronalCellsubtypes <- unlist(lapply(rownames(annot.table_), function(cell){
  class <- annot.table_[cell,]$Cell_class
  neuron <- annot.table_[cell,]$Neuron_cluster_ID
  if (neuron != ""){
    i <- neuron
  } else {
    i <- class
  }
  i
}))
length(neuronalCellsubtypes)
```

    ## [1] 59651

``` {.r}
## Using the gene counts of the individual cells, generate a ground truth gene expression profile
## for each of the cell-types

## major cell-types:
cellTypes <- majorCellTypes
cells <- rownames(annot.table_)
mat <- counts_[cells,]
mm <- model.matrix(~ 0 + factor(cellTypes))
colnames(mm) <- levels(factor(cellTypes))
gtGexpCellTypes <- t(t(as.matrix(mat)) %*% mm)
gtGexpCellTypes <- gtGexpCellTypes/rowSums(gtGexpCellTypes)
dim(gtGexpCellTypes)
```

    ## [1]   9 135

``` {.r}
## expanded neuronal cell-types:
cellTypes <- neuronalCellsubtypes
cells <- rownames(annot.table_)
mat <- counts_[cells,]
mm <- model.matrix(~ 0 + factor(cellTypes))
colnames(mm) <- levels(factor(cellTypes))
gtGexpNeuronalCellTypes <- t(t(as.matrix(mat)) %*% mm)
gtGexpNeuronalCellTypes <- gtGexpNeuronalCellTypes/rowSums(gtGexpNeuronalCellTypes)
dim(gtGexpNeuronalCellTypes)
```

    ## [1]  76 135

## 100um2 simulated pixels for 9 major cell-types

This generates hash tables of simulated spots for each bregma separately:

``` {.r}
## The function takes into account the cell-types in the "Cell_class" column of the metadata
## so make sure this is set to the cell-types of interest
annot.table_$Cell_class <- majorCellTypes

FN7_hash <- simulateBregmaSpots(annot.table_,
                                counts = counts_, # gene counts matrix
                                patch_size = 100) # size of the simulated pixels in um2 (units of the Centroids)
```

    ## [1] "0.26"
    ## [1] "0.21"
    ## [1] "0.16"
    ## [1] "0.11"
    ## [1] "0.06"
    ## [1] "0.01"
    ## [1] "-0.04"
    ## [1] "-0.09"
    ## [1] "-0.14"
    ## [1] "-0.19"
    ## [1] "-0.24"
    ## [1] "-0.29"

Tables of counts of the 9 major cell-types in each pixel for entire FN7 animal
(100um2):

``` {.r}
# table of number of cell types in each spot for entire FN7
# spot IDs for all spots in FN7
FN7_spotIDs <- unlist(lapply(hash::keys(FN7_hash), function(ix){
  # table to df
  rownames(FN7_hash[[ix]]$cellTypeTable)
}))
# combine all the cellTypeTables of each bregma in FN7 (counts of each cell type in each spot)
FN7_cellTypeTable <- lapply(hash::keys(FN7_hash), function(ix){
  # table to df
  as.data.frame.matrix(FN7_hash[[ix]]$cellTypeTable)
})
# combine into single df, and because some bregmas may be missing cell types,
# use rbindlist to keep all columns and add NAs to spots for cell types
# they are missing
FN7_cellTypeTable <- data.table::rbindlist(FN7_cellTypeTable, fill = TRUE)
# replace NAs with 0s
FN7_cellTypeTable[is.na(FN7_cellTypeTable)] <- 0
# spot IDs as row names
FN7_cellTypeTable <- as.matrix(FN7_cellTypeTable)
rownames(FN7_cellTypeTable) <- FN7_spotIDs
```

Generation of list that contains each simulated corpus of each bregma:

``` {.r}
simBregmasFN7 <- lapply(hash::keys(FN7_hash), function(ix){
  bregma <- buildBregmaCorpus(hashTable = FN7_hash, 
                                bregmaID = ix)
  print(bregma$sim)
  print(dim(bregma$gtSpotTopics))
  print(dim(bregma$gtCtGenes))
  bregma
})
```

    ## A 256x135 simple triplet matrix.
    ## [1] 256   9
    ## [1]   9 135
    ## A 256x135 simple triplet matrix.
    ## [1] 256   9
    ## [1]   9 135
    ## A 256x135 simple triplet matrix.
    ## [1] 256   9
    ## [1]   9 135
    ## A 256x135 simple triplet matrix.
    ## [1] 256   9
    ## [1]   9 135
    ## A 256x135 simple triplet matrix.
    ## [1] 256   9
    ## [1]   9 135
    ## A 256x135 simple triplet matrix.
    ## [1] 256   9
    ## [1]   9 135
    ## A 256x135 simple triplet matrix.
    ## [1] 256   9
    ## [1]   9 135
    ## A 256x135 simple triplet matrix.
    ## [1] 256   9
    ## [1]   9 135
    ## A 256x135 simple triplet matrix.
    ## [1] 256   9
    ## [1]   9 135
    ## A 256x135 simple triplet matrix.
    ## [1] 256   9
    ## [1]   9 135
    ## A 256x135 simple triplet matrix.
    ## [1] 256   8
    ## [1]   8 135
    ## A 256x135 simple triplet matrix.
    ## [1] 256   9
    ## [1]   9 135

``` {.r}
names(simBregmasFN7) <- hash::keys(FN7_hash)
```

Recall that the `$gtSpotTopics` of each simulated corpus are of primary interest
to compare the fitted models to for each bregma.

The `$gtCtGenes` for these corpuses are based only using the cells that were in
the given bregma and in simulated patches and so should be ignored if looking at
the entire animal.

Combine the simulated bregmas in the `simBregmaFN7` list to make a single corpus
for all the bregmas to fit a single model to. The `gtSpotTopics` can be combined
as well for a ground truth reference, but each `gtCtGenes` is built using just
the cells in each given bregma. So instead use the ground truth gene expression
profiles (ex: `gtGexpCellTypes`), which were average gene counts for cell types
across all cells of given type.

``` {.r}
# 1. sim
# combine the sim slam matrices to make the corpus for all spots across all bregmas
sim_N7 <- slam::as.simple_triplet_matrix(do.call(rbind, lapply(names(simBregmasFN7), function(ix){
  m <- as.matrix(simBregmasFN7[[ix]]$sim)
  m
})))
sim_N7
```

    ## A 3072x135 simple triplet matrix.

``` {.r}
# 2. gtSpotTopics
# each gtSpotTopics ref can have different numbers of cell types that are present in each bregma,
# at least for the neuro. So need a way to combine and have all 75 columns,
# and set 0 for patches that do not have any of one ct
gtSpotTopics_N7 <- lapply(names(simBregmasFN7), function(ix){
  simBregmasFN7[[ix]]$gtSpotTopics
})
names(gtSpotTopics_N7) <- names(simBregmasFN7)
gtSpotTopics_N7 <- data.table::rbindlist(gtSpotTopics_N7, fill = TRUE)
dim(gtSpotTopics_N7)
```

    ## [1] 3072    9

``` {.r}
gtSpotTopics_N7 <- as.matrix(gtSpotTopics_N7)
rownames(gtSpotTopics_N7) <- rownames(sim_N7)
# Cts not present in a bregma but columns added here have NA values
# replace NAs with 0's
gtSpotTopics_N7[is.na(gtSpotTopics_N7)] <- 0
gtSpotTopics_N7[1:3,]
```

    ##                           Astrocyte Endothelial Ependymal Excitatory Inhibitory
    ## -1971.778086_2745.402956 0.07142857   0.0000000         0 0.07142857  0.5000000
    ## -1971.778086_2845.402956 0.09090909   0.1818182         0 0.27272727  0.2727273
    ## -1971.778086_2945.402956 0.50000000   0.0625000         0 0.18750000  0.1875000
    ##                           Microglia OD Immature OD Mature Pericytes
    ## -1971.778086_2745.402956 0.07142857   0.1428571 0.1428571         0
    ## -1971.778086_2845.402956 0.00000000   0.0000000 0.1818182         0
    ## -1971.778086_2945.402956 0.00000000   0.0625000 0.0000000         0

``` {.r}
# 3. gtCtGenes
# the `gtCtGenes` is the beta of the average gene expression for each cell cluster,
# in this case using all of cells across the bregma in the animal
dim(gtGexpCellTypes)
```

    ## [1]   9 135

``` {.r}
# 4. cellCounts
# counts of cells in each simulated spot, but also has spot coordinates for easy plotting
cellCounts_N7 <- do.call(rbind, lapply(names(simBregmasFN7), function(ix){
  m <- simBregmasFN7[[ix]]$cellCounts
  m
}))
dim(cellCounts_N7)
```

    ## [1] 3072    3

``` {.r}
cellCounts_N7[1:10,]
```

    ##                                  x        y counts
    ## -1971.778086_2745.402956 -1971.778 2745.403     14
    ## -1971.778086_2845.402956 -1971.778 2845.403     11
    ## -1971.778086_2945.402956 -1971.778 2945.403     16
    ## -1971.778086_3045.402956 -1971.778 3045.403     17
    ## -1971.778086_3145.402956 -1971.778 3145.403     13
    ## -1971.778086_3245.402956 -1971.778 3245.403     19
    ## -1971.778086_3345.402956 -1971.778 3345.403     12
    ## -1971.778086_3445.402956 -1971.778 3445.403     13
    ## -1971.778086_3545.402956 -1971.778 3545.403     15
    ## -1971.778086_3645.402956 -1971.778 3645.403     21

``` {.r}
# 5. annotDf
# recall this is meta data data frame includes information for only the cells that were kept in simulated spots
annotDf_N7 <- do.call(rbind, lapply(names(simBregmasFN7), function(ix){
  m <- simBregmasFN7[[ix]]$annotDf
  m
}))
dim(annotDf_N7)
```

    ## [1] 49142     6

``` {.r}
annotDf_N7[1:10,]
```

    ##                                      Centroid_X Centroid_Y Bregma  Cell_class
    ## 8b2b13ab-f5ea-4127-8451-b984ee5189fc  -3453.744   2830.955  -0.04 Endothelial
    ## c9c9b0be-824b-4796-ae30-e1cdb7e8c2f5  -3431.094   2798.066  -0.04  Inhibitory
    ## 96691c30-a95a-40bc-9be1-fa743f2d462f  -3425.506   2808.064  -0.04  Inhibitory
    ## 51233fd6-29ef-4871-8922-5ac929dd4263  -3422.492   2759.921  -0.04 Endothelial
    ## 8be0499e-320f-41b4-bb8c-b6b807d15ef7  -3420.258   2826.501  -0.04 Endothelial
    ## da2659c2-5a97-41a5-a463-9aec5408a9f6  -3406.136   2813.243  -0.04  Inhibitory
    ## edcfafd3-d83a-4485-8c06-8315393044f5  -3404.233   2777.500  -0.04  Inhibitory
    ## 88d1b729-6000-40cc-b582-417573a9b64f  -3393.321   2821.759  -0.04  Excitatory
    ## 61e14c8a-772e-40a5-aea0-e40a4a9a4eab  -3392.659   2803.046  -0.04 Endothelial
    ## eec3145b-79cc-497c-b807-b2902554c3e5  -3392.615   2789.444  -0.04   Astrocyte
    ##                                      Neuron_cluster_ID                 patch_id
    ## 8b2b13ab-f5ea-4127-8451-b984ee5189fc                   -3471.778086_2745.402956
    ## c9c9b0be-824b-4796-ae30-e1cdb7e8c2f5               I-5 -3471.778086_2745.402956
    ## 96691c30-a95a-40bc-9be1-fa743f2d462f               I-6 -3471.778086_2745.402956
    ## 51233fd6-29ef-4871-8922-5ac929dd4263                   -3471.778086_2745.402956
    ## 8be0499e-320f-41b4-bb8c-b6b807d15ef7                   -3471.778086_2745.402956
    ## da2659c2-5a97-41a5-a463-9aec5408a9f6              I-20 -3471.778086_2745.402956
    ## edcfafd3-d83a-4485-8c06-8315393044f5               I-5 -3471.778086_2745.402956
    ## 88d1b729-6000-40cc-b582-417573a9b64f              E-20 -3471.778086_2745.402956
    ## 61e14c8a-772e-40a5-aea0-e40a4a9a4eab                   -3471.778086_2745.402956
    ## eec3145b-79cc-497c-b807-b2902554c3e5                   -3471.778086_2745.402956

``` {.r}
# construct the list similar to the output of `buildBregmaCorpus()`
simFN7 <- list(sim = sim_N7,
               gtSpotTopics = gtSpotTopics_N7,
               gtCtGenes = gtGexpCellTypes,
               cellCounts = cellCounts_N7,
               # classColors = classColors,
               annotDf = annotDf_N7)
```

`simFN7` is the simulated dataset for 9 major cell-types in 100um2 pixels across
all tissue sections of the animal.

`simFN7$sim`: A 3072x135 simple triplet matrix (pixels x genes count matrix).
Input for STdeconvolve

`simFN7$gtSpotTopics`: ground truth cell-type pixel proportions for 9 major
cell-types

`simFN7$gtCtGenes`: ground truth cell-type gene expression profiles

`simFN7$cellCounts`: number of cells in each simulated pixel

`simFN7$annotDf`: metadata table of cells kept in simulated pixels
