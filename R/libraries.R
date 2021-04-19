library(Matrix)
library(MERINGUE)
library(MUDAN) # some MERINGUE functions actually from MUDAN (?) cleanCounts, getDifferentialGenes, normalizeCounts, plotEmbedding
library(ggplot2)
library(scatterpie)
library(viridis) # for "scale_fill_viridis in `vizGeneCounts`
library(dynamicTreeCut) # for clustering topics
library(slam) # for formatting matrix for input into LDA 


# for figures but not for actual package
library(Seurat) # for building scRNAseq references
library(hash) # for organizing and processing MERFISH bregmas
library(gplots) # for `heatmap.2`
library(SPOTlight) # for comparison to SPOTlight
library(RCTD) # for comparison to RCTD
library(mltools) # for rmse
library(randomcoloR) # for distinct color palettes
library(ggpubr) # for stats on ggplots
library(jpeg) # for viz of BCL H&E stain jpegs

# for GSEA
library(liger)
library(GO.db)
library(AnnotationDbi)
library(msigdbr)
