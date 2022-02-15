The following code demonstrates how to create the BCL ST input data matrix for `STdeconvolve`

```{r}
library(STdeconvolve)
```

First load the BCL data. Note that this data object is stored on the [`devel` branch of `STdeconvolve`](https://github.com/JEFworks-Lab/STdeconvolve/tree/devel/data)

```{r}
data("BCL") # the 4 breast cancer sections
```

```{r}
## 4 layers, 1031 total spots.
bcl_pos <- BCL$pos[,1:2]
bcl_slice <- BCL$pos[,3]
names(bcl_slice) <- rownames(bcl_pos)
length(bcl_slice)
# clean up poor spots and genes
bclClean <- cleanCounts(counts = BCL$counts, # genes x spots mtx
                                 min.reads = 10,
                                 min.lib.size = 10,
                                 verbose=TRUE)

## merge the 4 sections together
bcl_paths <- list()
bcl_paths[[1]] <- as.matrix(t(BCL$counts[,names(bcl_slice[bcl_slice == 1])]))
bcl_paths[[2]] <- as.matrix(t(BCL$counts[,names(bcl_slice[bcl_slice == 2])]))
bcl_paths[[3]] <- as.matrix(t(BCL$counts[,names(bcl_slice[bcl_slice == 3])]))
bcl_paths[[4]] <- as.matrix(t(BCL$counts[,names(bcl_slice[bcl_slice == 4])]))
## lets also make sure we are using genes that are present across all slices
## because the final corpus will be all the spots across all the slices
## and each one will also need to contain the genes in the union
all_genes <- lapply(bcl_paths, function(p){
  print(dim(p))
  genes <- colnames(p)
  print(length(genes))
  genes
})
all_genes <- Reduce(intersect, all_genes)
length(all_genes)
## find shared set of OD genes across sections to try and identify common transcriptional patterns
bcl_ODgenes <- lapply(bcl_paths, function(p) {
  dat <- preprocess(p,
                   alignFile = NA,
                   extractPos = FALSE,
                   selected.genes = all_genes, # use the set of shared genes across all slices
                   nTopGenes = 5,
                   genes.to.remove = NA,
                   removeAbove = 0.95,
                   removeBelow = NA,
                   min.reads = 10,
                   min.lib.size = 10,
                   min.detected = 1,
                   ODgenes = TRUE, # select the overdispersed genes
                   od.genes.alpha = 0.05, # use a smaller alpha to use the union of the top OD genes in each section
                   gam.k = 5,
                   nTopOD = 100) # Select the top OD genes in each slice
  colnames(dat$corpus)
})
## because very few genes are OD across all sections, instead take the union and consider genes
## that are OD in any slice
unionBCLGenes <- Reduce(union, bcl_ODgenes)
length(unionBCLGenes)
```

```{r}
bclCorpus <- preprocess(t(as.matrix(BCL$counts)),
                       alignFile = NA,
                       extractPos = FALSE,
                       selected.genes = unionBCLGenes,
                       min.reads = 0, 
                       min.lib.size = 0, # because using specific set of genes, remove spots with 0 counts but keep others
                       ODgenes = FALSE)
bcl_pos <- BCL$pos[,1:2]
bcl_slice <- BCL$pos[,3]
names(bcl_slice) <- rownames(bcl_pos)
length(bcl_slice)
bclCorpus$pos <- bcl_pos[rownames(bclCorpus$corpus), ]
# filter spots in bcl_slice based on those kept in the corpus
bcl_slice_filt <- bcl_slice[names(bcl_slice) %in% rownames(bclCorpus$pos)]
length(bcl_slice_filt)
bclCorpus$slice <- bcl_slice_filt
dim(bclCorpus$corpus)
bclCorpus$slm
dim(bclCorpus$pos)
length(bclCorpus$slice)
```