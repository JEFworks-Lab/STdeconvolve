The following code demonstrates how to create the MOB ST input data matrix for STdeconvolve

```{r}
library(STdeconvolve)
```

```{r}
data("mOB") # data for MOB replicate #8 used primarily in the analyses
```

```{r}
## Remove poor pixels and genes and create a gene counts matrix 
## for the remaining 7365 genes and 260 pixels
mobClean <- cleanCounts(mOB$counts,
                        min.lib.size = 100,
                        max.lib.size = Inf,
                        min.reads = 100,
                        min.detected = 1,
                        verbose = TRUE,
                        plot=TRUE)

## Wrapper around `cleanCounts`, `restrictCorpus` and additional
## filtering options for feature selection of overdispersed genes
## (see Tutorial `additional_features` for more information)
mobCorpus <- preprocess(t(mOB$counts),
                       alignFile = NA,
                       extractPos = FALSE,
                       selected.genes = NA,
                       nTopGenes = NA,
                       genes.to.remove = NA,
                       removeAbove = NA,
                       removeBelow = NA,
                       min.reads = 100,
                       min.lib.size = 100,
                       min.detected = 1,
                       ODgenes = TRUE,
                       nTopOD = NA,
                       od.genes.alpha = 0.05,
                       gam.k = 5,
                       verbose = TRUE)
## add the positional information back into the corpus for the filtered pixels
mobCorpus$pos <- mOB$pos[rownames(mobCorpus$corpus), ]
```


