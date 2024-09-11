# Dependencies
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(dplyr)
  library(here)
  library(scater)
  library(scuttle)
  library(glue)
  library(qs)
  library(parallel)
  library(scran)
  library(BiocParallel)
  library(BiocNeighbors)
  library(BiocSingular)
  library(data.table)
  library(scales)
})

# read data
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))
subsets <- c("Myeloid", "T", "B", "NK")
sce_list <- lapply(subsets, function(subset) {
  qread(glue("{proj_dir}/data/stamp_5/processed/Lvl2/{subset}/lvl2_sce.qs"), nthreads = 8)
})
names(sce_list) <- subsets
# Remove PCA from obj
sce_list <- lapply(sce_list, function(sce) {
  reducedDim(sce, "PCA") <- NULL
  sce
})
sce <- do.call(cbind,sce_list)

reducedDim(sce, "PCA") <- NULL
reducedDim(sce,"UMAP") <- NULL
assay(sce,"logcounts") <- NULL
#Log normalize
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
sce <- logNormCounts(sce, BPPARAM = bp) 

## PCA
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
set.seed(101001)
sce <- fixedPCA(sce,BSPARAM=IrlbaParam(), subset.row = NULL) # approximate SVD with irlba
num_pcs_to_retain <- 8
percent.var <- attr(reducedDim(sce), "percentVar")
# Plot Elbow
data <- data.frame(PC = 1:length(percent.var), Variance = percent.var)
gg_var <- ggplot(data, aes(x = PC, y = Variance)) +
  geom_point() +
  xlab("PC") +
  ylab("Variance explained (%)") +
  geom_vline(xintercept = num_pcs_to_retain, color = "red") +
  theme_bw()
gg_var
# Plot PCA
reducedDim(sce, "PCA") <-  reducedDim(sce, "PCA")[,1:num_pcs_to_retain]
gg_pca <- plotPCA(sce, scattermore = TRUE) + ggtitle("PCA")
gg_pca

## UMAP
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
set.seed(123)
sce <- runUMAP(sce, dimred = "PCA", n_dimred = c(1:num_pcs_to_retain), n_threads = 8)
# Plot UMAP
gg_um <- plotReducedDim(sce, "UMAP", scattermore = TRUE) 
gg_um

# Save obj
qsave(sce, glue("{proj_dir}/data/stamp_5/processed/combined_sce.qs"), nthreads = 8)
