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
})

# Load bin functions
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

sce <- qread(glue("{proj_dir}/data/stamp_5/processed/final_obj.qs"))


# Log normalize
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
sce <- logNormCounts(sce, BPPARAM = bp) 

## PCA
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
set.seed(101001)
sce <- fixedPCA(sce,BSPARAM=IrlbaParam(), subset.row = NULL) # approximate SVD with irlba
num_pcs_to_retain <- 10
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
pal_lin <- Polychrome::createPalette(26,c("#99FFFF", "#FF99FF", "#FFFF99"))
names(pal_lin) <- unique(sce$Cell_Lineage)

plotReducedDim(sce, "UMAP", scattermore = TRUE, color_by = "Cell_Lineage") +
  scale_color_manual(values = pal_)


pal_full <- Polychrome::createPalette(26,c("#99FFFF", "#FF99FF", "#FFFF99"))
names(pal_full) <- unique(sce$Cell_Type)

plotReducedDim(sce, "UMAP", scattermore = TRUE, color_by = "Cell_Type") +
  scale_color_manual(values = pal_full)


qsave(sce, file = glue("{proj_dir}/data/stamp_5/processed/final_stamp5.qs"), nthreads = 8)








