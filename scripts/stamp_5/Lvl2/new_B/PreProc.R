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
# bin 
sub <- "B"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))
# load data
res_dir <- paste0(proj_dir, "/data/stamp_5/processed")
sce <- qread(glue("{res_dir}/lvl1_sce.qs"), nthreads = 8)
sce <- sce[,sce$lvl1 == sub]
#sce <- sce[,sample(colnames(sce),10000)]
sce
# Log normalize
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
sce <- logNormCounts(sce, BPPARAM = bp) 

# Find Variable Features
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
dec.var <- modelGeneVar(sce, BPPARAM = bp)
hvg <- getTopHVGs(dec.var,fdr.threshold = 1) # select hvg on fdr
dec.var$hvg <- "no" # Assign to dec.var column for plot
dec.var$hvg[rownames(dec.var) %in% hvg] <- "yes"
gg_hvg <- plot_hvg(dec.var = dec.var, sub = sample) # plot
gg_hvg

## PCA
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
set.seed(101001)
sce <- fixedPCA(sce,BSPARAM=IrlbaParam(), subset.row = hvg) # approximate SVD with irlba
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

## Combined plot
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
combined <- wrap_plots(gg_var, gg_hvg, gg_pca, gg_um, ncol = 2, nrow = 2) +
  plot_annotation(tag_levels = "A") + 
  plot_annotation(title = glue("Stamp 5 - {sub}"), subtitle = glue("N = {comma(ncol(sce))} cells"))
combined
# save plot
pdf(glue("{plt_dir}/stamp_5/Lvl2/{sub}/PreProc.pdf"), width = 12, height = 8)
combined
dev.off()

# Save obj
qsave(sce, glue("{proj_dir}/data/stamp_5/processed/Lvl2/{sub}/proc_sce.qs"), nthreads = 8)