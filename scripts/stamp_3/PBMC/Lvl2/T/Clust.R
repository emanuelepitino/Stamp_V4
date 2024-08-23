# Libraries
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
})

sub <- "T"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

res_dir <- paste0(proj_dir, "/data/stamp_3/processed/Lvl2")
sce <- qread(glue("{res_dir}/PreProcNew_{sub}.qs"))
sce

# sce <- sce[,1:1000]

# Build SNN graph
snn.gr <- buildSNNGraph(sce, type = "jaccard", BNPARAM=AnnoyParam(), use.dimred="PCA", BPPARAM = bp)

clusters <- igraph::cluster_louvain(snn.gr, resolution = 1)

ratio <- pairwiseModularity(snn.gr, clusters$membership, as.ratio=TRUE)

library(pheatmap)
pheatmap(log2(ratio+1), cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("gold", "navy"))(100))

# Assign labels
sce$label <- as.character(clusters$membership)
table(sce$label)

plotReducedDim(sce, "UMAP", colour_by ="label", point_size = 1, text_by = "label", scattermore = T)

# Save
qsave(sce, glue("{res_dir}/clust_sce_{sub}.qs"), nthreads = 8)
