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

# Misc and paths
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

# Load data
res_dir <- paste0(proj_dir, "/data/stamp_1/processed")
sce <- qread(glue("{res_dir}/proc_sce.qs"), nthreads = 8)

sce

# sce <- sce[,1:1000]

# Build SNN graph
snn.gr <- buildSNNGraph(sce, type = "jaccard", BNPARAM=AnnoyParam(), use.dimred="PCA", BPPARAM = bp)

# Run Leiden
clusters <- igraph::cluster_louvain(snn.gr, resolution = 0.8)

# Assign labels
sce$label <- as.character(clusters$membership)
table(sce$label)

plotReducedDim(sce, "UMAP", colour_by ="label", point_size = 1, text_by = "label", scattermore = T)

# Save
qsave(sce, glue("{res_dir}/clust_sce.qs"), nthreads = 8)
