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
sce <- HDF5Array::loadHDF5SummarizedExperiment(glue("{res_dir}/proc_sce"))

sce

# sce <- sce[,1:1000]

# Build SNN graph
snn.gr <- buildSNNGraph(sce, BNPARAM=AnnoyParam(), use.dimred="PCA", BPPARAM = bp)

# Run Louvain
clusters <- igraph::cluster_leiden(snn.gr, resolution_parameter = 0.0001)

# Assign labels
sce$label <- as.character(clusters$membership)
table(sce$label)

# Save
HDF5Array::saveHDF5SummarizedExperiment(sce, glue("{res_dir}/clust_sce"), replace = T)
