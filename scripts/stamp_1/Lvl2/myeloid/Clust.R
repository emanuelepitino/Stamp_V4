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
sub <- "Myeloid"
res_dir <- paste0(proj_dir, "/data/stamp_1/processed/Lvl2/",sub)
sce <- qread(glue("{res_dir}/proc_sce.qs"))
sce

# Annoy Algorithm
# Build SNN graph
snn.gr <- buildSNNGraph(sce, type = "jaccard", BNPARAM=AnnoyParam(ntrees = 200), use.dimred="PCA", BPPARAM = bp)
# Run Louvain
annoy <- igraph::cluster_louvain(snn.gr, resolution = 0.3)

# Assign labels
sce$label <- as.character(annoy$membership)

#table(sce$leiden)
table(sce$label)

plotReducedDim(sce, "UMAP", colour_by = "label", text_by = "label", point_size = 1, raster = F, scattermore = T)
# Save
qsave(sce, glue("{res_dir}/clust_lvl2_sce.qs"))
