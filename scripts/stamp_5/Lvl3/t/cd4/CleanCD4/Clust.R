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
sub <- "CD4"
res_dir <- paste0(proj_dir, "/data/stamp_5/processed/Lvl3/",sub)
sce <- qread(glue("{res_dir}/proc_sce.qs"))

sce

#sce <- sce[,1:50000]

# Annoy Algorithm
# Build SNN graph
snn.gr <- buildSNNGraph(sce, BNPARAM=AnnoyParam(ntrees = 200), use.dimred="PCA", BPPARAM = bp)
# Run Leiden
annoy <- igraph::cluster_leiden(snn.gr, resolution_parameter = 0.0003)


# Assign labels
#sce$leiden <- as.character(leiden)
sce$label <- as.character(annoy$membership)
#table(sce$leiden)
table(sce$label)

plotReducedDim(sce, "UMAP", colour_by = "label", text_by = "label", point_size = 1, rasterise = F, scattermore = T)

# Save
qsave(sce, glue("{res_dir}/clust_sce.qs"))



