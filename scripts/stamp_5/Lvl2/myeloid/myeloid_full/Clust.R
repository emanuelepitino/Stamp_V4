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
lin <- "Myeloid"
sub <- "Myeloid"
res_dir <- paste0(proj_dir, "/data/stamp_5/processed/Lvl2/",lin,"/",sub)
sce <- qread(glue("{res_dir}/proc_sce2.qs"))

sce

#sce <- sce[,1:100000]

# Annoy Algorithm
# Build SNN graph
snn.gr <- buildSNNGraph(sce, type = "jaccard", BNPARAM=AnnoyParam(ntrees = 200), use.dimred="PCA", BPPARAM = bp)
# Run Louvain
#annoy <- igraph::cluster_louvain(snn.gr, resolution = 0.5)
annoy2 <- igraph::cluster_louvain(snn.gr, resolution = 1)
# Assign labels
#sce$leiden <- as.character(leiden)
#sce$label <- as.character(annoy$membership)
sce$label2 <- as.character(annoy2$membership)

#table(sce$leiden)
table(sce$label2)

#plotReducedDim(sce, "UMAP", colour_by = "label", text_by = "label", point_size = 1, raster = F, scattermore = T)
plotReducedDim(sce, "UMAP", colour_by = "label2", text_by = "label2", point_size = 1, raster = F, scattermore = T)

# Save
qsave(sce, glue("{res_dir}/clust_lvl2_sce2.qs"))


