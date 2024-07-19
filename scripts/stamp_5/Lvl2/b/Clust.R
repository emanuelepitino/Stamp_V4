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
sub <- "B"
res_dir <- paste0(proj_dir, "/data/stamp_5/processed/Lvl2/",sub)
sce <- qread(glue("{res_dir}/proc_sce.qs"))

sce

#sce <- sce[,1:100000]

# Annoy Algorithm
# Build SNN graph
# Jaccard index
snn.gr <- buildSNNGraph(sce, type = "jaccard", BNPARAM=AnnoyParam(ntrees = 200), use.dimred="PCA", BPPARAM = bp)
# Run Leiden
# Louvain
louvain1 <- igraph::cluster_louvain(snn.gr, resolution = 1)
louvain05 <- igraph::cluster_louvain(snn.gr, resolution = 0.5)

#annoy <- igraph::cluster_leiden(snn.gr, resolution_parameter = 0.0001)
#annoy2 <- igraph::cluster_leiden(snn.gr, resolution_parameter = 0.0003)
#annoy3 <- igraph::cluster_leiden(snn.gr, resolution_parameter = 0.0005)
#annoy4 <- igraph::cluster_leiden(snn.gr, resolution_parameter = 0.0008)


# Assign labels
#sce$leiden <- as.character(leiden)
sce$label <- as.character(louvain1$membership)
sce$label2 <- as.character(louvain05$membership)
#sce$label3 <- as.character(annoy3$membership)
#sce$label4 <- as.character(annoy4$membership)

#table(sce$leiden)
table(sce$label)
table(sce$label2)

plotReducedDim(sce, "UMAP", colour_by = "label", text_by = "label", point_size = 1, raster = F, scattermore = T)
plotReducedDim(sce, "UMAP", colour_by = "label2", text_by = "label2", point_size = 1, raster = F, scattermore = T)
#plotReducedDim(sce, "UMAP", colour_by = "label3", text_by = "label3", point_size = 1, raster = F, scattermore = T)
#plotReducedDim(sce, "UMAP", colour_by = "label4", text_by = "label4", point_size = 1, raster = F, scattermore = T)

# Save
qsave(sce, glue("{res_dir}/clust_lvl2_sce.qs"))



