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
  library(scran)
  library(BiocParallel)
  library(BiocNeighbors)
  library(BiocSingular)
})

stamp <- "stamp_7b"
sample <- "iESC"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

res_dir <- glue("{proj_dir}/data/{stamp}/{sample}")
sce <- qread(glue("{res_dir}/PreProcNew.qs"), nthreads = 8)
#sce <- sce[,sample(colnames(sce),1000)]
sce 

# sce <- sce[,1:1000]

# Build SNN graph
snn.gr <- buildSNNGraph(sce, type = "jaccard", BNPARAM=AnnoyParam(), use.dimred="PCA", BPPARAM = bp)

# Run Louvain
clusters <- igraph::cluster_louvain(snn.gr, resolution = 0.8)

# Cluster Modularity score
ratio <- bluster::pairwiseModularity(snn.gr, clusters$membership, as.ratio=TRUE)
library(pheatmap)
pheatmap(log2(ratio+1), cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("gold", "navy"))(100))
# Assign labels
sce$label <- as.character(clusters$membership)
table(sce$label)

pal <- Polychrome::createPalette(50, c("#fc6060", "#74f774", "#7c7cfc"))
names(pal) <- unique(sce$label)

sce <- sce[, sample(ncol(sce))]
plotReducedDim(sce, "UMAP", colour_by ="label", point_size = 1, text_by = "label", scattermore = T) +
  scale_color_manual(values = pal)

# Save
qsave(sce, glue("{res_dir}/clust_sce.qs"), nthreads = 8)