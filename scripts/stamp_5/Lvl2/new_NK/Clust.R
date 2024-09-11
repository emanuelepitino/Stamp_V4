# dependencies
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
sub <- "NK"
scedir <- glue("{proj_dir}/data/stamp_5/processed/Lvl2/{sub}")
sce <- qread(glue("{scedir}/proc_sce.qs"), nthreads = 8)
#sce <- sce[,sample(colnames(sce),10000)]
sce
# Build SNN graph 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
snn.gr <- buildSNNGraph(sce, type = "jaccard",
                        BNPARAM=AnnoyParam(),  # annoy algorithm
                        use.dimred="PCA", BPPARAM = bp)
# Run Leiden
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
clusters <- igraph::cluster_louvain(snn.gr, resolution = 0.5)

# Viz
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sce$label <- as.character(clusters$membership) # assign labels
table(sce$label)
plotReducedDim(sce, "UMAP", colour_by ="label",text_by = "label", scattermore = T)

# Save
qsave(sce, glue("{scedir}/clust_sce.qs"), nthreads = 8)
