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
})

# Load bin functions
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

## Load data
myeloid <- qread(glue("{proj_dir}/data/stamp_5/processed/Myeloid/lvl2_sce.qs"), nthreads = 8)
B <- qread(glue("{proj_dir}/data/stamp_5/processed/B/lvl2_sce.qs"), nthreads = 8)
T <- qread(glue("{proj_dir}/data/stamp_5/processed/T/lvl2_sce.qs"), nthreads = 8)
NK <- qread(glue("{proj_dir}/data/stamp_5/processed/NK/proc_sce.qs"), nthreads = 8)
NK$lvl2 <- "NK"

sce <- cbind(myeloid,B,T,NK)
sce <- sce[,sce$lvl2 != "LowQ"]

sce$Cell_Lineage <- sce$lvl1
sce$Cell_Type <- sce$lvl2
sce$Cell_Lineage <- factor(sce$Cell_Lineage, levels = c("T","B","Myeloid","NK"))


sce$Cell_Lineage[sce$Cell_Lineage %in% c("Central Memory CD4","Effector Memory CD4","Naive CD4")] <- "Naive/Memory CD4"


# Lineages dot plot
feats_lin <- c("CD3E","CD3D","CD4","CD8A","CD8B",
               "CD19","CD79A","CD79B",
               "CD14","FCGR3A","CD68","CD86","NGK7","NCAM1","KLRC1","KLRD1","KLRF1")

dot_lin <- plotDots(sce, features = intersect(feats_lin, rownames(sce)), group = "Cell_Lineage", scale =TRUE, center = TRUE) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 20, color = "black")) &
  labs(x = "Cell Lineage")
dot_lin


# Cell Type dot plot

feats_type <- c("CD3E","CD3D","CD8B","CD8A","KLRK1","GZMK", "EOMES","LAG3","ZNF683","GZMH","KLRD1","FASLG","PDCD1","KLRB1","CD38","CCR6",
          "CCR4","CCL20","RORC","FOXP3","TNFRSF9","CCR7","IL7R","TCF7",
          "CD40LG","SELL","LEF1",
          "CD79A","CD79B","CD19","IGHM","TCL1A", "IGHG1","IGHG2","IGHG3","IGHG4","CD27","CD24","S100A9","CD14","FCGR3A","IL1B","SELL",
          "CD68","ITGAM","CD34","FLT3","CX3CR1","CLEC9A","CCR5",
            "S100A8","LYZ","CCR2","ITGAX","CD1C")
          
#sce$Cell_Type <- factor(sce$Cell_Type, levels = unique(sce$Cell_Type))

sce$Cell_Type <- factor(sce$Cell_Type, levels = c("Naive CD8","Cytotoxic CD8","Effector Memory CD8","NK", "T helper","Tregs", "Naive CD4", "Effector Memory CD4",
                            "Central Memory CD4",
                            "Memory B","Naive B","Unswitched B","Memory B ITGAX+",
                            "Classical Mono","Intermediate Mono","Non Classical Mono","type II Dendritic Cells"))

plotDots(sce, features = intersect(feats_type, rownames(sce)), group = "Cell_Type", scale =TRUE, center = TRUE) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 20, color = "black")) &
  labs(x = "Cell Lineage")




