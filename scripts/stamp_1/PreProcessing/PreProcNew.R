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

stamp <- "stamp_1"
sample <- "PBMCs"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

res_dir <- glue("{proj_dir}/data/{stamp}")
sce <- qread(glue("{res_dir}/qc_sce.qs"))
sce


sce <- sce[,sample(colnames(sce),250000)]
sce

counts <- t(counts(sce))

# normalize by total counts 
totalcounts <- Matrix::rowSums(counts)  
norm <- counts / totalcounts

# log1p transformation
lognorm <- log1p(norm)

assay(sce, "logcounts") <- t(lognorm)
assay(sce, "normcounts") <- t(norm)

# PC analysis
pc1 <- irlba::prcomp_irlba(sqrt(norm), n = 25)
# 1. Calculate the variance explained by each PC 
explained_variance <- pc1$sdev^2
# 2. Calculate the proportion of variance explained
proportion_variance_explained <- explained_variance / sum(explained_variance)
# 3. Create a data frame for plot
variance_df <- data.frame(
  PC = seq_along(proportion_variance_explained),
  ProportionVarianceExplained = proportion_variance_explained
)

# Run UMAP
um1 <- uwot::umap(pc1$x, n_neighbors = 40, spread = 1, min_dist = 0.1, metric = "cosine")
rownames(um1) <- rownames(norm)

par(mar = c(0,0,0,0))
plot(um1, pch = 16, cex = 0.1, col = "dodgerblue4")

reducedDim(sce,"PCA") <- pc1$x
reducedDim(sce,"UMAP") <- um1

# Save plots
pdf(glue("{plt_dir}/{stamp}/PreProc.pdf"), width = 7, height = 5)
plot(um1, pch = 16, cex = 0.1, col = "dodgerblue4")
dev.off()

# Save data
qsave(sce, file = glue("{res_dir}/PreProcNew.qs"))