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
  library(data.table)
})

dir <- glue("{here()}")
# Parameters and paths
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

# Load data
data_dir <- glue("{proj_dir}/data/stamp_5/raw/raw/PBMCs")
barcodes <- fread(glue("{data_dir}/cell_feature_matrix/barcodes.tsv.gz"), sep = "\t", nThread = 3, header = F)
features <- fread(glue("{data_dir}/cell_feature_matrix/features.tsv.gz"), sep = "\t", nThread = 3, header = F)
counts <- Matrix::readMM(glue("{data_dir}/cell_feature_matrix/matrix.mtx.gz"))
md <- data.table::fread(glue("{data_dir}/cells.csv"))

counts <- as(counts, "dgCMatrix")

colnames(counts) <- barcodes$V1
rownames(counts) <- features$V2

gs <- rownames(counts)
ncp <- grep("NegControlProbe", gs)
ncc <- grep("NegControlCodeword", gs)
uc <- grep("UnassignedCodeword", gs)

as <- list(counts=counts[-c(ncp, ncc,uc), ])
ae <- list(
  NegControlProbe=SingleCellExperiment(list(counts=counts[ncp, ])),
  NegControlCodeword=SingleCellExperiment(list(counts=counts[ncc,])),
  UnassignedCodeword = SingleCellExperiment(list(counts=counts[uc,])))

sce <- SingleCellExperiment(as, colData=md, altExps=ae)

colnames(colData(sce))

res_dir <- paste0(proj_dir, "/data/stamp_5/processed")
if (!dir.exists(res_dir)) {
  dir.create(res_dir)
}
qsave(sce, file = glue("{res_dir}/raw_sce.qs"))