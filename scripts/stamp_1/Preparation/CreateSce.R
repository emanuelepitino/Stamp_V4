# Create sce 
#Libraries
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(dplyr)
  library(patchwork)
  library(grid)
  library(ggpubr)
  library(here)
  library(scater)
  library(scuttle)
  library(glue)
  library(scran)
  library(patchwork)
  library(qs)
  library(data.table)
})

dir <- glue("{here()}/scripts")
# Parameters and paths
source(glue("{dir}/misc/paths.R"))
source(glue("{dir}/misc/BIN.R"))

data_dir <- glue("{proj_dir}/data/stamp_1/raw/raw_proc")
counts <- qread(glue("{data_dir}/counts_unfiltered.qs"), nthreads = 5)
negcounts <-  qread(glue("{data_dir}/negcounts_unfiltered.qs"), nthreads = 5)
falsecounts <-  qread(glue("{data_dir}/falsecounts_unfiltered.qs"), nthreads = 5)
md <-  qread(glue("{data_dir}/metadata_unfiltered.qs"), nthreads = 5)

sce <- SingleCellExperiment(
  assays = list(counts = t(counts)),
  colData = md
)

metadata(sce)$negcounts <- negcounts
metadata(sce)$falsecounts <- falsecounts

qsave(sce, file = glue("{proj_dir}/data/stamp_1/raw/raw_proc/raw_sce.qs"), nthreads = 8)






