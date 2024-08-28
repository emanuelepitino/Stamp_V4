# Libraries
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dplyr)
  library(here)
  library(glue)
  library(qs)
  library(AUCell)
})

dir <- glue("{here()}/scripts")
# Parameters and paths
source(glue("{dir}/misc/paths.R"))

# Load data
suppressMessages({
  sample <- "PBMCs_mix"
  stamp <- "stamp_8"
  source(glue("{dir}/misc/BIN.R"))
  data_dir <- glue("{proj_dir}/data/{stamp}/{sample}")
})
sce <- qread(glue("{data_dir}/qc_sce.qs"), nthreads = 8)
sce


# General signature of MCF7 & SKBR3 from Luciano
sign_general <- list('SKBR3' = c( 'ERBB2', 'DHRS2' ), 'MCF7' = c('EEF1A2', 'EPCAM'))

sce <- sce[,sample(colnames(sce),10000)]

cells_AUC <- AUCell_run(sce, geneSets = sign_general, BPPARAM = bp)
