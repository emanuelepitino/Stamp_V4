
# Unsupervised clustering - InSituType Stamp_7b ESCs

# Libraries
suppressPackageStartupMessages({
  library(SingleCellExperiment)
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
  library(data.table)
  library(InSituType)
})

# Data loading
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

stamp <- "stamp_7b"
sample <- "iESC"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

res_dir <- glue("{proj_dir}/data/{stamp}/{sample}")
sce <- qread(glue("{res_dir}/PreProcNew.qs"), nthreads = 8)
#sce <- sce[,sample(colnames(sce),1000)]

# Free memory
reducedDim(sce,"PCA") <- NULL
reducedDim(sce,"UMAP") <- NULL
assay(sce,"logcounts") <- NULL
assay(sce,"normcounts") <- NULL
gc()
sce 

#Run IST unsupervised clustering
ist <- function(sce, nk, gs=TRUE, pbs=NULL, bkg=TRUE) {
  # dependencies
  library(InSituType)
  library(SingleCellExperiment)
  # load counts
  mtx <- counts(sce[gs, ])
  mtx <- as(t(mtx), "dgCMatrix")
  # cohorting based on IF data
  j <- names(cd <- colData(sce))
  i <- grep("^Mean", j, value=TRUE)
  i <- setdiff(i, "Mean.G")
  i <- c("Area", "AspectRatio", i)
  coh <- fastCohorting(as.matrix(cd[i]))
  # background estimation
  neg <- grep("^neg", altExpNames(sce), value=TRUE)
  neg <- sce$nCount_negprobes/nrow(altExp(sce, neg))
  # update reference profiles
  pbs <- if (!is.null(pbs)) {
    bkg <- if (bkg) {
      rna <- sce$nCount_RNA
      rna*mean(neg)/mean(rna) 
    }
    updateReferenceProfiles(
      reference_profiles=pbs, counts=mtx,
      neg=neg, bg=bkg)$updated_profiles
  }
  # clustering
  insitutype(mtx, 
             reference_profiles=pbs,
             update_reference_profiles=FALSE,
             neg=neg, cohort=coh, n_clusts=nk)
  # Free memory
  rm(mtx, coh, pbs)
}

unsup <- ist(sce, 10:20)

##Save
dir <- glue("{proj_dir}/data/{stamp}/{sample}/Ist")
dir.create(dir, showWarnings = F, recursive = T)
qsave(unsup, file = glue("{dir}/unsup_10_20.qs"))

sessionInfo()