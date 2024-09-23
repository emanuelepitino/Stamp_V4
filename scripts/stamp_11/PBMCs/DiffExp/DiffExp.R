# dependencies
suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(SparseArray)
  library(SingleCellExperiment)
  library(qs)
  library(glue)
  library(here)
  library(AUCell)
  library(scran)
})
# data
dir <- glue("{here()}")
# Parameters and paths
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))
stamp <- "stamp_11"
sub <- "PBMCs"
sce <- qread(glue("{proj_dir}/data/{stamp}/{sub}/processed/clust_sce.qs"), nthreads = 8)

sce

mrk <- scoreMarkers(sce, groups = sce$label, BPPARAM = bp)

feat <- lapply(mrk, function(df) {
  as.data.frame(df) %>%
    arrange(desc(median.logFC.detected)) %>%
    head(10) %>%
    rownames()
})
feat <- unique(unlist(feat))

dot <- plotDots(sce, group = "label", features = feat, scale = T, center = T) 
dot

# save
pltdir <- glue("{plt_dir}/{stamp}/{sub}/")
pdf(glue("{pltdir}/DotPlotExtended.pdf"), height = 18)
dot
dev.off()
