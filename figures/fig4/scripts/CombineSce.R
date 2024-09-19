# dependencies
library(SingleCellExperiment)
library(qs)
library(glue)
library(here)
library(scales)

# paths
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))
stamps <- c("stamp_9","stamp_11","stamp_12")
sub <- c("MCF7","SKBR3")

# read sce objects
combinations <- expand.grid(stamps = stamps, sub = sub, stringsAsFactors = FALSE)
stamp_sub <- setNames(
  mapply(function(stamp, sub) {
    sce <- qread(glue("{proj_dir}/data/{stamp}/{sub}/processed/sce_filt.qs"))
    sce$stamp <- stamp
    sce$sub <- sub
    return(sce)
  }, combinations$stamps, combinations$sub, SIMPLIFY = FALSE),
  glue("{combinations$stamps}_{combinations$sub}")
)
# Extract coldata
cd <- lapply(stamp_sub, function(sce){
  as.data.frame(colData(sce)) %>% select(cell_area,nucleus_area,sum,detected, control_probe_counts, stamp, sub)
})
cd <- do.call(rbind,cd)

# save
dir <- glue("{proj_dir}/figures/fig4/obj")
dir.create(dir, showWarnings = F)
qsave(cd, file = glue("{dir}/cd.qs"))

sce <- do.call(cbind,stamp_sub)
qsave(sce, file = glue("{dir}/sce.qs"))