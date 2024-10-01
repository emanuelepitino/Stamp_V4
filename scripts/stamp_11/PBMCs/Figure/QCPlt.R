# dep
library(ggplot2)
library(qs)
library(here)
library(glue)


# Parameters and paths
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))
stamp <- "stamp_11"
sub <- "PBMCs"

# load data
sce_5k <- qread(glue("{proj_dir}/data/{stamp}/{sub}/processed/lvl1_sce.qs"), nthreads = 8)
sce_io <- qread(glue("{proj_dir}/data/stamp_5/processed/lvl1_sce.qs"))

# panels intersection
feat <- intersect(rownames(sce_io),rownames(sce_5k))
# subset
sce_5k <- sce_5k[feat,]
sce_io <- sce_io[feat,]

sce_5k <- addPerCellQCMetrics(sce_5k)
sce_io <- addPerCellQCMetrics(sce_5k)


cd_5k <- as.data.frame(colData(sce_5k)) %>% select(lvl1, detected, sum, cell_area, nucleus_area)
cd_io <- as.data.frame(colData(sce_io))

