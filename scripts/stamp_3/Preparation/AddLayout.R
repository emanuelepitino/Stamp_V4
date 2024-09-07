suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(HDF5Array)
  library(SparseArray)
  library(glue)
  library(here)
  library(qs)
  library(SingleCellExperiment)
})

dir <- glue("{here()}/scripts")
# Parameters and paths
source(glue("{dir}/misc/paths.R"))
source(glue("{dir}/misc/BIN.R"))

# # loading
dir <- glue("{proj_dir}/data/stamp_3/raw/raw_proc")

sce <- qread( glue("{proj_dir}/data/stamp_3/raw/raw_proc/raw_sce.qs"), nthreads = 8)

md <- as.data.frame(colData(sce))

gg_cells <- ggplot(md, aes(x = CenterX_global_px, y = CenterY_global_px)) + 
   ggrastr::rasterise(geom_point(shape = 16, size = 0.05), dpi = 800) +
  coord_equal()

# Add samples information 
samples <- list(
  PBMCs_250k = c(514:769),
  c100 = c(271:360),
  c250 = c(181:260),
  c500 = c(91:180),
  c1000 = c(1:90),
  c20k = c(361:369,378:386,395:403,412:420,429:437,446:454,463:471,480:488,497:505),
  n20k = c(371:377,388:394,405:411,422:428,439:445,456:462,473:479,490:496,507:513)
)

invisible(lapply(names(samples), function(name) {
  sce$sample[sce$fov %in% samples[[name]]] <<- name
}))

sce <- sce[,! is.na(sce$sample)]

qsave(sce, file = glue("{proj_dir}/data/stamp_3/raw/raw_proc/layout_sce.qs"), nthreads = 8)
