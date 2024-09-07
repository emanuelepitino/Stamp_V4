suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(HDF5Array)
  library(SparseArray)
  library(glue)
  library(here)
  library(SingleCellExperiment)
})

dir <- glue("{here()}/scripts")
# Parameters and paths
source(glue("{dir}/misc/paths.R"))
source(glue("{dir}/misc/BIN.R"))

# # loading
dir <- glue("{proj_dir}/data/stamp_2/raw/raw_proc")

sce <- qread( glue("{proj_dir}/data/stamp_2/raw/raw_proc/raw_sce.qs"), nthreads = 8)

md <- as.data.frame(colData(sce))

gg_cells <- ggplot(md, aes(x = CenterX_global_px, y = CenterY_global_px)) + 
   ggrastr::rasterise(geom_point(shape = 16, size = 0.05), dpi = 800) +
  coord_equal()

# Add samples information 
samples <- list(
  MIX = c(1:99),
  MCF7 = c(100:198),
  LnCAP = c(199:297),
  SKBR3 = c(298:396)
)

invisible(lapply(names(samples), function(name) {
  sce$sample[sce$fov %in% samples[[name]]] <<- name
}))

qsave(sce, file = glue("{proj_dir}/data/stamp_2/raw/raw_proc/layout_sce.qs"), nthreads = 8)





