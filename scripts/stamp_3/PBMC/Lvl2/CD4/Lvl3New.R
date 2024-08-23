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
  library(scales)
  library(ROGUE)
})

sub <- "CD4"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

res_dir <- paste0(proj_dir, "/data/stamp_3/processed/Lvl2")
sce <- qread(glue("{res_dir}/clust_sce_{sub}.qs"))

# Find markers
markers <- scoreMarkers(sce, as.factor(sce$label), BPPARAM = bp)

# Clean markers list
transform_marker <- function(marker_df, cluster_name) {
  marker_df <- as.data.frame(marker_df) %>%
    select(median.logFC.detected, self.detected) %>%
    arrange(desc(self.detected)) %>%
    #filter(median.logFC.cohen > 0.25) %>%
    mutate(gene = rownames(.), cluster = as.numeric(cluster_name))
  rownames(marker_df) <- NULL
  return(marker_df)
}

# Apply the function to each element of the markers list along with their names
markers <- mapply(transform_marker, markers, names(markers), SIMPLIFY = FALSE)

markers <- bind_rows(markers) %>%
  arrange(cluster,desc(median.logFC.detected), .by_group = TRUE) 

markers$median.logFC.detected <- round(markers$median.logFC.detected, 2)

top <- markers %>%
  group_by(cluster) %>%
  slice_head(n = 20)
feats <- unique(top$gene)


pal <- palette_general()
names(pal) <- unique(sce$label)
gg_clust1 <- create_plots2(sce, "label", feats, pal)

df <- as.data.frame(colData(sce))
gg_clust2 <- wrap_plots(
  plot_density(df, "sum", "label", pal, "Counts",500),
  plot_density(df, "detected", "label", pal, "Features",300),
  plot_density(df, "Area.um2", "label", pal, "Cell Area",500),
  ncol = 2, nrow = 2) +
  plot_annotation(tag_levels = "A")


sce$lvl2[sce$label %in% c(2,3,10)] <- "CD4"
sce$lvl2[sce$label %in% c(5,7,8,9)] <- "CD8"
sce$lvl2[sce$label %in% c(1,6,12,13)] <- "LowQ"
sce$lvl2[sce$label %in% c(4,11)] <- "NK"


# Find markers
markers <- scoreMarkers(sce, sce$lvl2, BPPARAM = bp)

# Clean markers list
transform_marker <- function(marker_df, cluster_name) {
  marker_df <- as.data.frame(marker_df) %>%
    select(median.logFC.detected, self.detected) %>%
    arrange(desc(self.detected)) %>%
    #filter(median.logFC.cohen > 0.25) %>%
    mutate(gene = rownames(.), cluster = as.character(cluster_name))
  rownames(marker_df) <- NULL
  return(marker_df)
}

# Apply the function to each element of the markers list along with their names
markers <- mapply(transform_marker, markers, names(markers), SIMPLIFY = FALSE)

markers <- bind_rows(markers) %>%
  arrange(cluster,desc(median.logFC.detected), .by_group = TRUE) 

markers$median.logFC.detected <- round(markers$median.logFC.detected, 2)

top <- markers %>%
  group_by(cluster) %>%
  slice_head(n = 10)
feats <- unique(top$gene)

pal <- palette_general()
names(pal) <- unique(sce$lvl2)

gg_anno1 <- create_plots2(sce, "lvl2", feats, pal)

df <- as.data.frame(colData(sce))
gg_anno2 <- wrap_plots(
  plot_density(df, "sum", "lvl2", pal, "Counts",500),
  plot_density(df, "detected", "lvl2", pal, "Features",300),
  plot_density(df, "Area.um2", "lvl2", pal, "Cell Area",500),
  ncol = 2, nrow = 2) +
  plot_annotation(tag_levels = "A")


wrap_plots(gg_anno1, gg_anno2, ncol = 1)

outdir <- glue("{plt_dir}/stamp_3/Lvl2")
if(!dir.exists(paste0(outdir))){
  dir.create(outdir)
}

pdf(paste0(outdir,"/lvl2_{sub}.pdf"), width = 15, height = 10)
gg_clust1
gg_clust2
gg_anno1
gg_anno2
dev.off()

qsave(sce, file = glue("{res_dir}/lvl2_sce_{sub}.qs"), nthreads = 8)
