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

sub <- "Whole PBMCs"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

res_dir <- paste0(proj_dir, "/data/stamp_3/processed")
sce <- qread(glue("{res_dir}/clust_sce.qs"))

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
  slice_head(n = 10)
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


sce$lvl1[sce$label %in% c(1,2,4)] <- "T"
sce$lvl1[sce$label == 3] <- "B"
sce$lvl1[sce$label == 5] <- "Myeloid"
sce$lvl1[sce$label == 6] <- "LowQ"


# Find markers
markers <- scoreMarkers(sce, sce$lvl1, BPPARAM = bp)

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
names(pal) <- unique(sce$lvl1)

gg_anno1 <- create_plots2(sce, "lvl1", feats, pal)

df <- as.data.frame(colData(sce))
gg_anno2 <- wrap_plots(
  plot_density(df, "sum", "lvl1", pal, "Counts",500),
  plot_density(df, "detected", "lvl1", pal, "Features",300),
  plot_density(df, "Area.um2", "lvl1", pal, "Cell Area",500),
  ncol = 2, nrow = 2) +
  plot_annotation(tag_levels = "A")


wrap_plots(gg_anno1, gg_anno2, ncol = 1)

outdir <- glue("{plt_dir}/stamp_3")
if(!dir.exists(paste0(outdir))){
  dir.create(outdir)
}

pdf(paste0(outdir,"/lvl1.pdf"), width = 15, height = 10)
gg_clust1
gg_clust2
gg_anno1
gg_anno2
dev.off()

qsave(sce, file = glue("{res_dir}/lvl1_sce.qs"), nthreads = 8)
