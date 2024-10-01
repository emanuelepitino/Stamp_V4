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
 })
# data
dir <- glue("{here()}")
# Parameters and paths
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))
stamp <- "stamp_11"
sub <- "PBMCs"
sce <- qread(glue("{proj_dir}/data/{stamp}/{sub}/processed/clust_sce.qs"), nthreads = 8)

# palette
pal <- palette_general()
names(pal) <- unique(sce$label)

# HSPC signature
ft <- c("KIT", "TRDC", "IL1R1", "SOX4", "TNFRSF18", "TYROBP", "TNFRSF4", "FCER1G", "IL2RA", "GATA3")

#sce <- sce[,sample(colnames(sce),10000)]
# run  aucell
out <- AUCell_run(sce, ft, BPPARAM = bp)
out <- out@assays@data@listData$AUC
# add score
sce$hspc <- out["geneSet", match(colnames(sce), colnames(out))]

# Plot df
df <- as.data.frame(colData(sce))
df <- df[sample(rownames(df)),] # shuffle
um <- reducedDim(sce,"UMAP") # umap coords
df <- merge(df, um, by = "row.names")
rownames(df) <- df$Row.names
df$Row.names <- NULL
# order clusters
df$label <- factor(df$label, levels = sort(unique(as.numeric(df$label)),decreasing = F))

# UMAP with score & clusters
um_score <- ggplot(df, aes( x = UMAP1, y = UMAP2, color = hspc)) +
  ggrastr::rasterise(geom_point(shape = 16, size = 0.5), dpi = 400) +
  scale_color_viridis_c() +
  theme_bw() +
  theme(panel.grid = element_blank())
#um_score

# Calculate the centroids for each label
centroids <- df %>%
  group_by(label) %>%
  summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))
# Create the plot
um_clust <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = label)) +
  ggrastr::rasterise(geom_point(shape = 16, size = 0.5), dpi = 400) +
  scale_color_manual(values = pal) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  guides(color = "none") +
  geom_label(data = centroids, aes(label = label), 
             color = "black", fill = "white", size = 4, fontface = "bold")

#um_clust
um <- wrap_plots(um_score,um_clust, ncol= 2)
#um

# Violin plot
violin <- ggplot(df, aes(x = label, y = hspc)) +
  ggrastr::rasterise(geom_beeswarm(aes(color = hspc), 
                corral = "random",  # Use 'wrap' to keep points closer to the center
                corral.width = 0.8,  # Narrow corral width to reduce spreading
                cex = 0.1), dpi = 400) +  # Adjust point size to reduce overlap
  geom_violin(alpha = 0.5, size = 0.2, width = 1) +  # Adjust the width of the violins
  scale_color_viridis_c() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
   theme_bw() +
   theme(panel.grid = element_blank()) + 
   labs(x = "Cluster", y = "HSPC score")

#violin
violin <- plotColData(sce, x = "label", y = "hspc", point_size = 0, color = "hspc")
# save plots
pltdir <- glue("{plt_dir}/{stamp}/{sub}/")
pdf(glue("{pltdir}/HSPC.pdf"), width = 18)
violin
um
dev.off()
# save obj
qsave(sce, file = glue("{proj_dir}/data/{stamp}/{sub}/processed/scored_sce.qs"), nthreads = 8)



