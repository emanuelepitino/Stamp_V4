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
sce <- qread(glue("{proj_dir}/data/stamp_4/processed/qc_sce.qs"), nthreads = 8)
# Square with 10 CTCs
sce <- sce[,sce$fov < 285]
# read markers
mkr <- qread(glue("{proj_dir}/data/pbmc_mcf7/mcf7_markers.qs"))
ft <- head(mkr[mkr$cluster == "MCF7", "gene"], 100)
pbmc <- head(mkr[mkr$cluster == "PBMC", "gene"], 5)

# subset
#rndm <- sample(colnames(sce), 50000)
#rndm <- c(rndm,ctcs)
#sce <- sce[,colnames(sce) %in% rndm]

# run  aucell
out <- AUCell_run(sce, ft, BPPARAM = bp)
out <- out@assays@data@listData$AUC

sce$score <- out["geneSet", match(colnames(sce), colnames(out))]

df <- as.data.frame(colData(sce))
df <- df %>% select(sum,score, Area.um2, Mean.PanCK)
df <- df[sample(rownames(df)),] # shuffle

# Filter the data to label points above the horizontal line and to the right of the vertical line
thr_pck <- 2500
thr_score <- 0.2
df_lab1 <- df[df$Mean.PanCK > thr_pck & df$score > thr_score, ]
df_lab1$lab <- as.character(paste0("2.",1:nrow(df_lab1)))

# Visualize
gg_ctc2 <- ggplot(df, aes(x = score, y = Mean.PanCK, color = sum, size = Area.um2)) + 
  ggrastr::rasterise(geom_point(shape = 16), dpi = 500) +
  scale_color_viridis_c(option = "H") +
  scale_size_continuous(breaks = c(250, 1000),range = c(0.1, 5)) +
  scale_color_viridis_c(option = "H") + 
  geom_hline(yintercept = thr_pck, color = "red", linetype = "dashed") +
  geom_vline(xintercept = thr_score, color = "red", linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(size = "Area.um2", color = "nCount", x = "MCF7 score", y = "PCK MFI") +
  geom_label_repel(
    data = subset(df_lab1, rownames(df_lab1) == "c_1_201_2776"), 
    aes(label = lab),
    fill = NA, 
    direction = "both",
    label.size = NA,
    segment.color = "black",
    segment.size = 1,
    size = 8, 
    box.padding = 0.5,
    point.padding = 0.3, 
    force = 1,
    force_pull = 0.5,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  theme(text = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 20, color = "black"))



pdf("/Users/emanuelepitino/Desktop/ctc_fig/CTC_2.40.pdf", width = 6, height = 4)
gg_ctc2
dev.off()


sce$cline <- "PBMCs"
sce$cline[colnames(sce) %in% rownames(df_lab1)] <- "CTCs"
# add lab info for metadata in napari
sce$lab <- df_lab1$lab[match(colnames(sce),rownames(df_lab1))]
sce$lab[is.na(sce$lab)] <- "none"

dir <- glue("{proj_dir}/data/stamp_4/processed/CTC2")
dir.create(dir, showWarnings = F)
qsave(sce, file = glue("{dir}/ctc2_sce.qs"), nthreads = 8)
cd <- as.data.frame(colData(sce))
qsave(cd, file = glue("{dir}/ctc2_cd.qs"), nthreads = 8)