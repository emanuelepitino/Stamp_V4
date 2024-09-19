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
stamp <- "stamp_6"
sample <- "CTC_1"
sce <- qread(glue("{proj_dir}/data/{stamp}/{sample}/processed/sce_filt.qs"), nthreads = 8)
sce

# read markers
mkr <- qread(glue("{proj_dir}/data/pbmc_mcf7/mcf7_markers.qs"))
ft <- head(mkr[mkr$cluster == "MCF7", "gene"], 20)
pbmc <- head(mkr[mkr$cluster == "PBMC", "gene"], 20)

#ft <- mkr$gene[mkr$cluster == "MCF7"]
#pbmc <- mkr$gene[mkr$cluster == "PBMC"]
#sce <- sce[,sample(colnames(sce),500000)]

# Intersect by ENSEMBL id
ft_x <- fread(glue("{data_dir}/cell_feature_matrix/features.tsv.gz"), sep = "\t", nThread = 3, header = F)
ft_c <- fread(glue("{data_dir}/../../features.tsv.gz"), sep = "\t", nThread = 3, header = F)

ft <- ft_c$V1[ft_c$V2 %in% ft]
ft <- ft_x$V2[ft_x$V1 %in% ft]

# run  aucell
out <- AUCell_run(sce, ft, BPPARAM = bp)
out <- out@assays@data@listData$AUC

sce$ctc <- out["geneSet", match(colnames(sce), colnames(out))]

# score PBMCs
ft <- ft_c$V1[ft_c$V2 %in% pbmc]
ft <- ft_x$V2[ft_x$V1 %in% ft]

# run  aucell
out <- AUCell_run(sce, ft, BPPARAM = bp)
out <- out@assays@data@listData$AUC

sce$pbmc <- out["geneSet", match(colnames(sce), colnames(out))]


df <- as.data.frame(colData(sce))
df <- df %>% select(sum,ctc,pbmc, cell_area, detected)
df <- df[sample(rownames(df)),] # shuffle

ggplot(df, aes(x = ctc, y = pbmc, color = log(sum))) + 
  geom_point(shape = 16, size = 2) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_viridis_c()

# Filter the data to label points above the horizontal line and to the right of the vertical line
thr_pck <- 2500
thr_score <- 0.12
df_lab1 <- df[df$Mean.PanCK > thr_pck & df$score > thr_score, ]

df_lab1$lab <- as.character(1:nrow(df_lab1))
# Visualize

ggplot(df, aes(x = log(score), y = log(pbmc), color = log(sum))) + 
  geom_point(shape = 16, size = 0.5)

 gg_ctc1 <- ggplot(df, aes(x = score, y = sum, color = sum, size = cell_area)) + 
  ggrastr::rasterise(geom_point(shape = 16), dpi = 500) +
  scale_color_viridis_c(option = "H") +
  geom_label_repel(data = df_lab1, aes(label = lab),
                   fill = NA, 
                   direction = "both",
                   label.size = NA,
                   segment.color = NA,
                   size = 5, 
                   box.padding = 0.5,
                   point.padding = 0.3, 
                   force = 1,
                   force_pull = 0.5,
                   max.overlaps = Inf,
                   min.segment.length = 0) +
  scale_size_continuous(range = c(0.1, 5)) +
  scale_color_viridis_c(option = "H") +
  geom_hline(yintercept = thr_pck, color = "red", linetype = "dashed") +
  geom_vline(xintercept = thr_score, color = "red", linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(size = "Area.um2", color = "nCount", x = "MCF7 score", y = "PCK MFI") +
  scale_y_continuous(labels = scientific_10) +
  scale_x_continuous(labels = scientific_10)

gg_ctc1


pltdir <- glue("{proj_dir}/figures/fig4/rds")
saveRDS(gg_ctc1, file = glue("{pltdir}/gg_ctc1.rds"))
saveRDS(df_lab1, file = glue("{pltdir}/df_lab1.rds"))

# save sce & cd
sce$cline <- "PBMCs"
sce$cline[colnames(sce) %in% rownames(df_lab1)] <- "CTCs"

dir <- glue("{proj_dir}/data/stamp_4/processed/CTC1")
dir.create(dir, showWarnings = F)
qsave(sce, file = glue("{dir}/ctc1_sce.qs"), nthreads = 8)
cd <- as.data.frame(colData(sce))
qsave(cd, file = glue("{dir}/ctc1_cd.qs"), nthreads = 8)