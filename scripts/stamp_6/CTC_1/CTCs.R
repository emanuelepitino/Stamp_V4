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
library(data.table)
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
ft_x <- fread(glue("{proj_dir}/data/{stamp}/CTC_1/raw/cell_feature_matrix/features.tsv.gz"), sep = "\t", nThread = 3, header = F)
ft_c <- fread(glue("{proj_dir}/data/{stamp}/features.tsv.gz"), sep = "\t", nThread = 3, header = F)

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

pbmc_thr <- 0.20
ctc_thr <- 0.60
ggplot(df, aes(x = ctc, y = pbmc, color = log(sum))) + 
  geom_point(shape = 16, size = 1.5) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_viridis_c() +
  geom_hline(yintercept = pbmc_thr, color = "red") +
  geom_vline(xintercept = ctc_thr, color = "red") 

cd <- df
cd$group <- "pbmc"
cd$group[cd$ctc > ctc_thr & cd$pbmc < pbmc_thr] <- "ctc"

sce_raw <- qread(paste0(proj_dir,"/data/stamp_6/CTC_1/processed/raw_sce.qs"))
rawcd <- as.data.frame(colData(sce_raw))
rawcd$group <- "pbmc"
rawcd$group[rawcd$cell_id %in% rownames(cd[cd$group == "ctc",])] <- "ctc"

rawcd <- rawcd %>% select(cell_id,group)
rownames(rawcd) <- NULL

write.csv(rawcd, file = glue("{proj_dir}/data/stamp_6/CTC_1/md.csv"), quote = F)

pltdir <- glue("{proj_dir}/figures/fig4X/rds")
saveRDS(gg_ctc1, file = glue("{pltdir}/gg_ctc1.rds"))
saveRDS(df_lab1, file = glue("{pltdir}/df_lab1.rds"))

# save sce & cd
sce$cline <- "PBMCs"
sce$cline[colnames(sce) %in% rownames(df_lab1)] <- "CTCs"

dir <- glue("{proj_dir}/data/stamp_6/processed/CTC1")
dir.create(dir, showWarnings = F)
qsave(sce, file = glue("{dir}/ctc1_sce.qs"), nthreads = 8)
cd <- as.data.frame(colData(sce))
qsave(cd, file = glue("{dir}/ctc1_cd.qs"), nthreads = 8)