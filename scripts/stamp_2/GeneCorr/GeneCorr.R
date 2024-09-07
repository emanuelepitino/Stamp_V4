suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dplyr)
  library(here)
  library(data.table)
  library(glue)
  library(qs)
})

dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

stamp <- "Stamp_2"
sample <- "MIX"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

outdir <- glue("{proj_dir}/data/{stamp}/processed/flex/flex_clines") # read flex data
flex <- qread(glue("{outdir}/proc_sce.qs"))

dir <- glue("{proj_dir}/data/{stamp}/{sample}/Ist")
cosmx <- qread(glue("{dir}/anno_sce.qs"), nthreads = 8)

# Take matrices
cs <- as(counts(cosmx),"dgCMatrix")
fl <- as(counts(flex), "dgCMatrix")

feat <- intersect(rownames(cs),rownames(fl)) # take features intersection
fl <- fl[feat, ] # subset 
cs <- cs[feat,] # subset

# Aggregate the flex  matrix
df <- as.data.frame(colData(flex))
colnames(fl) <- df[colnames(fl), "cluster"]
#fl <- t(rowsum(t(fl), group = colnames(fl)))
fl <- as(fl, "dgCMatrix")
# Aggregate the cosmx matrix
df <- as.data.frame(colData(cosmx))
colnames(cs) <- df[colnames(cs), "label"]
#cs <- t(rowsum(t(cs), group = colnames(cs)))
cs <- as(cs, "dgCMatrix")

# function to aggregate and normalize
aggr_norm <- function(mat){
  res <- do.call(cbind, lapply(unique(colnames(mat)), function(col_name) {
          sub <- mat[, colnames(mat) %in% col_name, drop = FALSE]
          agg <- rowSums(sub) / ncol(sub)
          agg <- rowSums(sub)
          agg
        }))
colnames(res) <- unique(colnames(mat))
return(res)
}

agg_fl <- as.data.frame(aggr_norm(fl))
agg_fl$gene <- rownames(agg_fl)
agg_cs <- as.data.frame(aggr_norm(cs))
agg_cs$gene <- rownames(agg_cs)


agg_cs <- melt(agg_cs, id.vars = "gene")
names(agg_cs)[names(agg_cs) == "value"] <- "cosmx"

agg_fl <- melt(agg_fl, id.vars = "gene")
names(agg_fl)[names(agg_fl) == "value"] <- "flex"



df <- merge(agg_cs, agg_fl, by = c("gene", "variable"), all = TRUE)

#df <- df %>% filter(flex > 1 & cosmx > 1)
# Plot
ggplot(df, aes(x = flex, y = cosmx)) +
  geom_point(shape = 16, size = 0.6, alpha = 0.8) +
  ggpubr::stat_cor(method = "spearman", label.x = -3, label.y = 2) + 
  scale_y_log10(labels = label_log()) +  #
  scale_x_log10(labels = label_log()) +
  geom_abline(slope = 1, x = 0, color = "red4", linetype = "dashed") +
  labs(title = "LnCAP") +
  theme_bw() +
  theme(
    text = element_text(size = 10, color = "black"), 
    axis.text = element_text(size = 10, color = "black"),
    panel.grid = element_blank()
  ) +
  facet_wrap(~variable)


corr_plot <- function(df,subset){
  
sub <- df[df$variable == subset,]

# Plot with labels for the top 10 genes
plt <- ggplot(sub, aes(x = flex, y = cosmx)) +
  ggrastr::rasterize(geom_point(shape = 16, size = 0.2, alpha = 0.8), dpi = 1000) +
  ggpubr::stat_cor(method = "spearman", label.x = -3, label.y = 10, size = 3) + 
  scale_y_log10() +
  scale_x_log10() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, color = "red4", linetype = "dashed") +
  labs(title = subset) +
  theme_bw() +
  theme(
    text = element_text(size = 10, color = "black"), 
    axis.text = element_text(size = 8, color = "black"),
    panel.grid = element_blank()
  ) 
return(plt)
}

gg_corr <- wrap_plots(
              corr_plot(df,"SKBR3"),
              corr_plot(df,"LnCAP"),
              corr_plot(df,"MCF7"),
              ncol = 3) +
                plot_layout(axis_titles = "collect")

pdf(glue("{proj_dir}/figures/fig1/raw/corr.pdf"), width = 12,height = 3)
gg_corr
dev.off()



