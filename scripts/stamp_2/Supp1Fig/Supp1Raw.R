library(glue)
library(here)
library(dplyr)
library(ggplot2)
library(here)
library(qs)

dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

scedir <- glue("{proj_dir}/data/stamp_2/processed/final/")
pooled <- qread(glue("{scedir}/sce_Pooled.qs"), nthreads = 8)
mix <- qread(glue("{scedir}/sce_MIX.qs"), nthreads = 8)

sce <- cbind(pooled,mix)
cd <- as.data.frame(colData(sce))
cd <- cd[sample(rownames(cd)),] # shuffle

pal <- Polychrome::createPalette(26,c("#99FFFF", "#FF99FF", "#FFFF99"))
names(pal) <- unique(sce$label)
## LAYOUT PLOT ##
################################################################################
# Calculate boundaries, center coordinates, and formatted row count for each group
cd_summary <- cd %>%
  group_by(sample) %>%
  summarise(
    xmin = min(CenterX_global_px),
    xmax = max(CenterX_global_px),
    ymin = min(CenterY_global_px),
    ymax = max(CenterY_global_px),
    center_x = (xmin + xmax) / 2,
    center_y = (ymin + ymax) / 2,
    row_count_label = formatC(n(), format = "d", big.mark = ".")
  )

# Plotting
gg_layout <- ggplot(cd_summary) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = sample), 
            fill = "transparent", color = "black") +
  geom_text(aes(x = center_x, y = center_y, label = sample), size = 7, color = "black") +
  geom_text(aes(x = center_x, y = center_y - 0.2 * (ymax - ymin), 
                label = paste("N =", row_count_label)), size = 4, color = "black") +
  coord_fixed() +
  labs(subtitle = "STAMP-C2 layout", x = "", y = "") +
  theme_bw() +
  theme(panel.grid = element_blank())

## CELL DENSITY PLOT ##
################################################################################
gg_cdens <- ggplot(cd, aes( x = CenterX_global_px, y = CenterY_global_px)) +
  stat_density2d(aes(fill = after_stat(level)), geom = "polygon", contour = T) +
  coord_equal() +
  scale_fill_viridis_c(option = "inferno", labels = scientific_10) +
  labs(fill = "", subtitle = "Cell density", x = "", y = "") +
  scale_x_continuous(labels = scientific_10) +
  scale_y_continuous(labels = scientific_10) +
  coord_fixed() +
  theme_bw() +
  theme(panel.grid = element_blank())

## Counts/Feat spatial plot ##
################################################################################
spatqc <- function(var,sub, pal){
  ggplot(cd, aes( x = CenterX_global_px, y = CenterY_global_px, color = log(!!sym(var)))) +
   ggrastr::rasterise(geom_point(shape = 16, size = 0.01), dpi = 500) +
   #scale_color_gradientn(colors = c("navy","red4")) +
    scale_color_viridis_c(option = "inferno") +
    #scale_color_gradientn(colors = c("navy","red")) +
   scale_x_continuous(labels = scientific_10) +
   scale_y_continuous(labels = scientific_10) +
   theme_bw() +
   coord_equal() +
   theme(panel.grid = element_blank()) + 
   labs(x = "", y = "", color = "", subtitle = sub) +
    theme(legend.position = "right")
}

gg_c <- spatqc("sum","nCount")
gg_f <- spatqc("detected","nFeature")
gg_a <- spatqc("Area.um2", "Area.um2")

gg_qcmet <- wrap_plots(gg_c,gg_f,gg_a, ncol = 3) + plot_layout(axis_titles = "collect")

## Annotation ##
################################################################################
gg_anno <- ggplot(cd, aes( x = CenterX_global_px, y = CenterY_global_px, color = label)) +
  ggrastr::rasterise(geom_point(shape = 16, size = 0.01), dpi = 500) +
  scale_x_continuous(labels = scientific_10) +
  scale_color_manual(values = pal) +
  scale_y_continuous(labels = scientific_10) +
  theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank()) + 
  labs(x = "", y = "", subtitle = "Annotation", color = "") +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 5)))

## InSituType Cluster Assignment Probability ##
################################################################################

gg_prob <- ggplot(cd, aes( x = CenterX_global_px, y = CenterY_global_px, color = label_prob)) +
  ggrastr::rasterise(geom_point(shape = 16, size = 0.01), dpi = 500) +
  scale_x_continuous(labels = scientific_10) +
  scale_color_gradientn(colors = c("navy","red")) +
  scale_y_continuous(labels = scientific_10) +
  theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank()) + 
  labs(x = "", y = "", color = "", subtitle = "InSituType cluster probability") +
  theme(legend.position = "right")



#wrap_plots(gg_c,gg_f,gg_a,gg_cdens,gg_prob, gg_anno, nrow = 2) 

## Correlation of gene expression CosMx - Flex ##
################################################################################
outdir <- glue("{proj_dir}/data/stamp_2/processed/flex/flex_clines") 
f <- qread(glue("{outdir}/proc_sce.qs"))
c <- sce

# take intersection
feat <- intersect(rownames(c),rownames(f))
# subset matrices
c <- c[feat,] 
f <- f[feat,]
# Take count
c <- as(counts(c), "dgCMatrix")
f <- as(counts(f), "dgCMatrix")
# rowsum & normalize
c <- rowSums(c) / ncol(c)
f <- rowSums(f) / ncol(f)
# convert
c <- as.data.frame(as.matrix(c))
f <- as.data.frame(as.matrix(f))
# rename
names(c)[1] <- "CosMx"
names(f)[1] <- "Flex"
# add gene names
c$gene <- rownames(c)
f$gene <- rownames(f)
df <- merge(c, f, by = "gene", all = TRUE)
# plot
gg_corr <- ggplot(df, aes(x = Flex, y = CosMx)) + 
  ggrastr::rasterize(geom_point(shape = 16, size = 0.9, alpha = 0.8), dpi = 300) +
  ggpubr::stat_cor(method = "spearman", label.x = -3, label.y = 2, size = 5) + 
  geom_smooth() +
  scale_y_log10(labels = scientific_10) +
  scale_x_log10(labels = scientific_10) +
  theme_bw() +
  theme(panel.grid =  element_blank()) +
  labs(x = "Mean nCount - Flex", y = "Mean nCount - CosMx")

gg_corr

### CORR SPLIT BY CELL TYPE
# Take matrices
cs <- as(counts(c),"dgCMatrix")
fl <- as(counts(f), "dgCMatrix")

feat <- intersect(rownames(cs),rownames(fl)) # take features intersection
fl <- fl[feat, ] # subset 
cs <- cs[feat,] # subset

# Aggregate the flex  matrix
df <- as.data.frame(colData(f))
colnames(fl) <- df[colnames(fl), "cluster"]
#fl <- t(rowsum(t(fl), group = colnames(fl)))
fl <- as(fl, "dgCMatrix")
# Aggregate the cosmx matrix
df <- as.data.frame(colData(c))
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


corr_plot <- function(df,subset){
  
  genes_skbr3 <- c("INHBB", "S100A8", "ARHGDIB", "ERBB2", "S100A4", "LGALS3BP", "S100A6", "S100A9", "DHRS2", "KRT7")
  genes_mcf7 <- c("CXCL12", "BASP1", "IL20", "MGP", "EPHA4", "TGFB1", "ESR1", "AREG", "LGALS3", "GATA3")
  genes_lncap <- c("VWF", "CSF2", "TNNT2", "DDR2", "SERPINB5", "EPHA7", "LIFR", "RPS4Y1", "TNFRSF19", "KLK3")
  
  if(subset == "MCF7"){ feat_color = genes_mcf7}
  if(subset == "LnCAP"){ feat_color = genes_lncap}
  if(subset == "SKBR3"){ feat_color = genes_skbr3}
  
  sub <- df[df$variable == subset,]
  # Plot 
  plt <- ggplot(sub, aes(x = flex, y = cosmx)) +
    geom_point(aes(color = ifelse(gene %in% feat_color,"red","black")),
                                  shape = 16, size = 1, alpha = 0.8) +
    ggpubr::stat_cor(method = "spearman", label.x = -3, label.y = 6, size = 3) + 
    scale_y_log10() +
    scale_x_log10() +
    geom_abline(slope = 1, intercept = 0, color = "red4", linetype = "dashed") +
    labs(title = subset, color = "Is Marker") +
    theme_bw() +
    theme(
      text = element_text(size = 10, color = "black"), 
      axis.text = element_text(size = 8, color = "black"),
      panel.grid = element_blank()
    ) +
    labs(x = "Mean nCount - Flex", y= "Mean nCount - CosMx")
  return(plt)
}

gg_corr <- wrap_plots(
  corr_plot(df,"SKBR3"),
  corr_plot(df,"LnCAP"),
  corr_plot(df,"MCF7"),
  ncol = 3) +
  plot_layout(axis_titles = "collect")

## SAVE
dir <- glue("{proj_dir}/figures/supp1/rds")
dir.create(dir, showWarnings = F, recursive = T)

saveRDS(gg_layout, glue("{dir}/layout.rds"))
saveRDS(gg_cdens, glue("{dir}/cdens.rds"))
saveRDS(gg_a, glue("{dir}/area.rds"))
saveRDS(gg_f, glue("{dir}/feat.rds"))
saveRDS(gg_c, glue("{dir}/counts.rds"))
saveRDS(gg_anno, glue("{dir}/anno.rds"))
saveRDS(gg_prob, glue("{dir}/prob.rds"))
saveRDS(gg_corr, glue("{dir}/corr.rds"))


corr_plot <- function(df, subset) {
  
  genes_skbr3 <- c("INHBB", "S100A8", "ARHGDIB", "ERBB2", "S100A4", "LGALS3BP", "S100A6", "S100A9", "DHRS2", "KRT7")
  genes_mcf7 <- c("CXCL12", "BASP1", "IL20", "MGP", "EPHA4", "TGFB1", "ESR1", "AREG", "LGALS3", "GATA3")
  genes_lncap <- c("VWF", "CSF2", "TNNT2", "DDR2", "SERPINB5", "EPHA7", "LIFR", "RPS4Y1", "TNFRSF19", "KLK3")
  if(subset == "MCF7") { feat_color = genes_mcf7 }
  if(subset == "LnCAP") { feat_color = genes_lncap }
  if(subset == "SKBR3") { feat_color = genes_skbr3 }
  
  sub <- df[df$variable == subset, ]
  # Create a new variable to indicate if a gene is a marker
  sub$Is_Marker <- ifelse(sub$gene %in% feat_color, "Yes", "No")
  
  # Separate the data
  sub_no <- sub[sub$Is_Marker == "No", ]
  sub_yes <- sub[sub$Is_Marker == "Yes", ]
  # Plot
  plt <- ggplot(sub, aes(x = flex, y = cosmx)) +
    ggpubr::stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) + 
    # Plot "No" points first
    geom_point(data = sub_no, aes(color = Is_Marker), shape = 16, size = 2, alpha = 0.8) +
    # Plot "Yes" points on top
    geom_point(data = sub_yes, aes(color = Is_Marker), shape = 16, size = 2, alpha = 0.8) +
    scale_y_log10() +
    scale_x_log10() +
    geom_abline(slope = 1, intercept = 0, color = "red4", linetype = "dashed") +
    labs(title = subset, color = "Is Marker") +
    theme_bw() +
    theme(
      text = element_text(size = 10, color = "black"), 
      axis.text = element_text(size = 8, color = "black"),
      panel.grid = element_blank(),
      legend.position = "none"  # Remove legend for combining later
    ) +
    labs(x = "Mean counts/gene - Flex", y= "Mean counts/gene - CosMx") +
    scale_color_manual(values = c("Yes" = "red", "No" = "grey"))
  
  return(plt)
}

wrap_plots(corr_plot(df, "SKBR3"),
           corr_plot(df, "LnCAP"),
           corr_plot(df, "MCF7"),
           ncol = 3) +
  plot_layout(ncol = 3, guides = "collect", axis_titles = "collect") &
  theme(legend.position = 'bottom')
