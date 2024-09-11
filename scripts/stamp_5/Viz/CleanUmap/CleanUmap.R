# Dependencies
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(dplyr)
  library(here)
  library(scater)
  library(scuttle)
  library(glue)
  library(qs)
  library(scran)
  library(scales)
})

# read data
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

sce <-  qread(glue("{proj_dir}/data/stamp_5/processed/combined_sce.qs"), nthreads = 8)
#sce <- sce[,sample(colnames(sce),1000)]
sce

# palette
val <- c(unique(as.character(sce$lvl1)),unique(sce$lvl2))
pal <- Polychrome::createPalette(26,c("#99FFFF", "#FF99FF", "#FFFF99"))
names(pal) <- val

um <- as.data.frame(reducedDim(sce,"UMAP"))

cd <- as.data.frame(colData(sce))
cd <- cbind(cd, um[rownames(cd), ])

cd <- cd[sample(rownames(cd)),] # shuffle


# UMAP cell lineages
gg_um_clineages <- ggplot(cd, aes(x = UMAP1, y = UMAP2, color = lvl1)) + 
  ggrastr::rasterise(geom_point(shape = 16, size = 0.1), dpi = 800) +
  scale_color_manual(values = pal) +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "top") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(color = "Cell Lineage", x= "x_px", y = "y_px") +
  guides(fill = guide_legend(nrow = 1))  +
  coord_equal()


# spatial plots
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
gg_ncount <- ggplot(cd, aes(x = x_centroid, y = y_centroid, color = log(sum))) + 
  ggrastr::rasterise(geom_point(shape = 16, size = 0.1), dpi = 800) + 
  scale_color_gradientn(colors = c("red4","navy")) +
  coord_equal() + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  scale_x_continuous(labels = scientific_10) +
  scale_y_continuous(labels = scientific_10) +
  labs(color = "nCount", x= "x_px", y = "y_px")

gg_nfeat <- ggplot(cd, aes(x = x_centroid, y = y_centroid, color = log(detected))) + 
  ggrastr::rasterise(geom_point(shape = 16, size = 0.1), dpi = 800) + 
  scale_color_gradientn(colors = c("red4","navy")) +
  coord_equal() + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  scale_x_continuous(labels = scientific_10) +
  scale_y_continuous(labels = scientific_10) +
  labs(color = "nFeature", x= "x_px", y = "y_px")

gg_dens <- ggplot(cd, aes(x = x_centroid, y = y_centroid)) +
  stat_density2d(aes( fill = after_stat(level)), geom = "polygon") +
  scale_fill_gradientn(colors = c("red4","navy")) +
  coord_equal() + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  scale_x_continuous(labels = scientific_10) +
  scale_y_continuous(labels = scientific_10) +
  labs(fill = "Density", x = "x_px", y = "y_px")

gg_ctype <- ggplot(cd, aes(x = x_centroid, y = y_centroid, color = lvl2)) + 
  ggrastr::rasterise(geom_point(shape = 16, size = 0.1), dpi = 800) + 
  scale_color_manual(values = pal) +
  coord_equal() + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  scale_x_continuous(labels = scientific_10) +
  scale_y_continuous(labels = scientific_10) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(color = "Cell Type", x= "x_px", y = "y_px")


# DotPlot
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
markers <- scoreMarkers(sce, sce$lvl2, BPPARAM = bp)
feat <- lapply(markers, function(df) {
  as.data.frame(df) %>%
    arrange(desc(median.logFC.detected)) %>%
    head(5) %>%
    rownames()
})
feat <- unique(unlist(feat))

#feat <- unique(c("LILRB4","CD14","CD86","CLEC12A","CSF1R","FCGR3A","CD1C","NCAM1","KLRC1","KLRD1","GNLY",
#         "GZMK","PRF1","GZMA","GZMB","IGHG1","IGHG3","IGHG4","CD80","CD19","CD79A",
#          "CD79B","CD8B","CD8A","CD3D","CD79B","CD28","CCR7","TCF7","CD4","CTLA4","CD3E","FOXP3","IL2RA"))

#feat <- c("CD4","CD8A","CD8B","NCAM1","FCGR3A","CD14","CD1C","IGHG1","GZMB","GZMK","CCR7","TCF7")

#sce <-  qread(glue("{proj_dir}/data/stamp_5/processed/combined_sce.qs"), nthreads = 8)
#sce$lvl2 <- factor(sce$lvl2, levels = c("pDCs","CD14+ Mono","CD16+ Mono","Plasmacells","CDC1+ DCs",
 #                                       "CD56dim CD16+ NK","CD56+ CD16- NK",
 #                                       "NK-T","GZMB+ CD8","GZMK+ CD8",
 #                                       "Naive CD8","Naive CD4","Th","Tregs",
 #                                       "Naive B","Class-Switched B"))
                                        
dotplot <- plotDots(sce, features = feat, group = "lvl2", scale =TRUE, center = TRUE) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


# Proportions
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
# full
df <- as.data.frame(table(as.character(colData(sce)$lvl1))) %>%
  mutate(pct = round(Freq / sum(Freq), 2)) %>%
  arrange(((pct)))  # Order by decreasing percentage

labels_with_pct <- paste0(df$Var1, " (", df$pct * 100, "%)")

full <- ggplot(df, aes(x = "", y = Freq, fill = factor(Var1, levels = Var1))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pal, labels = labels_with_pct) +
  labs(fill = "", x = "", y = "# Cells", subtitle = glue("Whole PBMCs; N = {comma(sum(df$Freq), big.mark = '.')}")) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  scale_y_continuous(labels = scientific_10) + 
  coord_flip() +
  guides(fill = guide_legend(nrow = 1)) 


# subsets
plot_stacked_bar <- function(sce, filter_value) {
  df <- as.data.frame(colData(sce)[sce$lvl1 == filter_value, ])
  df$lvl2 <- as.character(df$lvl2)
  
  df <- as.data.frame(table(df$lvl2)) %>%
    mutate(pct = round(Freq / sum(Freq), 2)) %>%
    arrange(pct)  # Order by decreasing percentage
  
  labels_with_pct <- paste0(df$Var1, " (", df$pct * 100, "%)")
  
  ggplot(df, aes(x = "", y = Freq, fill = factor(Var1, levels = Var1))) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = pal, labels = labels_with_pct) +
    labs(fill = "", x = "", y = "# Cells", subtitle = glue("{filter_value}; N = {comma(sum(df$Freq), big.mark = '.')}")) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "bottom") +
    scale_y_continuous(labels = scientific_10) +
    coord_flip() +
    guides(fill = guide_legend(nrow = 1)) 
}

gg_prop <- wrap_plots(
full,
plot_stacked_bar(sce,"T"),
plot_stacked_bar(sce,"Myeloid"),
plot_stacked_bar(sce,"NK"),
plot_stacked_bar(sce,"B"),
nrow = 5) + 
plot_layout(axis_titles = "collect")


# QC metrics density
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

df <- as.data.frame(colData(sce)) %>%
      select(lvl1,sum,detected,cell_area, nucleus_area)
df$lvl1 <- factor(df$lvl1, levels = c("Myeloid","NK","T","B"))

ggplot(df, aes(x = sum, fill = lvl1, color = lvl1)) +
  geom_density(alpha = 0.1, size = 0.5) +  # Increase line width with 'size'
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  theme_bw() +
  theme(panel.grid = element_blank())


create_boxplot <- function(df, var_name, title, ylab) {
  # Dynamically calculate the median for each sample and create labels
  labels <- setNames(
    lapply(unique(df$lvl1), function(cline) {
      median_value <- round(median(df[[var_name]][df$lvl1 == cline], na.rm = T), 0)
      glue("{median_value}")
    }),
    unique(df$lvl1)
  )
  
  # Create the plot
ggplot(df, aes(x = lvl1, y = !!sym(var_name), color = lvl1)) +
    geom_boxplot(aes(fill = lvl1), alpha = 0.1, outlier.shape = NA) + # Set alpha for fill only
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_color_manual(values = pal, labels = labels) +
    scale_fill_manual(values = pal) +
    labs(color = title) +
    labs(x = "Cell Lineage",y = ylab, color = "Median") +
    guides(fill = "none") +
    scale_y_continuous(labels = scientific_10) 
}

gg_qcmet <- wrap_plots(
create_boxplot(df, "sum","lvl1","nCount"),
create_boxplot(df, "detected","lvl1","nFeature"),
create_boxplot(df, "cell_area","lvl1","Cell Area"),
create_boxplot(df, "nucleus_area","lvl1","Nucleus Area"),
ncol = 2) +
plot_layout(axis_titles = "collect")





## FORMAT 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
library(extrafont)
#loadfonts(device = "pdf")
#extrafont::font_import()   

# plot theme
common_theme <- theme(
  axis.text = element_text(size = 10, color = "black", family = "Times New Roman"),     
  axis.title = element_text(size = 17, color = "black", family = "Times New Roman"),
  plot.subtitle =  element_text(size = 18, color = "black", family = "Times New Roman"),
  plot.title = element_text(size = 25, color = "black", family = "Times New Roman"),   
  legend.text = element_text(size = 10, color = "black", family = "Times New Roman"),                  
  legend.title = element_text(size = 13, color = "black", family = "Times New Roman"),
  panel.grid = element_blank()
)
# apply the same theme
um <- gg_um_clineages & common_theme
prop <- gg_prop & common_theme 
dotplot <- dotplot & common_theme
qcmet <- gg_qcmet & common_theme


A_B <- wrap_plots(um, dotplot) + plot_layout(widths = c(1,4))
C_D <- wrap_plots(prop, qcmet) + plot_layout(widths = c(1,1))
fig3 <- wrap_plots(A_B, C_D, ncol = 1) + plot_layout(heights = c(1,1.5))

pdf(glue("{proj_dir}/figures/fig3/fig3.pdf"), width = 20, height = 15)
fig3
dev.off()