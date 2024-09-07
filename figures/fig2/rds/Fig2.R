library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(grid)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

cnumb <- readRDS("gg_cnumb.rds")
corr <- readRDS("gg_corr.rds")
sum <- readRDS("gg_sum.rds")
feat <- readRDS("gg_feat.rds")

# load fonts
library(extrafont)
#loadfonts(device = "pdf")
#extrafont::font_import()   

# plot theme
common_theme <- theme(
  axis.text = element_text(size = 10, color = "black", family = "Times New Roman"),     
  axis.title = element_text(size = 17, color = "black", family = "Times New Roman"),   
  plot.title = element_text(size = 25, color = "black", family = "Times New Roman"),   
  legend.text = element_text(size = 12, color = "black", family = "Times New Roman"),                  
  legend.title = element_text(size = 13, color = "black", family = "Times New Roman"),
  panel.grid = element_blank()
)
# apply the same theme
cnumb <- cnumb & common_theme
corr <- corr & common_theme & labs(color = "Sample")
sum <- sum & common_theme 
feat <- feat & common_theme

A_B <- wrap_plots(cnumb, corr)
C_D <- wrap_plots(sum,feat) + plot_layout(axis_titles = "collect")

fig2 <- wrap_plots(A_B,C_D, ncol = 1)

dir <- glue("./../raw")
dir.create(dir, showWarnings = F)
pdf(file = glue("{dir}/fig2.pdf"), width = 15, height = 6)
fig2
dev.off()
