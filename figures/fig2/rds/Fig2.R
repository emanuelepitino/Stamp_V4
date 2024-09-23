library(ggplot2)
library(patchwork)
library(grid)
library(glue)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pos <-  readRDS("gg_pos.rds")
cnumb <- readRDS("gg_cnumb.rds")
sum <- readRDS("gg_sum.rds")
feat <- readRDS("gg_feat.rds")
area <- readRDS("gg_area.rds")
corr <- readRDS("gg_corr.rds")
gene_corr <- readRDS("gg_gene_corr.rds")
library(extrafont)
#loadfonts(device = "pdf")
#extrafont::font_import()   

# plot theme
common_theme <- theme(
  axis.text = element_text(size = 18, color = "black", family = "Times New Roman"),     
  axis.title = element_text(size = 22, color = "black", family = "Times New Roman"),   
  plot.title = element_text(size = 25, color = "black", family = "Times New Roman"),   
  legend.text = element_text(size = 18, color = "black", family = "Times New Roman"),                  
  legend.title = element_text(size = 22, color = "black", family = "Times New Roman"),
  panel.grid = element_blank()
)
# apply the same theme
pos <- pos & common_theme & theme(legend.position = "right")
cnumb <- cnumb & common_theme
corr <- corr & common_theme & labs(color = "Sample")
sum <- sum & common_theme 
feat <- feat & common_theme
area <- area & common_theme
gene_corr <- gene_corr & common_theme

a <- wrap_plots(pos,cnumb,corr, nrow = 1) + plot_layout(widths = c(2,2,2))
b <- wrap_plots(sum,feat,area, gene_corr, ncol = 2)

fig3 <- wrap_plots(a,b, ncol = 1) + plot_layout(heights = c(1,2))

pdf(file = glue("{dir}/../fig3.pdf"), width = 25, height = 18)
fig3
dev.off()
# Second panel
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
p2pos <- p2pos & common_theme 
p2cnumb <- p2cnumb & common_theme
p2corr <- p2corr & common_theme
p2sum <- p2sum & common_theme 
p2feat <- p2feat & common_theme

qcmet <- wrap_plots(p2sum,p2feat, nrow = 2) +
  plot_layout(axis_titles = "collect") 

p2 <- wrap_plots(p2pos,qcmet,p2corr, nrow = 1) +
  plot_layout(widths = c(3,1,2))

fig2 <- wrap_plots(p1,p2, nrow = 2) +
  plot_layout(heights = c(2,1))


dir <- glue("./../raw")
dir.create(dir, showWarnings = F)
pdf(file = glue("{dir}/fig2.pdf"), width = 25, height = 15)
fig2
dev.off()
