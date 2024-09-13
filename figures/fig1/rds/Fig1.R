library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(grid)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

cnumb <- readRDS("cnumb.rds")
bordeff <- readRDS("bord_eff.rds")
qcmet <- readRDS("qc.rds")
hm <- readRDS("hm.rds")
um <- readRDS("um.rds")

# load fonts
library(extrafont)

# plot theme
common_theme <- theme(
  axis.text = element_text(size = 25, color = "black", family = "Times New Roman"),    
  axis.title = element_text(size = 28, color = "black", family = "Times New Roman"),   
  plot.title = element_text(size = 22, color = "black", family = "Times New Roman"),
  legend.text = element_text(size = 20, color = "black", family = "Times New Roman"),  
  legend.title = element_text(size = 26, color = "black", family = "Times New Roman"),
  plot.subtitle = element_text(size = 28, color = "black", family = "Times New Roman"),
  panel.grid = element_blank())
# apply the same theme
cnumb <- cnumb & common_theme
bordeff <- bordeff & common_theme
qcmet <- qcmet & common_theme 
hm <- hm   & common_theme & theme(legend.title = element_blank(), axis.text = element_text(size = 28))
um <- um & common_theme & guides(color = guide_legend(override.aes = list(size = 4))) & 
  labs(color ="Cell Type") &
  theme(legend.position = "top")



B_C <- wrap_plots(qcmet,bordeff) + plot_layout(widths = c(1,3)) 
B_C_E <- wrap_plots(B_C,hm, ncol = 1) + plot_layout(heights = c(3,1))

D_F <- wrap_plots(um,cnumb, ncol = 1) + plot_layout(heights = c(1,3))

fig1 <- wrap_plots(B_C_E,D_F, ncol = 2) + plot_layout(widths = c(5,1))

pdf(file = "./../raw/fig1.pdf", width = 38, height = 25, bg = "transparent")
fig1
dev.off()