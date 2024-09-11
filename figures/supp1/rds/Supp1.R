library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(grid)
library(extrafont)

# plot theme
common_theme <- theme(
  axis.text = element_text(size = 10, color = "black", family = "Times New Roman"),     
  axis.title = element_text(size = 17, color = "black", family = "Times New Roman"),   
  plot.title = element_text(size = 25, color = "black", family = "Times New Roman"),   
  legend.text = element_text(size = 12, color = "black", family = "Times New Roman"),                  
  legend.title = element_text(size = 13, color = "black", family = "Times New Roman"),
  panel.grid = element_blank()
)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rds_files <- list.files(pattern = "\\.rds$")
names <- sub("\\.rds$", "", rds_files)
# read files and apply the common theme
list2env(setNames(lapply(rds_files, function(x) readRDS(x) & common_theme), names), envir = .GlobalEnv)

corr <- corr + labs(subtitle = "Pearson")
cdens <- cdens + scale_fill_gradientn(colors = c("white","red")) 

um_hm <- wrap_plots(hm,anno, nrow = 1) + plot_layout(widths = c(3,1))


qcmet <- wrap_plots(layout, counts, feat, area,cdens,corr, nrow = 2)

um <- um + theme(legend.position = "top")
hm <- hm & theme(axis.text.x = element_text(size = 16),
           axis.text.y = element_text(size = 10),
           legend.position = "top",
           legend.title = element_blank())

a <- wrap_plots(layout,um,hm, corr, nrow = 1)
b <- wrap_plots(counts, feat, area, prob, nrow = 1)

supp1 <- wrap_plots(a,b, nrow = 2)

dir <- "./../raw/"
dir.create(dir, showWarnings = F)
pdf(file = glue("{dir}/Supp1.pdf"), width = 22, height = 15)
supp1
dev.off()
