# dependencies
library(patchwork)
library(dplyr)
library(ggplot2)
library(glue)
library(extrafont)

# common theme
common_theme <- theme(
  axis.text = element_text(size = 22, color = "black", family = "Times New Roman"),     
  axis.title = element_text(size = 25, color = "black", family = "Times New Roman"),   
  plot.title = element_text(size = 25, color = "black", family = "Times New Roman"),   
  legend.text = element_text(size = 18, color = "black", family = "Times New Roman"),                  
  legend.title = element_text(size = 22, color = "black", family = "Times New Roman"),
  panel.grid = element_blank()
)
# load data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rds_files <- list.files(pattern = "\\.rds$")
names <- sub("\\.rds$", "", rds_files)
# read files and apply the common theme
list2env(setNames(lapply(rds_files, function(x) readRDS(x) & common_theme), names), envir = .GlobalEnv)

a <- wrap_plots(cnumb,count_feat_per_area, ncol = 2) + plot_layout(widths = c(1,2))
fig4 <- wrap_plots(a, gg_corr, ncol = 1) 

# save
dir <- "./../raw"
dir.create(dir, showWarnings = F)

pdf(glue("{dir}/fig4.pdf"), width = 25, height = 22)
wrap_plots(a, gg_corr, ncol = 1) 
dev.off()


