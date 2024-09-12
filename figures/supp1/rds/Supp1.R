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

corr <- corr + labs(subtitl