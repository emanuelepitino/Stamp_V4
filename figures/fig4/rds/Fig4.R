library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(grid)
library(scales)
library(ggtext)
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

# read files and apply the common theme
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # list files in dir
rds_files <- list.files(pattern = "\\.rds$")
names <- sub("\\.rds$", "", rds_files)
names <- names[!grepl("df", names)] # remove df_lab from those 
rds_files <- rds_files[!grepl("df",rds_files)]
list2env(setNames(lapply(rds_files, function(x) readRDS(x) & common_theme), names), envir = .GlobalEnv)
# read label df for the plot
df_lab1 <- readRDS("./df_lab1.rds")
df_lab2 <- readRDS("./df_lab2.rds")


gg_ctc1

gg_ctc2p2

# read cd to add cell numbers
cd1 <- qread(glue("{proj_dir}/data/stamp_4/processed/CTC1/ctc1_cd.qs"), nthreads = 8)
cd2 <- qread(glue("{proj_dir}/data/stamp_4/processed/CTC2/ctc2_cd.qs"), nthreads = 8)

# Define the function to add subtitle with percentages
add_ctc_subtitle <- function(plot, data, ctc_column = "cline", ctc_value = "CTCs") {
  # Calculate the values
  ctc_count <- nrow(data[data[[ctc_column]] == ctc_value,])
  total_count <- nrow(data)
  formatted_total <- format(total_count, big.mark = '.')
  percentage <- (ctc_count / total_count) * 100
    plot <- plot + 
    labs(subtitle = glue("CTC: {ctc_count} / {formatted_total} - **{format(round(percentage,4), scientific = FALSE)}%**")) +
    theme(plot.subtitle = element_markdown())
  return(plot)
}
# Apply the function to gg_ctc2p2 with cd2
gg_ctc2p2 <- add_ctc_subtitle(gg_ctc2p2, cd2)
# Apply the function to gg_ctc1 with cd1
gg_ctc1 <- add_ctc_subtitle(gg_ctc1, cd1)

gg_dot <- gg_dot & coord_flip()


ctc2 <- wrap_plots(gg_ctc2,gg_ctc2p2, nrow = 1) 
a <- wrap_plots(gg_dot, gg_ctc1, nrow = 1) 
fig4 <- wrap_plots(a, ctc2, nrow = 2)

dir <- "./../raw"
dir.create(dir, showWarnings = F)
pdf(glue("{dir}/fig4.pdf"), width = 22, height = 15)
fig4
dev.off()



