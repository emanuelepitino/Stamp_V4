library(ggplot2)
library(qs)
library(glue)
library(here)



dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

sce <- qread(glue("{proj_dir}/data/stamp_5/processed/Lvl2/annotated_sce.qs"))

cd <- as.data.frame(colData(sce))

cd$Cell_Lineage <- as.character(cd$Cell_Lineage)
cd$Cell_Lineage[cd$Cell_Lineage == "T"] <- "T/NK"
cd$Cell_Lineage <- factor(cd$Cell_Lineage, levels = c("B","T/NK","Myeloid"))

pal_lin <- Polychrome::createPalette(26,c("#99FFFF", "#FF99FF", "#FFFF99"))
names(pal_lin) <- unique(cd$Cell_Lineage)

create_boxplot <- function(df, var_name, title, ylab, log) {
  # Calculate medians
  medians <- df %>%
    group_by(Cell_Lineage) %>%
    summarise(median_value = median(!!sym(var_name), na.rm = TRUE))
  
  # Plot
  # Calculate the maximum value for each Cell_Lineage
  max_values <- df %>%
    group_by(Cell_Lineage) %>%
    summarize(max_value = max(!!sym(var_name)))
  
  # Adjust the plot
  plot <- ggplot(df, aes(x = Cell_Lineage, y = !!sym(var_name), color = Cell_Lineage)) +
    ggrastr::rasterise(geom_boxplot(aes(fill = Cell_Lineage), alpha = 0.5, outlier.size = 0.1), dpi = 500) +
    geom_text(data = medians, 
              aes(x = Cell_Lineage, y = max_values$max_value * 1.1, label = round(median_value, 0)), 
              vjust = 0.5, size = 5, color = "black") + 
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_color_manual(values = pal_lin) +
    scale_fill_manual(values = pal_lin) +
    labs(color = title, y = ylab, x = "", color = "Median", fill = "Cell Lineage") +
    guides(color = "none") +
    theme(text = element_text(size = 20, color = "black"),
          axis.text = element_text(size = 15, color = "black"),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "top") +
    coord_flip()
  
    
if(log == TRUE) {
  plot <- plot + scale_y_log10(limits = c(5, 600))
}
    return(plot)
}


cfa <- wrap_plots(create_boxplot(cd, "sum", "Cell_Lineage", "nCount", log = TRUE),
           create_boxplot(cd, "detected", "Cell_Lineage", "nFeature", log = TRUE),
           create_boxplot(cd, "cell_area", "Cell_Lineage", "Cell Area (um2)", log = FALSE),
           ncol = 1) + 
  plot_layout(axis_titles = "collect", guides = "collect") & 
  theme(legend.position = "top")

cfa

pdf("/Users/emanuelepitino/Desktop/fig_s11/cfa.pdf", height = 6, width = 5)
cfa
dev.off()








