suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(dplyr)
  library(here)
  library(scater)
  library(scuttle)
  library(glue)
  library(qs)
  library(parallel)
  library(scran)
  library(BiocParallel)
  library(BiocNeighbors)
  library(BiocSingular)
  library(data.table)
})

# Load bin functions
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

## Load data
myeloid <- qread(glue("{proj_dir}/data/stamp_5/processed/Myeloid/lvl2_sce.qs"), nthreads = 8)
B <- qread(glue("{proj_dir}/data/stamp_5/processed/B/lvl2_sce.qs"), nthreads = 8)
T <- qread(glue("{proj_dir}/data/stamp_5/processed/T/lvl2_sce.qs"), nthreads = 8)
NK <- qread(glue("{proj_dir}/data/stamp_5/processed/NK/proc_sce.qs"), nthreads = 8)
NK$lvl2 <- "NK"

sce <- cbind(myeloid,B,T,NK)
sce <- sce[,sce$lvl2 != "LowQ"]

sce$Cell_Lineage <- sce$lvl1
sce$Cell_Type <- sce$lvl2
sce$Cell_Lineage[sce$Cell_Lineage == "T"] <- "T/NK"
sce$Cell_Lineage <- factor(sce$Cell_Lineage, levels = c("T/NK","B","Myeloid"))


# Lineages dot plot
feats_lin <- c("CD3E","CD3D","CD4","CD8A","CD8B",
               "CD19","CD79A","CD79B",
               "CD14","FCGR3A","CD68","CD86","NGK7","NCAM1","KLRC1","KLRD1","KLRF1")

dot_lin <- plotDots(sce, features = intersect(feats_lin, rownames(sce)), group = "Cell_Lineage", scale =TRUE, center = TRUE) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 20, color = "black")) &
  labs(x = "Cell Lineage")
dot_lin


# Cell Type dot plot
sce$Cell_Type <- as.character(sce$Cell_Type)

sce$Cell_Type[sce$Cell_Type == "Central Memory CD4"] <- "CD4 TN/CM"
sce$Cell_Type[sce$Cell_Type == "Effector Memory CD4"] <- "CD4 TN/CM"
sce$Cell_Type[sce$Cell_Type == "Naive CD4"] <- "CD4 TN/CM"
sce$Cell_Type[sce$Cell_Type == "T helper"] <- "CD4 TEM"

sce$Cell_Type[sce$Cell_Type %in% c("Naive B","Unswitched Memory B")] <- "naive B"
sce$Cell_Type[sce$Cell_Type %in% c("Memory B","Memory B ITGAX+")] <- "B mem. / PC"

feats_type <- c("CD3E","CD3D","CD8B","CD8A","KLRK1","GZMK", "EOMES","LAG3","ZNF683","GZMH","KLRD1","FASLG","PDCD1","KLRB1","CD38",
                "FOXP3","TNFRSF9","IL2RA","CCR6",
          "CCR4","CCL20","RORC","CCR7","IL7R","TCF7",
          "CD40LG","SELL","LEF1","CD27",
          "CD79A","CD79B","CD19","IGHM","TCL1A", "IGHG1","IGHG2","IGHG3","IGHG4","CD24","S100A9","CD14","IL1B","SELL",
          "CD68","ITGAM","CD34","FLT3","CX3CR1","FCGR3A","CLEC9A","CCR5",
            "S100A8","LYZ","CCR2","ITGAX","CD1C")

sce$Cell_Type[sce$Cell_Type == "Classical Mono"] <- "class. mono."
sce$Cell_Type[sce$Cell_Type == "Non Classical Mono"] <- "non-class. mono."
sce$Cell_Type[sce$Cell_Type == "Intermediate Mono"] <- "int. mono."
sce$Cell_Type[sce$Cell_Type == "type II Dendritic Cells"] <- "DC"
sce$Cell_Type[sce$Cell_Type == "Class Switched B"] <- "B mem. / PC"
sce$Cell_Type[sce$Cell_Type == "Naive B"] <- "naive B"
sce$Cell_Type[sce$Cell_Type == "Cytotoxic CD8"] <- "CD8 TEM"
sce$Cell_Type[sce$Cell_Type == "Effector Memory CD8"] <- "CD8 TCM"
sce$Cell_Type[sce$Cell_Type == "Naive/Memory CD4"] <- "CD4 TN/CM"
sce$Cell_Type[sce$Cell_Type == "Naive CD8"] <- "CD8 TN"
sce$Cell_Type[sce$Cell_Type == "Tregs"] <- "Treg"



sce$Cell_Type <- factor(sce$Cell_Type, levels = c("naive CD8","CD8 TEM","CD8 TN","CD8 TCM","NK","Treg","CD4 TEM",
                                                  "CD4 TN/CM",
                                                  "naive B","B mem. / PC",
                                                  "class. mono.","int. mono.","non-class. mono.","DC"))


library(extrafont)
plotDots(sce, features = intersect(feats_type, rownames(sce)), group = "Cell_Type", scale = TRUE, center = TRUE) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, size = 14, color = "black", face = "italic", family = "Times New Roman"),
        axis.text.y = element_text(size = 20,color = "black", family = "Times New Roman"),
        axis.title = element_text(size = 20, color = "black", family = "Times New Roman")) &
  labs(x = "Cell Type", y = "Gene")


pdf("/Users/emanuelepitino/Desktop/dot/Dot.pdf", width= 14/2.54, height=6.2/2.54)
fs <- intersect(feats_type, rownames(sce))
plotDots(sce, features = fs, group = "Cell_Type", scale =TRUE, center = TRUE) +
  scale_color_gradient2("z-scaled\nmean expr.", low="blue4", mid="grey90", high="red4") +
  scale_size_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.5), range=c(0.1, 2)) +
  coord_flip() + 
  theme_minimal(6) + theme(
    aspect.ratio = length(unique(sce$Cell_Type))/length(fs),
    axis.text.y=element_text(color="black"),
    axis.text.x = element_text(angle = 45, hjust = 1, color="black"), 
        panel.grid.major = element_line(linewidth = 0.1, color = "lightgrey"),
        axis.title=element_blank(), 
        legend.key.size = unit(0.5, "lines")) 
dev.off()



 # proportions cell lineages,
sce$Cell_Lineage[sce$Cell_Type == "NK"] <- "T/NK"


pal_lin <- Polychrome::createPalette(26,c("#99FFFF", "#FF99FF", "#FFFF99"))
names(pal_lin) <- unique(sce$Cell_Lineage)

# Calculate cell counts and percentages
sce$Cell_Lineage <- as.character(sce$Cell_Lineage)
cells <- as.data.frame(table(sce$Cell_Lineage)) %>%
  dplyr::rename(Cluster = Var1) %>%
  mutate(pct = round(Freq / sum(Freq), 3)) %>%
  arrange(pct) %>%  # Arrange by decreasing percentage
  mutate(Cluster = factor(Cluster, levels = Cluster))  # Reorder factor levels

# Plot with ggplot2
prop_lin <- ggplot(cells, aes(x = "", y = pct, fill = Cluster)) +
  geom_bar(stat = "identity", width = 5) +
  scale_fill_manual(values = pal_lin) +
  geom_text(aes(label = paste0(pct * 100)), position = position_stack(vjust = 0.5), size = 5) +
  theme_minimal() +
  labs(x = "", y = "Proportion", fill = "", title = "",
       subtitle = glue("PBMCs: {number(sum(cells$Freq), big.mark = ',')} cells")) +
  theme(axis.text.x = element_blank(), axis.text = element_text(size = 15, color = "black"),
        text = element_text(size = 20, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        axis.ticks.x = element_blank(), legend.position = "bottom",
        panel.grid = element_blank(),
        plot.margin = unit(c(0.001, 0.1, 0.001, 0.1), "cm")) +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE, override.aes = list(size = 1)))


pal_full <- Polychrome::createPalette(26,c("#99FFFF", "#FF99FF", "#FFFF99"))
names(pal_full) <- unique(sce$Cell_Type)

# Proportion cell types
prop <- function(sce, subset){
  sub <- sce[, sce$Cell_Lineage == subset]
  sub$Cell_Type <- as.character(sub$Cell_Type)
  
  # Calculate cell counts and percentages
  cells <- as.data.frame(table(sub[["Cell_Type"]])) %>%
    dplyr::rename(Cluster = Var1) %>%
    mutate(pct = round(Freq / sum(Freq), 3)) %>%
    arrange(pct) %>%  # Arrange by decreasing percentage
    mutate(Cluster = factor(Cluster, levels = Cluster))  # Reorder factor levels
  
  # Plot with ggplot2
  p <- ggplot(cells, aes(x = "", y = pct, fill = Cluster)) +
    geom_bar(stat = "identity", width = 5) +
    scale_fill_manual(values = pal_full) +
    geom_text(aes(label = paste0(pct * 100)), position = position_stack(vjust = 0.5), size = 5) +
    theme_minimal() +
    labs(x = "", y = "Proportion", fill = "", title = "",
         subtitle = glue("{subset}: {number(sum(cells$Freq), big.mark = ',')} cells")) +
    theme(axis.text.x = element_blank(), axis.text = element_text(size = 15, color = "black"),
          text = element_text(size = 20, color = "black"),
          legend.text = element_text(size = 16, color = "black"),
          axis.ticks.x = element_blank(), legend.position = "bottom", 
          panel.grid = element_blank(),
          plot.margin = unit(c(0.001, 0.1, 0.001, 0.1),"cm")) +
    coord_flip() +
    guides(fill = guide_legend(reverse = TRUE, nrow = 1, override.aes = list(size = 1)))
  
  return(p)
}

pdf("/Users/emanuelepitino/Desktop/fig_s11/prop.pdf", width = 10, height = 8)
wrap_plots(prop_lin,prop(sce,"T/NK"),prop(sce,"Myeloid"),prop(sce,"B"), ncol = 1) +
  plot_layout(axis_titles = "collect")
dev.off()




qsave(sce, file = glue("{proj_dir}/data/stamp_5/processed/final_obj.qs"))





