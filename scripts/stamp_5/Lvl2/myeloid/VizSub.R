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

dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

# Myeloid ###################################################################### 
################################################################################ 
sub <- "Myeloid"
res_dir <- glue("{proj_dir}/data/stamp_5/processed/{sub}")
myeloid <- qread(glue("{res_dir}/lvl2_sce.qs"), nthreads = 8)

pal <- Polychrome::createPalette(26,c("#99FFFF", "#FF99FF", "#FFFF99"))
names(pal) <- unique(myeloid$lvl2)

feat <- c("S100A9","ITGAM","CD14","IL1B","FCGR3A","CX3CR1","CD68","ITGAX","CD1C","CCR2","CCR5")

plotDots(myeloid, features = intersect(feat, rownames(myeloid)), group = "lvl2", scale =TRUE, center = TRUE) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 20, color = "black")) &
  labs(x = "Cell Lineage", y = "Gene")
# B ############################################################################ 
################################################################################ 
sub <- "B"
res_dir <- glue("{proj_dir}/data/stamp_5/processed/{sub}")
B <- qread(glue("{res_dir}/lvl2_sce.qs"), nthreads = 8)

pal <- Polychrome::createPalette(26,c("#99FFFF", "#FF99FF", "#FFFF99"))
names(pal) <- unique(B$lvl2)


B$lvl2 <- factor(B$lvl2, levels = c("Memory B","Memory B ITGAX+","Unswitched Memory B","Naive B"))
feat <- c("CD19","CD79A","ITGAX","CD79B","CD27","CD38","IGHM","SELL","CCR7","TCL1A")
plotDots(B, features = intersect(feat, rownames(B)), group = "lvl2", scale =TRUE, center = TRUE) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 20, color = "black")) &
  labs(x = "Cell Lineage", y = "Gene")


# T ###################################################################### 
################################################################################ 
sub <- "T"
res_dir <- glue("{proj_dir}/data/stamp_5/processed/{sub}")
T <- qread(glue("{res_dir}/lvl2_sce.qs"), nthreads = 8)

res_dir <- glue("{proj_dir}/data/stamp_5/processed/NK")
NK <- qread(glue("{res_dir}/clust_sce.qs"), nthreads = 8)
NK$lvl2 <- NK$lvl1

T <- cbind(T,NK)

pal <- Polychrome::createPalette(26,c("#99FFFF", "#FF99FF", "#FFFF99"))
names(pal) <- unique(T$lvl2)

feat <- c("CCR5","GZMK","RORC","EOMES","LAG3","ZNF683","GZMH","KLRD1","PDCD1","FASLG","CD8B","CD8A","KLRK1","SELL","CCR7",
          "IL7R","KLRB1","CD38","FCGR3A","NCAM1","CCR6",
         "CCR4","CCL20","CD40LG","CD38","FOXP3","TNFRSF9",
         "LEF1","KLRB1")

T$lvl2 <- factor(T$lvl2, levels = c("Cytotoxic CD8", "Effector Memory CD8","Naive CD8","Naive CD4",
                                    "Central Memory CD4","Effector Memory CD4","NK","T helper","Tregs"))
plotDots(T, features = intersect(feat, rownames(T)), group = "lvl2", scale =TRUE, center = TRUE) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 20, color = "black")) &
  labs(x = "Cell Lineage", y = "Gene")


# Proportions ###################################################################### 
################################################################################ 
sce <- cbind(T,B,myeloid)
sce <- sce[,sce$lvl2 != "LowQ"]

sce$Cell_Lineage <- sce$lvl1
sce$Cell_Type <- sce$lvl2

pal_full <- Polychrome::createPalette(26,c("#99FFFF", "#FF99FF", "#FFFF99"))
names(pal_full) <- unique(sce$Cell_Type)
library(scales)
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
p <-  ggplot(cells, aes(x = "", y = pct, fill = Cluster)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = pal_full) +
    geom_text(aes(label = paste0(pct * 100, "%")), position = position_stack(vjust = 0.5), size = 3) +
    theme_minimal() +
    labs(x = "", y = "Proportion", fill = "", title = subset,
         subtitle = glue("N = {number(sum(cells$Freq), big.mark = '.')} cells")) +
    theme(axis.text.x = element_blank(), axis.text = element_text(size = 15, color = "black"),
          text = element_text(size = 20),
          axis.ticks.x = element_blank(), legend.position = "bottom") + 
    coord_flip()
return(p)
}

wrap_plots(prop(sce,"T"), prop(sce,"B"),prop(sce,"myeloid"), ncol = 1)


subset <- "myeloid"
sub <- sce[, sce$Cell_Lineage == subset]
sub$Cell_Type <- as.character(sub$Cell_Type)

# Calculate cell counts and percentages
cells <- as.data.frame(table(sub[["Cell_Type"]])) %>%
  dplyr::rename(Cluster = Var1) %>%
  mutate(pct = round(Freq / sum(Freq), 3)) %>%
  arrange(pct) %>%  # Arrange by decreasing percentage
  mutate(Cluster = factor(Cluster, levels = Cluster))  # Reorder factor levels

# Plot with ggplot2
ggplot(cells, aes(x = "", y = pct, fill = Cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pal_full) +
  geom_text(aes(label = paste0(pct * 100, "%")), position = position_stack(vjust = 0.5), size = 3) +
  theme_minimal() +
  labs(x = "", y = "Proportion", fill = "", title = subset,
       subtitle = glue("N = {number(sum(cells$Freq), big.mark = '.')} cells")) +
  theme(axis.text.x = element_blank(), axis.text = element_text(size = 15, color = "black"),
        text = element_text(size = 20),
        axis.ticks.x = element_blank(), legend.position = "bottom") + 
  coord_flip()



