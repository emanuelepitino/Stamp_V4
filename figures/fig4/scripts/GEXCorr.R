## Dependencies
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(qs)
  library(glue)
  library(here)
  library(scales)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
  library(dplyr)
})
## Data
cd <- qread("./../obj/cd.qs")
sce <- qread("./../obj/sce.qs")
source(glue("{here()}/scripts/misc/paths.R"))
source(glue("{here()}/scripts/misc/BIN.R"))
## Plot theme
# common theme
common_theme <- theme(text = element_text(color = "black", family = "Times New Roman", size = 20),
                      axis.text = element_text(color = "black", family = "Times New Roman", size = 18))

# Split by stamp and substamp, & lognorm counts
stamp <- unique(cd$stamp)
sub <- unique(cd$sub)
for(st in stamp){
  for(s in sub){
    tmp <- sce[,sce$stamp == st & sce$sub == s]
    tmp <- logNormCounts(tmp) # lognormalize
    tmp <- assay(tmp,"counts") 
    #tmp <- rowSums(tmp) / ncol(tmp)
    assign(paste0(st,"_",s),tmp)
  }
}

mcf7 <- data.frame(stamp_11 = stamp_11_MCF7, stamp_9 = stamp_9_MCF7, stamp_12= stamp_12_MCF7)
skbr3 <-  data.frame(stamp_11 = stamp_11_SKBR3, stamp_9 = stamp_9_SKBR3, stamp_12 = stamp_12_SKBR3)

df <- cbind(mcf7,skbr3)

df <- 
# Generate all unique pairwise combinations
  
corr_plot <- function(df){
pairwise_combinations <- combn(colnames(df), 2, simplify = FALSE)
# Create a list of ggplots for each pair
plot_list <- lapply(pairwise_combinations, function(pair) {
  ggplot(mcf7, aes_string(x = pair[1], y = pair[2])) +
    ggrastr::rasterise(geom_point(shape = 16, size = 0.5, color = "blue4"), dpi = 500) +
    geom_abline(slope = 1, color = "red4", linewidth = 0.3) +
    stat_cor(
      method = "pearson",
      label.x.npc = "left",
      label.y.npc = "top",
      size = 5
    ) +
    theme_bw() +
    theme(panel.grid = element_blank()) 
})
# Assign meaningful names to each plot in the list
names(plot_list) <- sapply(pairwise_combinations, function(pair) {
  paste(pair, collapse = "_vs_")
})
p <- wrap_plots(plot_list)
return(p)
}

gg_mcf7 <- corr_plot(mcf7) + 
  labs(subtitle = "MCF7") +
  theme(plot.subtitle = element_text(face = "bold")) &
  common_theme 

gg_skbr3 <- corr_plot(skbr3) + 
  labs(subtitle = "SKBR3") + 
  theme(plot.subtitle = element_text(face = "bold")) &
  common_theme

# Combine the plots vertically
gg_corr <- wrap_plots(gg_mcf7,gg_skbr3, ncol = 1)

dir <- glue("./../rds")
dir.create(dir,showWarnings = F)
saveRDS(gg_corr,file = glue("{dir}/gg_corr.rds"))