# dependencies
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dplyr)
  library(here)
  library(glue)
  library(qs)
})

# setup
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))
stamp <- "stamp_3"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

# read cd
dir <- glue("{proj_dir}/data/{stamp}/processed/combined")
fl <- list.files(dir, full.names = T)
cd <- lapply(fl, qread)

# filter
columns_to_keep <- c("sample", "sum", "detected", "fov", "Area.um2", "label")
cd <- lapply(cd, function(x) {
  x[, columns_to_keep, drop = FALSE]
})
cd <- do.call(rbind,cd)

#  pal
pal <- Polychrome::createPalette(26,c("#99FFFF", "#FF99FF", "#FFFF99"))
names(pal) <- unique(cd$label)

# CELL NUMBERS
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######
cd$sample <- factor(cd$sample, levels = c("c100","c250","c500","c1000","c20k","n20k"))
df <- as.data.frame(table(cd$sample, cd$label)) %>%
  group_by(Var1) %>%
  mutate(Proportion = round((Freq / sum(Freq)) * 100, 2))

gg_cnumb <- ggplot(df, aes(x = Var1, y = Proportion)) + 
  geom_col(aes(fill = Var2), alpha = 0.9) +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(text = element_text(size = 25, color = "black"),
        axis.text = element_text(size = 20, color = "black"),
        legend.position = "none",
        panel.grid = element_blank()) +
  labs(fill = "Cell Type", y = "Percentage") +
  labs(x = "Sample", fill = "Cluster") +
  geom_hline(yintercept = median(df$Proportion[df$Var2 == "SKBR3"]), color = "red", size = 1) +
  geom_hline(yintercept = (100-(median(df$Proportion[df$Var2 == "LnCAP"]))), color = "red", size = 1)
gg_cnumb
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######

# CORRELATION OF CELL NUMBERS
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######
psamp <- Polychrome::createPalette(10, c("#a6dba0", "#c2a5cf", "#f4a582"))
names(psamp) <- unique(cd$sample)

df <- as.data.frame(table(cd$sample)) %>%
      mutate(InputCells = c(100,225,500,1000,20000,20000))

correlation <- cor(df$InputCells, df$Freq, method = "pearson") # pearson
gg_corr <- ggplot(df, aes(x = InputCells, y = Freq, color = Var1)) +
  scale_color_manual(values = psamp) +
  geom_abline(slope = 1, color = "black", linetype = "dotted", alpha = 1) +
  geom_point(size = 5, alpha = 0.8) +
  labs(
    subtitle = paste("R²:", round(correlation,2)),
    x = "Input",
    y = "Recovered",
    color = "Sample"
  ) +
  scale_x_log10(labels = scientific_10) +
  scale_y_log10(labels = scientific_10) +
  theme_bw() +
  scale_y_log10(labels = scientific_10) +
  theme(panel.grid = element_blank()) 
gg_corr
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######

# GENE COUNTS & Features
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######
df <- cd
# plot function
create_boxplot <- function(df, var_name, title, ylab) {
  # Dynamically calculate the median for each sample and create labels
  labels <- setNames(
    lapply(unique(df$sample), function(sample) {
      median_value <- round(median(df[[var_name]][df$sample == sample]), 0)
      glue("{median_value}")
    }),
    unique(df$sample)
  )
  
  # Create the plot
  ggplot(df, aes(x = sample, y = !!sym(var_name), color = sample)) +
    geom_boxplot(aes(fill = sample), alpha = 0.1) + # Set alpha for fill only
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_color_manual(values = psamp, labels = labels) +
    scale_fill_manual(values = psamp) +
    labs(color = title) +
    labs(y = ylab, color = "Median") +
    guides(fill = "none") +
    scale_y_continuous(labels = scientific_10) 
}

gg_sum <- create_boxplot(df, "sum", "Sample", "nCount")
gg_feat <- create_boxplot(df, "detected", "Sample", "nFeature")

gg_sum
gg_feat
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######

# Save
dir <- glue("{proj_dir}/figures/fig2/rds")
dir.create(dir, showWarnings = F)
saveRDS(gg_cnumb, file = glue("{dir}/gg_cnumb.rds"))
saveRDS(gg_corr, file = glue("{dir}/gg_corr.rds"))
saveRDS(gg_sum, file = glue("{dir}/gg_sum.rds"))
saveRDS(gg_feat, file = glue("{dir}/gg_feat.rds"))