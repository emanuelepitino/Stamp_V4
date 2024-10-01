### Libraries
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dplyr)
  library(here)
  library(scater)
  library(scuttle)
  library(scran)
  library(glue)
  library(qs)
  library(scales)
  library(bluster)
})

### Paths
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))
sample <- "combined"
stamp <- "Stamp_7b"
### Load data

# Flex
base_dir <- glue("{proj_dir}/data/{stamp}/processed/flex/{sample}")
sce <- qread(glue("{base_dir}/integrated_sce.qs"), nthreads = 8)
sce

# CosMx
cosmx <- readRDS(glue("{proj_dir}/data/{stamp}/processed/irepan_data/stamp_palantir_subset.rds"))

mtx <- GetAssayData(cosmx, layer = "counts")
cd <- cosmx@meta.data

sce_cosmx <- SingleCellExperiment(
  assays = list(counts = mtx),
  colData = cd
)

rm(mtx)
rm(cosmx)
gc()

sce_cosmx$tech <- "cosmx"
sce$tech <- "flex"


# subset for intersection with cosmx panel
fs <- rownames(sce_cosmx)
fs <- intersect(fs, rownames(sce))
sce <- sce[fs,]
sce_cosmx <- sce_cosmx[fs,]

# take mtx 
mtx_fl <- counts(sce)
mtx_cs <- counts(sce_cosmx)
mtx <- cbind(mtx_fl,mtx_cs)

# take colData
cd_fl <- as.data.frame(colData(sce)) %>% select(sample,tech)
cd_cs <- as.data.frame(colData(sce_cosmx))  %>% select(sample,tech)

cd <- rbind(cd_fl,cd_cs)

# merge sce
sce <- SingleCellExperiment(
  assays = list(counts = mtx),
  colData = cd
)

# lognormalize
sce <- logNormCounts(sce)

# PCA
################ ################ ################ ################ ################ 
set.seed(101001)
sce <- fixedPCA(sce)

num_pcs_to_retain <- 20
percent.var <- attr(reducedDim(sce), "percentVar")

# Create a data frame for ggplot
data <- data.frame(PC = 1:length(percent.var), Variance = percent.var)
# Plot
gg_var <- ggplot(data, aes(x = PC, y = Variance)) +
  geom_point() +
  xlab("PC") +
  ylab("Variance explained (%)") +
  geom_vline(xintercept = num_pcs_to_retain, color = "red") +
  theme_bw()
gg_var

reducedDim(sce, "PCA") <-  reducedDim(sce, "PCA")[,1:num_pcs_to_retain]
wh(6,5)
gg_pca <- plotPCA(sce, scattermore = TRUE, point_size = 2, color_by = "tech") + ggtitle("PCA")
gg_pca
################ ################ ################ ################ ################ 

sce <- harmony::RunHarmony(
  sce,
  group.by.vars = "tech",
  dims.use = 1:20,
  verbose = TRUE,
  reduction.save = "HARMONY",
  ncores = 8
)

################ ################ ################ ################ ################ 

df <- as.data.frame(colData(sce)) # cd
pca <- reducedDim(sce,"PCA") # pc
harmony <- reducedDim(sce,"HARMONY") # pc

df <- cbind(df,pca) # merge
df <- df[sample(nrow(df)), ] # shuffle

pal <- Polychrome::createPalette(n_distinct(sce$tech), c("#FBB4AE", "#B3CDE3", "#CCEBC5"))
names(pal) <- unique(sce$tech)

df$sample <- factor(df$sample, levels = c("iESC_0h","iESC_6h","iESC_12h","iESC_24h","iESC_48h","iESC_72h","iESC_96h","iESC_120h"))

unint <- ggplot(df, aes(x = PC1, y = PC2, color = tech)) +
  geom_point(shape = 16, size = 0.5) +
  theme_bw() +
  scale_color_manual(values = pal) +
  theme(
    text = element_text(size = 25, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 18),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()) +
  coord_equal()

################ ################ ################ ################ ################ 

df <- as.data.frame(colData(sce)) # cd
harmony <- reducedDim(sce,"HARMONY") # pc

df <- cbind(df,harmony) # merge
df <- df[sample(nrow(df)), ] # shuffle

df$tp_misl <- df$sample
df$tp_misl[df$sample == "iESC_120h" & df$tech == "cosmx"] <- "iESC_48h"
df$tp_misl[df$sample == "iESC_96h" & df$tech == "cosmx"] <- "iESC_72h"
df$tp_misl[df$sample == "iESC_72h" & df$tech == "cosmx"] <- "iESC_96h"
df$tp_misl[df$sample == "iESC_48h" & df$tech == "cosmx"] <- "iESC_120h"

  
  
  
pal <- Polychrome::createPalette(n_distinct(sce$tech), c("#FBB4AE", "#B3CDE3", "#CCEBC5"))
names(pal) <- unique(sce$tech)

df$tp_misl <- factor(df$tp_misl, levels = c("iESC_0h","iESC_6h","iESC_12h","iESC_24h","iESC_48h","iESC_72h","iESC_96h","iESC_120h"))

df$label <- as.character(df$label)
int <- ggplot(df, aes(x = HARMONY_1, y = HARMONY_2, color = label)) +
  geom_point(shape = 16, size = 0.5) +
  theme_bw() +
 # scale_color_manual(values = pal) +
  theme(
    text = element_text(size = 25, color = "black"),
    axis.text = element_text(size = 20, color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 18),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()) +
  coord_equal() +
  facet_wrap(~tp_misl) +
  guides(color = guide_legend(override.aes = list(size = 6)))  


wrap_plots(unint,int, ncol = 2) 


# Cluster
################ ################ ################ ################ ################ 
## Build snn graph
snn <- buildSNNGraph(sce, type = "jaccard", use.dimred="HARMONY", BPPARAM = bp)
## Run louvain 
clusters <- igraph::cluster_louvain(snn, resolution = 0.2)

sce$label <- clusters$membership
df <- as.data.frame(colData(sce))
df$sample <- factor(df$sample,
                    levels = c("iESC_0h","iESC_6h","iESC_12h","iESC_24h","iESC_48h","iESC_72h","iESC_96h","iESC_120h"))

df$tp_misl <- df$sample
df$tp_misl[df$sample == "iESC_120h" & df$tech == "cosmx"] <- "iESC_48h"
df$tp_misl[df$sample == "iESC_96h" & df$tech == "cosmx"] <- "iESC_72h"
df$tp_misl[df$sample == "iESC_72h" & df$tech == "cosmx"] <- "iESC_96h"
df$tp_misl[df$sample == "iESC_48h" & df$tech == "cosmx"] <- "iESC_120h"


fs <- c("SOX2","PU5F1","NANOG","GATA3","HEY1","EOMES","KDR","SNAI2","PDGFRA",
        "DUSP6","FOXF1","SNAI1","CXCR4","APOA1","TGBI")

plotDots(sce, features = intersect(fs, rownames(sce)), group = "label", scale = T, center = T) + 
  coord_flip() 


df <- as.data.frame(table(df$label, df$tech,df$sample))

# Calculate the percentage and reorder by 'singlet' percentage
df <- df %>%
  group_by(Var1) %>%
  mutate(pct = round(Freq / sum(Freq), 2)) %>%
  ungroup()


# Create the plot
ggplot(df, aes(x = Var3, y = pct, fill = Var1)) +
  geom_col(position = "stack") + 
  labs(x = "Cluster", y = "Proportion", title = "Clusters proportion by tech") + 
  theme_bw() +
  theme(text = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 10, color = "black")) + 
  guides(color = guide_legend(override.aes = list(size = 3))) + 
  facet_wrap(~Var2)







