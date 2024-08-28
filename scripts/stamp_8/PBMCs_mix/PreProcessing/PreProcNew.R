suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dplyr)
  library(here)
  library(glue)
  library(qs)
  library(HDF5Array)
  library(BiocParallel)
  library(BiocNeighbors)
  library(BiocSingular)
  library(scater)
})

stamp <- "stamp_8"
sample <- "PBMCs_mix"
dir <- glue("{here()}")
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

res_dir <- glue("{proj_dir}/data/{stamp}/{sample}")
sce <- qread(glue("{res_dir}/qc_sce.qs"), nthreads = 8)
#sce <- HDF5Array::loadHDF5SummarizedExperiment(dir = glue("{res_dir}/qc_sce"))
sce



# logNorm
sce <- logNormCounts(sce)
head(# normalize by total counts 
norm <- t(counts(sce)) / Matrix::rowSums(t(counts(sce)))

# log1p transformation
lognorm <- log1p(norm)

assay(sce, "logcounts") <- t(lognorm)
#assay(sce, "normcounts") <- t(norm)

rm(counts)
rm(lognorm)

# PC analysis
set.seed(101000)
sce <- runPCA(sce, ncomponents=20, BSPARAM=RandomParam())

pca <- reducedDim(tst,"PCA")
var <- as.data.frame(attr(pca,"percentVar"))
var$PC <- 1:(nrow(var))
# 4. Plot the elbow plot 
elbow <- ggplot(var, aes(x =  PC , y = attr(pca, "percentVar"))) +
  geom_point() +
  geom_line() +
  xlab("Principal Component") +
  ylab("Proportion of Variance Explained") +
  ggtitle(glue("{stamp} - {sample}")) +
  geom_vline(xintercept = 6, linetype = "dashed", color = "red") + 
  theme_minimal()

# Keep a subset of principal components
pc1$x <- pc1$x[,1:6]
# Run UMAP
um1 <- uwot::umap(pc1$x, n_neighbors = 40, spread = 1, min_dist = 0.1, metric = "cosine", n_threads = 8)
rownames(um1) <- rownames(norm)

par(mar = c(0,0,0,0))
plot(um1, pch = 16, cex = 0.1, col = "dodgerblue4")

reducedDim(sce,"PCA") <- pc1$x
reducedDim(sce,"UMAP") <- um1

sce$sample <- factor(sce$sample, levels = c("iESC_0h","iESC_6h","iESC_12h","iESC_24h",
                                            "iESC_48h","iESC_72h","iESC_96h","iESC_120h"))
# Save plots
plots_dir <- glue("{plt_dir}/{stamp}/{sample}")
pdf(glue("{plots_dir}/PreProc.pdf"), width = 7, height = 5)
elbow
plot(um1, pch = 16, cex = 0.1, col = "dodgerblue4")
dev.off()

# Save data
qsave(sce, file = glue("{res_dir}/PreProcNew.qs"), nthreads = 8)