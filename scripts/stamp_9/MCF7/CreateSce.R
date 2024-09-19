# Libraries
library(SingleCellExperiment)
library(here)
library(glue)
library(data.table)
library(qs)
# paths
dir <- glue("{here()}")
# Parameters and paths
source(glue("{dir}/scripts/misc/paths.R"))
source(glue("{dir}/scripts/misc/BIN.R"))

stamp <- "stamp_9"
sub <- "MCF7"
# data
sub_dir <- list.dirs(glue("{proj_dir}/data/{stamp}/{sub}"), full.names = FALSE, recursive = FALSE)[1]
data_dir <- glue("{proj_dir}/data/{stamp}/{sub}/{sub_dir}")
barcodes <- fread(glue("{data_dir}/cell_feature_matrix/barcodes.tsv.gz"), sep = "\t", nThread = 8, header = F)
features <- fread(glue("{data_dir}/cell_feature_matrix/features.tsv.gz"), sep = "\t", nThread = 8, header = F)
counts <- Matrix::readMM(glue("{data_dir}/cell_feature_matrix/matrix.mtx.gz"))
md <- fread(glue("{data_dir}/cells.csv.gz"), sep = ",", nThread = 8, header = T)
# coerce
counts <- as(counts, "dgCMatrix")
#add row and colnames
colnames(counts) <- barcodes$V1
rownames(counts) <- features$V2
# split matrix
gs <- rownames(counts)
ncp <- grep("NegControlProbe", gs)
ncc <- grep("NegControlCodeword", gs)
uc <- grep("UnassignedCodeword", gs)
dep <- grep("Deprecated",gs)
int <- grep("Intergenic",gs)
# alternative experiment
as <- list(counts=counts[-c(ncp, ncc,uc,dep,int), ])
ae <- list(
  NegControlProbe=SingleCellExperiment(list(counts=counts[ncp, ])),
  NegControlCodeword=SingleCellExperiment(list(counts=counts[ncc,])),
  DeprecatedCodeword = SingleCellExperiment(list(counts = counts[dep,])),
  IntergenicRegion =  SingleCellExperiment(list(counts = counts[int,])),
  UnassignedCodeword = SingleCellExperiment(list(counts=counts[uc,])))
# create sce
sce <- SingleCellExperiment(as, colData=md, altExps=ae)
# save
res_dir <- glue("{proj_dir}/data/{stamp}/{sub}/processed")
dir.create(res_dir, showWarnings = F)
qsave(sce, file = glue("{res_dir}/raw_sce.qs"))