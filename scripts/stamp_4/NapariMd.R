# Create metadata file with CTC classification for Napari

# dependencies
library(here)
library(glue)
library(qs)

source(glue("{here()}/scripts/misc/paths.R"))
source(glue("{here()}/scripts/misc/BIN.R"))

dir <- glue("{proj_dir}/data/stamp_4/processed/")
cd1 <- qread(glue("{dir}/CTC1/ctc1_cd.qs"), nthreads = 8)
cd2 <- qread(glue("{dir}/CTC2/ctc2_cd.qs"), nthreads = 8)
cd <- rbind(cd1,cd2)

sce <- qread(glue("{proj_dir}/data/stamp_4/raw/raw_sce.qs"), nthreads = 8)
sce <- as.data.frame(colData(sce)) %>%
  select(cell_id)

sce$ctc <- cd$lab[match(sce$cell_id, cd$cell_id)]
sce$ctc[is.na(sce$ctc)] <- "none"

cd <- sce
rownames(cd) <- NULL # remove rownames
cd$cell_ID <- cd$cell_id
cd <- cd %>% relocate(cell_ID) # put cell_ID as first column
cd <- cd %>% select(cell_ID,ctc)
write.table(cd, file=glue("{proj_dir}/data/stamp_4/_metadata.csv"), 
            sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE) # save
