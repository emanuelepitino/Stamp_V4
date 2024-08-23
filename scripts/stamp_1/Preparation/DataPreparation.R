# Data preparation for Stamp1

# Libraries
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(dplyr)
  library(patchwork)
  library(grid)
  library(ggpubr)
  library(here)
  library(scater)
  library(scuttle)
  library(glue)
  library(scran)
  library(patchwork)
  library(qs)
  library(data.table)
})

dir <- glue("{here()}/scripts")
# Parameters and paths
source(glue("{dir}/misc/paths.R"))
source(glue("{dir}/misc/BIN.R"))

# load data
myflatfiledir <- glue("{proj_dir}/data/stamp_1/raw/PBMCs")


# Process flatfiles
### automatically get slide names:
slidenames <- dir(myflatfiledir)

### lists to collect the counts matrices and metadata, one per slide
countlist <- vector(mode='list', length=length(slidenames)) 
metadatalist <- vector(mode='list', length=length(slidenames)) 

for(i in 1:length(slidenames)){
  
  slidename <- slidenames[i] 
  
  msg <- paste0("Loading slide ", slidename, ", ", i, "/", length(slidenames), ".")
  message(msg)    
  # slide-specific files:
  thisslidesfiles <- dir(paste0(myflatfiledir, "/", slidename))
  
  # load in metadata:
  thisslidesmetadata <- thisslidesfiles[grepl("metadata\\_file", thisslidesfiles)]
  tempdatatable <- data.table::fread(paste0(myflatfiledir, "/", slidename, "/", thisslidesmetadata))
  
  # numeric slide ID 
  slide_ID_numeric <- tempdatatable[1,]$slide_ID 
  
  # load in counts as a data table:
  thisslidescounts <- thisslidesfiles[grepl("exprMat\\_file", thisslidesfiles)]
  
  countsfile <- paste0(myflatfiledir, "/", slidename, "/", thisslidescounts)
  nonzero_elements_perchunk <- 5*10**7
  ### Safely read in the dense (0-filled ) counts matrices in chunks.
  ### Note: the file is gzip compressed, so we don't know a priori the number of chunks needed.
  lastchunk <- FALSE 
  skiprows <- 0
  chunkid <- 1
  
  required_cols <- fread(countsfile, select=c("fov", "cell_ID"))
  stopifnot("columns 'fov' and 'cell_ID' are required, but not found in the counts file" = 
              all(c("cell_ID", "fov") %in% colnames(required_cols)))
  number_of_cells <- nrow(required_cols)
  
  number_of_cols <-  ncol(fread(countsfile, nrows = 2))
  number_of_chunks <- ceiling(number_of_cols * number_of_cells / (nonzero_elements_perchunk))
  chunk_size <- floor(number_of_cells / number_of_chunks)
  sub_counts_matrix <- vector(mode='list', length=number_of_chunks)
  
  pb <- txtProgressBar(min = 0, max = number_of_chunks, initial = 0, char = "=",
                       width = NA, title, label, style = 3, file = "")
  cellcount <- 0
  while(lastchunk==FALSE){
    read_header <- FALSE
    if(chunkid==1){
      read_header <- TRUE
    }
    
    countsdatatable <- data.table::fread(countsfile,
                                         nrows=chunk_size,
                                         skip=skiprows + (chunkid > 1),
                                         header=read_header)
    if(chunkid == 1){
      header <- colnames(countsdatatable)
    } else {
      colnames(countsdatatable) <- header
    }
    
    cellcount <- nrow(countsdatatable) + cellcount     
    if(cellcount == number_of_cells) lastchunk <- TRUE
    skiprows <- skiprows + chunk_size
    slide_fov_cell_counts <- paste0("c_",slide_ID_numeric, "_", countsdatatable$fov, "_", countsdatatable$cell_ID)
    sub_counts_matrix[[chunkid]] <- as(countsdatatable[,-c("fov", "cell_ID"),with=FALSE], "sparseMatrix") 
    rownames(sub_counts_matrix[[chunkid]]) <- slide_fov_cell_counts 
    setTxtProgressBar(pb, chunkid)
    chunkid <- chunkid + 1
  }
  
  close(pb)   
  
  countlist[[i]] <- do.call(rbind, sub_counts_matrix) 
  
  # ensure that cell-order in counts matches cell-order in metadata   
  slide_fov_cell_metadata <- paste0("c_",slide_ID_numeric, "_", tempdatatable$fov, "_", tempdatatable$cell_ID)
  countlist[[i]] <- countlist[[i]][match(slide_fov_cell_metadata, rownames(countlist[[i]])),] 
  metadatalist[[i]] <- tempdatatable 
  
  # track common genes and common metadata columns across slides
  if(i==1){
    sharedgenes <- colnames(countlist[[i]]) 
    sharedcolumns <- colnames(tempdatatable)
  }  else {
    sharedgenes <- intersect(sharedgenes, colnames(countlist[[i]]))
    sharedcolumns <- intersect(sharedcolumns, colnames(tempdatatable))
  }
  
}

# reduce to shared metadata columns and shared genes
for(i in 1:length(slidenames)){
  metadatalist[[i]] <- metadatalist[[i]][, ..sharedcolumns]
  countlist[[i]] <- countlist[[i]][, sharedgenes]
}

counts <- do.call(rbind, countlist)
metadata <- rbindlist(metadatalist)

# add to metadata: add a global non-slide-specific FOV ID:
metadata$FOV <- paste0("s", metadata$slide_ID, "f", metadata$fov)

# remove cell_ID metadata column, which only identifies cell within slides, not across slides:
metadata$cell_ID <- NULL

# isolate negative control matrices:
negcounts <- counts[, grepl("Negative", colnames(counts))]
falsecounts <- counts[, grepl("SystemControl", colnames(counts))]

# reduce counts matrix to only genes:
counts <- counts[, !grepl("Negative", colnames(counts)) & !grepl("SystemControl", colnames(counts))]


# Put all together
atomxdata <- list(counts = counts,
                  metadata = metadata,
                  negcounts = negcounts,
                  falsecounts = falsecounts)
str(atomxdata)

outdir <- paste0(myflatfiledir,"/../raw_proc")
if(!dir.exists(paste0(outdir))){
  dir.create(outdir)
}

qsave(atomxdata$counts, file = paste0(outdir, "/counts_unfiltered.qs"))
qsave(atomxdata$negcounts, file = paste0(outdir, "/negcounts_unfiltered.qs"))
qsave(atomxdata$falsecounts, file = paste0(outdir, "/falsecounts_unfiltered.qs"))
qsave(atomxdata$metadata, file = paste0(outdir, "/metadata_unfiltered.qs"))