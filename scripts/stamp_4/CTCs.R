

# # dependencies
 suppressPackageStartupMessages({
     library(dplyr)
     library(Matrix)
     library(HDF5Array)
     library(SparseArray)
     library(SingleCellExperiment)
   library(glue)
   library(qs)
   library(here)
 })
# 
# # loading
 dir <- glue("{here()}")
 # Parameters and paths
 source(glue("{dir}/scripts/misc/paths.R"))
 source(glue("{dir}/scripts/misc/BIN.R"))
 dir <- glue("{proj_dir}/data/stamp_6/CTC_1/raw")
 
 sce <- qread(glue("{proj_dir}/data/stamp_4/processed/qc_sce.qs"), nthreads = 8)

# Square with 10 CTCs
sub <- sce[,sce$fov > 285]
gs <- c("TFF1","SPTSSB","AREG","MDK","PSMD6","ACKR3","EEF1A2","USP32","GATA3")

gs <- c("EEF1A2", "KRT8","FASN","HSPB1","TRIM37", "CCND1","PSMD6","NELFCD", "TFF1","MACROD1","LMNA",
     "KRT19","APPBP2", "CTTN","CLDN7","STARD10", "FAM83H", "SPINT2", "GATA3","USP32","S100A14","XBP1",
     "EPCAM","FKBP10", "SREBF1", "THOC7","CD9","NME4","MDK","PARD6B")


area <- sub$Area.um2
counts <- sub$sum

gs <- intersect(rownames(sce), gs)
es <- counts(sub[gs, ])
score <- colSums(es)

df <- data.frame(score,
                 panck=asinh(sub$Mean.PanCK/150), 
                 cd45=asinh(sub$Mean.CD45/150),
                 counts = counts,
                 area = area)

idx <- sub$cell_id[tail(order(sub$Mean.PanCK), 5e4)] # plot 50k highest PanCK cells

# CTCs ids
ctcs <- c("c_1_413_3812","c_1_296_456","c_1_505_1852","c_1_349_3459","c_1_457_1061",
          "c_1_457_1074","c_1_503_1397","c_1_332_21","c_1_315_4399","c_1_552_319","c_1_555_168")
# Artifacts ids
artifacts <- c("c_1_435_2447","c_1_289_62","c_1_529_667","c_1_464_396","c_1_464_385","c_1_559_135","c_1_384_2638")

fd <- df[idx, ]
fd <- fd[order(fd$score), ]
fd_cs <- fd[rownames(fd) %in% cs, ]

# Combine CTCs and artifacts with a new column for label type
fd_labels <- rbind(
  data.frame(fd[rownames(fd) %in% ctcs, ], label_type = 'ctc'),
  data.frame(fd[rownames(fd) %in% artifacts, ], label_type = 'artifact')
)
ctc_10 <- wrap_plots(
ggplot(fd, aes(panck, cd45, col=log1p(score))) + 
  geom_point(shape=16, size=1.5) + 
  scale_color_gradientn("MCF7\nscore", colors=c("gold", "red", "navy")) +
  theme_bw() + 
  theme(panel.grid.minor=element_blank()) +
 # labs(x = "Mean PanCK", y = "Mean CD45") +
  geom_label_repel(data=fd_labels, aes(label=rownames(fd_labels)),
                  color = ifelse(fd_labels$label_type == 'ctc', 'darkgreen', 'grey'),
                  direction = "both",
                  max.overlaps = Inf,
                  min.segment.length = 0),

ggplot(fd, aes(panck, log1p(counts), col=log1p(score))) + 
  geom_point(shape=16, size=1.5) + 
  scale_color_gradientn("MCF7\nscore", colors=c("gold", "red", "navy")) +
  theme_bw() + 
  theme(panel.grid.minor=element_blank()) +
 # labs(x = "Mean PanCK", y = "Counts") +
  geom_label_repel(data=fd_labels, aes(label=rownames(fd_labels)),
                   color = ifelse(fd_labels$label_type == 'ctc', 'darkgreen', 'grey'),
                   direction = "both",
                   max.overlaps = Inf,
                   min.segment.length = 0),

ggplot(fd, aes(panck, area, col=log1p(score))) + 
  geom_point(shape=16, size=1.5) + 
  scale_color_gradientn("MCF7\nscore", colors=c("gold", "red", "navy")) +
  theme_bw() + 
  theme(panel.grid.minor=element_blank()) +
  #labs(x = "Mean PanCK", y = "Counts") +
  geom_label_repel(data=fd_labels, aes(label=rownames(fd_labels)),
                   color = ifelse(fd_labels$label_type == 'ctc', 'darkgreen', 'grey'),
                   direction = "both",
                   max.overlaps = Inf,
                   min.segment.length = 0,),

ggplot(fd, aes(panck, (counts/area), col=log1p(score))) + 
  geom_point(shape=16, size=1.5) + 
  scale_color_gradientn("MCF7\nscore", colors=c("gold", "red", "navy")) +
  theme_bw() + 
  theme(panel.grid.minor=element_blank()) +
  #labs(x = "Mean PanCK", y = "Counts") +
  geom_label_repel(data=fd_labels, aes(label=rownames(fd_labels)),
                   color = ifelse(fd_labels$label_type == 'ctc', 'darkgreen', 'grey'),
                   direction = "both",
                   max.overlaps = Inf,
                   min.segment.length = 0,),
ncol = 2)


pdf(glue("{plt_dir}/stamp_4/ctc_10.pdf"), width = 20, height = 20)
ctc_10
dev.off()




# Square with 20 CTCs
sub <- sce[,sce$fov < 285]
#gs <- c("TFF1","SPTSSB","AREG","MDK","PSMD6","ACKR3","EEF1A2","USP32","GATA3")

gs <- c("EEF1A2", "KRT8","FASN","HSPB1","TRIM37", "CCND1","PSMD6","NELFCD", "TFF1","MACROD1","LMNA",
        "KRT19","APPBP2", "CTTN","CLDN7","STARD10", "FAM83H", "SPINT2", "GATA3","USP32","S100A14","XBP1",
        "EPCAM","FKBP10", "SREBF1", "THOC7","CD9","NME4","MDK","PARD6B")


area <- sub$Area.um2
counts <- sub$sum

gs <- intersect(rownames(sce), gs)
es <- counts(sub[gs, ])
score <- colSums(es)

df <- data.frame(score,
                 panck=asinh(sub$Mean.PanCK/150), 
                 cd45=asinh(sub$Mean.CD45/150),
                 counts = counts,
                 area = area)

idx <- sub$cell_id[tail(order(sub$Mean.PanCK), 5e4)] # plot 50k highest PanCK cells

# CTCs ids
ctcs <- c("c_1_413_3812","c_1_296_456","c_1_505_1852","c_1_349_3459","c_1_457_1061",
          "c_1_457_1074","c_1_503_1397","c_1_332_21","c_1_315_4399","c_1_552_319","c_1_555_168")
# Artifacts ids
artifacts <- c("c_1_435_2447","c_1_289_62","c_1_529_667","c_1_464_396","c_1_464_385","c_1_559_135","c_1_384_2638")

fd <- df[idx, ]
fd <- fd[order(fd$score), ]
fd_cs <- fd[rownames(fd) %in% cs, ]

# Combine CTCs and artifacts with a new column for label type
fd_labels <- rbind(
  data.frame(fd[rownames(fd) %in% ctcs, ], label_type = 'ctc'),
  data.frame(fd[rownames(fd) %in% artifacts, ], label_type = 'artifact')
)
ctc_10 <- wrap_plots(
  ggplot(fd, aes(panck, cd45, col=log1p(score))) + 
    geom_point(shape=16, size=1.5) + 
    scale_color_gradientn("MCF7\nscore", colors=c("gold", "red", "navy")) +
    theme_bw() + 
    theme(panel.grid.minor=element_blank()) +
    # labs(x = "Mean PanCK", y = "Mean CD45") +
    geom_label_repel(data=fd_labels, aes(label=rownames(fd_labels)),
                     color = ifelse(fd_labels$label_type == 'ctc', 'darkgreen', 'grey'),
                     direction = "both",
                     max.overlaps = Inf,
                     min.segment.length = 0),
  
  ggplot(fd, aes(panck, log1p(counts), col=log1p(score))) + 
    geom_point(shape=16, size=1.5) + 
    scale_color_gradientn("MCF7\nscore", colors=c("gold", "red", "navy")) +
    theme_bw() + 
    theme(panel.grid.minor=element_blank()) +
    # labs(x = "Mean PanCK", y = "Counts") +
    geom_label_repel(data=fd_labels, aes(label=rownames(fd_labels)),
                     color = ifelse(fd_labels$label_type == 'ctc', 'darkgreen', 'grey'),
                     direction = "both",
                     max.overlaps = Inf,
                     min.segment.length = 0),
  
  ggplot(fd, aes(panck, area, col=log1p(score))) + 
    geom_point(shape=16, size=1.5) + 
    scale_color_gradientn("MCF7\nscore", colors=c("gold", "red", "navy")) +
    theme_bw() + 
    theme(panel.grid.minor=element_blank()) +
    #labs(x = "Mean PanCK", y = "Counts") +
    geom_label_repel(data=fd_labels, aes(label=rownames(fd_labels)),
                     color = ifelse(fd_labels$label_type == 'ctc', 'darkgreen', 'grey'),
                     direction = "both",
                     max.overlaps = Inf,
                     min.segment.length = 0,),
  
  ggplot(fd, aes(panck, (counts/area), col=log1p(score))) + 
    geom_point(shape=16, size=1.5) + 
    scale_color_gradientn("MCF7\nscore", colors=c("gold", "red", "navy")) +
    theme_bw() + 
    theme(panel.grid.minor=element_blank()) +
    #labs(x = "Mean PanCK", y = "Counts") +
    geom_label_repel(data=fd_labels, aes(label=rownames(fd_labels)),
                     color = ifelse(fd_labels$label_type == 'ctc', 'darkgreen', 'grey'),
                     direction = "both",
                     max.overlaps = Inf,
                     min.segment.length = 0,),
  ncol = 2)


pdf(glue("{plt_dir}/stamp_4/ctc_20.pdf"), width = 20, height = 20)
ctc_20
dev.off()
