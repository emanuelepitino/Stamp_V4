wh <- function(w,h){
    options(repr.plot.width = w, repr.plot.height = h)
}
###############################################################################
# 0. Palette general

library(Polychrome)

palette_general <- function(){
  
  set.seed(1234) # for reproducibility
  color_palette <- Polychrome::createPalette(60, c("#fc6060", "#74f774", "#7c7cfc"))
  return(unname(color_palette))
}
###############################################################################

###############################################################################
# 1. Palette for seurat cell_types_pred

library(Polychrome)


palette_cell_types <- function(){
  
  set.seed(1234) # for reproducibility
  color_palette <- Polychrome::createPalette(200, c("#fc6060", "#74f774", "#7c7cfc"))
  names(color_palette) <- NULL
  
  palette <- c("NKT cells" = color_palette[1], 
               "Naive CD4+ T cells" = color_palette[2],
               "Cytotoxic CD8+ T cells" = color_palette[3],
               "Central Memory CD4+ T cells" = color_palette[4],
               "Tr1" = color_palette[5],
               "Th" = color_palette[6],
               "Proliferative CD4+ T cells" = color_palette[7],
               "Treg" = color_palette[8],
               "Itgb1+ Memory CD4+ T cells" = color_palette[9],
               "Naive CD8+ T cells" = color_palette[10],
               "Activated Effector Memory CD8+ T cells" = color_palette[11],
               "Central Memory CD8+ T cells" = color_palette[12],
               "IFN+ Effector CD4+ T cells" = color_palette[13],
               "Cytotoxic Effector Memory CD4+ T cells" = color_palette[14],
               "IFNgR+ Memory CD8+ T cells" = color_palette[15],
               "dn T cells" = color_palette[16],
               "IFN+ Effector Memory CD8+ T cells" = color_palette[17],
               "Proliferative CD8+ T cells" = color_palette[18],
               "MAIT cells" = color_palette[19],
               "Il6st+ Naive CD4+ T cells" = color_palette[20])
  return(palette)
}
###############################################################################


###############################################################################
# 1.2 Full palette for cell types, including non T cells

library(Polychrome)


palette_cell_types_full <- function(){
  
  set.seed(1234) # for reproducibility
  color_palette <- Polychrome::createPalette(200, c("#fc6060", "#74f774", "#7c7cfc"))
  names(color_palette) <- NULL
  
  palette <- c("NKT cells" = color_palette[1], 
               "Naive CD4+ T cells" = color_palette[2],
               "Cytotoxic CD8+ T cells" = color_palette[3],
               "Central Memory CD4+ T cells" = color_palette[4],
               "Tr1" = color_palette[5],
               "Th" = color_palette[6],
               "Treg" = color_palette[8],
               "Itgb1+ Memory CD4+ T cells" = color_palette[9],
               "Naive CD8+ T cells" = color_palette[10],
               "Activated Effector Memory CD8+ T cells" = color_palette[11],
               "Central Memory CD8+ T cells" = color_palette[12],
               "Cytotoxic Effector Memory CD4+ T cells" = color_palette[14],
               "dn T cells" = color_palette[16],
               "IFN+ Effector Memory CD8+ T cells" = color_palette[17],
               "MAIT cells" = color_palette[19],
               "Il6st+ Naive CD4+ T cells" = color_palette[20],
               "Classical Monocytes Il1b+" = color_palette[21],
               "Classical Monocytes" = color_palette[22],
               "CD4+ T cell-Platelet complexes" = color_palette[23],
               "ILC" = color_palette[24],
               "Switched CD27- Memory B cells" = color_palette[25],
               "Non-Classical Monocytes" = color_palette[26],
               "IFNgR+ Memory CD8+ B cells" = color_palette[27],
               "Classical Monocytes CD163+" = color_palette[28],
               "Naive B cells" = color_palette[29],
               "Monocytes Csfr3+" = color_palette[30],
               "CD16+ high RB NK cells" = color_palette[31],
               "CD56dim CD16+ CD57- NK cells" = color_palette[32],
               "Classical Monocyte-platelet complexes" = color_palette[33],
               "Transitional B cells" = color_palette[32],
               "CD56dim CD16+ CD57+ NK cells" = color_palette[33],
               "DC4" = color_palette[34],
               "Unswitched Memory B cells" = color_palette[35],
               "XISThi B cells" = color_palette[36],
               "Intermediate Monocytes" = color_palette[37],
               "CD19hiCD11c+ Memory B cells" = color_palette[38],
               "DC2 Cd1c-" = color_palette[39],
               "Progenitor cells" = color_palette[40],
               "CD8+ T cell-Platelet complexes" = color_palette[41],
               "CD56bright NK cells" = color_palette[42]
               
  )
  return(palette)
}
###############################################################################
###############################################################################
# 1. Palette for clonal homeostasis plot

palette_homeostasis <- function(){
  
  library(Polychrome)
  set.seed(1234) # for reproducibility
  color_palette <- Polychrome::createPalette(32, c("#fc6060", "#74f774", "#7c7cfc"))
  names(color_palette) <- NULL
  
  palette <- color_palette[10:15]
  
  palette <- c("Expanded (180 < X <= 300)" = "#B2182B",
               "Medium (5 < X <= 20)" = "#90EE90" ,
               "Large (20 < X <= 80)" = "#FD8D3C",
               "Large_expanded (80 < X <= 180)" = "#D6604D",
               "Small (1 < X <= 5)" = "grey")
  return(palette)
}
###############################################################################

###############################################################################
# 1.B Continuous palette

palette_continuous <- function(){
  library(Polychrome)
  set.seed(1234) # for reproducibility
  color_palette <- colorRampPalette(c("#322783","#fbcc1c"))
  return(color_palette)
}
###############################################################################

###############################################################################
# 4. Function to calculate diversity within cell subtypes in scRepertoire

subtype_diversity = function(combined,subtype){
  
  
  diversities <- c("Shannon", "Chao")
  therapies <- c("JAKi","antiIL6")
  
  
  
  # 1. WEEK 0
  
  #Filter for subtype and Week 0 
  subtyped_combined <- lapply(combined, function(df) subset(df, cell_types_pred == subtype & week == "wk0"))
  
  ### ANTI IL 6 ####
  
  # 1. Calculate diversity in R-antiIL6-W0 
  antiIL6_R_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "R"))
  diversity_antiIL6_R_W0 <- clonalDiversity(antiIL6_R_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # 2. Calculate diversity in NR-antiIL6-W0 
  antiIL6_NR_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "NR"))
  diversity_antiIL6_NR_W0 <- clonalDiversity(antiIL6_NR_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # A. Merge these two tables 
  diversity_antiIL6_W0 <- bind_rows(diversity_antiIL6_NR_W0,diversity_antiIL6_R_W0)
  
  # B. put in long format    
  diversity_long_antiIL6_W0 <- gather(diversity_antiIL6_W0,index, value, colnames(diversity_antiIL6_W0[,1:5]))
  diversity_long_antiIL6_W0 <- separate(diversity_long_antiIL6_W0,col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  
  #### JAK INHIBITORS ####
  
  # 3. Calculate diversity in R-JAKi-W0 
  JAKi_R_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "R"))
  diversity_JAKi_R_W0 <- clonalDiversity(JAKi_R_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # 4. Calculate diversity in NR-JAKi-W0 
  JAKi_NR_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "NR"))
  diversity_JAKi_NR_W0 <- clonalDiversity(JAKi_NR_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # A. Merge these two tables 
  diversity_JAKi_W0 <- bind_rows(diversity_JAKi_NR_W0,diversity_JAKi_R_W0)
  # B. put in long format    
  diversity_long_JAKi_W0 <- gather(diversity_JAKi_W0,index, value, colnames(diversity_JAKi_W0[,1:5]))
  diversity_long_JAKi_W0 <- separate(diversity_long_JAKi_W0,
                                     col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  
  # 1. WEEK 12
  
  #Filter for subtype and Week 12 
  subtyped_combined <- lapply(combined, function(df) subset(df, cell_types_pred == subtype & week == "wk12"))
  
  ### ANTI IL 6 ####
  
  # 5.Calculate diversity in R-antiIL6-W12 
  antiIL6_R_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "R"))
  diversity_antiIL6_R_W12 <- clonalDiversity(antiIL6_R_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # 6. Calculate diversity in NR-antiIL6-W12 
  antiIL6_NR_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "NR"))
  diversity_antiIL6_NR_W12 <- clonalDiversity(antiIL6_NR_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # A. Merge these two tables 
  diversity_antiIL6_W12 <- bind_rows(diversity_antiIL6_NR_W12,diversity_antiIL6_R_W12)
  
  # B. put in long format    
  diversity_long_antiIL6_W12 <- gather(diversity_antiIL6_W12,index, value, colnames(diversity_antiIL6_W12[,1:5]))
  diversity_long_antiIL6_W12 <- separate(diversity_long_antiIL6_W12,col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  #### JAK INHIBITORS ####
  
  # 7. Calculate diversity in R-JAKi-W12 
  JAKi_R_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "R"))
  diversity_JAKi_R_W12 <- clonalDiversity(JAKi_R_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # 8. Calculate diversity in NR-JAKi-W12 
  JAKi_NR_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "NR"))
  diversity_JAKi_NR_W12 <- clonalDiversity(JAKi_NR_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # A. Merge these two tables 
  diversity_JAKi_W12 <- bind_rows(diversity_JAKi_NR_W12,diversity_JAKi_R_W12)
  
  # B. put in long format    
  diversity_long_JAKi_W12 <- gather(diversity_JAKi_W12,index, value, colnames(diversity_JAKi_W12[,1:5]))
  diversity_long_JAKi_W12 <- separate(diversity_long_JAKi_W12,
                                      col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  
  
  
  
  # put everything together
  full_diversity <- bind_rows(diversity_long_antiIL6_W0,
                              diversity_long_antiIL6_W12,
                              diversity_long_JAKi_W0,
                              diversity_long_JAKi_W12)
  
  full_diversity <- full_diversity %>% filter(index %in% diversities) #take diversities i care about
  
  
  # Calculate p values
  statistics <- full_diversity %>%
    group_by(therapy, index, week) %>%
    wilcox_test(data =.,
                formula = value~response,
                paired = FALSE) %>%
    adjust_pvalue(method = "fdr")
  
  
  # add xy_position and y_position to correctly put the p values on the x and y axis
  statistics <- statistics %>% 
    add_x_position(x = "response", dodge = 0.8)%>%
    add_y_position(scales = "free")
  
  #adjust the number of significant digits
  statistics$p.adj <- round(statistics$p.adj, digits = 5)
  
  
  function_output <- list(full_diversity, statistics, subtype)
  names(function_output) <- c("full_diversity","statistics","subtype")
  
  return(function_output)
  
}
###############################################################################

# 4. Function to calculate diversity within cell lineages in scRepertoire

cell_lineages_diversity = function(combined,cell_lineage){
  
  
  diversities <- c("Shannon", "Chao")
  therapies <- c("JAKi","antiIL6")
  
  
  
  # 1. WEEK 0
  
  #Filter for subtype and Week 0 
  subtyped_combined <- lapply(combined, function(df) subset(df, cell_lineages_pred == cell_lineage & week == "wk0"))
  
  ### ANTI IL 6 ####
  
  # 1. Calculate diversity in R-antiIL6-W0 
  antiIL6_R_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "R"))
  diversity_antiIL6_R_W0 <- clonalDiversity(antiIL6_R_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # 2. Calculate diversity in NR-antiIL6-W0 
  antiIL6_NR_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "NR"))
  diversity_antiIL6_NR_W0 <- clonalDiversity(antiIL6_NR_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # A. Merge these two tables 
  diversity_antiIL6_W0 <- bind_rows(diversity_antiIL6_NR_W0,diversity_antiIL6_R_W0)
  
  # B. put in long format    
  diversity_long_antiIL6_W0 <- gather(diversity_antiIL6_W0,index, value, colnames(diversity_antiIL6_W0[,1:5]))
  diversity_long_antiIL6_W0 <- separate(diversity_long_antiIL6_W0,col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  
  #### JAK INHIBITORS ####
  
  # 3. Calculate diversity in R-JAKi-W0 
  JAKi_R_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "R"))
  diversity_JAKi_R_W0 <- clonalDiversity(JAKi_R_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # 4. Calculate diversity in NR-JAKi-W0 
  JAKi_NR_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "NR"))
  diversity_JAKi_NR_W0 <- clonalDiversity(JAKi_NR_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # A. Merge these two tables 
  diversity_JAKi_W0 <- bind_rows(diversity_JAKi_NR_W0,diversity_JAKi_R_W0)
  # B. put in long format    
  diversity_long_JAKi_W0 <- gather(diversity_JAKi_W0,index, value, colnames(diversity_JAKi_W0[,1:5]))
  diversity_long_JAKi_W0 <- separate(diversity_long_JAKi_W0,
                                     col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  
  # 1. WEEK 12
  
  #Filter for subtype and Week 12 
  subtyped_combined <- lapply(combined, function(df) subset(df, cell_lineages_pred == cell_lineage & week == "wk12"))
  
  ### ANTI IL 6 ####
  
  # 5.Calculate diversity in R-antiIL6-W12 
  antiIL6_R_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "R"))
  diversity_antiIL6_R_W12 <- clonalDiversity(antiIL6_R_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # 6. Calculate diversity in NR-antiIL6-W12 
  antiIL6_NR_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "NR"))
  diversity_antiIL6_NR_W12 <- clonalDiversity(antiIL6_NR_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # A. Merge these two tables 
  diversity_antiIL6_W12 <- bind_rows(diversity_antiIL6_NR_W12,diversity_antiIL6_R_W12)
  
  # B. put in long format    
  diversity_long_antiIL6_W12 <- gather(diversity_antiIL6_W12,index, value, colnames(diversity_antiIL6_W12[,1:5]))
  diversity_long_antiIL6_W12 <- separate(diversity_long_antiIL6_W12,col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  #### JAK INHIBITORS ####
  
  # 7. Calculate diversity in R-JAKi-W12 
  JAKi_R_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "R"))
  diversity_JAKi_R_W12 <- clonalDiversity(JAKi_R_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # 8. Calculate diversity in NR-JAKi-W12 
  JAKi_NR_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "NR"))
  diversity_JAKi_NR_W12 <- clonalDiversity(JAKi_NR_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # A. Merge these two tables 
  diversity_JAKi_W12 <- bind_rows(diversity_JAKi_NR_W12,diversity_JAKi_R_W12)
  
  # B. put in long format    
  diversity_long_JAKi_W12 <- gather(diversity_JAKi_W12,index, value, colnames(diversity_JAKi_W12[,1:5]))
  diversity_long_JAKi_W12 <- separate(diversity_long_JAKi_W12,
                                      col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  
  
  
  
  # put everything together
  full_diversity <- bind_rows(diversity_long_antiIL6_W0,
                              diversity_long_antiIL6_W12,
                              diversity_long_JAKi_W0,
                              diversity_long_JAKi_W12)
  
  full_diversity <- full_diversity %>% filter(index %in% diversities) #take diversities i care about
  
  
  function_output <- list(full_diversity, cell_lineage)
  names(function_output) <- c("full_diversity","cell_lineage")
  
  return(function_output)
  
}

###############################################################################
#function to calculate celltype specific diversity in Responders

subtype_R_diversity = function(combined,subtype){
  
  
  diversities <- c("Shannon", "Chao")
  therapies <- c("JAKi","antiIL6")
  
  #filter for Responders
  subtyped_combined <- lapply(combined, function(df) subset(df, response == "R"))
  
  # 1. WEEK 0
  
  #Filter for subtype and Week 0 
  subtyped_combined <- lapply(combined, function(df) subset(df, cell_types_pred == subtype & week == "wk0"))
  
  ### ANTI IL 6 ####
  
  # 1. Calculate diversity in R-antiIL6-W0 
  antiIL6_R_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "R"))
  diversity_antiIL6_R_W0 <- clonalDiversity(antiIL6_R_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  
  # B. put in long format    
  diversity_long_antiIL6_W0 <- gather(diversity_antiIL6_R_W0,index, value, colnames(diversity_antiIL6_R_W0[,1:5]))
  diversity_long_antiIL6_W0 <- separate(diversity_long_antiIL6_W0,col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  
  #### JAK INHIBITORS ####
  
  # 3. Calculate diversity in R-JAKi-W0 
  JAKi_R_W0 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "R"))
  diversity_JAKi_R_W0 <- clonalDiversity(JAKi_R_W0, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # B. put in long format    
  diversity_long_JAKi_W0 <- gather(diversity_JAKi_R_W0,index, value, colnames(diversity_JAKi_R_W0[,1:5]))
  diversity_long_JAKi_W0 <- separate(diversity_long_JAKi_W0,
                                     col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  
  # 1. WEEK 12
  
  #Filter for subtype and Week 12 
  subtyped_combined <- lapply(combined, function(df) subset(df, cell_types_pred == subtype & week == "wk12"))
  
  ### ANTI IL 6 ####
  
  # 5.Calculate diversity in R-antiIL6-W12 
  antiIL6_R_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "antiIL6" & response == "R"))
  diversity_antiIL6_R_W12 <- clonalDiversity(antiIL6_R_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  
  # B. put in long format    
  diversity_long_antiIL6_W12 <- gather(diversity_antiIL6_R_W12,index, value, colnames(diversity_antiIL6_R_W12[,1:5]))
  diversity_long_antiIL6_W12 <- separate(diversity_long_antiIL6_W12,col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  #### JAK INHIBITORS ####
  
  # 7. Calculate diversity in R-JAKi-W12 
  JAKi_R_W12 <- lapply(subtyped_combined, function(df) subset(df, therapy == "JAKi" & response == "R"))
  diversity_JAKi_R_W12 <- clonalDiversity(JAKi_R_W12, cloneCall = "strict", chain = "both", exportTable = TRUE)
  
  # B. put in long format    
  diversity_long_JAKi_W12 <- gather(diversity_JAKi_R_W12,index, value, colnames(diversity_JAKi_R_W12[,1:5]))
  diversity_long_JAKi_W12 <- separate(diversity_long_JAKi_W12,
                                      col = "Group", sep = "_", into = c("disease","therapy","response","patient","week","timepoint"))
  
  # put everything together
  full_diversity <- bind_rows(diversity_long_antiIL6_W0,
                              diversity_long_antiIL6_W12,
                              diversity_long_JAKi_W0,
                              diversity_long_JAKi_W12)
  
  full_diversity <- full_diversity %>% filter(index %in% diversities) #take diversities i care about
  
  
  # Calculate p values
  statistics <- full_diversity %>%
    group_by(therapy, index, response) %>%
    wilcox_test(data =.,
                formula = value~week,
                paired = FALSE) %>%
    adjust_pvalue(method = "fdr")
  
  
  # add xy_position and y_position to correctly put the p values on the x and y axis
  statistics <- statistics %>% 
    add_x_position(x = "week", dodge = 0.8)%>%
    add_y_position(scales = "free")
  
  #adjust the number of significant digits
  statistics$p.adj <- round(statistics$p.adj, digits = 5)
  
  
  function_output <- list(full_diversity, statistics, subtype)
  names(function_output) <- c("full_diversity","statistics","subtype")
  
  return(function_output)
  
}

###############################################################################
# 5. Function to plot cell type specific diversity 

plot_celltypediversity_function <- function(diversity_output){
  
  therapies <- c("JAKi","antiIL6")
  
  diversity <- diversity_output[[1]]
  statistics <- diversity_output[[2]]
  subtype <- diversity_output[[3]]
  
  plot <- list()
  for (t in therapies) {
    
    df <- diversity %>% filter(therapy == t)
    df <- df %>% filter(index == "Shannon")
    stat <- statistics %>% filter(therapy == t)
    
    
    plot[[t]] <- ggplot(df, aes(x = response, y = value, fill = week)) +
      facet_wrap(~therapy, scales = "free") +
      geom_boxplot(alpha = 0, outlier.alpha = 0, aes(color = week), lwd = 1, position = "dodge") +
      scale_color_manual(values = c("lightblue","orchid4")) +
      geom_jitter(alpha= 0.8, size = 1.5, position = position_jitterdodge()) +
      theme(axis.text = element_text(size = 12),
            axis.text.x=element_text(angle=-45, vjust = 1),
            axis.text.x.bottom = element_text(size = 30),
            text = element_text(size = 30))
  }
  
  plot_wrapped <- patchwork::wrap_plots(plot)
  return(plot_wrapped)
}

###############################################################################

###############################################################################
# 5. Function to plot cell type specific diversity in responders across timepoints

plot_celltypediversity_R_function <- function(diversity_output){
  
  therapies <- c("JAKi","antiIL6")
  
  diversity <- diversity_output[[1]]
  statistics <- diversity_output[[2]]
  subtype <- diversity_output[[3]]
  
  plot <- list()
  for (t in therapies) {
    
    df <- diversity %>% filter(therapy == t)
    stat <- statistics %>% filter(therapy == t)
    
    
    plot[[t]] <- ggplot(df, aes(x = week, y = value)) +
      facet_wrap(~therapy~index, scales = "free") +
      geom_poi(aes(color = week), position = position_jitterdodge(0.1)) +
      geom_boxplot(alpha = 0.5, aes(fill = week), position = "dodge") +
      scale_color_manual(values = c("darkred","darkorchid4")) +
      scale_fill_manual(values = c("darkred", "darkorchid4")) +
      theme_pubr() +
      theme(axis.text = element_text(size = 12),
            axis.text.x=element_text(angle=-45, vjust = 1)) +
      stat_pvalue_manual(stat,label = "p.adj") +
      labs(title = t, subtitle = subtype) +
      theme(text = element_text(size = 20))  +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
  }
  
  plot_wrapped <- patchwork::wrap_plots(plot)
  return(plot_wrapped)
}

###############################################################################
# 6. Function to plot Symphony results

plotBasic = function(umap_labels,                  # metadata, with UMAP labels in UMAP1 and UMAP2 slots
                     title = NULL,                 # Plot title
                     color.by = NULL,              # metadata column name for coloring
                     facet.by = NULL,              # (optional) metadata column name for faceting
                     color.mapping = NULL,         # custom color mapping
                     legend.position = 'right') {  # Show cell type legend
  
  p = umap_labels %>%
    dplyr::sample_frac(1L) %>% # permute rows randomly
    ggplot(aes(x = UMAP1, y = UMAP2)) + 
    geom_point_rast(aes(col = get(color.by)), size = 1, stroke = 0.4, shape = 16)
  if (!is.null(color.mapping)) { p = p + scale_color_manual(values = color.mapping) }
  
  # Default formatting
  p = p + theme_bw() +
    labs(title = title, color = color.by) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position=legend.position) +
    theme(legend.text = element_text(size=8), legend.title=element_text(size=12)) + 
    guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = 'none')
  
  if(!is.null(facet.by)) {
    p = p + facet_wrap(~get(facet.by)) +
      theme(strip.text.x = element_text(size = 12)) }    
  return(p)
}

###############################################################################


###############################################################################
# 9. Function to export table summary for quality control

table_summary_qc <- function(seurat_obj, title){
  df <- list()
  
  libraries <- unique(seurat_obj$gem_id)
  for(library in libraries){
    seurat_fun <- subset(seurat_obj, gem_id == library)
    
    df[[library]] <- data.frame(
      "Number of cells" = length(Cells(seurat_fun)),
      "median umi per cell" = median(seurat_fun$nCount_RNA),
      "median gene per cell" = median(seurat_fun$nFeature_RNA),
      "mitochondrial fraction" = mean(seurat_fun$mt_pct))
    
    df[[library]]$gem_id <- library
    
  }
  
  df <- bind_rows(df)
  df <- cbind(df[,5],df[,-5])
  df <- rename(df, "gem_id" = "df[, 5]")
  tab <- gt(df, groupname_col = group_vars(df)) 
  tab <- tab_header(tab,title = title)
  
  return(tab)
}

# 11. Function to plot homeostasis plot with scRepertoire function 

homeostasis_plot <- function(seurat_obj){
  
  seurat <- seurat_obj
  
  #subset by therapy & response
  seurat_R_antiIL6 <-  subset(seurat, response == "R" & therapy == "antiIL6")
  seurat_NR_antiIL6 <- subset(seurat, response == "NR" & therapy == "antiIL6")
  
  seurat_R_JAKi <-  subset(seurat, response == "R" & therapy == "JAKi")
  seurat_NR_JAKi <- subset(seurat, response == "NR" & therapy == "JAKi")
  
  # Put in a list
  homeostasis_list <- list(seurat_NR_antiIL6,seurat_R_antiIL6,seurat_NR_JAKi,seurat_R_JAKi)
  names(homeostasis_list) <- c("seurat_NR_antiIL6","seurat_R_antiIL6","seurat_NR_JAKi","seurat_R_JAKi")
  
  
  #plot
  homeostasis_plot_list <- list()
  homeostasis_tables <- list()
  for(name in names(homeostasis_list)){
    df <- homeostasis_list[[name]]
    p <- occupiedscRepertoire(df, x.axis = "ident", 
                              facet.by = "week",
                              label = FALSE) +
      labs(title = name, x = "Cell types", y = "single cells") +
      coord_flip() +
      theme(text = element_text(size = 15)) +
      scale_fill_manual(values = palette)
    
    
    homeostasis_plot_list[[name]] <- p
    
    t <- occupiedscRepertoire(df, x.axis = "ident", 
                              facet.by = "week",
                              label = FALSE, 
                              exportTable = TRUE)
    homeostasis_tables[[name]] <- t
  }
  
  # remove redundant legend in plot[1]
  homeostasis_plot_list[[1]] <- homeostasis_plot_list[[1]] +
    guides(fill = "none")
  #wrap plots
  antiIL6 <- wrap_plots(homeostasis_plot_list[1:2], ncol = 2) +
    plot_annotation(title = "Clonal Homeostasis") +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", text= element_text(size = 15)) 
  
  JAKi <- wrap_plots(homeostasis_plot_list[3:4], ncol = 2) +
    plot_annotation(title = "Clonal Homeostasis") &
    guides(fill = "none")
  
  
  # Plot
  list_plot <- list(antiIL6, JAKi)
  plot <- patchwork::wrap_plots(list_plot, ncol = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  output <- list()
  output[[1]] <- plot
  output[[2]] <- homeostasis_tables
  return(output)
}


##### FUNCTION TO PRE PROCESS A SUBSET OF DATA REMOVING THE TCR GENES FROM THE CLUSTERING ######

process_sub <- function(t_seu, resolution){

## Select TCR genes to exclude from downstream analysis

# For TCR Alpha Chain (TRA) Genes
TRAV <- grep("^TRAV", rownames(t_seu), value = TRUE)
TRAJ <- grep("^TRAJ", rownames(t_seu), value = TRUE)
TRAC <- grep("^TRAC", rownames(t_seu), value = TRUE)

# For TCR Beta Chain (TRB) Genes
TRBV <- grep("^TRBV", rownames(t_seu), value = TRUE)
TRBD <- grep("^TRBD", rownames(t_seu), value = TRUE)
TRBJ <- grep("^TRBJ", rownames(t_seu), value = TRUE)
TRBC <- grep("^TRBC", rownames(t_seu), value = TRUE)

# For TCR Gamma Chain (TRG) Genes
TRGV <- grep("^TRGV", rownames(t_seu), value = TRUE)
TRGJ <- grep("^TRGJ", rownames(t_seu), value = TRUE)
TRGC <- grep("^TRGC", rownames(t_seu), value = TRUE)

# For TCR Delta Chain (TRD) Genes
TRDV <- grep("^TRDV", rownames(t_seu), value = TRUE)
TRDD <- grep("^TRDD", rownames(t_seu), value = TRUE)
TRDJ <- grep("^TRDJ", rownames(t_seu), value = TRUE)
TRDC <- grep("^TRDC", rownames(t_seu), value = TRUE)

TCR_genes <- c(TRAV, TRAJ, TRAC, TRBV, TRBD, TRBJ, TRBC, TRGV, TRGJ, TRGC, TRDV, TRDD, TRDJ, TRDC)

# Normalize data and identify highly variable genes
t_seu <- t_seu %>%
  NormalizeData(normalization.method = "LogNormalize") %>%
  FindVariableFeatures(nfeatures = 2000)

## Filter the HVGs
hvg_names <- VariableFeatures(t_seu)
hvg_names <- setdiff(hvg_names, TCR_genes)
VariableFeatures(t_seu) <- hvg_names


t_seu <- t_seu %>%
  ScaleData() %>%
  RunPCA()

## Integrate

t_seu <- RunHarmony(t_seu, group.by.vars = "gem_id")
t_seu <- RunUMAP(t_seu, dims = 1:35, reduction = "harmony")

t_seu <- t_seu %>%
    FindNeighbors(dims = 1:35, reduction = "harmony") %>%
    FindClusters(resolution = resolution)


return(t_seu)
}
###############################################################################

###############################################################################
'%!in%' <- function(x, table) {
  !(x %in% table)
}
###############################################################################
					   
###############################################################################
generate_plots <- function(tmp_seurat, assay_name, feature_name) {
	# Ensure the DefaultAssay is set to the one you're interested in
	DefaultAssay(tmp_seurat) <- assay_name
	
	# Adjusted to directly use the provided feature_name for plotting
	plot_list <- list()
	plot_list[['Pre']] <- FeaturePlot(tmp_seurat[, tmp_seurat$timepoint == "Pre"], features = feature_name, split.by = "response",
			   min.cutoff = "q05", max.cutoff = "q95", keep.scale = "all") &
	scale_colour_viridis_c() &
	theme(legend.position = "right") &
	plot_annotation(title = "Pre",
	theme = theme(plot.title = element_text(size = 30)))
	
	plot_list[['Post']] <- FeaturePlot(tmp_seurat[, tmp_seurat$timepoint == "Post"], features = feature_name, split.by = "response",
				min.cutoff = "q05", max.cutoff = "q95", keep.scale = "all") &
	scale_colour_viridis_c() &
	theme(legend.position = "right") &
	plot_annotation(title = "Post",
	theme = theme(plot.title = element_text(size = 30)))
	
	return(plot_list)
}
###############################################################################

###############################################################################
# Adjusted function to generate unique identifiers ensuring all are unique
generate_unique_ids <- function(n, seed = 123) {
library(stringr)
  potential_ids <- expand.grid(letters, letters, sprintf("%02d", 0:99), stringsAsFactors = FALSE)
  potential_ids <- apply(potential_ids, 1, paste, collapse = "")
  potential_ids <- sample(potential_ids, n) # Randomly sample without replacement
  
  if (length(unique(potential_ids)) < n) {
    stop("Not enough unique IDs can be generated with the current method.")
  }
  
  return(potential_ids)
}
###############################################################################
					   
###############################################################################
addAssayDataToMetadata <- function(seu, assayName) {
  # Extract assay data
  assayData <- GetAssayData(seu, assay = assayName, slot = "data")
  
  # Convert to dataframe and transpose
  assayData <- as.data.frame(t(assayData))
  
  # Remove the now redundant barcode column
  assayData$barcode <- NULL
  
  # Add the modified assay data to the Seurat object's metadata
  seu <- AddMetaData(object = seu, metadata = assayData)
  
  return(seu)
}
###############################################################################

