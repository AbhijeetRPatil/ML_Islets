## Load required libraries
library("Seurat")
library("stringi")
library("data.table")
library("ggpubr")
library("cowplot")
library("dplyr")
library("tidyverse")
#setwd("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/plots/MC/")
setwd("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/plots/revision/Fig5/")
source("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/progs/top_geneex_functions.R")
################################## functions ##############################################################

## Read data
local <- readRDS("/mnt/alvand/abhijeet/aab/apr_WO-T2D/objs/panc_raw_meta_filtered_sct_umap_WO-T2D_WO_RB_MT_04182022.rds")

## T1D vs Control
cond = "AAB"
disease_order_combined <- c("CTL", "AAb+", "MC-AAb+", "T1D")

## For T1DvsControl
local_cond <- subset(local, subset = disease_state != cond)

################################################## Beta cells #################################################################
## which cells
which_cell = "Beta"
## top 2
t1 <- '^HLA-A$'; t2 <- '^HLA-B$'; t3 <- '^HLA-C$'; t4 <- '^HLA-E$';
ylimit1 <- c(0,5); ylimit2 <- c(0,5); ylimit3 <- c(0,5); ylimit4 <- c(0,3);

##############################
## MC
cellbarcodes <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/beta_cellbarcodes.csv")
cellX <- cellbarcodes$X

## subset the MC cells
local_cond_filtered <- local[,colnames(local) %in% cellX]
##############################
## All combined
filename = "ALL_top4_MHCI_Beta.pdf"
p_combined <- expr_plot_4genes_combined(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
                            which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
pdf(filename, width = 15, height = 5)
p_combined
dev.off()
##############################################################################################################################

# ################################################## Immune cells #################################################################
# ## HLA-I
# ## which cells
# which_cell = "Immune"
# ## top 2
# t1 <- '^HLA-A$'; t2 <- '^HLA-B$'; t3 <- '^HLA-C$'; t4 <- '^HLA-E$';
# t1 <- '^HLA-A$'; t2 <- '^HLA-B$'; t3 <- '^HLA-E$'; t4 <- '^INS$';
# t1 <- '^HLA-A$'; t2 <- '^HLA-B$'; t3 <- '^HLA-E$'; t4 <- '^CD74$';
# t1 <- '^HLA-A$'; t2 <- '^HLA-B$'; t3 <- '^HLA-E$'; t4 <- '^CTSB$';
# # CD74, HLA-B, HLA-A, CTSB, HLA-E
# ylimit1 <- c(0,4); ylimit2 <- c(0,4); ylimit3 <- c(0,4); ylimit4 <- c(0,3);
# 
# ##############################
# ## MC
# cellbarcodes <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/immune_cellbarcodes.csv")
# cellX <- cellbarcodes$X
# 
# ## subset the MC cells
# local_cond_filtered <- local[,colnames(local) %in% cellX]
# ##############################
# ## All combined
# filename = "ALL_top4_AntigenProcesing_Immune_v4.pdf"
# p_combined <- expr_plot_4genes_combined(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
#                                         which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
# pdf(filename, width = 15, height = 5)
# p_combined
# dev.off()
# 
# ###############################
# ## HLA-I
# ## which cells
# which_cell = "Immune"
# ## top 2
# t1 <- '^HLA-DRA$'; t2 <- '^HLA-DPB1$'; 
# ylimit1 <- c(0,5); ylimit2 <- c(0,5); 
# 
# ##############################
# ## MC
# cellbarcodes <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/immune_cellbarcodes.csv")
# cellX <- cellbarcodes$X
# 
# ## subset the MC cells
# local_cond_filtered <- local[,colnames(local) %in% cellX]
# ##############################
# ## All combined
# filename = "ALL_top2_MHCII_Immune.pdf"
# p_combined <- expr_plot_2genes_combined(t1, t2, filename, local, local_cond_filtered, disease_order_combined,
#                                         which_cell, ylimit1, ylimit2)
# pdf(filename, width = 8, height = 5)
# p_combined
# dev.off()
# ##############################################################################################################################

################################################## Alpha cells #################################################################
## which cells
which_cell = "Alpha"
## top 2
t1 <- '^HLA-A$'; t2 <- '^HLA-B$'; t3 <- '^HLA-C$';  t4 <- '^HLA-E$';
ylimit1 <- c(0,6); ylimit2 <- c(0,6); ylimit3 <- c(0,5);  ylimit4 <- c(0,3); 

##############################
## MC
cellbarcodes <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/alpha_cellbarcodes.csv")
cellX <- cellbarcodes$X

## subset the MC cells
local_cond_filtered <- local[,colnames(local) %in% cellX]
##############################
## All combined
filename = "ALL_top4_MHCI_Alpha.pdf"
p_combined <- expr_plot_4genes_combined(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
                                                 which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
# p_combined <- expr_plot_4genes_combined_wolimits(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
#                                         which_cell)
pdf(filename, width = 15, height = 5)
p_combined
dev.off()
##############################################################################################################################

# ################################################## Endotheial cells #################################################################
# ## which cells
# which_cell = "Endothelial"
# ## top 2
# t1 <- '^HLA-A$'; t2 <- '^HLA-B$'; t3 <- '^HLA-C$'; 
# ylimit1 <- c(0,5); ylimit2 <- c(0,5); ylimit3 <- c(0,5); 
# ##############################
# ## MC
# cellbarcodes <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/endothelial_cellbarcodes.csv")
# cellX <- cellbarcodes$X
# 
# ## subset the MC cells
# local_cond_filtered <- local[,colnames(local) %in% cellX]
# ##############################
# ## All combined
# filename = "ALL_top4_MHCI_Endothelial.pdf"
# p_combined <- expr_plot_3genes_combined(t1, t2, t3, filename, local, local_cond_filtered, disease_order_combined,
#                                         which_cell, ylimit1, ylimit2, ylimit3)
# pdf(filename, width = 18, height = 5)
# p_combined
# dev.off()
# ########################################################################################################################################

# ################################################## Acinar cells #################################################################
# ## which cells
# which_cell = "Acinar"
# ## top 2
# t1 <- '^HLA-A$'; t2 <- '^HLA-B$'; t3 <- '^HLA-C$'; t4 <- '^HLA-E$'; t4 <- '^AZGP1$';
# ylimit1 <- c(0,4); ylimit2 <- c(0,4); ylimit3 <- c(0,4); ylimit4 <- c(0,4); ylimit5 <- c(0,2);
# 
# ##############################
# ## MC
# cellbarcodes <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/acinar_cellbarcodes.csv")
# cellX <- cellbarcodes$X
# 
# ## subset the MC cells
# local_cond_filtered <- local[,colnames(local) %in% cellX]
# ##############################
# ## All combined
# filename = "ALL_top4_MHCI_Acinar.pdf"
# p_combined <- expr_plot_5genes_combined(t1, t2, t3, t4, t5, filename, local, local_cond_filtered, disease_order_combined,
#                                         which_cell, ylimit1, ylimit2, ylimit3, ylimit4, ylimit5)
# pdf(filename, width = 18, height = 5)
# p_combined
# dev.off()
# 
# ##############################
# 
# ## HPAP092
# cellbarcodes <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/acinar_cellbarcodes.csv")
# cellX <- cellbarcodes$X
# cellX <- cellX[grepl("*-43$", cellX)]
# ## subset the MC cells
# local_cond_filtered <- local[,colnames(local) %in% cellX]
# ##############################
# ## All combined
# filename = "ALL_top4_MHCI_Acinar_HPAP092.pdf"
# p_combined <- expr_plot_5genes_combined(t1, t2, t3, t4, t5, filename, local, local_cond_filtered, disease_order_combined,
#                                         which_cell, ylimit1, ylimit2, ylimit3, ylimit4, ylimit5)
# pdf(filename, width = 18, height = 5)
# p_combined
# dev.off()
# 
# ##############################################################################################################################


# ################################################## Ductal cells #################################################################
# ## which cells
# which_cell = "Ductal"
# ## top 2
# t1 <- '^HLA-DRB1$'; t2 <- '^HLA-DQB1$'; 
# ylimit1 <- c(0,1); ylimit2 <- c(0,1); 
# 
# ##############################
# ## MC
# cellbarcodes <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/ductal_cellbarcodes.csv")
# cellX <- cellbarcodes$X
# 
# ## subset the MC cells
# local_cond_filtered <- local[,colnames(local) %in% cellX]
# ##############################
# ## All combined
# filename = "ALL_top2_MHCII_Ductal.pdf"
# p_combined <- expr_plot_2genes_combined(t1, t2, filename, local, local_cond_filtered, disease_order_combined,
#                                         which_cell, ylimit1, ylimit2)
# pdf(filename, width = 8, height = 5)
# p_combined
# dev.off()
# ##############################################################################################################################
# ################################################## Ductal cells (HLA-I) #################################################################
# ## which cells
# which_cell = "Ductal"
# ## top 2
# t1 <- '^HLA-A$'; t2 <- '^HLA-B$'; t3 <- '^HLA-C$'; t4 <- '^HLA-E$';
# t1 <- '^HLA-DRB5$'; t2 <- '^HLA-DQB1$'; t3 <- '^HLA-DRB1$'; t4 <- '^HLA-DRA$';
# t1 <- '^HLA-DRB5$'; t2 <- '^CELA3A$'; t3 <- '^PRSS1$'; t4 <- '^CELA3B$';
# 
# ylimit1 <- c(0,1); ylimit2 <- c(0,1); ylimit3 <- c(0,1); ylimit4 <- c(0,1);
# ylimit1 <- c(0,4); ylimit2 <- c(0,4); ylimit3 <- c(0,4); ylimit4 <- c(0,3);
# 
# ##############################
# ## MC
# cellbarcodes <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/ductal_cellbarcodes.csv")
# cellX <- cellbarcodes$X
# 
# ## subset the MC cells
# local_cond_filtered <- local[,colnames(local) %in% cellX]
# ##############################
# ## All combined
# filename = "ALL_top4_MHCI_Ductal.pdf"
# p_combined <- expr_plot_4genes_combined(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
#                                         which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
# pdf(filename, width = 15, height = 5)
# p_combined
# dev.off()
# #####################################################################################################################
################################################################################

## for 4 genes combined
expr_plot_4genes_combined <- function(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
                                      ylimit1, ylimit2, ylimit3, ylimit4, pal)
{
  local_cell <- local#subset(local, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE),
                grep(t3, rownames(local_cell), value = TRUE),
                grep(t4, rownames(local_cell), value = TRUE))
  local_cell_filtered <- local_cond_filtered
  
  p <- list()
  for (i in gene_var)
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df1 <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df1) <- c("Samples", "Type", "Gene")
    
    gene2<- FetchData(local_cell_filtered, vars = i)
    colnames(gene2) <- "gene"
    df2 <- data.frame(local_cell_filtered$disease_id, local_cell_filtered$disease_state, gene2$gene)
    colnames(df2) <- c("Samples", "Type", "Gene")
    df2$Samples <- str_replace(df2$Samples, "AAB_", "MC-AAB_")
    df2$Type <- str_replace(df2$Type, "AAB", "MC-AAB")
    df <- rbind(df1,df2)
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*', '^A.*', '^M.*', '^T.*'), 
                                                       replacement = disease_order_combined,
                                                       vectorize_all = FALSE), .after = "Samples") 
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",#color = 'Type', #
                       # palette = c("#00AFBB", "#E7B800", "#3CB371", "#FC4E07"), # "#FC4E07",
                       # palette = c("#00AFBB", "#E7B800", "#FC4E07", "#3CB371"), # "#FC4E07",                        #440154  482173
                       palette = pal,  
                       add=c("boxplot", "jitter"), add.params = list(fill="white"),
                       order = disease_order_combined,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #+ scale_x_discrete(labels=c("Control" = "Control", "AAB" = "MC-AAB",
    #"T1D" = "T1D")) + theme(legend.position = "none")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  p3 <- ggpar(p[[3]], ylim = ylimit3)
  p4 <- ggpar(p[[4]], ylim = ylimit4)
  
  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     p4 + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid(prow, ncol = 1, rel_heights = c(1, .2), legend_b) 
  return(p)
}
####################
## All combined
## MC
cellbarcodes1 <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/alpha_cellbarcodes.csv")
cellbarcodes2 <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/alpha_cellbarcodes.csv")
cellbarcodes3 <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/beta_cellbarcodes.csv")
cellbarcodes4 <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/ductal_cellbarcodes.csv")
cellbarcodes5 <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/delta_cellbarcodes.csv")
cellbarcodes6 <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/immune_cellbarcodes.csv")
cellbarcodes7 <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/endothelial_cellbarcodes.csv")
cellbarcodes8 <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/stellates_cellbarcodes.csv")

cellX <- c(cellbarcodes1$X, cellbarcodes2$X, cellbarcodes3$X, cellbarcodes4$X, 
           cellbarcodes5$X, cellbarcodes6$X, cellbarcodes7$X, cellbarcodes8$X)
##############################
t1 <- '^HLA-A$'; t2 <- '^HLA-B$'; t3 <- '^HLA-C$'; t4 <- '^HLA-E$';
ylimit1 <- c(0,4); ylimit2 <- c(0,4);  ylimit3 <- c(0,3); ylimit4 <- c(0,2.5); 

cellX1 <- cellX[grep(".+43$",cellX)] ## HPAP092
local_cond_filtered <- local[,colnames(local) %in% cellX1]
filename = "ALL_top4_MHCI_HPAP092_CELLTYPEAGGR_MC.pdf"
disease_order_combined <- c("CTL", "AAb+", "MC-AAb+\n(HPAP092)", "T1D")
p_combined <- expr_plot_4genes_combined(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
                                        ylimit1, ylimit2, ylimit3, ylimit4, pal = c('#21918c', '#440154', '#fde725'))
pdf(filename, width = 15, height = 5)
p_combined
dev.off()

cellX2 <- cellX[grep(".+50$",cellX)] ## HPAP107
local_cond_filtered <- local[,colnames(local) %in% cellX2]
filename = "ALL_top4_MHCI_HPAP107_CELLTYPEAGGR_MC.pdf"
disease_order_combined <- c("CTL", "AAb+", "MC-AAb+\n(HPAP107)", "T1D")
p_combined <- expr_plot_4genes_combined(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
                                        ylimit1, ylimit2, ylimit3, ylimit4, pal = c('#21918c', '#440154', '#fde725'))
pdf(filename, width = 15, height = 5)
p_combined
dev.off()

cellX_Comb <- c(cellX1, cellX2)
local_cond_filtered <- local[,colnames(local) %in% cellX_Comb]
filename = "ALL_top4_MHCI_HPAP092_HPAP107_CELLTYPEAGGR_MC.pdf"
disease_order_combined <- c("CTL", "AAb+", "MC-AAb+\n(HPAP092)\n(HPAP107)", "T1D")
p_combined <- expr_plot_4genes_combined(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
                                        ylimit1, ylimit2, ylimit3, ylimit4, pal = c('#21918c', '#440154', '#FC4E07', '#fde725'))
pdf(filename, width = 15, height = 5)
p_combined
dev.off()
