rm(list = ls())
#setwd("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/plots/HLAI/beta")
setwd("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/plots/revision/Fig6/")
## Load required libraries
library("Seurat")
library("stringi")
library("data.table")
library("ggpubr")
library("cowplot")
library("dplyr")
library("tidyverse")
library("ggplot2")
library("gridExtra")
## Read data
# local <- data.gen()
local <- readRDS("/mnt/alvand/abhijeet/aab/apr_WO-T2D/objs/panc_raw_meta_filtered_sct_umap_WO-T2D_WO_RB_MT_04182022.rds")

## For HPAP092
## for 4 genes combined
expr_plot_4genes_combined <- function(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
                                      ylimit1, ylimit2, ylimit3, ylimit4)
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
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       # palette = c("#00AFBB", "#E7B800", "#FC4E07"), # "#FC4E07", "#3CB371",
                       palette = c('#21918c', '#440154', '#fde725'),
                       # palette = c("#00AFBB", "#E7B800", "#FC4E07", "#3CB371"), # "#FC4E07",
                       
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
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
## MC
cellbarcodes <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/Overall_cellbarcodes.csv")
cellX <- cellbarcodes$X

cellX1 <- cellX[grep(".+43$",cellX)] ## HPAP092
# cellX2 <- cellX[grep(".+50$",cellX)] ## HPAP107
cellX <- c(cellX1)
## subset the MC cells
local_cond_filtered <- local[,colnames(local) %in% cellX]


##############################
## All combined
disease_order_combined <- c("CTL", "AAb+", "MC-AAb+\n(HPAP092)", "T1D")

t1 <- '^HLA-A$'; t2 <- '^HLA-B$'; t3 <- '^HLA-C$'; t4 <- '^HLA-E$'; 
ylimit1 <- c(0,5); ylimit2 <- c(0,5); ylimit3 <- c(0,5); ylimit4 <- c(0,5); 
filename = "HPAP092_MC_Allcells_d.pdf"
p_combined <- expr_plot_4genes_combined(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
                                        ylimit1, ylimit2, ylimit3, ylimit4)
pdf(filename, width = 15, height = 4)
p_combined
dev.off()

t1 <- '^HLA-DMA$'; t2 <- '^HLA-DQB1$'; t3 <- '^HLA-DRA$'; t4 <- '^HLA-DPA1$'; 
ylimit1 <- c(0,0.5); ylimit2 <- c(0,0.5); ylimit3 <- c(0,0.5); ylimit4 <- c(0,0.5); 
filename = "HPAP092_MC_Allcells_e.pdf"
p_combined <- expr_plot_4genes_combined(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
                                        ylimit1, ylimit2, ylimit3, ylimit4)
pdf(filename, width = 15, height = 4)
p_combined
dev.off()

t1 <- '^HLA-DRB1$'; t2 <- '^HLA-DQA2$'; t3 <- '^HLA-DPB1$'; t4 <- '^HLA-DRB5$';
ylimit1 <- c(0,0.5); ylimit2 <- c(0,0.5); ylimit3 <- c(0,0.5); ylimit4 <- c(0,0.5); 
filename = "HPAP092_MC_Allcells_f.pdf"
p_combined <- expr_plot_4genes_combined(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
                                        ylimit1, ylimit2, ylimit3, ylimit4)
pdf(filename, width = 15, height = 4)
p_combined
dev.off()

t1 <- '^INS$'; t2 <- '^IL32$'; t3 <- '^TNFAIP3$'; t4 <- '^LMO7$'; 
ylimit1 <- c(0,10); ylimit2 <- c(0,5); ylimit3 <- c(0,1); ylimit4 <- c(0,1); 
filename = "HPAP092_MC_Allcells_g.pdf"
p_combined <- expr_plot_4genes_combined(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
                                        ylimit1, ylimit2, ylimit3, ylimit4)
pdf(filename, width = 15, height = 4)
p_combined
dev.off()
