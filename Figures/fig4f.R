rm(list = ls())
#setwd("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/plots/HLAI/beta")
setwd('/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/plots/revision/Fig4/')
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

## T1D vs Control
cond = "AAB"
disease_order <- c("CTL", "AAB+", "T1D")
## For T1DvsControl
local_cond <- subset(local, subset = disease_state != cond)
# local_cond$disease_state <- droplevels(local_cond$disease_state)

##########################################################################################

## for 4 genes
expr_plot <- function(t1, t2, t3, t4, filename, which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
{
  local_cell <- subset(local_cond, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE),
                grep(t3, rownames(local_cell), value = TRUE),
                grep(t4, rownames(local_cell), value = TRUE))
  
  p <- list()
  for (i in gene_var) 
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df) <- c("Samples", "Type", "Gene")
    # table(df$local_beta.sample_id)
    
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*','^A.*', '^T.*'), 
                                                       replacement = disease_order,
                                                       vectorize_all = FALSE), .after = "Samples") 
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       #palette = c("#00AFBB", "#E7B800"), # "#FC4E07",
                       palette = c('#21918c', '#fde725'),
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #, xlab = "Groups")
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
                     # labels = c("A", "B", "C"),
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
  return(p)
}


which_cell = "Beta"

## for 4 genes
## Plot 1
filename = "T1DvsControl_top_HLAI_beta_v1.pdf"
t1 <- '^HLA-A$'; t2 <- '^HLA-B$'; t3 <- '^HLA-C$'; t4 <- '^HLA-E$';
ylimit1 <- c(0,5); ylimit2 <- c(0,4); ylimit3 <- c(0,4); ylimit4 <- c(0,2);
plot1 <- expr_plot(t1, t2, t3, t4, filename, which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
pdf(filename, width = 15, height = 5)
plot1
dev.off()

########################################################################################
## Below can be deleted

local_aab <- subset(local, subset = disease_state == cond)
df1 <- data.frame(table(local_aab$cell_type))


## All AAB cells
md <- local_aab@meta.data %>% as.data.table
df <- md[, .N, by = c("hpap_id")]
colnames(df) <- c("HPAP_ID", "Frequency")
# Outside bars
lab <- df$`Frequency`



# ## All AAB cells
# md <- local@meta.data %>% as.data.table
# df <- md[, .N, by = c("hpap_id")]
# colnames(df) <- c("HPAP_ID", "Frequency")
# # Outside bars
# lab <- df$`Frequency`

xtheme <- theme(plot.title = element_text(face = "bold" ,hjust = 0.5, size= 12)
                ,axis.text.y = element_text(face = "bold",angle = 0, size = 10, hjust=1)
                ,axis.title.y = element_text(face = "bold", size = rel(0.8))
                ,axis.text.x = element_text(face = "bold",angle = 30, size = 10, hjust = 1)
                ,axis.title.x = element_text(face = "bold", size = rel(0.8))
                ,axis.ticks.x=element_blank(), strip.text = element_text(size=12))

pdf("BarPlot_All_AAB_Cell_Counts_v1.pdf", height = 5, width = 5)
ggplot(data=df, aes(x=reorder(HPAP_ID, -Frequency), y=Frequency)) + 
  geom_bar(stat="identity", fill=rgb(0.8,0.1,0.1,0.6))+
  geom_text(aes(label=lab), vjust=-0.3, size=3.5)+ xlab("HPAP_ID") + ylab("Total number of cells in AAB donors") +
  theme_minimal() + xtheme
dev.off()

######################

md <- local_aab@meta.data %>% as.data.table
df <- md[, .N, by = c("hpap_id", "cell_type")]
df <- df[df$cell_type =="Beta",]
colnames(df) <- c("HPAP_ID", "Cell Type", "Frequency")
# Outside bars
lab <- df$`Frequency`

pdf("BarPlot_Beta_AAB_Cell_Counts_v1.pdf", height = 5, width = 5)
ggplot(data=df, aes(x=reorder(HPAP_ID, -Frequency), y=Frequency)) + 
  geom_bar(stat="identity", fill=rgb(0.8,0.1,0.1,0.6))+
  geom_text(aes(label=lab), vjust=-0.3, size=3.5)+ xlab("HPAP_ID") + ylab("Total number of Beta cells in AAB donors") +
  theme_minimal() + xtheme
dev.off()

##############################################################################################
## Percent AAB
local_aab <- subset(local, subset = disease_state == cond)
df1 <- data.frame(table(local_aab$cell_type))
## All AAB cells
md <- local_aab@meta.data %>% as.data.table
df1 <- md[, .N, by = c("hpap_id")]
colnames(df1) <- c("HPAP_ID", "Frequency")


md <- local_aab@meta.data %>% as.data.table
df <- md[, .N, by = c("hpap_id", "cell_type")]
df <- df[df$cell_type =="Beta",]
colnames(df) <- c("HPAP_ID", "Cell Type", "Frequency")

df$`Frequency` <- round(df$Frequency/df1$Frequency*100,2)

# Outside bars
lab <- df$`Frequency`


pdf("BarPlot_Beta_AAB_Cell_Perc_v1.pdf", height = 6, width = 6)
ggplot(data=df, aes(x=reorder(HPAP_ID, -Frequency), y=Frequency)) + 
  geom_bar(stat="identity", fill=rgb(0.8,0.1,0.1,0.6))+
  geom_text(aes(label=paste(lab,"%")), vjust=-0.3, size=3.5)+ xlab("HPAP ID") + ylab("Percentage of Beta Cells in AAB+ Donors") +
  theme_minimal() + xtheme
dev.off()

