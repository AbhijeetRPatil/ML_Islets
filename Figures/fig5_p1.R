## Load required libraries
library("Seurat")
library("stringi")
library("data.table")
library("ggpubr")
library("cowplot")
library("dplyr")
library("tidyverse")
setwd("./plots/revision/Fig5/")
source("./progs/top_geneex_functions.R")
################################## functions ##############################################################
xtheme <- theme(plot.title = element_text(face = "bold" ,hjust = 0.5, size= 12)
                ,axis.text.y = element_text(face = "bold",angle = 0, size = 10, hjust=1)
                ,axis.title.y = element_text(face = "bold", size = rel(0.8))
                ,axis.text.x = element_text(face = "bold",angle = 0, size = 10)
                ,axis.title.x = element_text(face = "bold", size = rel(0.8))
                ,axis.ticks.x=element_blank(), strip.text = element_text(size=12))
## Read data
local <- readRDS("/mnt/alvand/abhijeet/aab/apr_WO-T2D/objs/panc_raw_meta_filtered_sct_umap_WO-T2D_WO_RB_MT_04182022.rds")

## T1D vs Control
cond = "AAB"
disease_order_combined <- c("Control", "AAB", "MC-AAB", "T1D")
# 

################################################## Beta cells #################################################################
## MC
cellbarcodes <- read.csv("./res/third/10_200cv/beta_cellbarcodes.csv")
cellX <- cellbarcodes$X

## subset the MC cells
local_cond_filtered <- local[,colnames(local) %in% cellX]

df1 <- data.frame(table(local_cond_filtered$cell_type))

## All AAB cells
md <- local_cond_filtered@meta.data %>% as.data.table
df <- md[, .N, by = c("hpap_id")]
colnames(df) <- c("HPAP_ID", "Frequency")
# Outside bars
lab <- df$`Frequency`

pdf("BarPlot_AAB_Beta_Cell_MC_Cts.pdf", height = 4, width = 3)
ggplot(data=df, aes(x=reorder(HPAP_ID, -Frequency), y=Frequency)) + 
  geom_bar(stat="identity", fill=rgb(0.8,0.1,0.1,0.6))+ scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_text(aes(label=lab), vjust=-0.3, size=3.5)+ xlab("AAB+ HPAP Donors") + ylab("# of Beta cells classified as T1D") +
  theme_minimal() + xtheme
dev.off()

################################################## Ductal cells #################################################################
## MC
cellbarcodes <- read.csv("./res/third/10_200cv/ductal_cellbarcodes.csv")
cellX <- cellbarcodes$X

## subset the MC cells
local_cond_filtered <- local[,colnames(local) %in% cellX]

df1 <- data.frame(table(local_cond_filtered$cell_type))

## All AAB cells
md <- local_cond_filtered@meta.data %>% as.data.table
df <- md[, .N, by = c("hpap_id")]
colnames(df) <- c("HPAP_ID", "Frequency")
# Outside bars
lab <- df$`Frequency`

pdf("BarPlot_AAB_Ductal_Cell_MC_Cts.pdf", height = 4, width = 3)
ggplot(data=df, aes(x=reorder(HPAP_ID, -Frequency), y=Frequency)) + 
  geom_bar(stat="identity", fill=rgb(0.8,0.1,0.1,0.6))+ scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_text(aes(label=lab), vjust=-0.3, size=3.5)+ xlab("AAB+ HPAP Donors") + ylab("# of Ductal cells classified as T1D") +
  theme_minimal() + xtheme
dev.off()

################################################## Acinar cells #################################################################
## MC
cellbarcodes <- read.csv("./res/third/10_200cv/acinar_cellbarcodes.csv")
cellX <- cellbarcodes$X

## subset the MC cells
local_cond_filtered <- local[,colnames(local) %in% cellX]

df1 <- data.frame(table(local_cond_filtered$cell_type))

## All AAB cells
md <- local_cond_filtered@meta.data %>% as.data.table
df <- md[, .N, by = c("hpap_id")]
colnames(df) <- c("HPAP_ID", "Frequency")
# Outside bars
lab <- df$`Frequency`

pdf("BarPlot_AAB_Acinar_Cell_MC_Cts.pdf", height = 4, width = 4)
ggplot(data=df, aes(x=reorder(HPAP_ID, -Frequency), y=Frequency)) + 
  geom_bar(stat="identity", fill=rgb(0.8,0.1,0.1,0.6))+ scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_text(aes(label=lab), vjust=-0.3, size=3.5)+ xlab("AAB+ HPAP Donors") + ylab("# of Acinar cells classified as T1D") +
  theme_minimal() + xtheme
dev.off()

################################################## Alpha cells #################################################################
## MC
cellbarcodes <- read.csv("./res/third/10_200cv/alpha_cellbarcodes.csv")
cellX <- cellbarcodes$X

## subset the MC cells
local_cond_filtered <- local[,colnames(local) %in% cellX]

df1 <- data.frame(table(local_cond_filtered$cell_type))

## All AAB cells
md <- local_cond_filtered@meta.data %>% as.data.table
df <- md[, .N, by = c("hpap_id")]
colnames(df) <- c("HPAP_ID", "Frequency")
# Outside bars
lab <- df$`Frequency`

pdf("BarPlot_AAB_Alpha_Cell_MC_Cts.pdf", height = 4, width = 4)
ggplot(data=df, aes(x=reorder(HPAP_ID, -Frequency), y=Frequency)) + 
  geom_bar(stat="identity", fill=rgb(0.8,0.1,0.1,0.6))+ scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_text(aes(label=lab), vjust=-0.3, size=3.5)+ xlab("AAB+ HPAP Donors") + ylab("# of Alpha cells classified as T1D ") +
  theme_minimal() + xtheme
dev.off()

################################################## Immune cells #################################################################
## MC
cellbarcodes <- read.csv("./res/third/10_200cv/immune_cellbarcodes.csv")
cellX <- cellbarcodes$X

## subset the MC cells
local_cond_filtered <- local[,colnames(local) %in% cellX]

df1 <- data.frame(table(local_cond_filtered$cell_type))

## All AAB cells
md <- local_cond_filtered@meta.data %>% as.data.table
df <- md[, .N, by = c("hpap_id")]
colnames(df) <- c("HPAP_ID", "Frequency")
# Outside bars
lab <- df$`Frequency`

pdf("BarPlot_AAB_Immune_Cell_MC_Cts.pdf", height = 4, width = 2)
ggplot(data=df, aes(x=reorder(HPAP_ID, -Frequency), y=Frequency)) + 
  geom_bar(stat="identity", fill=rgb(0.8,0.1,0.1,0.6))+ scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_text(aes(label=lab), vjust=-0.3, size=3.5)+ xlab("AAB+ HPAP Donors") + ylab("# of Immune cells classified as T1D ") +
  theme_minimal() + xtheme
dev.off()

################################################## Endothelial cells #################################################################
## MC
cellbarcodes <- read.csv("./res/third/10_200cv/endothelial_cellbarcodes.csv")
cellX <- cellbarcodes$X

## subset the MC cells
local_cond_filtered <- local[,colnames(local) %in% cellX]

df1 <- data.frame(table(local_cond_filtered$cell_type))

## All AAB cells
md <- local_cond_filtered@meta.data %>% as.data.table
df <- md[, .N, by = c("hpap_id")]
colnames(df) <- c("HPAP_ID", "Frequency")
# Outside bars
lab <- df$`Frequency`

pdf("BarPlot_AAB_Endothelial_Cell_MC_Cts.pdf", height = 4, width = 2)
ggplot(data=df, aes(x=reorder(HPAP_ID, -Frequency), y=Frequency)) + 
  geom_bar(stat="identity", fill=rgb(0.8,0.1,0.1,0.6))+ scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_text(aes(label=lab), vjust=-0.3, size=3.5)+ xlab("AAB+ HPAP Donors") + ylab("# of Endothelial cells classified as T1D ") +
  theme_minimal() + xtheme
dev.off()

################################################## Stellates cells #################################################################
## MC
cellbarcodes <- read.csv("./res/third/10_200cv/stellates_cellbarcodes.csv")
cellX <- cellbarcodes$X

## subset the MC cells
local_cond_filtered <- local[,colnames(local) %in% cellX]

df1 <- data.frame(table(local_cond_filtered$cell_type))

## All AAB cells
md <- local_cond_filtered@meta.data %>% as.data.table
df <- md[, .N, by = c("hpap_id")]
colnames(df) <- c("HPAP_ID", "Frequency")
# Outside bars
lab <- df$`Frequency`

pdf("BarPlot_AAB_Stellates_Cell_MC_Cts.pdf", height = 4, width = 3)
ggplot(data=df, aes(x=reorder(HPAP_ID, -Frequency), y=Frequency)) + 
  geom_bar(stat="identity", fill=rgb(0.8,0.1,0.1,0.6))+ scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_text(aes(label=lab), vjust=-0.3, size=3.5)+ xlab("AAB+ HPAP Donors") + ylab("# of Stellates cells classified as T1D ") +
  theme_minimal() + xtheme
dev.off()

################################################## Delta cells #################################################################
## MC
cellbarcodes <- read.csv("./res/third/10_200cv/delta_cellbarcodes.csv")
cellX <- cellbarcodes$X

## subset the MC cells
local_cond_filtered <- local[,colnames(local) %in% cellX]

df1 <- data.frame(table(local_cond_filtered$cell_type))

## All AAB cells
md <- local_cond_filtered@meta.data %>% as.data.table
df <- md[, .N, by = c("hpap_id")]
colnames(df) <- c("HPAP_ID", "Frequency")
# Outside bars
lab <- df$`Frequency`

pdf("BarPlot_AAB_Delta_Cell_MC_Cts.pdf", height = 4, width = 2)
ggplot(data=df, aes(x=reorder(HPAP_ID, -Frequency), y=Frequency)) + 
  geom_bar(stat="identity", fill=rgb(0.8,0.1,0.1,0.6))+ scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_text(aes(label=lab), vjust=-0.3, size=3.5)+ xlab("AAB+ HPAP Donors") + ylab("# of Delta cells classified as T1D ") +
  theme_minimal() + xtheme
dev.off()


