rm(list=ls())
# setwd('/revision/Fig1/')
library('Seurat')
library("Seurat")
library("dplyr")
library("Matrix")
library("ggplot2")
library("viridis")
library("ggprism")
library("viridisLite")
library("patchwork")
dat <- readRDS('panc.rds')
dat$disease_id <- gsub(' ', '', dat$disease_id)
dat$disease_id <- gsub('_', '', dat$disease_id)
dat$disease_state <- factor(x = dat$disease_state, levels = c('Control', 'AAB', 'T1D'))
dat$disease_state <- gsub("Control", "CTL", dat$disease_state)
dat$disease_state <- gsub("AAB", "AAb+", dat$disease_state)
dat$disease_state <- factor(x = dat$disease_state, levels = c('CTL', 'AAb+', 'T1D'))

## Stacked bar plot
df <- data.frame("Cell_Type"=dat@meta.data$cell_type, "Condition"=dat@meta.data$disease_state)
df<-df[df$Cell_Type!="Unknown",]
# df$Condition <- gsub("Control", "CTL", df$Condition)
df$Cell_Type <- gsub("Stellates_Mesenchymal", "Stellates", df$Cell_Type)
df$Cell_Type <- gsub("PP_Gamma", "PP", df$Cell_Type)
df$Condition <- factor(df$Condition, levels = c('Control', 'AAB', 'T1D'))


df <- df %>% group_by(Condition, Cell_Type) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(Mean = mean(Nb)) %>%
  mutate(SE = sd(Nb)/sqrt(length(Nb))) %>%
  # mutate(SE = sd(Nb)/sqrt(Nb)) %>%
  mutate(SD = sd(Nb)) %>%
  mutate(Percent = Nb/C*100) 
# %>%  mutate(SD = sd(Percent))

cols = c('#21918c', '#440154', '#fde725')
## Theme for plotting
xtheme <- theme_bw()+ theme(plot.title = element_text(face = "bold" ,hjust = 0.5, size= 10)
                            ,axis.text.y = element_text(face = "bold",angle = 0, size = 10, hjust=1)
                            ,axis.title.y = element_text(face = "bold", size = rel(0.8))
                            ,axis.text.x = element_text(face = "bold",angle = 0, size = 10)
                            ,axis.title.x = element_text(face = "bold", size = rel(0.8))
                            ,axis.ticks.x=element_blank(), strip.text = element_text(size=10))

pdf('./celltype_grouped_barchart_revised.pdf', height = 8, width = 15)
ggplot(df, aes(x=as.factor(Cell_Type), y=Nb, fill=Condition)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymax = Nb + SE, ymin = ifelse(Nb - SE < 0, 0, Nb - SE)), width=0.2,position=position_dodge((0.9)))+
  # scale_fill_viridis(discrete = T) + 
  scale_fill_manual(values = cols) +
  xlab("") + ylab("# of Cells") + xtheme +
  theme(legend.position='top', axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()


## Stacked bar plot
df <- data.frame("Disease_ID"=dat@meta.data$disease_id, "Condition"=dat@meta.data$disease_state)

df <- df %>% group_by(Condition, Disease_ID) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(Mean = mean(Nb)) %>%
  mutate(SE = sd(Nb)/sqrt(length(Nb))) %>%
  mutate(SD = sd(Nb)) %>%
  mutate(Percent = Nb/C*100) 

# library(stringr)
# class(df)
# df[str_sort(df$Disease_ID, numeric = TRUE),]
# df[df$Disease_ID[order(nchar(df$Disease_ID), df$Disease_ID)],]

library(gtools)
tmp <- mixedsort(df$Disease_ID)

df <- df %>% slice(match(tmp, Disease_ID))
df$Disease_ID <- factor(df$Disease_ID, levels = tmp)
df$Disease_ID <- factor(df$Disease_ID, levels = df$Disease_ID)
## Theme for plotting
xtheme <- theme_bw()+ theme(plot.title = element_text(face = "bold" ,hjust = 0.5, size= 10)
                            ,axis.text.y = element_text(face = "bold",angle = 0, size = 10, hjust=1)
                            ,axis.title.y = element_text(face = "bold", size = rel(0.8))
                            ,axis.text.x = element_text(face = "bold",angle = 0, size = 10)
                            ,axis.title.x = element_text(face = "bold", size = rel(0.8))
                            ,axis.ticks.x=element_blank(), strip.text = element_text(size=10))

pdf('./donor_distribution_barchart_revised.pdf', height = 6, width = 15)
ggplot(df, aes(x=Disease_ID, y=Percent, fill=Condition)) +
  # ggplot(df, aes(x=as.factor(Condition), y=Percent, fill=Disease_ID)) +
  
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  # geom_errorbar(aes(ymax = Nb + SE, ymin = ifelse(Nb - SE < 0, 0, Nb - SE)), width=0.2,position=position_dodge((0.9)))+
  # scale_fill_viridis(discrete = T) + 
  scale_fill_manual(values = cols) +
  xlab("") +  ylab("Percent of Cells in corresponding condition") + 
  xtheme +
  theme(legend.position='top', axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

##################################################################################
## Read data
dat <- readRDS("panc.rds")
## Create response variable
df <- dat@meta.data %>%
  mutate(class_label = recode(disease_state,
                              "Control" = 0,
                              "AAB" = 1,
                              "T1D" = 2))
dat[["class_label_num"]] <- as.numeric(df$class_label)

########################################################################################################################
## Stacked bar plot
## Read data
library("Seurat")
library("dplyr")
library("Matrix")
library("ggplot2")
library("viridis")
library("ggprism")
library("viridisLite")
library("patchwork")

dat <- readRDS("panc.rds")
## Create response variable
df <- dat@meta.data %>%
  mutate(class_label = recode(disease_state,
                              "Control" = 0,
                              "AAB" = 1,
                              "T1D" = 2))
dat[["class_label_num"]] <- as.numeric(df$class_label)

df <- data.frame("Cell_Type"=dat@meta.data$cell_type, "Condition"=dat@meta.data$disease_state)
df <- df %>% group_by(Condition, Cell_Type) %>%
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(Percent = Nb/C*100)

df<-df[df$Cell_Type!="Unknown",]
# df$Condition <- gsub("Control", "CTL", df$Condition)
df$Cell_Type <- gsub("Stellates_Mesenchymal", "Stellates", df$Cell_Type)
df$Cell_Type <- gsub("PP_Gamma", "PP", df$Cell_Type)
df$Condition <- factor(df$Condition, levels = c('Control', 'AAB', 'T1D'))
## Theme for plotting
xtheme <- theme_bw()+ theme(plot.title = element_text(face = "bold" ,hjust = 0.5, size= 10)
                            ,axis.text.y = element_text(face = "bold",angle = 0, size = 10, hjust=1)
                            ,axis.title.y = element_text(face = "bold", size = rel(0.8))
                            ,axis.text.x = element_text(face = "bold",angle = 0, size = 10)
                            ,axis.title.x = element_text(face = "bold", size = rel(0.8))
                            ,axis.ticks.x=element_blank(), strip.text = element_text(size=10))


pdf('celltype_stackbarchart_celltypes_condition_filtered.pdf', height = 4, width = 5)
ggplot(df, aes(fill=Condition, y=Percent, x=Cell_Type)) +
  geom_bar(position="fill", stat="identity") + 
  # scale_fill_viridis(discrete = T) + 
  scale_fill_manual(values = cols) +
  xlab("") + xtheme +
  theme(legend.position='top', axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

##################################################################################
## Remove unknown cells
dat <- subset(dat, cell_type !='Unknown')
dat$cell_type <- droplevels(dat$cell_type)

## Plot umap disease_state
pdf("./UMAP_Markers_Disease_state.pdf", width=15, height=40)
FeaturePlot(dat, features = c("PRSS1", "GCG", "INS", "SST", "KRT19", 
                              "VWF", #"GHRL", "NCF2", 
                              "PPY", "COL1A1"), label = T, repel = T, split.by = 'disease_state',
            pt.size = 0.001, ncol = 5) & theme(axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()

pdf("./UMAP_Markers_Disease_state_new.pdf", width=8, height=11)
FeaturePlot(dat, features = c("PRSS1", "GCG", "INS", "SST"), label = T, repel = T, split.by = 'disease_state',
            pt.size = 0.001) & theme(axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()


# pdf("./UMAP_Markers_Disease_state_new_v1.pdf", width=8, height=11)
# FeaturePlot(dat, features = c("PRSS1", "GCG", "INS", "SST"), label = T, repel = T, split.by = 'disease_state',
#             pt.size = NULL) & theme(axis.title.x = element_blank(), axis.title.y = element_blank())
# dev.off()


# pdf("./UMAP_Markers_Disease_state_v2.pdf", width=15, height=40)
# patchwork::wrap_plots(FeaturePlot(dat, features= c("PRSS1", "GCG", "INS", "SST", "KRT19",
#                                   "VWF",
#                                   #"GHRL", "NCF2",
#                                   "PPY", "COL1A1"), split.by = "disease_state", combine=FALSE,
#                                   label = T, repel = T))
# dev.off()



dat$disease_id <- factor(dat$disease_id, levels = tmp)
pdf('./UMAP_SampleIDs.pdf', width=10, height=6)
DimPlot(dat, group.by = "disease_id", reduction = "umap")
dev.off()

pdf('./UMAP_HPAPIDs.pdf', width=10, height=6)
DimPlot(dat, group.by = "hpap_id", reduction = "umap")
dev.off()

pdf('./UMAP_Disease_State.pdf', width=7, height=6)
DimPlot(dat, group.by = "disease_state", reduction = "umap", cols = c('#21918c', '#440154', '#fde725')) + xtheme #, cols= c("#e30800", "#f56505", "#dec400", "#006630", "#0223c7","#5b02c7", "#00b0e6", "#c40080", "#02f00a", "#7d3301", "#000000"))
dev.off()

pdf('./UMAP_Disease_State_v1.pdf', width=7, height=6)
DimPlot(dat, group.by = "disease_state", reduction = "umap", cols = c('#21918c', '#440154', '#fde725'), raster = F) + xtheme #, cols= c("#e30800", "#f56505", "#dec400", "#006630", "#0223c7","#5b02c7", "#00b0e6", "#c40080", "#02f00a", "#7d3301", "#000000"))
dev.off()

# '#21918c', '#440154', '#FC4E07', '#fde725'
pdf('./Topgenes_diseasestate.pdf', width=20, height=5)
VlnPlot(dat, features = c("HLA-A", "HLA-B", "HLA-C", 'HLA-E'), ncol = 4, group.by = "disease_state",
        pt.size = 0.001, combine = T,  cols = c('#21918c', '#440154', '#fde725')) & theme(legend.position = 'none',
                                                                                          axis.title.x = element_blank())
dev.off()

pdf('./Topgenes_diseasestate_v1.pdf', width=20, height=5)
VlnPlot(dat, features = c("HLA-A", "HLA-B", "HLA-C", 'HLA-E'), ncol = 4, group.by = "disease_state",
        pt.size = 0, combine = T,  cols = c('#21918c', '#440154', '#fde725')) & theme(legend.position = 'none',
                                                                                      axis.title.x = element_blank())
dev.off()


################################################################################
## Create stacked bar chart of cell types and aab+ donors
## Create response variable
df <- dat@meta.data %>% 
  mutate(class_label = recode(disease_state, 
                              "Control" = 0, 
                              "AAB" = 1,
                              "T1D" = 2)) 
dat[["class_label_num"]] <- as.numeric(df$class_label)

########################################################################################################################
## Stacked bar plot 
df <- data.frame("Cell_Type"=dat@meta.data$cell_type, "Condition"=dat@meta.data$disease_state, "DonorID"=dat@meta.data$hpap_id) 
df <- df[df$Condition=='AAB',]
df <- df[df$Cell_Type!="Unknown",]
df$Cell_Type <- gsub("Stellates_Mesenchymal", "Stellates", df$Cell_Type)
df$Cell_Type <- gsub("PP_Gamma", "PP", df$Cell_Type)

df <- df %>% group_by(DonorID, Cell_Type) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(Percent = Nb/C*100) 

# df$Condition <- gsub("Control", "CTL", df$Condition)

## Theme for plotting
xtheme <- theme_bw()+ theme(plot.title = element_text(face = "bold" ,hjust = 0.5, size= 10)
                            ,axis.text.y = element_text(face = "bold",angle = 0, size = 10, hjust=1)
                            ,axis.title.y = element_text(face = "bold", size = rel(0.8))
                            ,axis.text.x = element_text(face = "bold",angle = 0, size = 10)
                            ,axis.title.x = element_text(face = "bold", size = rel(0.8))
                            ,axis.ticks.x=element_blank(), strip.text = element_text(size=10)) 


# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# cols = gg_color_hue(8)

cols <-  c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#999999", "#000000", 'red')

pdf('./aab_celltype_donorids_stackbarchart.pdf', height = 5, width = 7)
ggplot(df, aes(fill=Cell_Type, y=Percent, x=DonorID)) + 
  geom_bar(position="fill", stat="identity") + #scale_fill_viridis(discrete = T) +
  # scale_fill_viridis_c(option = "magma") +
  # scale_fill_manual(values = cols)+
  # scale_fill_discrete(palette=cols) +
  scale_fill_manual(values =cols) +
  xlab("") + xtheme + 
  theme(legend.position='top', axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()


pdf('./aab_celltype_grouped_barchart_revised.pdf', height = 8, width = 15)
ggplot(df, aes(x=as.factor(DonorID), y=Nb, fill=Cell_Type)) +
  geom_bar(position=position_dodge2(width = 0.9, preserve = "single"), stat="identity", colour='black') +
  # geom_errorbar(aes(ymax = Nb + SE, ymin = ifelse(Nb - SE < 0, 0, Nb - SE)), width=0.2,position=position_dodge((0.9)))+
  # scale_fill_viridis(discrete = T) + 
  xlab("") + ylab("# of Cells") + xtheme +
  scale_fill_manual(values =cols) +
  theme(legend.position='top', axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()


pdf('./aab_distribution_barchart_revised.pdf', height = 6, width = 15)
ggplot(df, aes(x=DonorID, y=Percent, fill=Cell_Type)) +
  # ggplot(df, aes(x=as.factor(Condition), y=Percent, fill=Disease_ID)) +
  geom_bar(position=position_dodge(width = 0.9, preserve = "single"), stat="identity", colour='black') +
  # geom_errorbar(aes(ymax = Nb + SE, ymin = ifelse(Nb - SE < 0, 0, Nb - SE)), width=0.2,position=position_dodge((0.9)))+
  # scale_fill_viridis(discrete = T) + 
  scale_fill_manual(values =cols) +
  xlab("") +  ylab("Percent of cells in each donor") + 
  xtheme +
  theme(legend.position='top', axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()
##################################
## CXCL8
##################################

cols = c('#21918c', '#440154', '#fde725')
cols <-  c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#999999", "#000000", 'red')

rm(list=ls())
setwd('/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/plots/revision/Fig1/')
library('Seurat')
library("Seurat")
library("dplyr")
library("Matrix")
library("ggplot2")
library("viridis")
library("ggprism")
library("viridisLite")
library("patchwork")
dat <- readRDS('panc.rds')

dat <- subset(dat, cell_type !='Unknown')
# dat_test <- subset(dat, disease_state =='T1D')

dat$cell_type <- droplevels(dat$cell_type)
dat$disease_state <- factor(x = dat$disease_state, levels = c('Control', 'AAB', 'T1D'))

# Idents(dat) <- 'cell_type'

# VlnPlot(dat, "CXCL8", group.by = "hpap_id", split.by = "disease_state")
# VlnPlot(dat, "CXCL8", group.by = "disease_state", split.by = "hpap_id")

## custom donor colors
df <- table(dat$disease_id, dat$hpap_id)
df <- as.data.frame(table(dat$disease_id, dat$hpap_id))
df <- df[df$Freq>0,]

a <- '#440154'
c <- '#21918c'
t <- '#fde725'
cols <- c(t,c,t,a,c,c,t,a,t,c,c,c,c,a,c,c,c,a,c,a,c,a,a,c,c,c,t,c,c,c,t,c,c,t,a,c,c,c,c,c,t,t,a,c,c,c,c,c,c,a)
pdf("./CXCL8_Vln_v1.pdf", height = 5, width = 10)
VlnPlot(dat, "CXCL8", group.by = "disease_state", split.by = "hpap_id", pt.size = 0) + scale_fill_manual(values = cols) # + theme(legend.position = 'none')
dev.off()

pdf("./CXCL8_Vln_v2.pdf", height = 5, width = 10)
VlnPlot(dat, "CXCL8", group.by = "disease_state", split.by = "hpap_id", pt.size = 0) + scale_fill_manual(values = cols) + theme(legend.position = 'none') # 
dev.off()


pdf("./CXCL8_Vln_v3.pdf", height = 5, width = 20)
VlnPlot(dat, "CXCL8", group.by = "disease_state", split.by = "hpap_id", pt.size = 0) + scale_fill_manual(values = cols) + theme(legend.position = 'none') # 
dev.off()

cols <-  c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#999999", "#000000", 'red')
pdf("./CXCL8_Vln_Celltype.pdf", height = 5, width = 20)
VlnPlot(dat, "CXCL8", group.by = "disease_state", split.by = "cell_type", pt.size = 0) + scale_fill_manual(values = cols) 
dev.off()


cols <-  c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#999999", "#000000", 'red')
pdf("./CXCL8_Vln_Celltype_v1.pdf", height = 5, width = 15)
VlnPlot(dat, "CXCL8", group.by = "disease_state", split.by = "cell_type", pt.size = 0) + scale_fill_manual(values = cols)  + theme(legend.position = 'none') # 
dev.off()


# VlnPlot(dat, "CXCL8", group.by = "hpap_id", split.by = "disease_state", cols = c('#21918c', '#440154', '#fde725'))
# VlnPlot(dat, "CXCL8", group.by = "disease_state", split.by = "hpap_id") + scale_fill_manual(values = c(rep('#21918c',31),rep('#440154',10),rep('#fde725',9)))

# VlnPlot(dat, "CXCL8", group.by = "hpap_id", split.by = "disease_state", cols = cols, split.plot = T)
# VlnPlot(dat, "CXCL8", group.by = "disease_state", cols = cols)

###############################
## Pseudobulk plot
## For HPAP092
## for 4 genes combined
library("Seurat")
library("stringi")
library("data.table")
library("ggpubr")
library("cowplot")
library("dplyr")
library("tidyverse")
library("ggplot2")
library("gridExtra")
local <- readRDS('panc.rds')
local <- subset(dat, cell_type !='Unknown')
# dat_test <- subset(dat, disease_state =='T1D')

local$cell_type <- droplevels(dat$cell_type)
local$disease_state <- factor(x = local$disease_state, levels = c('Control', 'AAB', 'T1D'))

expr_plot_4genes_combined_allcells <- function(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
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
disease_order_combined <- c("CTL", "AAb+", "MC-AAb+\n(HPAP092)", "T1D")

t1 <- '^CXCL8$'; t2 <- '^IL32$'; t3 <- '^HLA-C$'; t4 <- '^HLA-E$'; 
ylimit1 <- c(0,3); ylimit2 <- c(0,4); ylimit3 <- c(0,5); ylimit4 <- c(0,5); 
filename = "HPAP092_MC_CXCL8.pdf"
p_combined <- expr_plot_4genes_combined_allcells(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
                                                 ylimit1, ylimit2, ylimit3, ylimit4)
pdf(filename, width = 15, height = 4)
p_combined
dev.off()

## MC cell types
expr_plot_4genes_combined_cell <- function(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
                                           ylimit1, ylimit2, ylimit3, ylimit4, which_cell)
{
  local_cell <- subset(local, subset = cell_type == which_cell)
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

which_cell <- 'Ductal'
## MC
cellbarcodes <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/ductal_cellbarcodes.csv")
cellX <- cellbarcodes$X

cellX1 <- cellX[grep(".+43$",cellX)] ## HPAP092
# cellX2 <- cellX[grep(".+50$",cellX)] ## HPAP107
cellX <- c(cellX1)
## subset the MC cells
local_cond_filtered <- local[,colnames(local) %in% cellX]
disease_order_combined <- c("CTL", "AAb+", "MC-AAb+\n(HPAP092)", "T1D")

t1 <- '^CXCL8$'; t2 <- '^IL32$'; t3 <- '^HLA-C$'; t4 <- '^HLA-E$'; 
ylimit1 <- c(0,3); ylimit2 <- c(0,4); ylimit3 <- c(0,5); ylimit4 <- c(0,5); 
filename = "HPAP092_MC_CXCL8_Ductal.pdf"
p_combined <- expr_plot_4genes_combined_cell(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
                                             ylimit1, ylimit2, ylimit3, ylimit4, which_cell)
pdf(filename, width = 15, height = 4)
p_combined
dev.off()

## Average Expression and Percent Expressed
Avg_Exp_Pct_Exp_3grps <- function(sc_file, gene_queryGV1)
{
  ## Subset the cells for percent expressed evaluation
  ## table(sc_file$Type)
  sc_file_Control <- subset(sc_file, subset = disease_state=="Control")
  sc_file_AAB <- subset(sc_file, subset = disease_state=="AAB")
  sc_file_T1D <- subset(sc_file, subset = disease_state == "T1D")
  
  ## Subset the cells for percent expressed evaluation
  
  ## Control
  tmp.case <- DotPlot(sc_file_Control, features=gene_queryGV1, group.by=c("cell_type"), scale = F)
  tmp.case <- tmp.case$data
  tmp.case <- tmp.case[,c(1,2,4)]
  # tmp.case <- round(tmp.case[,c(1,2)],2)
  # rownames(tmp.case) <- tmp.case$id
  colnames(tmp.case) <- c(paste0("Average Expression"," of ", gene_queryGV1, " in ", "CTL"),
                          paste0("Percent Expressed", " in ", "CTL"),
                          "Cell Types")
  
  ## AAB
  tmp.case1 <- DotPlot(sc_file_AAB, features=gene_queryGV1, group.by=c("cell_type"), scale = F)
  tmp.case1 <- tmp.case1$data
  tmp.case1 <- tmp.case1[,c(1,2,4)]
  # tmp.case <- round(tmp.case[,c(1,2)],2)
  # rownames(tmp.case1) <- tmp.case1$id
  colnames(tmp.case1) <- c(paste0("Average Expression"," of ", gene_queryGV1, " in ", "AAb+"),
                           paste0("Percent Expressed", " in ", "AAb+"),
                           "Cell Types")
  
  ## T1D
  tmp.ctl <- DotPlot(sc_file_T1D, features=gene_queryGV1, group.by=c("cell_type"), scale = F)
  tmp.ctl <- tmp.ctl$data
  tmp.ctl <- tmp.ctl[,c(1,2,4)]
  # rownames(tmp.ctl) <- tmp.ctl$id
  colnames(tmp.ctl) <- c(paste0("Average Expression"," of ", gene_queryGV1, " in ", "T1D"),
                         paste0("Percent Expressed", " in ",  "T1D"),
                         "Cell Types")
  
  df <- merge(tmp.case, tmp.case1, by="Cell Types")
  df <- merge(df, tmp.ctl, by="Cell Types")
  
  rownames(df) <- df$`Cell Types`
  df$`Cell Types` <- NULL
  # df[,c(2,3,4,5)]
  df <- round(df, 2)
  # print(df, row.names = F)
  df <- cbind('Cell Types'=rownames(df), df)
  
  cell_cts <- data.frame(table(sc_file$cell_type))
  colnames(cell_cts) <- c("Cell Types", "Number of cells")
  
  # df$`Cell Types` <- rownames(df)
  print(df, row.names = F)
  write.csv(df, paste0("./",gene_queryGV1, "_avgexp", ".csv"))
  return(df)
}

## CXCL8
Avg_Exp_Pct_Exp_3grps(sc_file = dat, gene_queryGV1 = "CXCL8")

################################################################################
## Validation cohort
rm(list=ls())

setwd('/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/plots/revision/Fig1/')
library("Seurat")
library("stringi")
library("data.table")
library("ggpubr")
library("cowplot")
library("dplyr")
library("tidyverse")
library("ggplot2")
library("gridExtra")
local <- readRDS('/mnt/alvand/abhijeet/aab/apr_WO-T2D/objs/4_donors/results/12.04.23/panc_WO_RB_MT_validation_n4_12042023.rds')
local$disease_state <- local$grp
local$disease_id <- local$sample_ids
## T1D vs Control
cond = "AAB"
disease_order <- c("CTL", "AAB+", "T1D")
## For T1DvsControl
local_cond <- subset(local, subset = disease_state != cond)
# local_cond$disease_state <- droplevels(local_cond$disease_state)

##########################################################################################

## for 4 genes
expr_plot_cell <- function(t1, t2, t3, t4, filename, which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
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
filename = "Validation_T1DvsControl_HLAI_beta.pdf"
t1 <- '^HLA-A$'; t2 <- '^HLA-B$'; t3 <- '^HLA-C$'; t4 <- '^HLA-E$';
ylimit1 <- c(0,4); ylimit2 <- c(0,4); ylimit3 <- c(0,4); ylimit4 <- c(0,2);
plot1 <- expr_plot_cell(t1, t2, t3, t4, filename, which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
pdf(filename, width = 15, height = 5)
plot1
dev.off()

## for 4 genes
expr_plot_all <- function(t1, t2, t3, t4, filename, which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
{
  local_cell <- local_cond
  # local_cell <- subset(local_cond, subset = cell_type == which_cell)
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
    
    # p[[i]] <- ggboxplot(df, x = "Type", y = "Gene", fill = "Type",#color = "Type",
    #                    #palette = c("#00AFBB", "#E7B800"), # "#FC4E07",
    #                    palette = c('#21918c', '#fde725'),
    #                    add=c("jitter"),#add.params = list(fill="white"),
    #                    order = disease_order,
    #                    shape = "Type", #size = 0.1,
    #                    ylab = i, xlab = FALSE) #, xlab = "Groups")
    
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


filename = "Validation_T1DvsControl_HLAI_All_Test.pdf"
t1 <- '^HLA-A$'; t2 <- '^HLA-B$'; t3 <- '^HLA-C$'; t4 <- '^HLA-E$';
ylimit1 <- c(0,4); ylimit2 <- c(0,4); ylimit3 <- c(0,2); ylimit4 <- c(0,2);
plot1 <- expr_plot_all(t1, t2, t3, t4, filename, which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
pdf(filename, width = 15, height = 5)
plot1
dev.off()

t1 <- '^HLA-DMA$'; t2 <- '^HLA-DQB1$'; t3 <- '^HLA-DRA$'; t4 <- '^HLA-DPA1$'; 
ylimit1 <- c(0,0.5); ylimit2 <- c(0,0.5); ylimit3 <- c(0,0.5); ylimit4 <- c(0,0.5); 
filename = "Validation_T1DvsControl_HLAII_p1.pdf"
plot1 <- expr_plot_all(t1, t2, t3, t4, filename, which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
pdf(filename, width = 15, height = 5)
plot1
dev.off()

t1 <- '^HLA-DRB1$'; t2 <- '^HLA-DQA2$'; t3 <- '^HLA-DPB1$'; t4 <- '^HLA-DRB5$';
ylimit1 <- c(0,0.5); ylimit2 <- c(0,0.5); ylimit3 <- c(0,0.5); ylimit4 <- c(0,0.5); 
filename = "Validation_T1DvsControl_HLAII_p2.pdf"
plot1 <- expr_plot_all(t1, t2, t3, t4, filename, which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
pdf(filename, width = 15, height = 5)
plot1
dev.off()

t1 <- '^INS$'; t2 <- '^IL32$'; t3 <- '^TNFAIP3$'; t4 <- '^LMO7$'; 
ylimit1 <- c(0,10); ylimit2 <- c(0,5); ylimit3 <- c(0,1); ylimit4 <- c(0,1); 
filename = "Validation_T1DvsControl_Others.pdf"
plot1 <- expr_plot_all(t1, t2, t3, t4, filename, which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
pdf(filename, width = 15, height = 5)
plot1
dev.off()

t1 <- '^CXCL8$'; t2 <- '^HLA-F$'; t3 <- '^HLA-G$'; t4 <- '^LMO7$'; 
ylimit1 <- c(0,1); ylimit2 <- c(0,1); ylimit3 <- c(0,0.5); ylimit4 <- c(0,0.5); 
filename = "Validation_T1DvsControl_Others2.pdf"
plot1 <- expr_plot_all(t1, t2, t3, t4, filename, which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
pdf(filename, width = 15, height = 5)
plot1
dev.off()


which_cell = "Ductal"

t1 <- '^CXCL8$'; t2 <- '^HLA-F$'; t3 <- '^HLA-G$'; t4 <- '^LMO7$'; 
ylimit1 <- c(0,1); ylimit2 <- c(0,1); ylimit3 <- c(0,0.5); ylimit4 <- c(0,2); 
filename = "Validation_T1DvsControl_Others_Ductal.pdf"
plot1 <- expr_plot_cell(t1, t2, t3, t4, filename, which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
pdf(filename, width = 15, height = 5)
plot1
dev.off()

which_cell = "Immune"

t1 <- '^CXCL8$'; t2 <- '^HLA-F$'; t3 <- '^HLA-G$'; t4 <- '^LMO7$'; 
ylimit1 <- c(0,1); ylimit2 <- c(0,1); ylimit3 <- c(0,0.5); ylimit4 <- c(0,0.5); 
filename = "Validation_T1DvsControl_Others_Immune.pdf"
plot1 <- expr_plot_cell(t1, t2, t3, t4, filename, which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
pdf(filename, width = 15, height = 5)
plot1
dev.off()


## Average Expression and Percent Expressed
Avg_Exp_Pct_Exp_2grps <- function(sc_file, gene_queryGV1)
{
  ## Subset the cells for percent expressed evaluation
  ## table(sc_file$Type)
  sc_file_Control <- subset(sc_file, subset = disease_state=="CTL")
  sc_file_T1D <- subset(sc_file, subset = disease_state == "T1D")
  
  ## Subset the cells for percent expressed evaluation
  
  ## Control
  tmp.case <- DotPlot(sc_file_Control, features=gene_queryGV1, group.by=c("cell_type"), scale = F)
  tmp.case <- tmp.case$data
  tmp.case <- tmp.case[,c(1,2,4)]
  # tmp.case <- round(tmp.case[,c(1,2)],2)
  # rownames(tmp.case) <- tmp.case$id
  colnames(tmp.case) <- c(paste0("Average Expression"," of ", gene_queryGV1, " in ", "CTL"),
                          paste0("Percent Expressed", " in ", "CTL"),
                          "Cell Types")
  
  ## T1D
  tmp.ctl <- DotPlot(sc_file_T1D, features=gene_queryGV1, group.by=c("cell_type"), scale = F)
  tmp.ctl <- tmp.ctl$data
  tmp.ctl <- tmp.ctl[,c(1,2,4)]
  # rownames(tmp.ctl) <- tmp.ctl$id
  colnames(tmp.ctl) <- c(paste0("Average Expression"," of ", gene_queryGV1, " in ", "T1D"),
                         paste0("Percent Expressed", " in ",  "T1D"),
                         "Cell Types")
  
  df <- merge(tmp.case, tmp.ctl, by="Cell Types")
  
  rownames(df) <- df$`Cell Types`
  df$`Cell Types` <- NULL
  # df[,c(2,3,4,5)]
  df <- round(df, 2)
  # print(df, row.names = F)
  df <- cbind('Cell Types'=rownames(df), df)
  
  cell_cts <- data.frame(table(sc_file$cell_type))
  colnames(cell_cts) <- c("Cell Types", "Number of cells")
  
  # df$`Cell Types` <- rownames(df)
  print(df, row.names = F)
  write.csv(df, paste0("./",gene_queryGV1, "_avgexp_validation", ".csv"))
  return(df)
}

## CXCL8
Avg_Exp_Pct_Exp_2grps(sc_file = local, gene_queryGV1 = "CXCL8")
#####################################################################
## Overall classifiers
Classifier <- c('T1D-CTL', 'T1D-AAb+', 'AAb+-CTL')
Frequency <- c(100,59,50)
df <- as.data.frame(cbind(Classifier, Frequency))
df$Classifier <- factor(x = df$Classifier, levels = c('T1D-CTL', 'T1D-AAb+', 'AAb+-CTL'))

xtheme <- theme(plot.title = element_text(face = "bold" ,hjust = 0.5, size= 12)
                ,axis.text.y = element_text(face = "bold",angle = 0, size = 10, hjust=1)
                ,axis.title.y = element_text(face = "bold", size = rel(0.8))
                ,axis.text.x = element_text(face = "bold",angle = 30, size = 10, hjust = 1)
                ,axis.title.x = element_text(face = "bold", size = rel(0.8))
                ,axis.ticks.x=element_blank(), strip.text = element_text(size=12))

pdf("./OverallClassifier_CXCL8_SelectionFrequency.pdf", height = 4, width = 3)
ggplot(df, aes(fill=Classifier, y=as.numeric(df$Frequency), x=Classifier)) + geom_col(position=position_dodge2(preserve="single")) +
  #  geom_bar(position="dodge", stat="identity")
  ylab("CXCL8 Gene Selection Frequency") +
  # scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw()+ xtheme
dev.off()


#####################################################################
## CellType classifiers
Classifier <- rep(c('T1D-CTL', 'T1D-AAb+', 'AAb+-CTL'),8)
Frequency <- c(100,2,0,8,76,8,29,8,
               79,0,0,0,28,3,6,21,
               0,62,13,0,32,28,46,69)
Cell_Type <- rep(c('Acinar', 'Alpha', 'Beta', 'Delta', 'Ductal', 'Endothelial', 'Immune', 'Stellates'),3)
df <- as.data.frame(cbind(Cell_Type, Classifier, Frequency))
df$Classifier <- factor(x = df$Classifier, levels = c('T1D-CTL', 'T1D-AAb+', 'AAb+-CTL'))
df$Cell_Type <- factor(df$Cell_Type,                                    # Change ordering manually
                       levels = c('Acinar', 'Alpha', 'Beta', 'Delta', 'Ductal', 'Endothelial', 'Immune', 'Stellates'))

xtheme <- theme(plot.title = element_text(face = "bold" ,hjust = 0.5, size= 12)
                ,axis.text.y = element_text(face = "bold",angle = 0, size = 10, hjust=1)
                ,axis.title.y = element_text(face = "bold", size = rel(0.8))
                ,axis.text.x = element_text(face = "bold",angle = 30, size = 10, hjust = 1)
                ,axis.title.x = element_text(face = "bold", size = rel(0.8))
                ,axis.ticks.x=element_blank(), strip.text = element_text(size=12))

pdf("./CellTypeClassifier_CXCL8_SelectionFrequency.pdf", height = 4, width = 8)
ggplot(df, aes(fill=Classifier, y=as.numeric(df$Frequency), x=Cell_Type)) + geom_col(position=position_dodge2(preserve="single")) +
  #  geom_bar(position="dodge", stat="identity")
  ylab("CXCL8 Gene Selection Frequency") +
  # scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw()+ xtheme
dev.off()

################################################################################
## UMAP for donors HPAp092 and HPAp107
donor1 <- colnames(dat)[dat$hpap_id=='HPAP092']
donor2 <- colnames(dat)[dat$hpap_id=='HPAP107']
# new_obj <- subset(dat, cells = donor1)

pdf('./UMAP_Disease_State_HPAP092_Higlight.pdf', width=7, height=6)
DimPlot(object = dat, cells.highlight = donor1, cols.highlight = c('#440154'), cols = "gray", order = TRUE, label = T, repel = T)
dev.off()

pdf('./UMAP_Disease_State_HPAP107_Higlight.pdf', width=7, height=6)
DimPlot(object = dat, cells.highlight = donor2, cols.highlight = c('#440154'), cols = "gray", order = TRUE, label = T, repel = T)
dev.off()


pdf('./UMAP_Disease_State_HPAP092_HPAP107_Higlight.pdf', width=7, height=6)
DimPlot(object = dat, cells.highlight = c(donor1,donor2), cols.highlight = c('#440154'), cols = "gray", order = TRUE, label = T, repel = T)
dev.off()

################################################################################
## Module scores
library('RColorBrewer')
library('patchwork')
dat$disease_state <- factor(x = dat$disease_state, levels = c('Control', 'AAB', 'T1D'))

## HLA-I Module score
genes_HLA_I <- c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-E')
dat <- AddModuleScore(dat,
                      features = list(genes_HLA_I),
                      name="HLAI_genes")

# Plot scores
pdf('./Modulescore_HLA-I_Combined.pdf', width=7, height=6)
FeaturePlot(dat,
            features = "HLAI_genes1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

# Plot scores
pdf('./Modulescore_HLA-I_Condition.pdf', width=20, height=6)
# FeaturePlot(dat,
#             features = "HLAI_genes1", split.by = 'disease_state', label = TRUE, repel = TRUE) +
#   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
patchwork::wrap_plots(FeaturePlot(dat, features= 'HLAI_genes1', split.by = "disease_state", combine=FALSE,
                                  label = T, repel = T)) &
  theme_minimal() & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

## HLA_II module score
genes_HLA_II <- c('HLA-DMA', 'HLA-DQB1', 'HLA-DRA', 'HLA-DPA1', 'HLA-DRB1', 'HLA-DQA2', 'HLA-DPB1', 'HLA-DRB5')
dat <- AddModuleScore(dat,
                      features = list(genes_HLA_II),
                      name="HLAII_genes")
# Plot scores
pdf('./Modulescore_HLA-II_Combined.pdf', width=7, height=6)
FeaturePlot(dat,
            features = "HLAII_genes1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

# Plot scores
pdf('./Modulescore_HLA-II_Condition.pdf', width=20, height=6)
# FeaturePlot(dat,
#             features = "HLAII_genes1", split.by = 'disease_state', label = TRUE, repel = TRUE) +
scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
patchwork::wrap_plots(FeaturePlot(dat, features= 'HLAII_genes1', split.by = "disease_state", combine=FALSE,
                                  label = T, repel = T)) &
  theme_minimal() & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

################################################################################
gc()
rm(list=ls())
options(future.globals.maxSize = 1000000 * 1024^2)
library(Seurat)
library(readxl)
# pseudobulk the counts based on donor-condition-celltype
dat <- readRDS("panc.rds")

# ## All comparisons
# pseudo <- AggregateExpression(dat, assays = "RNA", slot = "counts", return.seurat = T, group.by = c("disease_state", "hpap_id"))
# cts <- pseudo@assays$RNA@counts

## T1DvsAAB comparison
dat_T1DvsAAB <- subset(dat, subset = disease_state != "Control")
pseudo_T1DvsAAB <- AggregateExpression(dat_T1DvsAAB, assays = "RNA", slot = "counts", return.seurat = T, group.by = c("disease_state", "hpap_id"))
Idents(pseudo_T1DvsAAB)
cts <- pseudo_T1DvsAAB@assays$RNA@counts
head(cts)
dim(cts)
features <- readxl::read_xlsx('./t1dvsaab_features_imp.xlsx')
features <- features$Genes
cts <- cts[rownames(cts) %in% features,]
# dim(cts)
# saveRDS(cts, './cts_t1dvsaab_features_imp.rds')


##################
setwd('C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/Revised.submission/tmp/')
library(ggplot2)
library(DESeq2)
library(ggrepel)
cts <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/Revised.submission/tmp/cts_t1dvsaab_features_imp.rds")
# c1 <- colnames(cts)
# names(c1)<- NULL
names(colnames(cts)) <- NULL
rownames(cts)
colnames(cts)
metadata <- as.data.frame(colnames(cts))
colnames(metadata) <- 'Samples'
metadata$group_id <- sub("^([^_]+)_.*$", "\\1", metadata$Samples)
metadata$hpap_id <- sub("^[^_]*_", "", metadata$Samples)

colnames(metadata)
head(metadata)
metadata$Samples <- as.factor(metadata$Samples)
metadata$group_id <- as.factor(metadata$group_id)
rownames(metadata) <- metadata$Samples

# Create DESeq2 object        
dds <- DESeqDataSetFromMatrix(cts, 
                              colData = metadata, 
                              design = ~ group_id)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
# z <- DESeq2::plotPCA(rld, ntop = nrow(rld), intgroup = "group_id")
# z + geom_label(aes(label=metadata$Samples))

z <- DESeq2::plotPCA(rld, ntop = 500, intgroup = "group_id")
# z + geom_label(aes(label=metadata$Samples))
nudge <- position_nudge(y = 1)

pdf('./rld_pb.pdf', height = 8, width = 8)
z + geom_text_repel(aes(label = metadata$hpap_id), position = nudge) + theme_bw()
dev.off()


# vsd <- vst(dds, blind=FALSE)
# z <- DESeq2::plotPCA(vsd, ntop = 500, intgroup = "group_id")
# # z + geom_label(aes(label=metadata$Samples))
# nudge <- position_nudge(y = 1)
# 
# pdf('./vsd_test.pdf', height = 8, width = 8)
# z + geom_text(aes(label = metadata$hpap_id), position = nudge)
# dev.off()

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
rownames(metadata) <- metadata$Samples
# Plot heatmap
library(pheatmap)

pdf('./heatmap_pb.pdf', height = 8, width = 12)
pheatmap(rld_cor, annotation = metadata[, c("group_id"), drop=F])
dev.off()

#################
## Test
#
# ## HLA-I Module score
# genes_HLA_I_Test <- c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G')
# dat <- AddModuleScore(dat,
#                       features = list(genes_HLA_I_Test),
#                       name="HLAI_genes_Test")
# 
# # Plot scores
# pdf('./Modulescore_HLA-I_Combined_Test.pdf', width=7, height=6)
# FeaturePlot(dat,
#             features = "HLAI_genes_Test1", label = TRUE, repel = TRUE) +
#   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
# dev.off()
# 
# dat$disease_state <- factor(x = dat$disease_state, levels = c('Control', 'AAB', 'T1D'))
#
# # Plot scores
# pdf('./Modulescore_HLA-I_Condition_Test2.pdf', width=15, height=6)
# FeaturePlot(dat,
#             features = "HLAI_genes_Test1", split.by = 'disease_state', label = TRUE, repel = TRUE, 
#             # combine = F,
#             keep.scale=NULL) +
#   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
# dev.off()
# 
# pdf('./Modulescore_HLA-I_Condition_Test2.pdf', width=20, height=6)
# patchwork::wrap_plots(FeaturePlot(dat, features= 'HLAI_genes_Test1', split.by = "disease_state", combine=FALSE,
#                                   label = T, repel = T)) &
#   theme_minimal() & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
# dev.off()

#######################################################################################################################
rm(list=ls())

setwd('/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/plots/revision/Fig1/')
library('Seurat')
library("Seurat")
library("dplyr")
library("Matrix")
library("ggplot2")
library("viridis")
library("ggprism")
library("viridisLite")
library("patchwork")
dat <- readRDS('panc.rds')
# subcells <- sample(Cells(dat), size=1000, replace=F)
# sub_dat <- subset(x = dat, cells = subcells)

# dat <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/Demo_T1D.rds")

dat$disease_id <- gsub(' ', '', dat$disease_id)
dat$disease_id <- gsub('_', '', dat$disease_id)
dat$disease_state <- factor(x = dat$disease_state, levels = c('Control', 'AAB', 'T1D'))
dat$disease_state <- gsub("Control", "CTL", dat$disease_state)
dat$disease_state <- gsub("AAB", "AAb+", dat$disease_state)
dat$disease_state <- factor(x = dat$disease_state, levels = c('CTL', 'AAb+', 'T1D'))

## Remove unknown cells
dat <- subset(dat, cell_type !='Unknown')
dat$cell_type <- droplevels(dat$cell_type)

## MC
cellbarcodes <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/Overall_cellbarcodes.csv")
cellX <- cellbarcodes$X

# cellX1 <- cellX[grep(".+43$",cellX)] ## HPAP092
# cellX2 <- cellX[grep(".+50$",cellX)] ## HPAP107
# cellX <- c(cellX1)
## subset the MC cells
dat_cond_filtered <- dat[,colnames(dat) %in% cellX]

features <- readxl::read_xlsx('./t1dvsaab_features_imp.xlsx')
features <- features$Genes

# tmp.metadata <- dat@meta.data
# tmp.metadata$disease_state <- tmp.metadata[rownames(tmp.metadata) %in% colnames(dat_cond_filtered),]
# tmp.metadata$disease_state <- as.character(tmp.metadata$disease_state)
# tmp.metadata$disease_state[rownames(tmp.metadata) %in% colnames(dat_cond_filtered)] <- 'MC-AAb+'

dat$disease_state <- as.character(dat$disease_state)
dat$disease_state[colnames(dat) %in% colnames(dat_cond_filtered)] <- 'MC-AAb+'
dat$disease_state <- factor(x = dat$disease_state, levels = c('CTL', 'AAb+', 'MC-AAb+', 'T1D'))
# dat$disease_state[colnames(dat) %in% colnames(dat_cond_filtered)] <- 'MC-AAb+ (HPAP092)'
# dat$disease_state <- factor(x = dat$disease_state, levels = c('CTL', 'AAb+', 'MC-AAb+ (HPAP092)', 'T1D'))

# disease_order_combined <- c("CTL", "AAb+", "MC-AAb+\n(HPAP092)", "T1D")
xtheme <- theme_bw()+ theme(plot.title = element_text(face = "bold" ,hjust = 0.5, size= 10)
                            ,axis.text.y = element_text(face = "bold",angle = 0, size = 10, hjust=1)
                            ,axis.title.y = element_text(face = "bold", size = rel(0.8))
                            ,axis.text.x = element_text(face = "bold",angle = 0, size = 10)
                            ,axis.title.x = element_text(face = "bold", size = rel(0.8))
                            ,axis.ticks.x=element_blank(), strip.text = element_text(size=10))
pdf('./UMAP_Disease_State_4groups.pdf', width=7, height=6)
# pdf('./UMAP_Disease_State_4groups_HPAP092.pdf', width=7, height=6)
DimPlot(dat, group.by = "disease_state", reduction = "umap", cols = c('#21918c', '#440154', '#FC4E07', '#fde725')) + xtheme #, cols= c("#e30800", "#f56505", "#dec400", "#006630", "#0223c7","#5b02c7", "#00b0e6", "#c40080", "#02f00a", "#7d3301", "#000000"))
dev.off()
######################################################################################





################################################################################
## validation dataset
gc()
rm(list=ls())
options(future.globals.maxSize = 1000000 * 1024^2)
library(Seurat)
library(readxl)
# pseudobulk the counts based on donor-condition-celltype
dat <- readRDS("panc.rds")
dat$disease_id <- gsub(' ', '', dat$disease_id)
dat$disease_id <- gsub('_', '', dat$disease_id)
dat$disease_state <- factor(x = dat$disease_state, levels = c('Control', 'AAB', 'T1D'))
dat$disease_state <- gsub("Control", "CTL", dat$disease_state)
dat$disease_state <- gsub("AAB", "AAb+", dat$disease_state)
dat$disease_state <- factor(x = dat$disease_state, levels = c('CTL', 'AAb+', 'T1D'))
# ## All comparisons
# pseudo <- AggregateExpression(dat, assays = "RNA", slot = "counts", return.seurat = T, group.by = c("disease_state", "hpap_id"))
# cts <- pseudo@assays$RNA@counts

## T1DvsCTL comparison
dat_T1DvsCTL <- subset(dat, subset = disease_state != "AAb+")
pseudo_T1DvsCTL <- AggregateExpression(dat_T1DvsCTL, assays = "RNA", slot = "counts", return.seurat = T, group.by = c("disease_state", "hpap_id"))
Idents(pseudo_T1DvsCTL)
cts <- pseudo_T1DvsCTL@assays$RNA@counts
head(cts)
dim(cts)
features <- readxl::read_xlsx('./t1dvsctl_features_imp.xlsx')
features <- features$Genes
cts <- cts[rownames(cts) %in% features,]

library(ggplot2)
library(DESeq2)
library(ggrepel)
# c1 <- colnames(cts)
# names(c1)<- NULL
names(colnames(cts)) <- NULL
rownames(cts)
colnames(cts)
metadata <- as.data.frame(colnames(cts))
colnames(metadata) <- 'Samples'
metadata$group_id <- sub("^([^_]+)_.*$", "\\1", metadata$Samples)
metadata$hpap_id <- sub("^[^_]*_", "", metadata$Samples)

colnames(metadata)
head(metadata)
metadata$Samples <- as.factor(metadata$Samples)
metadata$group_id <- as.factor(metadata$group_id)
rownames(metadata) <- metadata$Samples

# Create DESeq2 object        
dds <- DESeqDataSetFromMatrix(cts, 
                              colData = metadata, 
                              design = ~ group_id)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
# z <- DESeq2::plotPCA(rld, ntop = nrow(rld), intgroup = "group_id")
# z + geom_label(aes(label=metadata$Samples))

z <- DESeq2::plotPCA(rld, ntop = 500, intgroup = "group_id")
# z + geom_label(aes(label=metadata$Samples))
nudge <- position_nudge(y = 1)

pdf('./rld_pb.pdf', height = 8, width = 8)
z + geom_text_repel(aes(label = metadata$hpap_id), position = nudge) + theme_bw()
dev.off()


# vsd <- vst(dds, blind=FALSE)
# z <- DESeq2::plotPCA(vsd, ntop = 500, intgroup = "group_id")
# # z + geom_label(aes(label=metadata$Samples))
# nudge <- position_nudge(y = 1)
# 
# pdf('./vsd_test.pdf', height = 8, width = 8)
# z + geom_text(aes(label = metadata$hpap_id), position = nudge)
# dev.off()

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
rownames(metadata) <- metadata$Samples
# Plot heatmap
library(pheatmap)

pdf('./heatmap_pb.pdf', height = 8, width = 12)
pheatmap(rld_cor, annotation = metadata[, c("group_id"), drop=F])
dev.off()

#################
## validation dataset
gc()
rm(list=ls())
options(future.globals.maxSize = 1000000 * 1024^2)

library("readxl")
library("Seurat")
library("stringi")
library("data.table")
library("ggpubr")
library("cowplot")
library("dplyr")
library("tidyverse")
library("ggplot2")
library("gridExtra")
# pseudobulk the counts based on donor-condition-celltype
dat <- readRDS("panc.rds")
dat$disease_id <- gsub(' ', '', dat$disease_id)
dat$disease_id <- gsub('_', '', dat$disease_id)
dat$disease_state <- factor(x = dat$disease_state, levels = c('Control', 'AAB', 'T1D'))
dat$disease_state <- gsub("Control", "CTL", dat$disease_state)
dat$disease_state <- gsub("AAB", "AAb+", dat$disease_state)
dat$disease_state <- factor(x = dat$disease_state, levels = c('CTL', 'AAb+', 'T1D'))
# ## All comparisons
# pseudo <- AggregateExpression(dat, assays = "RNA", slot = "counts", return.seurat = T, group.by = c("disease_state", "hpap_id"))
# cts <- pseudo@assays$RNA@counts

## T1DvsCTL comparison
dat_T1DvsCTL <- subset(dat, subset = disease_state != "AAb+")
dat_T1DvsCTL <- subset(dat_T1DvsCTL, subset = cell_type != "Ductal")
pseudo_T1DvsCTL <- AggregateExpression(dat_T1DvsCTL, assays = "RNA", slot = "counts", return.seurat = T, group.by = c("disease_state", "hpap_id"))
# pseudo_T1DvsCTL <- AggregateExpression(dat_T1DvsCTL, assays = "SCT", slot = "data", return.seurat = T, group.by = c("disease_state", "hpap_id"))
Idents(pseudo_T1DvsCTL)
cts <- pseudo_T1DvsCTL@assays$RNA@counts
# cts <- pseudo_T1DvsCTL@assays$SCT@data
head(cts)
dim(cts)
# features <- readxl::read_xlsx('./t1dvsctl_features_imp.xlsx')
features <- readxl::read_xlsx('./t1dvsctl_features_imp_alpha.xlsx')
# features <- readxl::read_xlsx('./t1dvsctl_features_imp_beta.xlsx')
# features <- readxl::read_xlsx('./t1dvsctl_features_imp_ductal.xlsx')
features <- features$Genes
cts <- cts[rownames(cts) %in% features,]

##
local <- readRDS('/mnt/alvand/abhijeet/aab/apr_WO-T2D/objs/4_donors/results/12.04.23/panc_WO_RB_MT_validation_n4_12042023.rds')
local$disease_state <- local$grp
local$disease_id <- local$sample_ids
## T1D vs Control
cond = "AAB"
## For T1DvsControl
local_cond <- subset(local, subset = disease_state != cond)
local_cond <- subset(local_cond, subset = cell_type != "Ductal")
pseudo_local_cond_T1DvsCTL <- AggregateExpression(local_cond, assays = "RNA", slot = "counts", return.seurat = T, group.by = c("disease_id"))
# pseudo_local_cond_T1DvsCTL <- AggregateExpression(local_cond, assays = "SCT", slot = "data", return.seurat = T, group.by = c("disease_id"))

Idents(pseudo_local_cond_T1DvsCTL)
cts_local_cond <- pseudo_local_cond_T1DvsCTL@assays$RNA@counts
# cts_local_cond <- pseudo_local_cond_T1DvsCTL@assays$SCT@data
head(cts_local_cond)
dim(cts_local_cond)
cts_local_cond <- cts_local_cond[rownames(cts_local_cond) %in% features,]
