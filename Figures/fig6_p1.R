setwd("./revision/Fig6/")
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
local <- readRDS("panc.rds")

########################################################
# Starts here
########################################################
###############################################################################################################
## Fig 6a, 6b, 6c
######### Bar plot showing misclassification in AAB donors overallcellbarcodes

xtheme <- theme(plot.title = element_text(face = "bold" ,hjust = 0.5, size= 12)
                ,axis.text.y = element_text(face = "bold",angle = 0, size = 10, hjust=1)
                ,axis.title.y = element_text(face = "bold", size = rel(0.8))
                ,axis.text.x = element_text(face = "bold",angle = 30, size = 10, hjust = 1)
                ,axis.title.x = element_text(face = "bold", size = rel(0.8))
                ,axis.ticks.x=element_blank(), strip.text = element_text(size=12))

## MC
cellbarcodes <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/third/10_200cv/Overall_cellbarcodes.csv")
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

pdf("./BarPlot_AAB_Overall_MC_Cts_v1.pdf", height = 4, width = 3)
ggplot(data=df, aes(x=reorder(HPAP_ID, -Frequency), y=Frequency)) +
  geom_bar(stat="identity", fill=rgb(0.8,0.1,0.1,0.6))+ scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_text(aes(label=lab), vjust=-0.3, size=3.5)+ xlab("AAB+ HPAP Donors") + ylab("# of AAB cells classified as T1D") +
  theme_minimal() + xtheme
dev.off()

library(readr)
xgb_t1dvsctrl <- read_csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/csv/t1dvsctrl/xgb_t1dvsctrl.csv")
xgb_t1dvsctrl$...1<- NULL
df1 <- data.frame(xgb_t1dvsctrl[grep("HLA", xgb_t1dvsctrl$xgb),])
Classifier=rep("T1D-CTL",nrow(df1))
df1 <- cbind(df1,Classifier)


xgb_t1dvsaab <- read_csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/csv/t1dvsaab/xgb_t1dvsaab.csv")
xgb_t1dvsaab$...1<- NULL
df2 <- data.frame(xgb_t1dvsaab[grep("HLA", xgb_t1dvsaab$xgb),])
Classifier=rep("T1D-AAB",nrow(df2))
df2 <- cbind(df2,Classifier)


xgb_aabvsctrl <- read_csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/csv/aabvsctrl/xgb_aabvsctrl.csv")
xgb_aabvsctrl$...1<- NULL
df3 <- data.frame(xgb_aabvsctrl[grep("HLA", xgb_aabvsctrl$xgb),])
Classifier=rep("AAB-CTL",nrow(df3))
df3 <- cbind(df3,Classifier)

df <- rbind(df1,df2,df3)
colnames(df) <- c("HLA", "Frequency", "Classifier")

df_HLA_I <- df[df$HLA %in% c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G"),]
colnames(df_HLA_I) <- c("HLA", "Frequency", "Classifier")
df_HLA_I[order(df_HLA_I$Frequency, decreasing = T),]

df_HLA_II <- df[!df$HLA %in% c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G"),]
colnames(df_HLA_II) <- c("HLA", "Frequency", "Classifier")
df_HLA_II[order(df_HLA_II$Frequency, decreasing = T),]

df <- rbind(df_HLA_I, df_HLA_II)
df$HLA <- factor(df$HLA,                                    # Change ordering manually
                 levels = c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
                            "HLA-DMA", "HLA-DQB1", "HLA-DRA", "HLA-DPA1", "HLA-DRB1", "HLA-DQA2", "HLA-DPB1", "HLA-DRB5"))
# library(wesanderson)
# names(wes_palettes)
# Grouped
df$Classifier <- factor(x = df$Classifier, levels = c('T1D-CTL', 'T1D-AAB', 'AAB-CTL'))

pdf("./BarPlot_OverallClassifier_HLA_v1.pdf", height = 4, width = 8)
ggplot(df, aes(fill=Classifier, y=Frequency, x=HLA)) + geom_col(position=position_dodge2(preserve="single")) +
  #  geom_bar(position="dodge", stat="identity") #+
  ylab("Selection Frequency") +
  # scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw()+ xtheme
dev.off()

xgb_t1dvsctrl <- read_csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/csv/t1dvsctrl/xgb_t1dvsctrl.csv")
xgb_t1dvsctrl$...1<- NULL
df1 <- data.frame(xgb_t1dvsctrl[xgb_t1dvsctrl$xgb %in% c("INS", "TNFAIP3", "LMO7", "IL32"),])
Classifier=rep("T1D-CTL",nrow(df1))
df1 <- cbind(df1,Classifier)

xgb_t1dvsaab <- read_csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/csv/t1dvsaab/xgb_t1dvsaab.csv")
xgb_t1dvsaab$...1<- NULL
df2 <- data.frame(xgb_t1dvsaab[xgb_t1dvsaab$xgb %in% c("INS", "TNFAIP3", "LMO7", "IL32"),])
Classifier=rep("T1D-AAB",nrow(df2))
df2 <- cbind(df2,Classifier)


xgb_aabvsctrl <- read_csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/csv/aabvsctrl/xgb_aabvsctrl.csv")
xgb_aabvsctrl$...1<- NULL
df3 <- data.frame(xgb_aabvsctrl[xgb_aabvsctrl$xgb %in% c("INS", "TNFAIP3", "LMO7", "IL32"),])
Classifier=rep("AAB-CTL",nrow(df3))
df3 <- cbind(df3,Classifier)

df <- rbind(df1,df2,df3)
colnames(df) <- c("Non-HLA", "Frequency", "Classifier")
df$`Non-HLA` <- factor(df$`Non-HLA`,
                       levels = c("INS", "IL32", "TNFAIP3", "LMO7"))
df$Classifier <- factor(x = df$Classifier, levels = c('T1D-CTL', 'T1D-AAB', 'AAB-CTL'))

pdf("./BarPlot_OverallClassifier_Non_HLA_v2.pdf", height = 4, width = 4)
ggplot(df, aes(fill=Classifier, y=Frequency, x=`Non-HLA`)) + geom_col(position=position_dodge2(preserve="single")) +
  #  geom_bar(position="dodge", stat="identity")
  ylab("Selection Frequency") +
  # scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw()+ xtheme
dev.off()

