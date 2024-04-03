library("splitstackshape")
library("dplyr")
library("ggplot2")
rm(list = ls())
## Load all functions required
source("./fea-plot-lib.R")
## Set the path of cell type 
#setwd("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/plots/t1dvsctrl/beta/")
setwd("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/plots/revision/Fig4/")
#############################################################################
#############################################################################
### KEGG
#############################################################################
## Without FDR corrected
KEGG_t1dvsctrl <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/csv/t1dvsctrl/beta/KEGG.tsv", sep = "\t")
KEGG_t1dvsctrl <- concat.split.multiple(data = KEGG_t1dvsctrl, split.cols = c("Term"), seps = ":")
colnames(KEGG_t1dvsctrl)[length(colnames(KEGG_t1dvsctrl))] = "KEGG_Pathways"
## Pass to function
df <- KEGG_t1dvsctrl
file_name <- "KEGG_t1dvsctrl_PVAL_005.pdf"
plot_title <- "Beta- T1D vs CTRL (KEGG)  P-value < 0.05"
pdf(file_name, width = 8, height = 8)
KEGG_plot(df, plot_title)
dev.off()
# 
#############################################################################
## with FDR corrected
KEGG_t1dvsctrl <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/csv/t1dvsctrl/beta/KEGG.tsv", sep = "\t")
KEGG_t1dvsctrl <- KEGG_t1dvsctrl %>%
  filter(FDR<0.05)
KEGG_t1dvsctrl <- concat.split.multiple(data = KEGG_t1dvsctrl, split.cols = c("Term"), seps = ":")
colnames(KEGG_t1dvsctrl)[length(colnames(KEGG_t1dvsctrl))] = "KEGG_Pathways"

## Pass to function
df <- KEGG_t1dvsctrl
file_name <- "KEGG_t1dvsctrl_FDR_005.pdf"
plot_title <- "Beta- T1D vs CTRL (KEGG)  FDR < 0.05"
pdf(file_name, width = 8, height = 8)
KEGG_plot_FDR(df, plot_title)
dev.off()
#############################################################################
#############################################################################
### REACTOME
#############################################################################
## Without FDR corrected
REACTOME_t1dvsctrl <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/csv/t1dvsctrl/beta/REACTOME.tsv", sep = "\t")
REACTOME_t1dvsctrl <- concat.split.multiple(data = REACTOME_t1dvsctrl, split.cols = c("Term"), seps = "~")
colnames(REACTOME_t1dvsctrl)[length(colnames(REACTOME_t1dvsctrl))] = "REACTOME_Pathways"
## Pass to function
df <- REACTOME_t1dvsctrl
file_name <- "REACTOME_t1dvsctrl_PVAL_005.pdf"
plot_title <- "Beta- T1D vs CTRL (REACTOME)  P-value < 0.05"
pdf(file_name, width = 20, height = 10)
REACTOME_plot(df, plot_title)
dev.off()
#############################################################################
## with FDR corrected
REACTOME_t1dvsctrl <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/csv/t1dvsctrl/beta/REACTOME.tsv", sep = "\t")
REACTOME_t1dvsctrl <- REACTOME_t1dvsctrl %>%
  filter(FDR<0.05)
REACTOME_t1dvsctrl <- concat.split.multiple(data = REACTOME_t1dvsctrl, split.cols = c("Term"), seps = "~")
colnames(REACTOME_t1dvsctrl)[length(colnames(REACTOME_t1dvsctrl))] = "REACTOME_Pathways"
## Pass to function
df <- REACTOME_t1dvsctrl
file_name <- "REACTOME_t1dvsctrl_FDR_005.pdf"
plot_title <- "Beta- T1D vs CTRL (REACTOME)  FDR < 0.05"
pdf(file_name, width = 20, height = 10)
REACTOME_plot_FDR(df, plot_title)
dev.off()
#################################################################################
#################################################################################
### GO BP
#################################################################################
## Without FDR corrected
GO_BP_t1dvsctrl <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/csv/t1dvsctrl/beta/GO_BP.tsv", sep="\t")
GO_BP_t1dvsctrl <- concat.split.multiple(data = GO_BP_t1dvsctrl, split.cols = c("Term"), seps = "~")
colnames(GO_BP_t1dvsctrl)[length(colnames(GO_BP_t1dvsctrl))] = "GO_Term"
## Pass to function
df <- GO_BP_t1dvsctrl
file_name <- "GO_BP_t1dvsctrl_PVAL_005.pdf"
plot_title <- "Beta- T1D vs CTRL (GO BP)  P-value < 0.05"
pdf(file_name, width = 18, height = 12)
GO_plot(df, plot_title)
dev.off()
#################################################################################
## With FDR corrected
GO_BP_t1dvsctrl <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/csv/t1dvsctrl/beta/GO_BP.tsv", sep="\t")
GO_BP_t1dvsctrl <- GO_BP_t1dvsctrl %>%
  filter(FDR<0.05)
GO_BP_t1dvsctrl <- concat.split.multiple(data = GO_BP_t1dvsctrl, split.cols = c("Term"), seps = "~")
colnames(GO_BP_t1dvsctrl)[length(colnames(GO_BP_t1dvsctrl))] = "GO_Term"
## Pass to function
df <- GO_BP_t1dvsctrl
file_name <- "GO_BP_t1dvsctrl_FDR_005.pdf"
plot_title <- "Beta- T1D vs CTRL (GO BP) FDR < 0.05"
pdf(file_name, width = 15, height = 10)
GO_plot_FDR(df, plot_title)
dev.off()
#################################################################################
#################################################################################
### GO CC
#################################################################################
## Without FDR corrected
GO_CC_t1dvsctrl <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/csv/t1dvsctrl/beta/GO_CC.tsv", sep="\t")
GO_CC_t1dvsctrl <- concat.split.multiple(data = GO_CC_t1dvsctrl, split.cols = c("Term"), seps = "~")
colnames(GO_CC_t1dvsctrl)[length(colnames(GO_CC_t1dvsctrl))] = "GO_Term"
## Pass to function
df <- GO_CC_t1dvsctrl
file_name <- "GO_CC_t1dvsctrl_PVAL_005.pdf"
plot_title <- "Beta- T1D vs CTRL (GO CC)  P-value < 0.05"

pdf(file_name, width = 18, height = 12)
GO_plot(df, plot_title)
dev.off()
#################################################################################
## With FDR corrected
GO_CC_t1dvsctrl <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/csv/t1dvsctrl/beta/GO_CC.tsv", sep="\t")
GO_CC_t1dvsctrl <- GO_CC_t1dvsctrl %>%
  filter(FDR<0.05)
GO_CC_t1dvsctrl <- concat.split.multiple(data = GO_CC_t1dvsctrl, split.cols = c("Term"), seps = "~")
colnames(GO_CC_t1dvsctrl)[length(colnames(GO_CC_t1dvsctrl))] = "GO_Term"
## Pass to function
df <- GO_CC_t1dvsctrl
file_name <- "GO_CC_t1dvsctrl_FDR_005.pdf"
plot_title <- "Beta- T1D vs CTRL (GO CC) FDR < 0.05"
PVALue <- "FDR"
pdf(file_name, width = 15, height = 10)
GO_plot_FDR(df, plot_title)
dev.off()
#################################################################################
#################################################################################
### GO MF
#################################################################################
## Without FDR corrected
GO_MF_t1dvsctrl <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/csv/t1dvsctrl/beta/GO_MF.tsv", sep="\t")
GO_MF_t1dvsctrl <- concat.split.multiple(data = GO_MF_t1dvsctrl, split.cols = c("Term"), seps = "~")
colnames(GO_MF_t1dvsctrl)[length(colnames(GO_MF_t1dvsctrl))] = "GO_Term"
## Pass to function
df <- GO_MF_t1dvsctrl
file_name <- "GO_MF_t1dvsctrl_PVAL_005.pdf"
plot_title <- "Beta- T1D vs CTRL (GO MF)  P-value < 0.05"

pdf(file_name, width = 18, height = 12)
GO_plot(df, plot_title)
dev.off()
#################################################################################
## With FDR corrected
GO_MF_t1dvsctrl <- read.csv("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/csv/t1dvsctrl/beta/GO_MF.tsv", sep="\t")
GO_MF_t1dvsctrl <- GO_MF_t1dvsctrl %>%
  filter(FDR<0.05)
GO_MF_t1dvsctrl <- concat.split.multiple(data = GO_MF_t1dvsctrl, split.cols = c("Term"), seps = "~")
colnames(GO_MF_t1dvsctrl)[length(colnames(GO_MF_t1dvsctrl))] = "GO_Term"
## Pass to function
df <- GO_MF_t1dvsctrl
file_name <- "GO_MF_t1dvsctrl_FDR_005.pdf"
plot_title <- "Beta- T1D vs CTRL (GO MF) FDR < 0.05"
pdf(file_name, width = 15, height = 10)
GO_plot_FDR(df, plot_title)
dev.off()
#################################################################################
