rm(list=ls())
options(future.globals.maxSize = 1000000 * 1024^2)
library(Seurat)
# pseudobulk the counts based on donor-condition-celltype
dat <- readRDS('/mnt/alvand/abhijeet/aab/apr_WO-T2D/objs/4_donors/results/12.04.23/panc_WO_RB_MT_validation_n4_12042023.rds')
dat$disease_state <- dat$grp
dat$disease_id <- dat$sample_ids
## T1D vs Control
dat_T1DvsCTL <- subset(dat, subset = disease_state != "AAB")
pseudo_T1DvsCTL <- AggregateExpression(dat_T1DvsCTL, assays = "RNA", slot = "counts", return.seurat = T, group.by = c("disease_state", "disease_id", "cell_type"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_T1DvsCTL))

# the metadata for the pseudobulk object is missing, so we need to add it back
pseudo_T1DvsCTL$cell_type <- sapply(strsplit(Cells(pseudo_T1DvsCTL), split = "_"), "[", 4)
pseudo_T1DvsCTL$sample_id <- sapply(strsplit(Cells(pseudo_T1DvsCTL), split = "_"), "[", 3)
pseudo_T1DvsCTL$disease_state <- sapply(strsplit(Cells(pseudo_T1DvsCTL), split = "_"), "[", 2)
# pseudo_T1DvsCTL$disease_state <- sapply(strsplit(Cells(pseudo_T1DvsCTL), split = "_"), "[", 1)
pseudo_T1DvsCTL$cell_type.disease_state <- paste(pseudo_T1DvsCTL$cell_type, pseudo_T1DvsCTL$disease_state, sep = "_")

# Idents(pseudo_T1DvsCTL) <- "disease_state"
# PB_DE.Allcells <- FindMarkers(object = pseudo_T1DvsCTL,
#                               ident.1 = "T1D",
#                               ident.2 = "Control",
#                               test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
# saveRDS(PB_DE.Allcells, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/PB_DE_T1DvsCTL_Allcells.rds")
# 
# 
Idents(pseudo_T1DvsCTL) <- "cell_type.disease_state"
# 
PB_DE.Alpha <- FindMarkers(object = pseudo_T1DvsCTL,
                           ident.1 = "Alpha_T1D",
                           ident.2 = "Alpha_CTL",
                           test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, 
                           min.cells.group = 0,
                           verbose = FALSE)
# saveRDS(PB_DE.Alpha, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/PB_DE_T1DvsCTL_Alpha.rds")
write.csv(PB_DE.Alpha, './PB_DE_T1DvsCTL_valid_Alpha.csv')

PB_DE.Acinar <- FindMarkers(object = pseudo_T1DvsCTL,
                           ident.1 = "Acinar_T1D",
                           ident.2 = "Acinar_CTL",
                           test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, 
                           min.cells.group = 0,
                           verbose = FALSE)
# saveRDS(PB_DE.Acinar, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/PB_DE_T1DvsCTL_Acinar.rds")
write.csv(PB_DE.Acinar, './PB_DE_T1DvsCTL_valid_Acinar.csv')

PB_DE.Ductal <- FindMarkers(object = pseudo_T1DvsCTL,
                           ident.1 = "Ductal_T1D",
                           ident.2 = "Ductal_CTL",
                           test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, 
                           min.cells.group = 0,
                           verbose = FALSE)
# saveRDS(PB_DE.Ductal, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/PB_DE_T1DvsCTL_Ductal.rds")
write.csv(PB_DE.Ductal, './PB_DE_T1DvsCTL_valid_Ductal.csv')

PB_DE.Immune <- FindMarkers(object = pseudo_T1DvsCTL,
                           ident.1 = "Immune_T1D",
                           ident.2 = "Immune_CTL",
                           test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, 
                           min.cells.group = 0,
                           verbose = FALSE)
# saveRDS(PB_DE.Immune, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/PB_DE_T1DvsCTL_Immune.rds")
write.csv(PB_DE.Immune, './PB_DE_T1DvsCTL_valid_Immune.csv')

PB_DE.Beta <- FindMarkers(object = pseudo_T1DvsCTL, 
                          ident.1 = "Beta_T1D", 
                          ident.2 = "Beta_CTL",
                          test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, 
                          min.cells.group = 0,
                          verbose = FALSE)
write.csv(PB_DE.Beta, './PB_DE_T1DvsCTL_valid_Beta.csv')
# saveRDS(PB_DE.Beta, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/PB_DE_T1DvsCTL_valid_Beta.rds")







rm(list=ls())
options(future.globals.maxSize = 1000000 * 1024^2)
library(Seurat)
# pseudobulk the counts based on donor-condition-celltype
dat <- readRDS('/mnt/alvand/abhijeet/aab/apr_WO-T2D/objs/4_donors/results/12.04.23/panc_WO_RB_MT_validation_n4_12042023.rds')
dat$disease_state <- dat$grp
dat$disease_id <- dat$sample_ids
## T1D vs Control
dat_T1DvsCTL <- subset(dat, subset = disease_state != "AAB")
dat_T1DvsCTL_beta <- subset(dat_T1DvsCTL, subset = cell_type == "Beta")
# saveRDS(dat_T1DvsCTL_beta, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/T1DvsCTL_beta.rds")
dat_T1DvsCTL_beta$cell_type <- droplevels(dat_T1DvsCTL_beta$cell_type)
dat_T1DvsCTL_beta$cell_type.disease_state <- paste(dat_T1DvsCTL_beta$cell_type, dat_T1DvsCTL_beta$disease_state, sep = "_")
Idents(dat_T1DvsCTL_beta) <- "cell_type.disease_state"
dat_T1DvsCTL_beta.DE <- FindMarkers(dat_T1DvsCTL_beta, ident.1 = "Beta_T1D", ident.2 = "Beta_CTL", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
write.csv(dat_T1DvsCTL_beta.DE, './DE_T1DvsCTL_valid_Beta.csv')