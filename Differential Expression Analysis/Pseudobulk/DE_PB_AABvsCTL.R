rm(list=ls())
options(future.globals.maxSize = 1000000 * 1024^2)
library(Seurat)
# pseudobulk the counts based on donor-condition-celltype
dat <- readRDS("/mnt/alvand/abhijeet/aab/apr_WO-T2D/objs/panc_raw_meta_filtered_sct_umap_WO-T2D_WO_RB_MT_04182022.rds")
dat_AABvsCTL <- subset(dat, subset = disease_state != "T1D")
pseudo_AABvsCTL <- AggregateExpression(dat_AABvsCTL, assays = "RNA", slot = "counts", return.seurat = T, group.by = c("disease_state", "hpap_id", "cell_type"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_AABvsCTL))

# the metadata for the pseudobulk object is missing, so we need to add it back
pseudo_AABvsCTL$cell_type <- sapply(strsplit(Cells(pseudo_AABvsCTL), split = "_"), "[", 3)
pseudo_AABvsCTL$sample_id <- sapply(strsplit(Cells(pseudo_AABvsCTL), split = "_"), "[", 2)
pseudo_AABvsCTL$disease_state <- sapply(strsplit(Cells(pseudo_AABvsCTL), split = "_"), "[", 1)
pseudo_AABvsCTL$cell_type.disease_state <- paste(pseudo_AABvsCTL$cell_type, pseudo_AABvsCTL$disease_state, sep = "_")

Idents(pseudo_AABvsCTL) <- "disease_state"
PB_DE.Allcells <- FindMarkers(object = pseudo_AABvsCTL,
                         ident.1 = "AAB",
                         ident.2 = "Control",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Allcells, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/PB_DE_AABvsCTL_Allcells.rds")


Idents(pseudo_AABvsCTL) <- "cell_type.disease_state"

PB_DE.Alpha <- FindMarkers(object = pseudo_AABvsCTL, 
                         ident.1 = "Alpha_AAB", 
                         ident.2 = "Alpha_Control",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Alpha, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/PB_DE_AABvsCTL_Alpha.rds")

PB_DE.Beta <- FindMarkers(object = pseudo_AABvsCTL, 
                         ident.1 = "Beta_AAB", 
                         ident.2 = "Beta_Control",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Beta, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/PB_DE_AABvsCTL_Beta.rds")

PB_DE.Acinar <- FindMarkers(object = pseudo_AABvsCTL, 
                         ident.1 = "Acinar_AAB", 
                         ident.2 = "Acinar_Control",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Acinar, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/PB_DE_AABvsCTL_Acinar.rds")

PB_DE.Ductal <- FindMarkers(object = pseudo_AABvsCTL, 
                         ident.1 = "Ductal_AAB", 
                         ident.2 = "Ductal_Control",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Ductal, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/PB_DE_AABvsCTL_Ductal.rds")

PB_DE.Delta <- FindMarkers(object = pseudo_AABvsCTL, 
                         ident.1 = "Delta_AAB", 
                         ident.2 = "Delta_Control",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Delta, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/PB_DE_AABvsCTL_Delta.rds")

PB_DE.Endothelial <- FindMarkers(object = pseudo_AABvsCTL, 
                         ident.1 = "Endothelial_AAB", 
                         ident.2 = "Endothelial_Control",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Endothelial, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/PB_DE_AABvsCTL_Endothelial.rds")

PB_DE.Stellates <- FindMarkers(object = pseudo_AABvsCTL, 
                         ident.1 = "Stellates_AAB", 
                         ident.2 = "Stellates_Control",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Stellates, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/PB_DE_AABvsCTL_Stellates.rds")

PB_DE.Immune<- FindMarkers(object = pseudo_AABvsCTL, 
                               ident.1 = "Immune_AAB", 
                               ident.2 = "Immune_Control",
                               test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Immune, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/PB_DE_AABvsCTL_Immune.rds")
#dat_AABvsAAB <- subset(dat, subset = disease_state != "CTL")
#saveRDS(dat_AABvsAAB, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/AABvsAAB_Overall.rds")
#dat_AABvsCTL <- subset(dat, subset = disease_state != "AAB")
#saveRDS(dat_AABvsCTL, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/AABvsCTL_Overall.rds")
