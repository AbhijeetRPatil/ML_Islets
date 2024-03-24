rm(list=ls())
options(future.globals.maxSize = 1000000 * 1024^2)
library(Seurat)
# pseudobulk the counts based on donor-condition-celltype
dat <- readRDS("panc.rds")
dat_T1DvsCTL <- subset(dat, subset = disease_state != "AAB")
pseudo_T1DvsCTL <- AggregateExpression(dat_T1DvsCTL, assays = "RNA", slot = "counts", return.seurat = T, group.by = c("disease_state", "hpap_id", "cell_type"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_T1DvsCTL))

# the metadata for the pseudobulk object is missing, so we need to add it back
pseudo_T1DvsCTL$cell_type <- sapply(strsplit(Cells(pseudo_T1DvsCTL), split = "_"), "[", 3)
pseudo_T1DvsCTL$sample_id <- sapply(strsplit(Cells(pseudo_T1DvsCTL), split = "_"), "[", 2)
pseudo_T1DvsCTL$disease_state <- sapply(strsplit(Cells(pseudo_T1DvsCTL), split = "_"), "[", 1)
pseudo_T1DvsCTL$cell_type.disease_state <- paste(pseudo_T1DvsCTL$cell_type, pseudo_T1DvsCTL$disease_state, sep = "_")

Idents(pseudo_T1DvsCTL) <- "disease_state"
PB_DE.Allcells <- FindMarkers(object = pseudo_T1DvsCTL,
                         ident.1 = "T1D",
                         ident.2 = "Control",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Allcells, "PB_DE_T1DvsCTL_Allcells.rds")


Idents(pseudo_T1DvsCTL) <- "cell_type.disease_state"

PB_DE.Alpha <- FindMarkers(object = pseudo_T1DvsCTL, 
                         ident.1 = "Alpha_T1D", 
                         ident.2 = "Alpha_Control",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Alpha, "PB_DE_T1DvsCTL_Alpha.rds")

PB_DE.Beta <- FindMarkers(object = pseudo_T1DvsCTL, 
                         ident.1 = "Beta_T1D", 
                         ident.2 = "Beta_Control",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Beta, "PB_DE_T1DvsCTL_Beta.rds")

PB_DE.Acinar <- FindMarkers(object = pseudo_T1DvsCTL, 
                         ident.1 = "Acinar_T1D", 
                         ident.2 = "Acinar_Control",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Acinar, "PB_DE_T1DvsCTL_Acinar.rds")

PB_DE.Ductal <- FindMarkers(object = pseudo_T1DvsCTL, 
                         ident.1 = "Ductal_T1D", 
                         ident.2 = "Ductal_Control",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Ductal, "PB_DE_T1DvsCTL_Ductal.rds")

PB_DE.Delta <- FindMarkers(object = pseudo_T1DvsCTL, 
                         ident.1 = "Delta_T1D", 
                         ident.2 = "Delta_Control",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Delta, "PB_DE_T1DvsCTL_Delta.rds")

PB_DE.Endothelial <- FindMarkers(object = pseudo_T1DvsCTL, 
                         ident.1 = "Endothelial_T1D", 
                         ident.2 = "Endothelial_Control",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Endothelial, "PB_DE_T1DvsCTL_Endothelial.rds")

PB_DE.Stellates <- FindMarkers(object = pseudo_T1DvsCTL, 
                         ident.1 = "Stellates_T1D", 
                         ident.2 = "Stellates_Control",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Stellates, "PB_DE_T1DvsCTL_Stellates.rds")

PB_DE.Immune<- FindMarkers(object = pseudo_T1DvsCTL, 
                               ident.1 = "Immune_T1D", 
                               ident.2 = "Immune_Control",
                               test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Immune, "PB_DE_T1DvsCTL_Immune.rds")
#dat_T1DvsAAB <- subset(dat, subset = disease_state != "CTL")
#saveRDS(dat_T1DvsAAB, "T1DvsAAB_Overall.rds")
#dat_AABvsCTL <- subset(dat, subset = disease_state != "T1D")
#saveRDS(dat_AABvsCTL, "AABvsCTL_Overall.rds")
