rm(list=ls())
options(future.globals.maxSize = 1000000 * 1024^2)
library(Seurat)
# pseudobulk the counts based on donor-condition-celltype
dat <- readRDS("panc.rds")
dat_T1DvsAAB <- subset(dat, subset = disease_state != "Control")
pseudo_T1DvsAAB <- AggregateExpression(dat_T1DvsAAB, assays = "RNA", slot = "counts", return.seurat = T, group.by = c("disease_state", "hpap_id", "cell_type"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_T1DvsAAB))

# the metadata for the pseudobulk object is missing, so we need to add it back
pseudo_T1DvsAAB$cell_type <- sapply(strsplit(Cells(pseudo_T1DvsAAB), split = "_"), "[", 3)
pseudo_T1DvsAAB$sample_id <- sapply(strsplit(Cells(pseudo_T1DvsAAB), split = "_"), "[", 2)
pseudo_T1DvsAAB$disease_state <- sapply(strsplit(Cells(pseudo_T1DvsAAB), split = "_"), "[", 1)
pseudo_T1DvsAAB$cell_type.disease_state <- paste(pseudo_T1DvsAAB$cell_type, pseudo_T1DvsAAB$disease_state, sep = "_")

Idents(pseudo_T1DvsAAB) <- "disease_state"
PB_DE.Allcells <- FindMarkers(object = pseudo_T1DvsAAB,
                         ident.1 = "T1D",
                         ident.2 = "AAB",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Allcells, "PB_DE_T1DvsAAB_Allcells.rds")


Idents(pseudo_T1DvsAAB) <- "cell_type.disease_state"

PB_DE.Alpha <- FindMarkers(object = pseudo_T1DvsAAB, 
                         ident.1 = "Alpha_T1D", 
                         ident.2 = "Alpha_AAB",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Alpha, "PB_DE_T1DvsAAB_Alpha.rds")

PB_DE.Beta <- FindMarkers(object = pseudo_T1DvsAAB, 
                         ident.1 = "Beta_T1D", 
                         ident.2 = "Beta_AAB",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Beta, "PB_DE_T1DvsAAB_Beta.rds")

PB_DE.Acinar <- FindMarkers(object = pseudo_T1DvsAAB, 
                         ident.1 = "Acinar_T1D", 
                         ident.2 = "Acinar_AAB",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Acinar, "PB_DE_T1DvsAAB_Acinar.rds")

PB_DE.Ductal <- FindMarkers(object = pseudo_T1DvsAAB, 
                         ident.1 = "Ductal_T1D", 
                         ident.2 = "Ductal_AAB",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Ductal, "PB_DE_T1DvsAAB_Ductal.rds")

PB_DE.Delta <- FindMarkers(object = pseudo_T1DvsAAB, 
                         ident.1 = "Delta_T1D", 
                         ident.2 = "Delta_AAB",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Delta, "PB_DE_T1DvsAAB_Delta.rds")

PB_DE.Endothelial <- FindMarkers(object = pseudo_T1DvsAAB, 
                         ident.1 = "Endothelial_T1D", 
                         ident.2 = "Endothelial_AAB",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Endothelial, "PB_DE_T1DvsAAB_Endothelial.rds")

PB_DE.Stellates <- FindMarkers(object = pseudo_T1DvsAAB, 
                         ident.1 = "Stellates_T1D", 
                         ident.2 = "Stellates_AAB",
                         test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Stellates, "PB_DE_T1DvsAAB_Stellates.rds")

PB_DE.Immune<- FindMarkers(object = pseudo_T1DvsAAB, 
                               ident.1 = "Immune_T1D", 
                               ident.2 = "Immune_AAB",
                               test.use = "DESeq2", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(PB_DE.Immune, "PB_DE_T1DvsAAB_Immune.rds")
