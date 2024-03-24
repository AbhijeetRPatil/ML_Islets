options(future.globals.maxSize = 1000000 * 1024^2)
library("Seurat")
dat <- readRDS("/mnt/alvand/abhijeet/aab/apr_WO-T2D/objs/panc_raw_meta_filtered_sct_umap_WO-T2D_WO_RB_MT_04182022.rds")
dat_T1DvsCTL <- subset(dat, subset = disease_state != "AAB")
# saveRDS(dat_T1DvsCTL, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/T1DvsCTL_Overall.rds")
Idents(dat_T1DvsCTL) <- "disease_state"
dat_T1DvsCTL.DE <- FindMarkers(dat_T1DvsCTL, ident.1 = "T1D", ident.2 = "Control", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsCTL.DE, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/DE_T1DvsCTL.rds")

## Alpha
dat_T1DvsCTL_alpha <- subset(dat_T1DvsCTL, subset = cell_type == "Alpha")
# saveRDS(dat_T1DvsCTL_alpha, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/T1DvsCTL_alpha.rds")
dat_T1DvsCTL_alpha$cell_type <- droplevels(dat_T1DvsCTL_alpha$cell_type)
dat_T1DvsCTL_alpha$cell_type.disease_state <- paste(dat_T1DvsCTL_alpha$cell_type, dat_T1DvsCTL_alpha$disease_state, sep = "_")
Idents(dat_T1DvsCTL_alpha) <- "cell_type.disease_state"
dat_T1DvsCTL_alpha.DE <- FindMarkers(dat_T1DvsCTL_alpha, ident.1 = "Alpha_T1D", ident.2 = "Alpha_Control", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsCTL_alpha.DE, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/DE_T1DvsCTL_alpha.rds")

## Beta
dat_T1DvsCTL_beta <- subset(dat_T1DvsCTL, subset = cell_type == "Beta")
# saveRDS(dat_T1DvsCTL_beta, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/T1DvsCTL_beta.rds")
dat_T1DvsCTL_beta$cell_type <- droplevels(dat_T1DvsCTL_beta$cell_type)
dat_T1DvsCTL_beta$cell_type.disease_state <- paste(dat_T1DvsCTL_beta$cell_type, dat_T1DvsCTL_beta$disease_state, sep = "_")
Idents(dat_T1DvsCTL_beta) <- "cell_type.disease_state"
dat_T1DvsCTL_beta.DE <- FindMarkers(dat_T1DvsCTL_beta, ident.1 = "Beta_T1D", ident.2 = "Beta_Control", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsCTL_beta.DE, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/DE_T1DvsCTL_beta.rds")

## Acinar
dat_T1DvsCTL_acinar <- subset(dat_T1DvsCTL, subset = cell_type == "Acinar")
# saveRDS(dat_T1DvsCTL_acinar, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/T1DvsCTL_acinar.rds")
dat_T1DvsCTL_acinar$cell_type <- droplevels(dat_T1DvsCTL_acinar$cell_type)
dat_T1DvsCTL_acinar$cell_type.disease_state <- paste(dat_T1DvsCTL_acinar$cell_type, dat_T1DvsCTL_acinar$disease_state, sep = "_")
Idents(dat_T1DvsCTL_acinar) <- "cell_type.disease_state"
dat_T1DvsCTL_acinar.DE <- FindMarkers(dat_T1DvsCTL_acinar, ident.1 = "Acinar_T1D", ident.2 = "Acinar_Control",  min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsCTL_acinar.DE, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/DE_T1DvsCTL_acinar.rds")

## Ductal
dat_T1DvsCTL_ductal <- subset(dat_T1DvsCTL, subset = cell_type == "Ductal")
# saveRDS(dat_T1DvsCTL_ductal, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/T1DvsCTL_ductal.rds")
dat_T1DvsCTL_ductal$cell_type <- droplevels(dat_T1DvsCTL_ductal$cell_type)
dat_T1DvsCTL_ductal$cell_type.disease_state <- paste(dat_T1DvsCTL_ductal$cell_type, dat_T1DvsCTL_ductal$disease_state, sep = "_")
Idents(dat_T1DvsCTL_ductal) <- "cell_type.disease_state"
dat_T1DvsCTL_ductal.DE <- FindMarkers(dat_T1DvsCTL_ductal, ident.1 = "Ductal_T1D", ident.2 = "Ductal_Control",  min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsCTL_ductal.DE, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/DE_T1DvsCTL_ductal.rds")


## Immune
dat_T1DvsCTL_immune <- subset(dat_T1DvsCTL, subset = cell_type == "Immune")
# saveRDS(dat_T1DvsCTL_immune, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/T1DvsCTL_immune.rds")
dat_T1DvsCTL_immune$cell_type <- droplevels(dat_T1DvsCTL_immune$cell_type)
dat_T1DvsCTL_immune$cell_type.disease_state <- paste(dat_T1DvsCTL_immune$cell_type, dat_T1DvsCTL_immune$disease_state, sep = "_")
Idents(dat_T1DvsCTL_immune) <- "cell_type.disease_state"
dat_T1DvsCTL_immune.DE <- FindMarkers(dat_T1DvsCTL_immune, ident.1 = "Immune_T1D", ident.2 = "Immune_Control",  min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsCTL_immune.DE, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/DE_T1DvsCTL_immune.rds")

## Delta
dat_T1DvsCTL_delta <- subset(dat_T1DvsCTL, subset = cell_type == "Delta")
# saveRDS(dat_T1DvsCTL_delta, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/T1DvsCTL_delta.rds")
dat_T1DvsCTL_delta$cell_type <- droplevels(dat_T1DvsCTL_delta$cell_type)
dat_T1DvsCTL_delta$cell_type.disease_state <- paste(dat_T1DvsCTL_delta$cell_type, dat_T1DvsCTL_delta$disease_state, sep = "_")
Idents(dat_T1DvsCTL_delta) <- "cell_type.disease_state"
dat_T1DvsCTL_delta.DE <- FindMarkers(dat_T1DvsCTL_delta, ident.1 = "Delta_T1D", ident.2 = "Delta_Control",  min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsCTL_delta.DE, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/DE_T1DvsCTL_delta.rds")

## Endothelial
dat_T1DvsCTL_endothelial <- subset(dat_T1DvsCTL, subset = cell_type == "Endothelial")
# saveRDS(dat_T1DvsCTL_endothelial, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/T1DvsCTL_endothelial.rds")
dat_T1DvsCTL_endothelial$cell_type <- droplevels(dat_T1DvsCTL_endothelial$cell_type)
dat_T1DvsCTL_endothelial$cell_type.disease_state <- paste(dat_T1DvsCTL_endothelial$cell_type, dat_T1DvsCTL_endothelial$disease_state, sep = "_")
Idents(dat_T1DvsCTL_endothelial) <- "cell_type.disease_state"
dat_T1DvsCTL_endothelial.DE <- FindMarkers(dat_T1DvsCTL_endothelial, ident.1 = "Endothelial_T1D", ident.2 = "Endothelial_Control", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsCTL_endothelial.DE, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/DE_T1DvsCTL_endothelial.rds")

## Stellates_Mesenchymal
dat_T1DvsCTL_stellates <- subset(dat_T1DvsCTL, subset = cell_type == "Stellates_Mesenchymal")
# saveRDS(dat_T1DvsCTL_stellates, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/T1DvsCTL_stellates.rds")
dat_T1DvsCTL_stellates$cell_type <- droplevels(dat_T1DvsCTL_stellates$cell_type)
dat_T1DvsCTL_stellates$cell_type.disease_state <- paste(dat_T1DvsCTL_stellates$cell_type, dat_T1DvsCTL_stellates$disease_state, sep = "_")
Idents(dat_T1DvsCTL_stellates) <- "cell_type.disease_state"
dat_T1DvsCTL_stellates.DE <- FindMarkers(dat_T1DvsCTL_stellates, ident.1 = "Stellates_Mesenchymal_T1D", ident.2 = "Stellates_Mesenchymal_Control", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsCTL_stellates.DE, "/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/DE/DE_T1DvsCTL_stellates.rds")
