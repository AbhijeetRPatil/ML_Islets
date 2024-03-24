options(future.globals.maxSize = 1000000 * 1024^2)
library("Seurat")
dat <- readRDS("panc.rds")
dat_AABvsCTL <- subset(dat, subset = disease_state != "T1D")
# saveRDS(dat_AABvsCTL, "AABvsCTL_Overall.rds")
Idents(dat_AABvsCTL) <- "disease_state"
dat_AABvsCTL.DE <- FindMarkers(dat_AABvsCTL, ident.1 = "AAB", ident.2 = "Control", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_AABvsCTL.DE, "DE_AABvsCTL.rds")

## Alpha
dat_AABvsCTL_alpha <- subset(dat_AABvsCTL, subset = cell_type == "Alpha")
# saveRDS(dat_AABvsCTL_alpha, "AABvsCTL_alpha.rds")
dat_AABvsCTL_alpha$cell_type <- droplevels(dat_AABvsCTL_alpha$cell_type)
dat_AABvsCTL_alpha$cell_type.disease_state <- paste(dat_AABvsCTL_alpha$cell_type, dat_AABvsCTL_alpha$disease_state, sep = "_")
Idents(dat_AABvsCTL_alpha) <- "cell_type.disease_state"
dat_AABvsCTL_alpha.DE <- FindMarkers(dat_AABvsCTL_alpha, ident.1 = "Alpha_AAB", ident.2 = "Alpha_Control", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_AABvsCTL_alpha.DE, "DE_AABvsCTL_alpha.rds")

## Beta
dat_AABvsCTL_beta <- subset(dat_AABvsCTL, subset = cell_type == "Beta")
# saveRDS(dat_AABvsCTL_beta, "AABvsCTL_beta.rds")
dat_AABvsCTL_beta$cell_type <- droplevels(dat_AABvsCTL_beta$cell_type)
dat_AABvsCTL_beta$cell_type.disease_state <- paste(dat_AABvsCTL_beta$cell_type, dat_AABvsCTL_beta$disease_state, sep = "_")
Idents(dat_AABvsCTL_beta) <- "cell_type.disease_state"
dat_AABvsCTL_beta.DE <- FindMarkers(dat_AABvsCTL_beta, ident.1 = "Beta_AAB", ident.2 = "Beta_Control", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_AABvsCTL_beta.DE, "DE_AABvsCTL_beta.rds")

## Acinar
dat_AABvsCTL_acinar <- subset(dat_AABvsCTL, subset = cell_type == "Acinar")
# saveRDS(dat_AABvsCTL_acinar, "AABvsCTL_acinar.rds")
dat_AABvsCTL_acinar$cell_type <- droplevels(dat_AABvsCTL_acinar$cell_type)
dat_AABvsCTL_acinar$cell_type.disease_state <- paste(dat_AABvsCTL_acinar$cell_type, dat_AABvsCTL_acinar$disease_state, sep = "_")
Idents(dat_AABvsCTL_acinar) <- "cell_type.disease_state"
dat_AABvsCTL_acinar.DE <- FindMarkers(dat_AABvsCTL_acinar, ident.1 = "Acinar_AAB", ident.2 = "Acinar_Control",  min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_AABvsCTL_acinar.DE, "DE_AABvsCTL_acinar.rds")

## Ductal
dat_AABvsCTL_ductal <- subset(dat_AABvsCTL, subset = cell_type == "Ductal")
# saveRDS(dat_AABvsCTL_ductal, "AABvsCTL_ductal.rds")
dat_AABvsCTL_ductal$cell_type <- droplevels(dat_AABvsCTL_ductal$cell_type)
dat_AABvsCTL_ductal$cell_type.disease_state <- paste(dat_AABvsCTL_ductal$cell_type, dat_AABvsCTL_ductal$disease_state, sep = "_")
Idents(dat_AABvsCTL_ductal) <- "cell_type.disease_state"
dat_AABvsCTL_ductal.DE <- FindMarkers(dat_AABvsCTL_ductal, ident.1 = "Ductal_AAB", ident.2 = "Ductal_Control",  min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_AABvsCTL_ductal.DE, "DE_AABvsCTL_ductal.rds")


## Immune
dat_AABvsCTL_immune <- subset(dat_AABvsCTL, subset = cell_type == "Immune")
# saveRDS(dat_AABvsCTL_immune, "AABvsCTL_immune.rds")
dat_AABvsCTL_immune$cell_type <- droplevels(dat_AABvsCTL_immune$cell_type)
dat_AABvsCTL_immune$cell_type.disease_state <- paste(dat_AABvsCTL_immune$cell_type, dat_AABvsCTL_immune$disease_state, sep = "_")
Idents(dat_AABvsCTL_immune) <- "cell_type.disease_state"
dat_AABvsCTL_immune.DE <- FindMarkers(dat_AABvsCTL_immune, ident.1 = "Immune_AAB", ident.2 = "Immune_Control",  min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_AABvsCTL_immune.DE, "DE_AABvsCTL_immune.rds")

## Delta
dat_AABvsCTL_delta <- subset(dat_AABvsCTL, subset = cell_type == "Delta")
# saveRDS(dat_AABvsCTL_delta, "AABvsCTL_delta.rds")
dat_AABvsCTL_delta$cell_type <- droplevels(dat_AABvsCTL_delta$cell_type)
dat_AABvsCTL_delta$cell_type.disease_state <- paste(dat_AABvsCTL_delta$cell_type, dat_AABvsCTL_delta$disease_state, sep = "_")
Idents(dat_AABvsCTL_delta) <- "cell_type.disease_state"
dat_AABvsCTL_delta.DE <- FindMarkers(dat_AABvsCTL_delta, ident.1 = "Delta_AAB", ident.2 = "Delta_Control",  min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_AABvsCTL_delta.DE, "DE_AABvsCTL_delta.rds")

## Endothelial
dat_AABvsCTL_endothelial <- subset(dat_AABvsCTL, subset = cell_type == "Endothelial")
# saveRDS(dat_AABvsCTL_endothelial, "AABvsCTL_endothelial.rds")
dat_AABvsCTL_endothelial$cell_type <- droplevels(dat_AABvsCTL_endothelial$cell_type)
dat_AABvsCTL_endothelial$cell_type.disease_state <- paste(dat_AABvsCTL_endothelial$cell_type, dat_AABvsCTL_endothelial$disease_state, sep = "_")
Idents(dat_AABvsCTL_endothelial) <- "cell_type.disease_state"
dat_AABvsCTL_endothelial.DE <- FindMarkers(dat_AABvsCTL_endothelial, ident.1 = "Endothelial_AAB", ident.2 = "Endothelial_Control", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_AABvsCTL_endothelial.DE, "DE_AABvsCTL_endothelial.rds")

## Stellates_Mesenchymal
dat_AABvsCTL_stellates <- subset(dat_AABvsCTL, subset = cell_type == "Stellates_Mesenchymal")
# saveRDS(dat_AABvsCTL_stellates, "AABvsCTL_stellates.rds")
dat_AABvsCTL_stellates$cell_type <- droplevels(dat_AABvsCTL_stellates$cell_type)
dat_AABvsCTL_stellates$cell_type.disease_state <- paste(dat_AABvsCTL_stellates$cell_type, dat_AABvsCTL_stellates$disease_state, sep = "_")
Idents(dat_AABvsCTL_stellates) <- "cell_type.disease_state"
dat_AABvsCTL_stellates.DE <- FindMarkers(dat_AABvsCTL_stellates, ident.1 = "Stellates_Mesenchymal_AAB", ident.2 = "Stellates_Mesenchymal_Control", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_AABvsCTL_stellates.DE, "DE_AABvsCTL_stellates.rds")
