options(future.globals.maxSize = 1000000 * 1024^2)
library("Seurat")
dat <- readRDS("panc.rds")
dat_T1DvsAAB <- subset(dat, subset = disease_state != "Control")
# saveRDS(dat_T1DvsAAB, "T1DvsAAB_Overall.rds")
Idents(dat_T1DvsAAB) <- "disease_state"
dat_T1DvsAAB.DE <- FindMarkers(dat_T1DvsAAB, ident.1 = "T1D", ident.2 = "AAB", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsAAB.DE, "DE_T1DvsAAB.rds")

## Alpha
dat_T1DvsAAB_alpha <- subset(dat_T1DvsAAB, subset = cell_type == "Alpha")
# saveRDS(dat_T1DvsAAB_alpha, "T1DvsAAB_alpha.rds")
dat_T1DvsAAB_alpha$cell_type <- droplevels(dat_T1DvsAAB_alpha$cell_type)
dat_T1DvsAAB_alpha$cell_type.disease_state <- paste(dat_T1DvsAAB_alpha$cell_type, dat_T1DvsAAB_alpha$disease_state, sep = "_")
Idents(dat_T1DvsAAB_alpha) <- "cell_type.disease_state"
dat_T1DvsAAB_alpha.DE <- FindMarkers(dat_T1DvsAAB_alpha, ident.1 = "Alpha_T1D", ident.2 = "Alpha_AAB", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsAAB_alpha.DE, "DE_T1DvsAAB_alpha.rds")

## Beta
dat_T1DvsAAB_beta <- subset(dat_T1DvsAAB, subset = cell_type == "Beta")
# saveRDS(dat_T1DvsAAB_beta, "T1DvsAAB_beta.rds")
dat_T1DvsAAB_beta$cell_type <- droplevels(dat_T1DvsAAB_beta$cell_type)
dat_T1DvsAAB_beta$cell_type.disease_state <- paste(dat_T1DvsAAB_beta$cell_type, dat_T1DvsAAB_beta$disease_state, sep = "_")
Idents(dat_T1DvsAAB_beta) <- "cell_type.disease_state"
dat_T1DvsAAB_beta.DE <- FindMarkers(dat_T1DvsAAB_beta, ident.1 = "Beta_T1D", ident.2 = "Beta_AAB", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsAAB_beta.DE, "DE_T1DvsAAB_beta.rds")

## Acinar
dat_T1DvsAAB_acinar <- subset(dat_T1DvsAAB, subset = cell_type == "Acinar")
# saveRDS(dat_T1DvsAAB_acinar, "T1DvsAAB_acinar.rds")
dat_T1DvsAAB_acinar$cell_type <- droplevels(dat_T1DvsAAB_acinar$cell_type)
dat_T1DvsAAB_acinar$cell_type.disease_state <- paste(dat_T1DvsAAB_acinar$cell_type, dat_T1DvsAAB_acinar$disease_state, sep = "_")
Idents(dat_T1DvsAAB_acinar) <- "cell_type.disease_state"
dat_T1DvsAAB_acinar.DE <- FindMarkers(dat_T1DvsAAB_acinar, ident.1 = "Acinar_T1D", ident.2 = "Acinar_AAB",  min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsAAB_acinar.DE, "DE_T1DvsAAB_acinar.rds")

## Ductal
dat_T1DvsAAB_ductal <- subset(dat_T1DvsAAB, subset = cell_type == "Ductal")
# saveRDS(dat_T1DvsAAB_ductal, "T1DvsAAB_ductal.rds")
dat_T1DvsAAB_ductal$cell_type <- droplevels(dat_T1DvsAAB_ductal$cell_type)
dat_T1DvsAAB_ductal$cell_type.disease_state <- paste(dat_T1DvsAAB_ductal$cell_type, dat_T1DvsAAB_ductal$disease_state, sep = "_")
Idents(dat_T1DvsAAB_ductal) <- "cell_type.disease_state"
dat_T1DvsAAB_ductal.DE <- FindMarkers(dat_T1DvsAAB_ductal, ident.1 = "Ductal_T1D", ident.2 = "Ductal_AAB",  min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsAAB_ductal.DE, "DE_T1DvsAAB_ductal.rds")


## Immune
dat_T1DvsAAB_immune <- subset(dat_T1DvsAAB, subset = cell_type == "Immune")
# saveRDS(dat_T1DvsAAB_immune, "T1DvsAAB_immune.rds")
dat_T1DvsAAB_immune$cell_type <- droplevels(dat_T1DvsAAB_immune$cell_type)
dat_T1DvsAAB_immune$cell_type.disease_state <- paste(dat_T1DvsAAB_immune$cell_type, dat_T1DvsAAB_immune$disease_state, sep = "_")
Idents(dat_T1DvsAAB_immune) <- "cell_type.disease_state"
dat_T1DvsAAB_immune.DE <- FindMarkers(dat_T1DvsAAB_immune, ident.1 = "Immune_T1D", ident.2 = "Immune_AAB",  min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsAAB_immune.DE, "DE_T1DvsAAB_immune.rds")

## Delta
dat_T1DvsAAB_delta <- subset(dat_T1DvsAAB, subset = cell_type == "Delta")
# saveRDS(dat_T1DvsAAB_delta, "T1DvsAAB_delta.rds")
dat_T1DvsAAB_delta$cell_type <- droplevels(dat_T1DvsAAB_delta$cell_type)
dat_T1DvsAAB_delta$cell_type.disease_state <- paste(dat_T1DvsAAB_delta$cell_type, dat_T1DvsAAB_delta$disease_state, sep = "_")
Idents(dat_T1DvsAAB_delta) <- "cell_type.disease_state"
dat_T1DvsAAB_delta.DE <- FindMarkers(dat_T1DvsAAB_delta, ident.1 = "Delta_T1D", ident.2 = "Delta_AAB",  min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsAAB_delta.DE, "DE_T1DvsAAB_delta.rds")

## Endothelial
dat_T1DvsAAB_endothelial <- subset(dat_T1DvsAAB, subset = cell_type == "Endothelial")
# saveRDS(dat_T1DvsAAB_endothelial, "T1DvsAAB_endothelial.rds")
dat_T1DvsAAB_endothelial$cell_type <- droplevels(dat_T1DvsAAB_endothelial$cell_type)
dat_T1DvsAAB_endothelial$cell_type.disease_state <- paste(dat_T1DvsAAB_endothelial$cell_type, dat_T1DvsAAB_endothelial$disease_state, sep = "_")
Idents(dat_T1DvsAAB_endothelial) <- "cell_type.disease_state"
dat_T1DvsAAB_endothelial.DE <- FindMarkers(dat_T1DvsAAB_endothelial, ident.1 = "Endothelial_T1D", ident.2 = "Endothelial_AAB", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsAAB_endothelial.DE, "DE_T1DvsAAB_endothelial.rds")

## Stellates_Mesenchymal
dat_T1DvsAAB_stellates <- subset(dat_T1DvsAAB, subset = cell_type == "Stellates_Mesenchymal")
# saveRDS(dat_T1DvsAAB_stellates, "T1DvsAAB_stellates.rds")
dat_T1DvsAAB_stellates$cell_type <- droplevels(dat_T1DvsAAB_stellates$cell_type)
dat_T1DvsAAB_stellates$cell_type.disease_state <- paste(dat_T1DvsAAB_stellates$cell_type, dat_T1DvsAAB_stellates$disease_state, sep = "_")
Idents(dat_T1DvsAAB_stellates) <- "cell_type.disease_state"
dat_T1DvsAAB_stellates.DE <- FindMarkers(dat_T1DvsAAB_stellates, ident.1 = "Stellates_Mesenchymal_T1D", ident.2 = "Stellates_Mesenchymal_AAB", min.pct = 0, logfc.threshold = 0, verbose = FALSE)
saveRDS(dat_T1DvsAAB_stellates.DE, "DE_T1DvsAAB_stellates.rds")
