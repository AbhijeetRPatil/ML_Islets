data.gen <- function()
{
  ## Read data
  dat <- readRDS("panc.rds")
  ## Use RNA assay
  DefaultAssay(dat) <- "RNA"
  dat <- NormalizeData(dat)

  dat@meta.data$DiseaseState <- as.factor(dat@meta.data$disease_state)
  
  ## Create response variable
  df <- dat@meta.data %>% 
    mutate(class_label = recode(disease_state, 
                                "Control" = 0, 
                                "AAB" = 1,
                                "T1D" = 2)) 
  dat[["class_label_num"]] <- as.numeric(df$class_label)
  
  ## Split data
  dat_T1DvsAAB <- subset(dat, subset = disease_state != "Control")
  
  ## Alpha
  # if(celltypeX=="Alpha")
  dat_T1DvsAAB_alpha <- subset(dat_T1DvsAAB, subset = cell_type == "Alpha")
  ## Beta
  dat_T1DvsAAB_beta <- subset(dat_T1DvsAAB, subset = cell_type == "Beta")
  ## Acinar
  dat_T1DvsAAB_acinar <- subset(dat_T1DvsAAB, subset = cell_type == "Acinar")
  ## Ductal
  dat_T1DvsAAB_ductal <- subset(dat_T1DvsAAB, subset = cell_type == "Ductal")
  ## Immune
  dat_T1DvsAAB_immune <- subset(dat_T1DvsAAB, subset = cell_type == "Immune")
  ## Delta
  dat_T1DvsAAB_delta <- subset(dat_T1DvsAAB, subset = cell_type == "Delta")
  ## Endothelial
  dat_T1DvsAAB_endothelial <- subset(dat_T1DvsAAB, subset = cell_type == "Endothelial")
  ## Stellates_Mesenchymal
  dat_T1DvsAAB_stellates <- subset(dat_T1DvsAAB, subset = cell_type == "Stellates_Mesenchymal")
  
  ## Get the data out of slot
  
  ## Alpha
  data_T1DvsAAB_alpha <-GetAssayData(dat_T1DvsAAB_alpha, assay="RNA",slot="data")
  data_T1DvsAAB_alpha <- t(data_T1DvsAAB_alpha)
  y_alpha <- dat_T1DvsAAB_alpha$class_label_num
  y_alpha <- ifelse(y_alpha == 1, 0, y_alpha)
  y_alpha <- ifelse(y_alpha == 2, 1, y_alpha)
  
  alpha <- list(data_T1DvsAAB_alpha, y_alpha)
  
  ## Beta
  data_T1DvsAAB_beta <-GetAssayData(dat_T1DvsAAB_beta, assay="RNA",slot="data")
  data_T1DvsAAB_beta <- t(data_T1DvsAAB_beta)
  y_beta <- dat_T1DvsAAB_beta$class_label_num
  y_beta <- ifelse(y_beta == 1, 0, y_beta)
  y_beta <- ifelse(y_beta == 2, 1, y_beta)
  
  beta <- list(data_T1DvsAAB_beta, y_beta)
  
  ## Acinar
  data_T1DvsAAB_acinar <-GetAssayData(dat_T1DvsAAB_acinar, assay="RNA",slot="data")
  data_T1DvsAAB_acinar <- t(data_T1DvsAAB_acinar)
  y_acinar <- dat_T1DvsAAB_acinar$class_label_num
  y_acinar <- ifelse(y_acinar == 1, 0, y_acinar)
  y_acinar <- ifelse(y_acinar == 2, 1, y_acinar)
  
  acinar <- list(data_T1DvsAAB_acinar, y_acinar)
  
  ## Ductal
  data_T1DvsAAB_ductal <-GetAssayData(dat_T1DvsAAB_ductal, assay="RNA",slot="data")
  data_T1DvsAAB_ductal <- t(data_T1DvsAAB_ductal)
  y_ductal <- dat_T1DvsAAB_ductal$class_label_num
  y_ductal <- ifelse(y_ductal == 1, 0, y_ductal)
  y_ductal <- ifelse(y_ductal == 2, 1, y_ductal)
  
  ductal <- list(data_T1DvsAAB_ductal, y_ductal)
  
  ## Immune
  data_T1DvsAAB_immune <-GetAssayData(dat_T1DvsAAB_immune, assay="RNA",slot="data")
  data_T1DvsAAB_immune <- t(data_T1DvsAAB_immune)
  y_immune <- dat_T1DvsAAB_immune$class_label_num
  y_immune <- ifelse(y_immune == 1, 0, y_immune)
  y_immune <- ifelse(y_immune == 2, 1, y_immune)
  
  immune <- list(data_T1DvsAAB_immune, y_immune)
  
  ## Delta
  data_T1DvsAAB_delta <-GetAssayData(dat_T1DvsAAB_delta, assay="RNA",slot="data")
  data_T1DvsAAB_delta <- t(data_T1DvsAAB_delta)
  y_delta <- dat_T1DvsAAB_delta$class_label_num
  y_delta <- ifelse(y_delta == 1, 0, y_delta)
  y_delta <- ifelse(y_delta == 2, 1, y_delta)
  
  delta <- list(data_T1DvsAAB_delta, y_delta)
  
  ## Endothelial
  data_T1DvsAAB_endothelial <-GetAssayData(dat_T1DvsAAB_endothelial, assay="RNA",slot="data")
  data_T1DvsAAB_endothelial <- t(data_T1DvsAAB_endothelial)
  y_endothelial <- dat_T1DvsAAB_endothelial$class_label_num
  y_endothelial <- ifelse(y_endothelial == 1, 0, y_endothelial)
  y_endothelial <- ifelse(y_endothelial == 2, 1, y_endothelial)
  
  endothelial <- list(data_T1DvsAAB_endothelial, y_endothelial)
  
  ## Stellates
  data_T1DvsAAB_stellates <-GetAssayData(dat_T1DvsAAB_stellates, assay="RNA",slot="data")
  data_T1DvsAAB_stellates <- t(data_T1DvsAAB_stellates)
  y_stellates <- dat_T1DvsAAB_stellates$class_label_num
  y_stellates <- ifelse(y_stellates == 1, 0, y_stellates)
  y_stellates <- ifelse(y_stellates == 2, 1, y_stellates)
  
  stellates <- list(data_T1DvsAAB_stellates, y_stellates)
  
  return(list("alpha"=alpha, "beta"=beta, "acinar"=acinar, "ductal"=ductal, 
              "immune"=immune, "delta"=delta, "endothelial"=endothelial,
              "stellates"=stellates))
}
