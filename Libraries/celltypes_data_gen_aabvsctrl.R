data.gen <- function()
{
  ## Read data
  dat <- readRDS("panc.rds")
  dat@meta.data$DiseaseState <- as.factor(dat@meta.data$disease_state)
  
  ## Create response variable
  df <- dat@meta.data %>% 
    mutate(class_label = recode(disease_state, 
                                "Control" = 0, 
                                "AAB" = 1,
                                "T1D" = 2)) 
  dat[["class_label_num"]] <- as.numeric(df$class_label)
  
  ## Split data
  dat_AABvsCtrl <- subset(dat, subset = disease_state != "T1D")
  
  ## Alpha
  # if(celltypeX=="Alpha")
  dat_AABvsCtrl_alpha <- subset(dat_AABvsCtrl, subset = cell_type == "Alpha")
  ## Beta
  dat_AABvsCtrl_beta <- subset(dat_AABvsCtrl, subset = cell_type == "Beta")
  ## Acinar
  dat_AABvsCtrl_acinar <- subset(dat_AABvsCtrl, subset = cell_type == "Acinar")
  ## Ductal
  dat_AABvsCtrl_ductal <- subset(dat_AABvsCtrl, subset = cell_type == "Ductal")
  ## Immune
  dat_AABvsCtrl_immune <- subset(dat_AABvsCtrl, subset = cell_type == "Immune")
  ## Delta
  dat_AABvsCtrl_delta <- subset(dat_AABvsCtrl, subset = cell_type == "Delta")
  ## Endothelial
  dat_AABvsCtrl_endothelial <- subset(dat_AABvsCtrl, subset = cell_type == "Endothelial")
  ## Stellates_Mesenchymal
  dat_AABvsCtrl_stellates <- subset(dat_AABvsCtrl, subset = cell_type == "Stellates_Mesenchymal")
  
  ## Get the data out of slot
  
  ## Alpha
  data_AABvsCtrl_alpha <-GetAssayData(dat_AABvsCtrl_alpha, assay="SCT",slot="data")
  data_AABvsCtrl_alpha <- t(data_AABvsCtrl_alpha)
  y_alpha <- dat_AABvsCtrl_alpha$class_label_num
  
  alpha <- list(data_AABvsCtrl_alpha, y_alpha)
  
  ## Beta
  data_AABvsCtrl_beta <-GetAssayData(dat_AABvsCtrl_beta, assay="SCT",slot="data")
  data_AABvsCtrl_beta <- t(data_AABvsCtrl_beta)
  y_beta <- dat_AABvsCtrl_beta$class_label_num
  
  beta <- list(data_AABvsCtrl_beta, y_beta)
  
  ## Acinar
  data_AABvsCtrl_acinar <-GetAssayData(dat_AABvsCtrl_acinar, assay="SCT",slot="data")
  data_AABvsCtrl_acinar <- t(data_AABvsCtrl_acinar)
  y_acinar <- dat_AABvsCtrl_acinar$class_label_num
  
  acinar <- list(data_AABvsCtrl_acinar, y_acinar)
  
  ## Ductal
  data_AABvsCtrl_ductal <-GetAssayData(dat_AABvsCtrl_ductal, assay="SCT",slot="data")
  data_AABvsCtrl_ductal <- t(data_AABvsCtrl_ductal)
  y_ductal <- dat_AABvsCtrl_ductal$class_label_num
  
  ductal <- list(data_AABvsCtrl_ductal, y_ductal)
  
  ## Immune
  data_AABvsCtrl_immune <-GetAssayData(dat_AABvsCtrl_immune, assay="SCT",slot="data")
  data_AABvsCtrl_immune <- t(data_AABvsCtrl_immune)
  y_immune <- dat_AABvsCtrl_immune$class_label_num
  
  immune <- list(data_AABvsCtrl_immune, y_immune)
  
  ## Delta
  data_AABvsCtrl_delta <-GetAssayData(dat_AABvsCtrl_delta, assay="SCT",slot="data")
  data_AABvsCtrl_delta <- t(data_AABvsCtrl_delta)
  y_delta <- dat_AABvsCtrl_delta$class_label_num
  
  delta <- list(data_AABvsCtrl_delta, y_delta)
  
  ## Endothelial
  data_AABvsCtrl_endothelial <-GetAssayData(dat_AABvsCtrl_endothelial, assay="SCT",slot="data")
  data_AABvsCtrl_endothelial <- t(data_AABvsCtrl_endothelial)
  y_endothelial <- dat_AABvsCtrl_endothelial$class_label_num
  
  endothelial <- list(data_AABvsCtrl_endothelial, y_endothelial)
  
  ## Stellates
  data_AABvsCtrl_stellates <-GetAssayData(dat_AABvsCtrl_stellates, assay="SCT",slot="data")
  data_AABvsCtrl_stellates <- t(data_AABvsCtrl_stellates)
  y_stellates <- dat_AABvsCtrl_stellates$class_label_num
  
  stellates <- list(data_AABvsCtrl_stellates, y_stellates)
  
  return(list("alpha"=alpha, "beta"=beta, "acinar"=acinar, "ductal"=ductal, 
              "immune"=immune, "delta"=delta, "endothelial"=endothelial,
              "stellates"=stellates))
}
