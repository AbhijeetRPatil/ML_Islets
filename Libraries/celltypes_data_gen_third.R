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
  
  ## Alpha
  alpha <- subset(dat, subset = cell_type == "Alpha")
  ## Beta
  beta <- subset(dat, subset = cell_type == "Beta")
  ## Acinar
  acinar <- subset(dat, subset = cell_type == "Acinar")
  ## Ductal
  ductal <- subset(dat, subset = cell_type == "Ductal")
  ## Immune
  immune <- subset(dat, subset = cell_type == "Immune")
  ## Delta
  delta <- subset(dat, subset = cell_type == "Delta")
  ## Endothelial
  endothelial <- subset(dat, subset = cell_type == "Endothelial")
  ## Stellates_Mesenchymal
  stellates <- subset(dat, subset = cell_type == "Stellates_Mesenchymal")
  
  
  return(list("alpha"=alpha, "beta"=beta, "acinar"=acinar, "ductal"=ductal, 
              "immune"=immune, "delta"=delta, "endothelial"=endothelial,
              "stellates"=stellates))
}