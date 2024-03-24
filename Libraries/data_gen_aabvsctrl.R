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
  
  data <-GetAssayData(dat_AABvsCtrl, assay="SCT",slot="data")
  data <- t(data)
  y <- dat_AABvsCtrl$class_label_num
  
  return(list(data, y))
}