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
  dat_T1DvsAAB <- subset(dat, subset = disease_state != "Control")
  
  data <-GetAssayData(dat_T1DvsAAB, assay="SCT",slot="data")
  data <- t(data)
  y <- dat_T1DvsAAB$class_label_num
  y <- ifelse(y == 1, 0, y)
  y <- ifelse(y == 2, 1, y)
  
  return(list(data, y))
}