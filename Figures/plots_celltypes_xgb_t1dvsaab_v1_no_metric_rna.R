rm(list=ls())
library("xgboost")  # the main algorithm
library("caret")    # for the confusionmatrix() function (also needs e1071 package)
library("dplyr")    # for some data preperation
library("Ckmeans.1d.dp") # for xgb.ggplot.importance
library("Matrix")
library("Seurat")
library("parallel")
library("pROC")
library("viridis")
library("ggpubr")
## current working directory 
setwd("C:/Users/abhij/Desktop/Spring 2022/ML paper/LOOCV_New/")

memory.limit(size=10000000)
xtheme <- theme_bw()+ theme(plot.title = element_text(face = "bold" ,hjust = 0.5, size= 10)
                            ,axis.text.y = element_text(face = "bold",angle = 0, size = 10, hjust=1)
                            ,axis.title.y = element_text(face = "bold", size = rel(0.8))
                            ,axis.text.x = element_text(face = "bold",angle = 0, size = 10)
                            ,axis.title.x = element_text(face = "bold", size = rel(0.8))
                            ,axis.ticks.x=element_blank(), strip.text = element_text(size=10))

gen.data <- function(obj)
{
  names(obj[[3]]) <- c("HPAP038", "HPAP043", "HPAP045", "HPAP049", "HPAP050", "HPAP072", 
                       "HPAP024", "HPAP092", "HPAP107", "HPAP029")
  
  names(obj[[2]]) <- c("HPAP038", "HPAP043", "HPAP045", "HPAP049", "HPAP050", "HPAP072",
                       "HPAP024", "HPAP092", "HPAP107", "HPAP029")
  df <- list()
  obj[[3]] <- obj[[3]][1:10]
  obj[[2]] <- obj[[2]][1:10]
  obj[[1]] <- obj[[1]][1:10]
  for (i in names(obj[[2]])) 
  {
    # print(i)
    print(obj[[3]][[i]])
    df[[i]]$Accuracy <- obj[[3]][[i]]
    Freq <- obj[[2]][[i]]
    df[[i]]$T1D <- round(sum(Freq)/length(Freq),5)
    df[[i]]$AAB <- round(1-sum(Freq)/length(Freq),5)

    # df[[i]]$Freq <- round(df[[i]]$Freq/sum(df[[i]]$Freq)*100,1)
    # df[[i]]$Name <- rep(i,nrow(df[[i]]))
  }

  final_df <- do.call(rbind, df)
  final_df <- as.data.frame(final_df)
  final_df$T1D <- unlist(final_df$T1D)
  final_df$AAB <- unlist(final_df$AAB)
  final_df$Accuracy <- NULL
  final_df$Donor_Ids <- rownames(final_df)
  rownames(final_df) <- NULL
  library(dplyr)
  final_df <- final_df %>%
    relocate(Donor_Ids)
  library(reshape2)
  final_df = melt(final_df, id = c("Donor_Ids"))
  final_df <- as.data.frame(final_df)
  colnames(final_df)[2] <- "Condition"
  colnames(final_df)[1] <- "Name"
  colnames(final_df)[3] <- "Freq"
  final_df$Condition <- relevel(final_df$Condition, ref = 'AAB')
  return(final_df)
}


alpha_xgb_t1dvsaab_04162023 <- readRDS("./cv_AAB_preds/t1dvsaab/metric/cv/alpha_xgb_t1dvsaab_03162024_cv_rna.rds")
alpha.df <- gen.data(alpha_xgb_t1dvsaab_04162023)

beta_xgb_t1dvsaab_04162023 <- readRDS("./cv_AAB_preds/t1dvsaab/metric/cv/beta_xgb_t1dvsaab_03162024_cv_rna.rds")
beta.df <- gen.data(beta_xgb_t1dvsaab_04162023)

acinar_xgb_t1dvsaab_04162023 <- readRDS("./cv_AAB_preds/t1dvsaab/metric/cv/acinar_xgb_t1dvsaab_03162024_cv_rna.rds")
acinar.df <- gen.data(acinar_xgb_t1dvsaab_04162023)

delta_xgb_t1dvsaab_04162023 <- readRDS("./cv_AAB_preds/t1dvsaab/metric/cv/delta_xgb_t1dvsaab_03162024_cv_rna.rds")
delta.df <- gen.data(delta_xgb_t1dvsaab_04162023)

ductal_xgb_t1dvsaab_04162023 <- readRDS("./cv_AAB_preds/t1dvsaab/metric/cv/ductal_xgb_t1dvsaab_03162024_cv_rna.rds")
ductal.df <- gen.data(ductal_xgb_t1dvsaab_04162023)

endothelial_xgb_t1dvsaab_04162023 <- readRDS("./cv_AAB_preds/t1dvsaab/metric/cv/endothelial_xgb_t1dvsaab_03162024_cv_rna.rds")
endothelial.df <- gen.data(endothelial_xgb_t1dvsaab_04162023)

immune_xgb_t1dvsaab_04162023 <- readRDS("./cv_AAB_preds/t1dvsaab/metric/cv/immune_xgb_t1dvsaab_03162024_cv_rna.rds")
immune.df <- gen.data(immune_xgb_t1dvsaab_04162023)

stellates_xgb_t1dvsaab_04162023 <- readRDS("./cv_AAB_preds/t1dvsaab/metric/cv/stellates_xgb_t1dvsaab_03162024_cv_rna.rds")
stellates.df <- gen.data(stellates_xgb_t1dvsaab_04162023)

# alpha.df <- alpha.df %>% group_by(Name, Var1) %>%
#   summarise(Nb = n()) %>%
#   mutate(C = sum(Nb)) %>%
#   mutate(Percent=Nb/C*100)
plot.stack <- function(df, cell)
{
  p1 <- ggplot(df, aes(fill=Condition, y=Freq, x=Name)) + ggtitle(cell) + geom_bar(position = "fill", stat = "identity") + 
    scale_fill_viridis(discrete = T)+ xlab("") + ylab("Percentage of AAb+ cells") + xtheme + 
    theme(legend.position='top', axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))  
  return(p1)
}
acinar.plot <- plot.stack(acinar.df, "Acinar")
alpha.plot <- plot.stack(alpha.df, "Alpha")
beta.plot <- plot.stack(beta.df, "Beta")
delta.plot <- plot.stack(delta.df, "Delta")
ductal.plot <- plot.stack(ductal.df, "Ductal")
endothelial.plot <- plot.stack(endothelial.df, "Endothelial")
immune.plot <- plot.stack(immune.df, "Immune")
stellates.plot <- plot.stack(stellates.df, "Stellates")
# pdf(file = "./stacked_AAB_Preds_Individualdonors.pdf", width = 10, height = 10)
# CombinePlots(list(acinar.plot, alpha.plot, beta.plot,
#                   delta.plot, ductal.plot, 
#                   endothelial.plot, immune.plot, stellates.plot))
# dev.off()

pdf(NULL)
res <- ggarrange(acinar.plot, alpha.plot, beta.plot, delta.plot, ductal.plot, endothelial.plot, 
                 immune.plot, stellates.plot, nrow = 4, ncol = 2, common.legend = T, legend = "bottom")
dev.off()
pdf(file = "./stacked_AAB_Preds_Individualdonors_rna_new.pdf", width = 8, height = 10)
res
dev.off()
