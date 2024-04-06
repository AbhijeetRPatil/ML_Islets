##################################################################################
rm(list = ls())
gc()
library(ggplot2)
library(vioplot)
library(ggpubr)
library(gridExtra)

#######################################################################################
setwd("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/")
#######################################################################################
plot_box <- function(ML1, ML2, ML3, ML4, title)
{
  level_order <- c('XGBoost', 'SVM-Linear', 'SVM-Radial', 'Naive Bayes')

  Metrics_num <-  c(round(ML1[,1]*100,2), ML2[,1], ML3[,1], ML4[,1])
  Metrics <- c(rep("Accuracy",130))
  
  df <- data.frame(cbind(Metrics_num, Metrics))
  Condition <- c(rep("XGBoost",100), rep("SVM-Linear",10), rep("SVM-Radial",10), rep("Naive Bayes",10))
  
  df <- cbind(df,Condition)
  df$Metrics_num <- as.numeric(df$Metrics_num)
  
  p <- ggplot(df, aes(x= factor(Condition, level = level_order) , y = Metrics_num, fill = Metrics)) +
    geom_boxplot(outlier.size=0.1, lwd=0.1, fatten=1) + #lwd=0.1, 
    # scale_fill_manual(terrain.colors(3))+
    # scale_fill_brewer(palette="BuPu") +
  
    # scale_fill_discrete()+
    # facet_wrap(~factor(Condition, level = level_order), scale="free_x") + #"free" "fixed" "free_y"
    ggtitle(title) + xlab("") + ylab("Accuracy") + theme_bw() 
  
  xtheme <- theme(plot.title = element_text(hjust = 0.5, size= 12)
                  ,axis.text.y = element_text(angle = 0, size = 12, hjust=1)
                  ,axis.title.y = element_text( size = 12)
                  ,axis.text.x = element_text(angle = 30, size = 12, hjust = 1)
                  ,axis.title.x = element_text( size = rel(0.7))
                  ,axis.ticks.x=element_blank(), strip.text = element_text(size=12), legend.position = "none")
  p <- p + xtheme
  return(p)
}
#########################################
## Overall
## T1D vs Ctrl
xgb_t1dvsctrl <- readRDS("./ML_xgb/xgb_t1dvsctrl.rds")
xgb_t1dvsctrl <- xgb_t1dvsctrl[[2]]

svm_lin_t1dvsctrl <- readRDS("./ML_svm/svm_t1dvsctrl_lin.rds")
svm_lin_t1dvsctrl <- svm_lin_t1dvsctrl[[1]]

svm_rad_t1dvsctrl <- readRDS("./ML_svm/svm_t1dvsctrl_rad.rds")
svm_rad_t1dvsctrl <- svm_rad_t1dvsctrl[[1]]

naiveBayes_t1dvsctrl <- readRDS("./ML_naiveBayes/naiveBayes_t1dvsctrl.rds")
naiveBayes_t1dvsctrl <- naiveBayes_t1dvsctrl[[1]]

p1 <- plot_box(ML1=xgb_t1dvsctrl, ML2= svm_lin_t1dvsctrl, ML3= svm_rad_t1dvsctrl, ML4= naiveBayes_t1dvsctrl, title = "All Cells")

###########################################
## Acinar
acinar_xgb_t1dvsctrl <- readRDS("./ML_xgb/acinar_xgb_t1dvsctrl.rds")
acinar_xgb_t1dvsctrl <- acinar_xgb_t1dvsctrl[[2]]

acinar_svm_lin_t1dvsctrl <- readRDS("./ML_svm/acinar_svm_lin_t1dvsctl.rds")
acinar_svm_lin_t1dvsctrl <- acinar_svm_lin_t1dvsctrl[[1]]

acinar_svm_rad_t1dvsctrl <- readRDS("./ML_svm/acinar_svm_rad_t1dvsctl.rds")
acinar_svm_rad_t1dvsctrl <- acinar_svm_rad_t1dvsctrl[[1]]

acinar_naiveBayes_t1dvsctrl <- readRDS("./ML_naiveBayes/acinar_naiveBayes_t1dvsctl.rds")
acinar_naiveBayes_t1dvsctrl <- acinar_naiveBayes_t1dvsctrl[[1]]

p2 <- plot_box(ML1=acinar_xgb_t1dvsctrl, ML2= acinar_svm_lin_t1dvsctrl, ML3 = acinar_svm_rad_t1dvsctrl, ML4= acinar_naiveBayes_t1dvsctrl, title = "Acinar cells")
############################################
## Alpha
alpha_xgb_t1dvsctrl <- readRDS("./ML_xgb/alpha_xgb_t1dvsctrl.rds")
alpha_xgb_t1dvsctrl <- alpha_xgb_t1dvsctrl[[2]]

alpha_svm_lin_t1dvsctrl <- readRDS("./ML_svm/alpha_svm_lin_t1dvsctl.rds")
alpha_svm_lin_t1dvsctrl <- alpha_svm_lin_t1dvsctrl[[1]]

alpha_svm_rad_t1dvsctrl <- readRDS("./ML_svm/alpha_svm_rad_t1dvsctl.rds")
alpha_svm_rad_t1dvsctrl <- alpha_svm_rad_t1dvsctrl[[1]]

alpha_naiveBayes_t1dvsctrl <- readRDS("./ML_naiveBayes/alpha_naiveBayes_t1dvsctl.rds")
alpha_naiveBayes_t1dvsctrl <- alpha_naiveBayes_t1dvsctrl[[1]]

p3 <- plot_box(ML1=alpha_xgb_t1dvsctrl, ML2= alpha_svm_lin_t1dvsctrl, ML3 = alpha_svm_rad_t1dvsctrl, ML4= alpha_naiveBayes_t1dvsctrl, title = "Alpha cells")
############################################
## Beta
beta_xgb_t1dvsctrl <- readRDS("./ML_xgb/beta_xgb_t1dvsctrl.rds")
beta_xgb_t1dvsctrl <- beta_xgb_t1dvsctrl[[2]]

beta_svm_lin_t1dvsctrl <- readRDS("./ML_svm/beta_svm_lin_t1dvsctl.rds")
beta_svm_lin_t1dvsctrl <- beta_svm_lin_t1dvsctrl[[1]]

beta_svm_rad_t1dvsctrl <- readRDS("./ML_svm/beta_svm_rad_t1dvsctl.rds")
beta_svm_rad_t1dvsctrl <- beta_svm_rad_t1dvsctrl[[1]]

beta_naiveBayes_t1dvsctrl <- readRDS("./ML_naiveBayes/beta_naiveBayes_t1dvsctl.rds")
beta_naiveBayes_t1dvsctrl <- beta_naiveBayes_t1dvsctrl[[1]]

p4 <- plot_box(ML1=beta_xgb_t1dvsctrl, ML2= beta_svm_lin_t1dvsctrl, ML3 = beta_svm_rad_t1dvsctrl, ML4= beta_naiveBayes_t1dvsctrl, title = "Beta Cells")
############################################
## Delta
delta_xgb_t1dvsctrl <- readRDS("./ML_xgb/delta_xgb_t1dvsctrl.rds")
delta_xgb_t1dvsctrl <- delta_xgb_t1dvsctrl[[2]]

delta_svm_lin_t1dvsctrl <- readRDS("./ML_svm/delta_svm_lin_t1dvsctl.rds")
delta_svm_lin_t1dvsctrl <- delta_svm_lin_t1dvsctrl[[1]]

delta_svm_rad_t1dvsctrl <- readRDS("./ML_svm/delta_svm_rad_t1dvsctl.rds")
delta_svm_rad_t1dvsctrl <- delta_svm_rad_t1dvsctrl[[1]]

delta_naiveBayes_t1dvsctrl <- readRDS("./ML_naiveBayes/delta_naiveBayes_t1dvsctl.rds")
delta_naiveBayes_t1dvsctrl <- delta_naiveBayes_t1dvsctrl[[1]]

p5 <- plot_box(ML1=delta_xgb_t1dvsctrl, ML2= delta_svm_lin_t1dvsctrl, ML3 = delta_svm_rad_t1dvsctrl, ML4= delta_naiveBayes_t1dvsctrl, title = "Delta Cells")
############################################
## Ductal
ductal_xgb_t1dvsctrl <- readRDS("./ML_xgb/ductal_xgb_t1dvsctrl.rds")
ductal_xgb_t1dvsctrl <- ductal_xgb_t1dvsctrl[[2]]

ductal_svm_lin_t1dvsctrl <- readRDS("./ML_svm/ductal_svm_lin_t1dvsctl.rds")
ductal_svm_lin_t1dvsctrl <- ductal_svm_lin_t1dvsctrl[[1]]

ductal_svm_rad_t1dvsctrl <- readRDS("./ML_svm/ductal_svm_rad_t1dvsctl.rds")
ductal_svm_rad_t1dvsctrl <- ductal_svm_rad_t1dvsctrl[[1]]

ductal_naiveBayes_t1dvsctrl <- readRDS("./ML_naiveBayes/ductal_naiveBayes_t1dvsctl.rds")
ductal_naiveBayes_t1dvsctrl <- ductal_naiveBayes_t1dvsctrl[[1]]

p6 <- plot_box(ML1=ductal_xgb_t1dvsctrl, ML2= ductal_svm_lin_t1dvsctrl, ML3 = ductal_svm_rad_t1dvsctrl, ML4= ductal_naiveBayes_t1dvsctrl, title = "Ductal Cells")
############################################
## Immune
immune_xgb_t1dvsctrl <- readRDS("./ML_xgb/immune_xgb_t1dvsctrl.rds")
immune_xgb_t1dvsctrl <- immune_xgb_t1dvsctrl[[2]]

immune_svm_lin_t1dvsctrl <- readRDS("./ML_svm/immune_svm_lin_t1dvsctl.rds")
immune_svm_lin_t1dvsctrl <- immune_svm_lin_t1dvsctrl[[1]]

immune_svm_rad_t1dvsctrl <- readRDS("./ML_svm/immune_svm_rad_t1dvsctl.rds")
immune_svm_rad_t1dvsctrl <- immune_svm_rad_t1dvsctrl[[1]]

immune_naiveBayes_t1dvsctrl <- readRDS("./ML_naiveBayes/immune_naiveBayes_t1dvsctl.rds")
immune_naiveBayes_t1dvsctrl <- immune_naiveBayes_t1dvsctrl[[1]]

p7 <- plot_box(ML1=immune_xgb_t1dvsctrl, ML2= immune_svm_lin_t1dvsctrl, ML3 = svm_rad_t1dvsctrl, ML4= immune_naiveBayes_t1dvsctrl, title = "Immune Cells")
############################################
## Endothelial
endothelial_xgb_t1dvsctrl <- readRDS("./ML_xgb/endothelial_xgb_t1dvsctrl.rds")
endothelial_xgb_t1dvsctrl <- endothelial_xgb_t1dvsctrl[[2]]

endothelial_svm_lin_t1dvsctrl <- readRDS("./ML_svm/endothelial_svm_lin_t1dvsctl.rds")
endothelial_svm_lin_t1dvsctrl <- endothelial_svm_lin_t1dvsctrl[[1]]

endothelial_svm_rad_t1dvsctrl <- readRDS("./ML_svm/endothelial_svm_rad_t1dvsctl.rds")
endothelial_svm_rad_t1dvsctrl <- endothelial_svm_rad_t1dvsctrl[[1]]

endothelial_naiveBayes_t1dvsctrl <- readRDS("./ML_naiveBayes/endothelial_naiveBayes_t1dvsctl.rds")
endothelial_naiveBayes_t1dvsctrl <- endothelial_naiveBayes_t1dvsctrl[[1]]

p8 <- plot_box(ML1=endothelial_xgb_t1dvsctrl, ML2= endothelial_svm_lin_t1dvsctrl, ML3 = svm_rad_t1dvsctrl, ML4= endothelial_naiveBayes_t1dvsctrl, title = "Endothelial Cells")
############################################
## Stellates
stellates_xgb_t1dvsctrl <- readRDS("./ML_xgb/stellates_xgb_t1dvsctrl.rds")
stellates_xgb_t1dvsctrl <- stellates_xgb_t1dvsctrl[[2]]

stellates_svm_lin_t1dvsctrl <- readRDS("./ML_svm/stellates_svm_lin_t1dvsctl.rds")
stellates_svm_lin_t1dvsctrl <- stellates_svm_lin_t1dvsctrl[[1]]

stellates_svm_rad_t1dvsctrl <- readRDS("./ML_svm/stellates_svm_rad_t1dvsctl.rds")
stellates_svm_rad_t1dvsctrl <- stellates_svm_rad_t1dvsctrl[[1]]

stellates_naiveBayes_t1dvsctrl <- readRDS("./ML_naiveBayes/stellates_naiveBayes_t1dvsctl.rds")
stellates_naiveBayes_t1dvsctrl <- stellates_naiveBayes_t1dvsctrl[[1]]

p9 <- plot_box(ML1=stellates_xgb_t1dvsctrl, ML2= stellates_svm_lin_t1dvsctrl, ML3 = svm_rad_t1dvsctrl, ML4= stellates_naiveBayes_t1dvsctrl, title = "Stellates Cells")

ggarrange(p1, p2, p3, p4, p5, p6, p7)
dev.copy(pdf,"./revised.figures.drafts/Fig1/ML_models_Acc_Box.pdf")
# ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9)
# dev.copy(pdf,"./revised.figures.drafts/Fig1/ML_models_Acc_Box_v1.pdf")
dev.off()
dev.off()
