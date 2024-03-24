rm(list = ls())
gc()
library(ggplot2)
library(vioplot)
library(ggpubr)
library(gridExtra)
library(dplyr)
#######################################################################################
### Overall
setwd("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/ML_svm/")

## T1D vs Ctrl
# svm_rad_t1dvsctrl <- readRDS("beta_svm_rad_t1dvsctrl.rds")
# svm_rad_t1dvsctrl <- svm_rad_t1dvsctrl[[2]]
svm_rad_t1dvsctrl <- readRDS("svm_t1dvsctrl_rad.rds")
svm_rad_t1dvsctrl <- svm_rad_t1dvsctrl[[1]]
##################################################################################################
### Plotting

Accuracy <- cbind(svm_rad_t1dvsctrl[,1])
Sensitivity <- cbind(svm_rad_t1dvsctrl[,2])
Specificity <- cbind(svm_rad_t1dvsctrl[,3])
Metrics_num <- c(Accuracy, Sensitivity, Specificity)
Metrics <- c(rep("Accuracy",length(Accuracy)), rep("Sensitivity",length(Sensitivity)), rep("Specificity",length(Specificity)))

df <- data.frame(cbind(Metrics_num, Metrics))
Condition <- c(rep("T1D-CTL",length(Accuracy)+length(Sensitivity)+length(Specificity)))

df <- cbind(df,Condition)
df$Metrics_num <- as.numeric(df$Metrics_num)

df_overall_svm_rad = df %>% group_by(Condition, Metrics)  %>%
  summarise(Mean = round(mean(Metrics_num),2),
            SD = round(sd(Metrics_num),2),
            .groups = 'drop')

library(reshape2) 
df_overall_svm_rad <- reshape2::melt(df_overall_svm_rad) 
df_overall_svm_rad <- dcast(df_overall_svm_rad, Condition ~...) 
write.csv(df_overall_svm_rad, "../revised.figures.drafts/df_overall_svm_rad.csv", row.names = F)
stargazer::stargazer(df_overall_svm_rad, summary = F, rownames = F, type = 'text', out = '../revised.figures.drafts/df_overall_svm_rad.txt')

#################################################################################################
#################################################################################################
## Cell type
# rm(list = ls())
# gc()

## T1D vs CTL
acinar_svm_rad_t1dvsctrl <- readRDS("acinar_svm_rad_t1dvsctl.rds")
acinar_svm_rad_t1dvsctrl <- acinar_svm_rad_t1dvsctrl[[1]]

alpha_svm_rad_t1dvsctrl <- readRDS("alpha_svm_rad_t1dvsctl.rds")
alpha_svm_rad_t1dvsctrl <- alpha_svm_rad_t1dvsctrl[[1]]

beta_svm_rad_t1dvsctrl <- readRDS("beta_svm_rad_t1dvsctl.rds")
beta_svm_rad_t1dvsctrl <- beta_svm_rad_t1dvsctrl[[1]]

ductal_svm_rad_t1dvsctrl <- readRDS("ductal_svm_rad_t1dvsctl.rds")
ductal_svm_rad_t1dvsctrl <- ductal_svm_rad_t1dvsctrl[[1]]

delta_svm_rad_t1dvsctrl <- readRDS("delta_svm_rad_t1dvsctl.rds")
delta_svm_rad_t1dvsctrl <- delta_svm_rad_t1dvsctrl[[1]]

immune_svm_rad_t1dvsctrl <- readRDS("immune_svm_rad_t1dvsctl.rds")
immune_svm_rad_t1dvsctrl <- immune_svm_rad_t1dvsctrl[[1]]

# stellates_svm_rad_t1dvsctrl <- readRDS("stellates_svm_rad_t1dvsctl.rds")
# stellates_svm_rad_t1dvsctrl <- stellates_svm_rad_t1dvsctrl[[1]]

level_order_1 <- c('Acinar', 'Alpha', 
                   'Beta', 'Ductal',
                   'Delta', 'Immune')

Accuracy <- cbind(acinar_svm_rad_t1dvsctrl[,1], alpha_svm_rad_t1dvsctrl[,1], 
                  beta_svm_rad_t1dvsctrl[,1], ductal_svm_rad_t1dvsctrl[,1],
                  delta_svm_rad_t1dvsctrl[,1], immune_svm_rad_t1dvsctrl[,1])
Sensitivity <- cbind(acinar_svm_rad_t1dvsctrl[,2], alpha_svm_rad_t1dvsctrl[,2], 
                     beta_svm_rad_t1dvsctrl[,2], ductal_svm_rad_t1dvsctrl[,2],
                     delta_svm_rad_t1dvsctrl[,2], immune_svm_rad_t1dvsctrl[,2])
Specificity <- cbind(acinar_svm_rad_t1dvsctrl[,3], alpha_svm_rad_t1dvsctrl[,3], 
                     beta_svm_rad_t1dvsctrl[,3], ductal_svm_rad_t1dvsctrl[,3],
                     delta_svm_rad_t1dvsctrl[,3], immune_svm_rad_t1dvsctrl[,3])
Metrics_num <- c(Accuracy, Sensitivity, Specificity)
Metrics <- c(rep("Accuracy",60), rep("Sensitivity",60), rep("Specificity",60))

df1 <- data.frame(cbind(Metrics_num, Metrics))
Condition <- c(rep("Acinar",10), rep("Alpha",10), rep("Beta",10), rep("Ductal",10), rep("Delta",10), rep("Immune",10),
               rep("Acinar",10), rep("Alpha",10), rep("Beta",10), rep("Ductal",10), rep("Delta",10), rep("Immune",10),
               rep("Acinar",10), rep("Alpha",10), rep("Beta",10), rep("Ductal",10), rep("Delta",10), rep("Immune",10))

df1 <- cbind(df1,Condition)
df1$Metrics_num <- as.numeric(df1$Metrics_num)

df_overall_svm_rad_t1dvsctl = df1 %>% group_by(Condition, Metrics)  %>%
  summarise(Mean = round(mean(Metrics_num),2),
            SD = round(sd(Metrics_num),2),
            .groups = 'drop')

library(reshape2) 
df_overall_svm_rad_t1dvsctl <- reshape2::melt(df_overall_svm_rad_t1dvsctl) 
df_overall_svm_rad_t1dvsctl <- dcast(df_overall_svm_rad_t1dvsctl, Condition ~...) 
write.csv(df_overall_svm_rad_t1dvsctl, "../revised.figures.drafts/df_overall_svm_rad_t1dvsctl.csv", row.names = F)
stargazer::stargazer(df_overall_svm_rad_t1dvsctl, summary = F, rownames = F, type = 'text', out = '../revised.figures.drafts/df_overall_svm_rad_t1dvsctl.txt')

##################################################################################################

