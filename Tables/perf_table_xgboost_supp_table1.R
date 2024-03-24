rm(list = ls())
gc()
library(ggplot2)
library(vioplot)
library(ggpubr)
library(gridExtra)
library(dplyr)
#######################################################################################
### Overall
setwd("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/ML_xgb/")

## T1D vs Ctrl
xgb_t1dvsctrl <- readRDS("xgb_t1dvsctrl.rds")
xgb_t1dvsctrl <- xgb_t1dvsctrl[[2]]

# tmp <- xgb_t1dvsctrl[[2]]
# mean(tmp[,3])
# $ : num [1:3] 0.991 0.997 0.969
# $ : num [1:3] 0.000414 0.000352 0.001964
# round(0.000414,4)*100

## T1D vs AAB
xgb_t1dvsaab <- readRDS("xgb_t1dvsaab.rds")
xgb_t1dvsaab <- xgb_t1dvsaab[[2]]
# $ : num [1:3] 0.992 0.994 0.989
# $ : num [1:3] 0.00074 0.000848 0.001201

## AAB vs Ctrl
xgb_aabvsctrl <- readRDS("xgb_aabvsctrl.rds")
xgb_aabvsctrl <- xgb_aabvsctrl[[2]]
# $ : num [1:3] 0.963 0.987 0.884
# $ : num [1:3] 0.00226 0.00133 0.00875

level_order <- c('T1D-CTL', 'T1D-AAB', 'AAB-CTL')

##################################################################################################
### Plotting

Accuracy <- cbind(xgb_t1dvsctrl[,1], xgb_t1dvsaab[,1], xgb_aabvsctrl[,1])
Sensitivity <- cbind(xgb_t1dvsctrl[,2], xgb_t1dvsaab[,2], xgb_aabvsctrl[,2])
Specificity <- cbind(xgb_t1dvsctrl[,3], xgb_t1dvsaab[,3], xgb_aabvsctrl[,3])
Metrics_num <- c(Accuracy, Sensitivity, Specificity)
# Metrics <- c(rep("Accuracy",100), rep("Sensitivity",100), rep("Specificity",100),
#              rep("Accuracy",100), rep("Sensitivity",100), rep("Specificity",100),
#              rep("Accuracy",100), rep("Sensitivity",100), rep("Specificity",100))

Metrics <- c(rep("Accuracy",100), rep("Accuracy",100), rep("Accuracy",100),
             rep("Sensitivity",100), rep("Sensitivity",100), rep("Sensitivity",100),
             rep("Specificity",100), rep("Specificity",100), rep("Specificity",100))

df <- data.frame(cbind(Metrics_num, Metrics))
Condition <- c(rep("T1D-CTL",100), rep("T1D-AAB",100), rep("AAB-CTL",100),
               rep("T1D-CTL",100), rep("T1D-AAB",100), rep("AAB-CTL",100),
               rep("T1D-CTL",100), rep("T1D-AAB",100), rep("AAB-CTL",100))

df <- cbind(df,Condition)
df$Metrics_num <- as.numeric(df$Metrics_num)


df_overall_xgb = df %>% group_by(Condition, Metrics)  %>%
  summarise(Mean = round(mean(Metrics_num),4)*100,
            SD = round(sd(Metrics_num),4)*100,
            .groups = 'drop')

library(reshape2) 
df_overall_xgb <- reshape2::melt(df_overall_xgb) 
df_overall_xgb <- dcast(df_overall_xgb, Condition ~...) 
write.csv(df_overall_xgb, "../revised.figures.drafts/df_overall_xgb.csv", row.names = F)
stargazer::stargazer(df_overall_xgb, summary = F, rownames = F, type = 'text', out = '../revised.figures.drafts/df_overall_xgb.txt')

#################################################################################################
## Cell type
# rm(list = ls())
# gc()

## T1D vs CTL
acinar_xgb_t1dvsctrl <- readRDS("acinar_xgb_t1dvsctrl.rds")
acinar_xgb_t1dvsctrl <- acinar_xgb_t1dvsctrl[[2]]

alpha_xgb_t1dvsctrl <- readRDS("alpha_xgb_t1dvsctrl.rds")
alpha_xgb_t1dvsctrl <- alpha_xgb_t1dvsctrl[[2]]

beta_xgb_t1dvsctrl <- readRDS("beta_xgb_t1dvsctrl.rds")
beta_xgb_t1dvsctrl <- beta_xgb_t1dvsctrl[[2]]

ductal_xgb_t1dvsctrl <- readRDS("ductal_xgb_t1dvsctrl.rds")
ductal_xgb_t1dvsctrl <- ductal_xgb_t1dvsctrl[[2]]

delta_xgb_t1dvsctrl <- readRDS("delta_xgb_t1dvsctrl.rds")
delta_xgb_t1dvsctrl <- delta_xgb_t1dvsctrl[[2]]

immune_xgb_t1dvsctrl <- readRDS("immune_xgb_t1dvsctrl.rds")
immune_xgb_t1dvsctrl <- immune_xgb_t1dvsctrl[[2]]

level_order_1 <- c('Acinar', 'Alpha', 
                   'Beta', 'Ductal',
                   'Delta', 'Immune')

Accuracy <- cbind(acinar_xgb_t1dvsctrl[,1], alpha_xgb_t1dvsctrl[,1], 
                  beta_xgb_t1dvsctrl[,1], ductal_xgb_t1dvsctrl[,1],
                  delta_xgb_t1dvsctrl[,1], immune_xgb_t1dvsctrl[,1])
Sensitivity <- cbind(acinar_xgb_t1dvsctrl[,2], alpha_xgb_t1dvsctrl[,2], 
                     beta_xgb_t1dvsctrl[,2], ductal_xgb_t1dvsctrl[,2],
                     delta_xgb_t1dvsctrl[,2], immune_xgb_t1dvsctrl[,2])
Specificity <- cbind(acinar_xgb_t1dvsctrl[,3], alpha_xgb_t1dvsctrl[,3], 
                     beta_xgb_t1dvsctrl[,3], ductal_xgb_t1dvsctrl[,3],
                     delta_xgb_t1dvsctrl[,3], immune_xgb_t1dvsctrl[,3])
Metrics_num <- c(Accuracy, Sensitivity, Specificity)
Metrics <- c(rep("Accuracy",600), rep("Sensitivity",600), rep("Specificity",600))

df1 <- data.frame(cbind(Metrics_num, Metrics))
Condition <- c(rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100),
               rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100),
               rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100))

df1 <- cbind(df1,Condition)
df1$Metrics_num <- as.numeric(df1$Metrics_num)

df_overall_xgb_t1dvsctl = df1 %>% group_by(Condition, Metrics)  %>%
  summarise(Mean = round(mean(Metrics_num),4)*100,
            SD = round(sd(Metrics_num),4)*100,
            .groups = 'drop')

library(reshape2) 
df_overall_xgb_t1dvsctl <- reshape2::melt(df_overall_xgb_t1dvsctl) 
df_overall_xgb_t1dvsctl <- dcast(df_overall_xgb_t1dvsctl, Condition ~...) 
write.csv(df_overall_xgb_t1dvsctl, "../revised.figures.drafts/df_overall_xgb_t1dvsctl.csv", row.names = F)
stargazer::stargazer(df_overall_xgb_t1dvsctl, summary = F, rownames = F, type = 'text', out = '../revised.figures.drafts/df_overall_xgb_t1dvsctl.txt')

##################################################################################################

## T1D vs AAB
acinar_xgb_t1dvsaab <- readRDS("acinar_xgb_t1dvsaab.rds")
acinar_xgb_t1dvsaab <- acinar_xgb_t1dvsaab[[2]]

alpha_xgb_t1dvsaab <- readRDS("alpha_xgb_t1dvsaab.rds")
alpha_xgb_t1dvsaab <- alpha_xgb_t1dvsaab[[2]]

beta_xgb_t1dvsaab <- readRDS("beta_xgb_t1dvsaab.rds")
beta_xgb_t1dvsaab <- beta_xgb_t1dvsaab[[2]]

ductal_xgb_t1dvsaab <- readRDS("ductal_xgb_t1dvsaab.rds")
ductal_xgb_t1dvsaab <- ductal_xgb_t1dvsaab[[2]]

delta_xgb_t1dvsaab <- readRDS("delta_xgb_t1dvsaab.rds")
delta_xgb_t1dvsaab <- delta_xgb_t1dvsaab[[2]]

immune_xgb_t1dvsaab <- readRDS("immune_xgb_t1dvsaab.rds")
immune_xgb_t1dvsaab <- immune_xgb_t1dvsaab[[2]]

level_order_2 <- c('Acinar', 'Alpha', 
                   'Beta', 'Ductal',
                   'Delta', 'Immune')

Accuracy <- cbind(acinar_xgb_t1dvsaab[,1], alpha_xgb_t1dvsaab[,1], 
                  beta_xgb_t1dvsaab[,1], ductal_xgb_t1dvsaab[,1],
                  delta_xgb_t1dvsaab[,1], immune_xgb_t1dvsaab[,1])
Sensitivity <- cbind(acinar_xgb_t1dvsaab[,2], alpha_xgb_t1dvsaab[,2], 
                     beta_xgb_t1dvsaab[,2], ductal_xgb_t1dvsaab[,2],
                     delta_xgb_t1dvsaab[,2], immune_xgb_t1dvsaab[,2])
Specificity <- cbind(acinar_xgb_t1dvsaab[,3], alpha_xgb_t1dvsaab[,3], 
                     beta_xgb_t1dvsaab[,3], ductal_xgb_t1dvsaab[,3],
                     delta_xgb_t1dvsaab[,3], immune_xgb_t1dvsaab[,3])
Metrics_num <- c(Accuracy, Sensitivity, Specificity)
Metrics <- c(rep("Accuracy",600), rep("Sensitivity",600), rep("Specificity",600))

df2 <- data.frame(cbind(Metrics_num, Metrics))
Condition <- c(rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100),
               rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100),
               rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100))

df2 <- cbind(df2,Condition)
df2$Metrics_num <- as.numeric(df2$Metrics_num)


df_overall_xgb_t1dvsaab = df2 %>% group_by(Condition, Metrics)  %>%
  summarise(Mean = round(mean(Metrics_num),4)*100,
            SD = round(sd(Metrics_num),4)*100,
            .groups = 'drop')

library(reshape2) 
df_overall_xgb_t1dvsaab <- reshape2::melt(df_overall_xgb_t1dvsaab) 
df_overall_xgb_t1dvsaab <- dcast(df_overall_xgb_t1dvsaab, Condition ~...) 
write.csv(df_overall_xgb_t1dvsaab, "../revised.figures.drafts/df_overall_xgb_t1dvsaab.csv", row.names = F)
stargazer::stargazer(df_overall_xgb_t1dvsaab, summary = F, rownames = F, type = 'text', out = '../revised.figures.drafts/df_overall_xgb_t1dvsaab.txt')

#########################################################################################################

## AAB vs CTL
acinar_xgb_aabvsctrl <- readRDS("acinar_xgb_aabvsctrl.rds")
acinar_xgb_aabvsctrl <- acinar_xgb_aabvsctrl[[2]]

alpha_xgb_aabvsctrl <- readRDS("alpha_xgb_aabvsctrl.rds")
alpha_xgb_aabvsctrl <- alpha_xgb_aabvsctrl[[2]]

beta_xgb_aabvsctrl <- readRDS("beta_xgb_aabvsctrl.rds")
beta_xgb_aabvsctrl <- beta_xgb_aabvsctrl[[2]]

ductal_xgb_aabvsctrl <- readRDS("ductal_xgb_aabvsctrl.rds")
ductal_xgb_aabvsctrl <- ductal_xgb_aabvsctrl[[2]]

delta_xgb_aabvsctrl <- readRDS("delta_xgb_aabvsctrl.rds")
delta_xgb_aabvsctrl <- delta_xgb_aabvsctrl[[2]]

immune_xgb_aabvsctrl <- readRDS("immune_xgb_aabvsctrl.rds")
immune_xgb_aabvsctrl <- immune_xgb_aabvsctrl[[2]]

level_order_3 <- c('Acinar', 'Alpha', 
                   'Beta', 'Ductal',
                   'Delta', 'Immune')

Accuracy <- cbind(acinar_xgb_aabvsctrl[,1], alpha_xgb_aabvsctrl[,1], 
                  beta_xgb_aabvsctrl[,1], ductal_xgb_aabvsctrl[,1],
                  delta_xgb_aabvsctrl[,1], immune_xgb_aabvsctrl[,1])
Sensitivity <- cbind(acinar_xgb_aabvsctrl[,2], alpha_xgb_aabvsctrl[,2], 
                     beta_xgb_aabvsctrl[,2], ductal_xgb_aabvsctrl[,2],
                     delta_xgb_aabvsctrl[,2], immune_xgb_aabvsctrl[,2])
Specificity <- cbind(acinar_xgb_aabvsctrl[,3], alpha_xgb_aabvsctrl[,3], 
                     beta_xgb_aabvsctrl[,3], ductal_xgb_aabvsctrl[,3],
                     delta_xgb_aabvsctrl[,3], immune_xgb_aabvsctrl[,3])
Metrics_num <- c(Accuracy, Sensitivity, Specificity)
Metrics <- c(rep("Accuracy",600), rep("Sensitivity",600), rep("Specificity",600))

df3 <- data.frame(cbind(Metrics_num, Metrics))
Condition <- c(rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100),
               rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100),
               rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100))

df3 <- cbind(df3,Condition)
df3$Metrics_num <- as.numeric(df3$Metrics_num)


df_overall_xgb_aabvsctl = df3 %>% group_by(Condition, Metrics)  %>%
  summarise(Mean = round(mean(Metrics_num),4)*100,
            SD = round(sd(Metrics_num),4)*100,
            .groups = 'drop')

library(reshape2) 
df_overall_xgb_aabvsctl <- reshape2::melt(df_overall_xgb_aabvsctl) 
df_overall_xgb_aabvsctl <- dcast(df_overall_xgb_aabvsctl, Condition ~...) 
write.csv(df_overall_xgb_aabvsctl, "../revised.figures.drafts/df_overall_xgb_aabvsctl.csv", row.names = F)
stargazer::stargazer(df_overall_xgb_aabvsctl, summary = F, rownames = F, type = 'text', out = '../revised.figures.drafts/df_overall_xgb_aabvsctl.txt')
##################################################################################

