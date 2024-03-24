##################################################################################
rm(list = ls())
gc()
library(ggplot2)
library(vioplot)
library(ggpubr)
library(gridExtra)

#######################################################################################
### Overall
setwd("/mnt/alvand/abhijeet/aab/apr_WO-T2D/ml/res/")

## T1D vs Ctrl
xgb_t1dvsctrl <- readRDS("xgb_t1dvsctrl.rds")
xgb_t1dvsctrl <- xgb_t1dvsctrl[[2]]

## T1D vs AAB
xgb_t1dvsaab <- readRDS("xgb_t1dvsaab.rds")
xgb_t1dvsaab <- xgb_t1dvsaab[[2]]

## AAB vs Ctrl
xgb_aabvsctrl <- readRDS("xgb_aabvsctrl.rds")
xgb_aabvsctrl <- xgb_aabvsctrl[[2]]

level_order <- c('T1D-CTL', 'T1D-AAb+', 'AAb+-CTL')

##################################################################################################
### Plotting

Accuracy <- cbind(xgb_t1dvsctrl[,1], xgb_t1dvsaab[,1], xgb_aabvsctrl[,1])*100
Sensitivity <- cbind(xgb_t1dvsctrl[,2], xgb_t1dvsaab[,2], xgb_aabvsctrl[,2])*100
Specificity <- cbind(xgb_t1dvsctrl[,3], xgb_t1dvsaab[,3], xgb_aabvsctrl[,3])*100
Metrics_num <- c(Accuracy, Sensitivity, Specificity)
Metrics <- c(rep("Accuracy",100), rep("Accuracy",100), rep("Accuracy",100),
             rep("Sensitivity",100), rep("Sensitivity",100), rep("Sensitivity",100),
             rep("Specificity",100), rep("Specificity",100), rep("Specificity",100))

df <- data.frame(cbind(Metrics_num, Metrics))
Condition <- c(rep("T1D-CTL",100), rep("T1D-AAb+",100), rep("AAb+-CTL",100),
               rep("T1D-CTL",100), rep("T1D-AAb+",100), rep("AAb+-CTL",100),
               rep("T1D-CTL",100), rep("T1D-AAb+",100), rep("AAb+-CTL",100))

df <- cbind(df,Condition)
df$Metrics_num <- as.numeric(df$Metrics_num)

p <- ggplot(df, aes(x= factor(Condition, level = level_order) , y = Metrics_num, fill = Metrics)) +
  geom_boxplot(outlier.size=0.1, lwd=0.1, fatten=1)+ #lwd=0.1, 
  facet_wrap(~factor(Condition, level = level_order), scale="free_x") + #"free" "fixed" "free_y"
  ggtitle("All Cells") + xlab("") + ylab("") + theme_bw() 


xtheme <- theme(plot.title = element_text(hjust = 0.5, size= 12)
                ,axis.text.y = element_text(angle = 0, size = 8, hjust=1)
                ,axis.title.y = element_text( size = rel(0.7))
                ,axis.text.x = element_blank()
                ,axis.title.x = element_text( size = rel(0.7))
                ,axis.ticks.x=element_blank(), strip.text = element_text(size=12))
legtheme <- theme(legend.title = element_text( size = 12)
                  ,legend.text=element_text( size=12)
                  ,legend.key.size = unit(1.5, "cm")
                  ,legend.key.width = unit(1.5,"cm")
                  , legend.position = "right"
                  ,legend.box = "vertical")

p + legtheme 
legend <- get_legend(p)
p <- p + xtheme
# pdf('../plots/Cell_Metrics_Overall.pdf', width=8, height = 8)
# p
# dev.off()
# dev.off()

# ggarrange(p, ncol=1, nrow=1, common.legend = TRUE, legend="right")
# dev.copy(pdf,"Metrics_1.pdf")
# dev.off()

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
                  delta_xgb_t1dvsctrl[,1], immune_xgb_t1dvsctrl[,1])*100
Sensitivity <- cbind(acinar_xgb_t1dvsctrl[,2], alpha_xgb_t1dvsctrl[,2], 
                     beta_xgb_t1dvsctrl[,2], ductal_xgb_t1dvsctrl[,2],
                     delta_xgb_t1dvsctrl[,2], immune_xgb_t1dvsctrl[,2])*100
Specificity <- cbind(acinar_xgb_t1dvsctrl[,3], alpha_xgb_t1dvsctrl[,3], 
                     beta_xgb_t1dvsctrl[,3], ductal_xgb_t1dvsctrl[,3],
                     delta_xgb_t1dvsctrl[,3], immune_xgb_t1dvsctrl[,3])*100
Metrics_num <- c(Accuracy, Sensitivity, Specificity)
Metrics <- c(rep("Accuracy",600), rep("Sensitivity",600), rep("Specificity",600))

df1 <- data.frame(cbind(Metrics_num, Metrics))
Condition <- c(rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100),
               rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100),
               rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100))

df1 <- cbind(df1,Condition)
df1$Metrics_num <- as.numeric(df1$Metrics_num)

p_1 <- ggplot(df1, aes(x= factor(Condition, level = level_order_1) , y = Metrics_num, fill = Metrics)) +
  geom_boxplot(outlier.size=0.1, lwd=0.1, fatten=1)+ 
  facet_wrap(~factor(Condition, level = level_order_1), scale="free_x") + #"free" "fixed" "free_y"
  ggtitle("T1D-CTL") + xlab("") + ylab("") + theme_bw() #Major Cell Types in 

p_1 + legtheme 
legend <- get_legend(p_1)
p_1 <- p_1 + xtheme
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
                  delta_xgb_t1dvsaab[,1], immune_xgb_t1dvsaab[,1])*100
Sensitivity <- cbind(acinar_xgb_t1dvsaab[,2], alpha_xgb_t1dvsaab[,2], 
                     beta_xgb_t1dvsaab[,2], ductal_xgb_t1dvsaab[,2],
                     delta_xgb_t1dvsaab[,2], immune_xgb_t1dvsaab[,2])*100
Specificity <- cbind(acinar_xgb_t1dvsaab[,3], alpha_xgb_t1dvsaab[,3], 
                     beta_xgb_t1dvsaab[,3], ductal_xgb_t1dvsaab[,3],
                     delta_xgb_t1dvsaab[,3], immune_xgb_t1dvsaab[,3])*100
Metrics_num <- c(Accuracy, Sensitivity, Specificity)
Metrics <- c(rep("Accuracy",600), rep("Sensitivity",600), rep("Specificity",600))


df2 <- data.frame(cbind(Metrics_num, Metrics))
Condition <- c(rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100),
               rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100),
               rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100))

df2 <- cbind(df2,Condition)
df2$Metrics_num <- as.numeric(df2$Metrics_num)

p_2 <- ggplot(df2, aes(x= factor(Condition, level = level_order_2) , y = Metrics_num, fill = Metrics)) +
  geom_boxplot(outlier.size=0.1, lwd=0.1, fatten=1)+ 
  facet_wrap(~factor(Condition, level = level_order_2), scale="free_x") + #"free" "fixed" "free_y"
  ggtitle("T1D-AAb+") + xlab("") + ylab("") + theme_bw() #Major Cell Types in 

p_2 + legtheme 
legend <- get_legend(p_2)
p_2 <- p_2 + xtheme

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
                  delta_xgb_aabvsctrl[,1], immune_xgb_aabvsctrl[,1])*100
Sensitivity <- cbind(acinar_xgb_aabvsctrl[,2], alpha_xgb_aabvsctrl[,2], 
                     beta_xgb_aabvsctrl[,2], ductal_xgb_aabvsctrl[,2],
                     delta_xgb_aabvsctrl[,2], immune_xgb_aabvsctrl[,2])*100
Specificity <- cbind(acinar_xgb_aabvsctrl[,3], alpha_xgb_aabvsctrl[,3], 
                     beta_xgb_aabvsctrl[,3], ductal_xgb_aabvsctrl[,3],
                     delta_xgb_aabvsctrl[,3], immune_xgb_aabvsctrl[,3])*100
Metrics_num <- c(Accuracy, Sensitivity, Specificity)
Metrics <- c(rep("Accuracy",600), rep("Sensitivity",600), rep("Specificity",600))


df3 <- data.frame(cbind(Metrics_num, Metrics))
Condition <- c(rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100),
               rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100),
               rep("Acinar",100), rep("Alpha",100), rep("Beta",100), rep("Ductal",100), rep("Delta",100), rep("Immune",100))

df3 <- cbind(df3,Condition)
df3$Metrics_num <- as.numeric(df3$Metrics_num)

p_3 <- ggplot(df3, aes(x= factor(Condition, level = level_order_3) , y = Metrics_num, fill = Metrics)) +
  geom_boxplot(outlier.size=0.1, lwd=0.1, fatten=1)+ 
  facet_wrap(~factor(Condition, level = level_order_3), scale="free_x") + #"free" "fixed" "free_y"
  ggtitle("AAb+-CTL") + xlab("") + ylab("") + theme_bw() # Major Cell Types in 

p_3 + legtheme 
legend <- get_legend(p_3)
p_3 <- p_3 + xtheme
ggarrange(p, p_1, p_2, p_3, common.legend = T, legend="bottom") # ncol=2, nrow=2, 
dev.copy(pdf,"../plots/Cell_Metrics_Final_v3_box.pdf")
dev.off()
dev.off()