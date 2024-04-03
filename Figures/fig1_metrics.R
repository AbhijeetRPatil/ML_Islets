rm(list = ls())
setwd("./")
library("Seurat")
library("dplyr")
library("Matrix")
library("ggplot2")
library("viridis")
library("ggprism")
library("viridisLite")
library("patchwork")
################################################################################################################
## Before filtering plot
panc_raw_meta_04182022 <- readRDS("panc.rds")

# Create Data
## All- AAB Control T1D
Prop <- c(table(panc_raw_meta_04182022$disease_state)[[1]], 
          table(panc_raw_meta_04182022$disease_state)[[2]] , 
          table(panc_raw_meta_04182022$disease_state)[[3]])

# number of colors in the palette
no_of_colors <- length(Prop)

# options represent the color types, there are altogether 8 options.
palette <- viridis_pal(option = "D")(no_of_colors)

# hex color codes
palette
# represents colors in a pie chart manner
pdf('pie_chart.pdf', height = 5, width = 5)
pie(Prop, labels = c(paste0("AAB \n ncells= ", Prop[1], "\n ndonors= 10"),
                     paste0("CTL \n ncells= ", Prop[2], "\n ndonors= 31"),
                     paste0("\n T1D \n ncells= ", Prop[3], "\n ndonors= 9")),
    col = palette)
dev.off()

#################
## For different groups
## all cells
bf_all <- table(panc_raw_meta_04182022@meta.data$disease_id)

## aab cells
bf_aab <- bf_all[grep("AAB", names(bf_all))]
names(bf_aab) <- gsub("_ ", "", names(bf_aab))
bf_aab <- as.data.frame(bf_aab)


## ctl cells
bf_ctl <- bf_all[grep("Control", names(bf_all))]
names(bf_ctl) <- gsub("_ ", "", names(bf_ctl))
bf_ctl <- as.data.frame(bf_ctl)
bf_ctl$Var1 <- gsub("Control", "CTL", bf_ctl$Var1)

## t1d cells
bf_t1d <- bf_all[grep("T1D", names(bf_all))]
names(bf_t1d) <- gsub("_ ", "", names(bf_t1d))
bf_t1d <- as.data.frame(bf_t1d)

## Plotting

## AAB
pie_labels <- paste0(bf_aab$Var1, "  ", round(100 * bf_aab$Freq/sum(bf_aab$Freq), 2), "%")
pdf('pie_chart_aab.pdf', height = 5, width = 5)
pie(bf_aab$Freq, labels = pie_labels, col =  hcl.colors(length(bf_aab$Var1), "Purples"))  #BluYl Dynamic
dev.off()

## CTL
pie_labels <- paste0(bf_ctl$Var1, "  ", round(100 * bf_ctl$Freq/sum(bf_ctl$Freq), 2), "%")
pdf('pie_chart_ctl.pdf', height = 13, width = 13)
pie(bf_ctl$Freq, labels = pie_labels, col =  hcl.colors(length(bf_ctl$Var1), "TealGrn")) #BluYl
dev.off()

## T1D
pie_labels <- paste0(bf_t1d$Var1, "  ", round(100 * bf_t1d$Freq/sum(bf_t1d$Freq), 2), "%")
pdf('pie_chart_t1d.pdf', height = 5, width = 5)
pie(bf_t1d$Freq, labels = pie_labels, col =  c("#FFEA00", "#FCF55F", "#FADA5E",
                                               "#FAFA33", "#F4BB44", "#FBEC5D",
                                               "#FFFF00", "#FFFAA0", "#FFE5B4")) #BluYl
dev.off()

###############################################################################################################
## Median Gene Per Sample
df_mgps <- data.frame(panc_raw_meta_04182022@meta.data$nFeature_RNA, panc_raw_meta_04182022@meta.data$disease_state) 
colnames(df_mgps) <- c("Median_Gene_Number", "Condition")
pdf('boxplot_mediangenenumber.pdf', height = 5, width = 5)
ggplot(df_mgps, aes(Condition, Median_Gene_Number)) + geom_boxplot(aes(fill = Condition), 
                                                                   width=0.5, outlier.size = 0.2) +
  scale_fill_viridis(discrete = TRUE) + xlab("") + theme(legend.position = "none")
dev.off()
#######################################################################################################################

########################################################################################################################
## After filtering
########################################################################################################################
## Read data
dat <- readRDS("/mnt/alvand/abhijeet/aab/apr_WO-T2D/objs/panc_raw_meta_filtered_sct_umap_WO-T2D_WO_RB_MT_04182022.rds")
## Create response variable
df <- dat@meta.data %>% 
  mutate(class_label = recode(disease_state, 
                              "Control" = 0, 
                              "AAB" = 1,
                              "T1D" = 2)) 
dat[["class_label_num"]] <- as.numeric(df$class_label)

########################################################################################################################
## Stacked bar plot 
df <- data.frame("Cell_Type"=dat@meta.data$cell_type, "Condition"=dat@meta.data$disease_state) 
df <- df %>% group_by(Condition, Cell_Type) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(Percent = Nb/C*100) 

df<-df[df$Cell_Type!="Unknown",]
# df$Condition <- gsub("Control", "CTL", df$Condition)
df$Cell_Type <- gsub("Stellates_Mesenchymal", "Stellates", df$Cell_Type)
df$Cell_Type <- gsub("PP_Gamma", "PP", df$Cell_Type)

## Theme for plotting
xtheme <- theme_bw()+ theme(plot.title = element_text(face = "bold" ,hjust = 0.5, size= 10)
                            ,axis.text.y = element_text(face = "bold",angle = 0, size = 10, hjust=1)
                            ,axis.title.y = element_text(face = "bold", size = rel(0.8))
                            ,axis.text.x = element_text(face = "bold",angle = 0, size = 10)
                            ,axis.title.x = element_text(face = "bold", size = rel(0.8))
                            ,axis.ticks.x=element_blank(), strip.text = element_text(size=10)) 


pdf('celltype_stackbarchart_celltypes_condition_filtered.pdf', height = 4, width = 5)
ggplot(df, aes(fill=Condition, y=Percent, x=Cell_Type)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_viridis(discrete = T) + xlab("") + xtheme + 
  theme(legend.position='top', axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()
###############################################################################################################
## Median Gene Per Sample filtered
df_mgps <- data.frame(dat@meta.data$nFeature_RNA, dat@meta.data$disease_state) 
colnames(df_mgps) <- c("Median_Gene_Number", "Condition")
pdf('boxplot_mediangenenumber_filtered.pdf', height = 5, width = 5)
ggplot(df_mgps, aes(Condition, Median_Gene_Number)) + geom_boxplot(aes(fill = Condition), 
                                                                   width=0.5, outlier.size = 0.2) +
  scale_fill_viridis(discrete = TRUE) + xlab("") + theme(legend.position = "none")
dev.off()
#######################################################################################################################
## Pie chart

# Create Data
## AAB Control T1D T2D
Prop <- c(table(dat$disease_state)[[1]], 
          table(dat$disease_state)[[2]], 
          table(dat$disease_state)[[3]])

# number of colors in the palette
no_of_colors <- length(Prop)

# options represent the color types, there are altogether 8 options.
palette <- viridis_pal(option = "D")(no_of_colors)

# hex color codes
palette
# represents colors in a pie chart manner
pdf('pie_chart_filtered.pdf', height = 5, width = 5)
pie(Prop, labels = c(paste0("AAB \n ncells= ", Prop[1], "\n ndonors= 10"),
                     paste0("CTL \n ncells= ", Prop[2], "\n ndonors= 31"),
                     paste0("\n T1D \n ncells= ", Prop[3], "\n ndonors= 9")),
    col = palette)
dev.off()

############
## Metrics for all samples
## all cells
bf_all <- table(dat$disease_id)

## aab cells
bf_aab <- bf_all[grep("AAB", names(bf_all))]
names(bf_aab) <- gsub("_ ", "", names(bf_aab))
bf_aab <- as.data.frame(bf_aab)


## ctl cells
bf_ctl <- bf_all[grep("Control", names(bf_all))]
names(bf_ctl) <- gsub("_ ", "", names(bf_ctl))
bf_ctl <- as.data.frame(bf_ctl)
bf_ctl$Var1 <- gsub("Control", "CTL", bf_ctl$Var1)

## t1d cells
bf_t1d <- bf_all[grep("T1D", names(bf_all))]
names(bf_t1d) <- gsub("_ ", "", names(bf_t1d))
bf_t1d <- as.data.frame(bf_t1d)


## AAB
pie_labels <- paste0(bf_aab$Var1, "  ", round(100 * bf_aab$Freq/sum(bf_aab$Freq), 2), "%")
pdf('pie_chart_aab_filtered.pdf', height = 5, width = 5)
pie(bf_aab$Freq, labels = pie_labels, col =  hcl.colors(length(bf_aab$Var1), "Purples"))  #BluYl Dynamic
dev.off()

## CTL
pie_labels <- paste0(bf_ctl$Var1, "  ", round(100 * bf_ctl$Freq/sum(bf_ctl$Freq), 2), "%")
pdf('pie_chart_ctl_filtered.pdf', height = 13, width = 13)
pie(bf_ctl$Freq, labels = pie_labels, col =  hcl.colors(length(bf_ctl$Var1), "TealGrn")) #BluYl
dev.off()

## T1D
pie_labels <- paste0(bf_t1d$Var1, "  ", round(100 * bf_t1d$Freq/sum(bf_t1d$Freq), 2), "%")
pdf('pie_chart_t1d_filtered.pdf', height = 5, width = 5)
pie(bf_t1d$Freq, labels = pie_labels, col =  c("#FFEA00", "#FCF55F", "#FADA5E",
                                               "#FAFA33", "#F4BB44", "#FBEC5D",
                                               "#FFFF00", "#FFFAA0", "#FFE5B4")) #BluYl
dev.off()
