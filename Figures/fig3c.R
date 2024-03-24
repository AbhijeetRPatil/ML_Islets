rm(list = ls())
setwd("C:/Users/abhij/Desktop/RRA_New/rankedlists/")

################################################################################################################################
################################################################
## Acinar
################################################################
Acinar <- read.csv("./acinar_xgb_t1dvsctrl_RRA_GL.csv")
Acinar <- Acinar$Genes
Acinar <- as.character(unique(Acinar))
################################################################
## Alpha
################################################################
Alpha <- read.csv("./alpha_xgb_t1dvsctrl_RRA_GL.csv")
Alpha <- Alpha$Genes
Alpha <- as.character(unique(Alpha))
################################################################
## Beta
################################################################
Beta <- read.csv("./beta_xgb_t1dvsctrl_RRA_GL.csv")
Beta <- Beta$Genes
Beta <- as.character(unique(Beta))
################################################################
## Ductal
################################################################
Ductal <- read.csv("ductal_xgb_t1dvsctrl_RRA_GL.csv")
Ductal <- Ductal$Genes
Ductal <- as.character(unique(Ductal))
################################################################
## Delta
################################################################
Delta <- read.csv("delta_xgb_t1dvsctrl_RRA_GL.csv")
Delta <- Delta$Genes
Delta <- as.character(unique(Delta))
################################################################
## Immune
################################################################
Immune <- read.csv("immune_xgb_t1dvsctrl_RRA_GL.csv")
Immune <- Immune$Genes
Immune <- as.character(unique(Immune))
################################################################
## Endothelial
################################################################
Endothelial <- read.csv("endothelial_xgb_t1dvsctrl_RRA_GL.csv")
Endothelial <- Endothelial$Genes
Endothelial <- as.character(unique(Endothelial))
################################################################
## Stellates
################################################################
Stellates <- read.csv("stellates_xgb_t1dvsctrl_RRA_GL.csv")
Stellates <- Stellates$Genes
Stellates <- as.character(unique(Stellates))
################################################################################################################################
################################################################################################################################
#########################################################################################################
### Create dataframe

na.pad <- function(x,len){
  x[1:len]
}

makePaddedDataFrame <- function(l,...){
  maxlen <- max(sapply(l,length))
  data.frame(lapply(l,na.pad,len=maxlen),...)
}
#########################################################################################################
## Gene Integration
#########################################################################################################

## rlfs
group1_rlfs_genelist1 <- list(Acinar = Acinar,
                              Alpha = Alpha, 
                              Beta = Beta,
                              Ductal = Ductal, 
                              Delta = Delta,
                              Immune = Immune,
                              Endothelial = Endothelial,
                              Stellates = Stellates)

group1_rlfs <- makePaddedDataFrame(group1_rlfs_genelist1)

# write.csv(group1_rlfs,'group1_rlfs.csv')
# group1_rlfs <- read.csv("group1_rlfs.csv")

group1_rlfs <- sapply(group1_rlfs, as.character)
group1_rlfs[is.na(group1_rlfs)] <- ""
group1_rlfs <- as.data.frame(group1_rlfs)
# write.csv(group1_rlfs,'group1_rlfs.csv')

group1_rlfs_ll = c(Acinar = Acinar,
                   Alpha = Alpha, 
                   Beta = Beta,
                   Ductal = Ductal, 
                   Delta = Delta,
                   Immune = Immune,
                   Endothelial = Endothelial,
                   Stellates = Stellates)

group1_rlfs_N1 = length(names(table(group1_rlfs_ll)))

library("RobustRankAggreg")
library("dplyr")
group1_rlfs_list <- aggregateRanks(glist = group1_rlfs_genelist1, N = group1_rlfs_N1)

group1_rlfs_list_fin <- filter(group1_rlfs_list, Score<0.05)
write.csv(group1_rlfs_list_fin, "T1DvsCTL_CellType_RRA.csv")
####################################################################################################################

RRA_Combined <- group1_rlfs_list_fin$Name
group1_rlfs_genelist1 <- list(Acinar = Acinar,
                              Alpha = Alpha, 
                              Beta = Beta,
                              Ductal = Ductal, 
                              Delta = Delta,
                              Immune = Immune,
                              Endothelial = Endothelial,
                              Stellates = Stellates,
                              RRA_Combined = RRA_Combined)

## UpSetR Plots
common_rlfs <- Reduce(union, group1_rlfs_genelist1)
data.names = c("Acinar" , "Alpha", "Beta", 
               "Ductal", "Delta", "Immune",
               "Endothelial", "Stellates",
               "RRA_Combined")

mat0 <- matrix(0, nrow = length(common_rlfs), ncol = length(data.names))
colnames(mat0) <- data.names
rownames(mat0) <- common_rlfs
mat0[] <- do.call(cbind, lapply(group1_rlfs_genelist1, function(i) + (rownames(mat0) %in% i)))
# write.csv(as.data.frame(mat0), "rlfs_group1.csv")
mat0 <- as.data.frame(mat0)

#########################################################################################################################

library("UpSetR")
library("ComplexUpset")
library("ggplot2")

### UpSetR

# tiff(file="Venn_T1D_CTL_Upset.png",
#      width=8, height=6, units="in", res=600)
# upset(mat0, sets = data.names, nintersects=20, mb.ratio = c(0.7, 0.3), order.by = "freq", mainbar.y.label = "Gene Intersections",
#       sets.x.label = "Significant Genes")
# dev.off()

#####################################################################
# upset(mat0,
#       c("Acinar", "Alpha", "Beta", "Ductal", "Delta", "Immune", "Endothelial", "Stellates", "RRA_Combined"),
#       width_ratio = 0.1,
#       intersections = list(
#         c("Acinar", "Alpha"),
#         c("Acinar", "Alpha", "Beta")
#       ))
# dev.off()

# Reduce(intersect, list(group1_rlfs_genelist1$Acinar, group1_rlfs_genelist1$Alpha, group1_rlfs_genelist1$Beta, 
#                        group1_rlfs_genelist1$Ductal, group1_rlfs_genelist1$Delta, group1_rlfs_genelist1$Immune,
#                        group1_rlfs_genelist1$Endothelial, group1_rlfs_genelist1$Stellates,
#                        group1_rlfs_genelist1$RRA_Combined))

pdf(file="Venn_T1D_CTL_Upset.pdf", width=10, height=6)
upset(
  mat0,
  # sets = data.names,
  c("Acinar", "Alpha", "Beta", "Ductal", "Delta", "Immune", "Endothelial", "Stellates", "RRA_Combined"),
  n_intersections=20,
  intersections='all',
  name = "",
  queries=list(
    upset_query(set='Acinar', fill='#009392'),
    upset_query(set='Alpha', fill='#39B185'),
    upset_query(set='Beta', fill='#9CCB86'),
    upset_query(set='Ductal', fill='#E9E29C'),
    upset_query(set='Delta', fill='#EEB479'),
    upset_query(set='Immune', fill='#CF597E'),
    upset_query(set='Endothelial', fill='#DF098E'),
    upset_query(set='Stellates', fill='#aa66cd'),
    upset_query(set='RRA_Combined', fill='#FDCE30')
  ),
  base_annotations=list(
    'Intersection size'=(
      intersection_size(
        bar_number_threshold=1,  # show all numbers on top of bars
        width=0.5,   # reduce width of the bars
      )
      # 'Intersection size'=(intersection_size(text_mapping=aes(label=paste0(round(
      #   !!get_size_mode('exclusive_intersection')/!!get_size_mode('inclusive_union') * 100
      # ), '%'))))
      
      # add some space on the top of the bars
      + scale_y_continuous(expand=expansion(mult=c(0, 0.1)))
      + theme(
        # hide grid lines
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # show axis lines
        axis.line=element_line(colour='black')
      )
    )
  ),
  stripes=upset_stripes(
    geom=geom_segment(size=12),  # make the stripes larger
    colors=c('grey95', 'white')
  ),
  # to prevent connectors from getting the colorured
  # use `fill` instead of `color`, together with `shape='circle filled'`
  matrix=intersection_matrix(
    geom=geom_point(
      shape='circle filled',
      size=3.5,
      stroke=0.45
    )
  ),
  set_sizes=(
    upset_set_size(geom=geom_bar(width=0.5))
    + theme(
      axis.line.x=element_line(colour='black'),
      axis.ticks.x=element_line()
    )
  )
  ,
  sort_sets='descending',
  sort_intersections='descending'
)
dev.off()



# nsets=5,
# sets= c("Acinar", "Alpha", "Beta", "Ductal","RRA_Combined"),
# sets= c(c("Acinar", "Alpha", "Beta", "Ductal"), c("RRA_Combined", "Beta"), c("RRA_Combined", "Ductal")),
# nintersects=20,
# intersections=list("Acinar", "Alpha", "Beta", "Ductal", "Delta", "Immune", "Endothelial", "Stellates",
#                    c("Acinar", "Alpha"),
#                    # c("Acinar", "Alpha", "Beta", "Ductal", "Delta", "Immune", "Endothelial", "Stellates"),
#                    c("Acinar", "Alpha", "Beta"),
#                    c("Acinar", "Alpha", "Beta", "Ductal"),
#                    c("Acinar", "Alpha", "Beta", "Ductal", "Delta"),
#                    c("Acinar", "Alpha", "Beta", "Ductal", "Delta", "Immune"),
#                    c("Acinar", "Alpha", "Beta", "Ductal", "Delta", "Immune", "Endothelial"),
#                    c("Acinar", "Alpha", "Beta", "Ductal", "Delta", "Immune", "Endothelial", "Stellates"),
#                    c("Acinar", "Alpha", "Beta", "Ductal", "Delta", "Immune", "Endothelial", "Stellates", "RRA_Combined")
#                    ),
# group_by = 'sets',