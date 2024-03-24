rm(list=ls())
library("ggplot2")
library("stringr")
setwd("C:/Users/abhij/Desktop/RRA_New")

#########################
## All cells
KEGG_G1_All <- read.csv("./rankedlists/aabvsctrl/All_KEGG.csv")

KEGG_G1_All$Term <- KEGG_G1_All$Description
KEGG_G1_All$Term <- gsub("\\ -.*","",KEGG_G1_All$Term)
colnames(KEGG_G1_All)[2] <- "KEGGPathways"

## KEGG Combined plot
KEGG_G1_All <- cbind("var"=rep(" All Cells",nrow(KEGG_G1_All)), KEGG_G1_All)

############################################################################################################################

## Alpha cells
KEGG_G1_Alpha <- read.csv("./rankedlists/aabvsctrl/alpha_KEGG.csv")
KEGG_G1_Alpha$Term <- KEGG_G1_Alpha$Description
KEGG_G1_Alpha$Term <- gsub("\\ -.*","",KEGG_G1_Alpha$Term)
colnames(KEGG_G1_Alpha)[2] <- "KEGGPathways"

## KEGG Combined plot
KEGG_G1_Alpha <- cbind("var"=rep("Alpha",nrow(KEGG_G1_Alpha)), KEGG_G1_Alpha)

############################################################################################################################

## Ductal cells
KEGG_G1_Ductal <- read.csv("./rankedlists/aabvsctrl/ductal_KEGG.csv")
KEGG_G1_Ductal$Term <- KEGG_G1_Ductal$Description
KEGG_G1_Ductal$Term <- gsub("\\ -.*","",KEGG_G1_Ductal$Term)
colnames(KEGG_G1_Ductal)[2] <- "KEGGPathways"

## KEGG Combined plot
KEGG_G1_Ductal <- cbind("var"=rep("Ductal",nrow(KEGG_G1_Ductal)), KEGG_G1_Ductal)

############################################################################################################################
## Immune cells

KEGG_G1_Immune <- read.csv("./rankedlists/aabvsctrl/immune_KEGG.csv")
KEGG_G1_Immune$Term <- KEGG_G1_Immune$Description
KEGG_G1_Immune$Term <- gsub("\\ -.*","",KEGG_G1_Immune$Term)
colnames(KEGG_G1_Immune)[2] <- "KEGGPathways"

## KEGG Combined plot
KEGG_G1_Immune <- cbind("var"=rep("Immune",nrow(KEGG_G1_Immune)), KEGG_G1_Immune)

############################################################################################################################

## Beta cells
KEGG_G1_Beta <- read.csv("./rankedlists/aabvsctrl/beta_KEGG.csv")
KEGG_G1_Beta$Term <- KEGG_G1_Beta$Description
KEGG_G1_Beta$Term <- gsub("\\ -.*","",KEGG_G1_Beta$Term)
colnames(KEGG_G1_Beta)[2] <- "KEGGPathways"

## KEGG Combined plot
KEGG_G1_Beta <- cbind("var"=rep("Beta",nrow(KEGG_G1_Beta)), KEGG_G1_Beta)

############################################################################################################################

## Endothelial cells
KEGG_G1_Endothelial <- read.csv("./rankedlists/aabvsctrl/endothelial_KEGG.csv")
KEGG_G1_Endothelial$Term <- KEGG_G1_Endothelial$Description
KEGG_G1_Endothelial$Term <- gsub("\\ -.*","",KEGG_G1_Endothelial$Term)
colnames(KEGG_G1_Endothelial)[2] <- "KEGGPathways"

## KEGG Combined plot
KEGG_G1_Endothelial <- cbind("var"=rep("Endothelial",nrow(KEGG_G1_Endothelial)), KEGG_G1_Endothelial)

############################################################################################################################

## Acinar cells
KEGG_G1_Acinar <- read.csv("./rankedlists/aabvsctrl/acinar_KEGG.csv")
KEGG_G1_Acinar$Term <- KEGG_G1_Acinar$Description
KEGG_G1_Acinar$Term <- gsub("\\ -.*","",KEGG_G1_Acinar$Term)
colnames(KEGG_G1_Acinar)[2] <- "KEGGPathways"

## KEGG Combined plot
KEGG_G1_Acinar <- cbind("var"=rep("Acinar",nrow(KEGG_G1_Acinar)), KEGG_G1_Acinar)

############################################################################################################################

## Stellates cells
KEGG_G1_Stellates <- read.csv("./rankedlists/aabvsctrl/stellates_KEGG.csv")
KEGG_G1_Stellates$Term <- KEGG_G1_Stellates$Description
KEGG_G1_Stellates$Term <- gsub("\\ -.*","",KEGG_G1_Stellates$Term)
colnames(KEGG_G1_Stellates)[2] <- "KEGGPathways"

## KEGG Combined plot
KEGG_G1_Stellates <- cbind("var"=rep("Stellates",nrow(KEGG_G1_Stellates)), KEGG_G1_Stellates)

############################################################################################################################

## Delta cells
KEGG_G1_Delta <- read.csv("./rankedlists/aabvsctrl/delta_KEGG.csv")
KEGG_G1_Delta$Term <- KEGG_G1_Delta$Description
KEGG_G1_Delta$Term <- gsub("\\ -.*","",KEGG_G1_Delta$Term)
colnames(KEGG_G1_Delta)[2] <- "KEGGPathways"

## KEGG Combined plot
KEGG_G1_Delta <- cbind("var"=rep("Delta",nrow(KEGG_G1_Delta)), KEGG_G1_Delta)

############################################################################################################################

KEGG_All <- rbind(KEGG_G1_All, KEGG_G1_Acinar, KEGG_G1_Alpha, KEGG_G1_Beta, KEGG_G1_Delta, KEGG_G1_Ductal, 
                  KEGG_G1_Endothelial, KEGG_G1_Immune, KEGG_G1_Stellates)

kegg_subplot <- function(dat, keggtitle)
{
  ggplot(data = dat, aes(x = var, y = reorder(KEGGPathways,setSize), size = setSize, color = -log10(p.adjust)))+
    geom_point()+
    ggtitle(keggtitle)+
    scale_color_gradientn(colors = c("grey", "red"))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust=0.5, color="black", size=9, face="bold"), 
          axis.text.x = element_text(angle=45, hjust=1, color="black", size=9, face="bold"),
          axis.text.y = element_text(color="black", size=9, face="bold"),
          axis.title.x=element_blank(), axis.title.y=element_blank())
}
p8 <- kegg_subplot(KEGG_All, "AAB vs CTL- KEGG Pathways (FDR < 0.05)")
pdf("./rankedlists/aabvsctrl/KEGG_AABvsCTL_ALLFDR.pdf", width = 8, height = 12)
p8
dev.off()
