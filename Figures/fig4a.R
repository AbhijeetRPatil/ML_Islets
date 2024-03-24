library("ggplot2")
library("stringr")
setwd("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/")

#########################
## All cells
KEGG_G1_All <- read.csv("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/t1dvsctrl/overall/KEGG.tsv", sep = "\t")
Term <- str_split_fixed(KEGG_G1_All$Term, ":", 2)
colnames(Term) <- c("Term1", "Term2")
KEGG_G1_All$Term <- Term[,2]
KEGG_G1_All$Term <- gsub("\\ -.*","",KEGG_G1_All$Term)
KEGG_G1_All <- KEGG_G1_All[KEGG_G1_All$FDR<0.05,]
colnames(KEGG_G1_All)[2] <- "KEGGPathways"

KEGG_G2_All <- read.csv("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/t1dvsaab/overall/KEGG.tsv", sep = "\t")
Term <- str_split_fixed(KEGG_G2_All$Term, ":", 2)
colnames(Term) <- c("Term1", "Term2")
KEGG_G2_All$Term <- Term[,2]
KEGG_G2_All$Term <- gsub("\\ -.*","",KEGG_G2_All$Term)
KEGG_G2_All <- KEGG_G2_All[KEGG_G2_All$FDR<0.05,]
colnames(KEGG_G2_All)[2] <- "KEGGPathways"


KEGG_G3_All <- read.csv("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/aabvsctrl/overall/KEGG.tsv", sep = "\t")
Term <- str_split_fixed(KEGG_G3_All$Term, ":", 2)
colnames(Term) <- c("Term1", "Term2")
KEGG_G3_All$Term <- Term[,2]
KEGG_G3_All$Term <- gsub("\\ -.*","",KEGG_G3_All$Term)
KEGG_G3_All <- KEGG_G3_All[KEGG_G3_All$FDR<0.05,]
colnames(KEGG_G3_All)[2] <- "KEGGPathways"


## KEGG Combined plot
KEGG_G1_All <- KEGG_G1_All[c(1:10),]
KEGG_G1_All <- cbind("var"=rep("T1D-Ctrl",nrow(KEGG_G1_All)), KEGG_G1_All)

KEGG_G2_All <- KEGG_G2_All[c(1:10),]
KEGG_G2_All <- cbind("var"=rep("T1D-AAB",nrow(KEGG_G2_All)), KEGG_G2_All)

KEGG_G3_All <- KEGG_G3_All[c(1:10),]
KEGG_G3_All <- cbind("var"=rep("AAB-Ctrl",nrow(KEGG_G3_All)), KEGG_G3_All)

KEGG_All <- rbind(KEGG_G1_All, KEGG_G2_All, KEGG_G3_All)

kegg_subplot <- function(dat, keggtitle)
{
  ggplot(data = dat, aes(x = reorder(var, Count), y = reorder(KEGGPathways, -Count), size = Count, color = -log10(FDR)))+
    geom_point()+
    ggtitle(keggtitle)+
    scale_color_gradientn(colors = rainbow(5))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, color="black", size=10, face="bold"), 
          axis.text.x = element_text(color="black", size=10, face="bold"),
          axis.text.y = element_text(color="black", size=10, face="bold"),
          axis.title.x=element_blank(), axis.title.y=element_blank())
}
p8 <- kegg_subplot(KEGG_All, "All cells- KEGG Pathways (FDR < 0.05)")
pdf("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/GO and KEGG Plots/KEGG_All.pdf", width = 8, height = 5)
p8
dev.off()

############################################################################################################################

## Alpha cells
KEGG_G1_All <- read.csv("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/t1dvsctrl/alpha/KEGG.tsv", sep = "\t")
Term <- str_split_fixed(KEGG_G1_All$Term, ":", 2)
colnames(Term) <- c("Term1", "Term2")
KEGG_G1_All$Term <- Term[,2]
KEGG_G1_All$Term <- gsub("\\ -.*","",KEGG_G1_All$Term)
KEGG_G1_All <- KEGG_G1_All[KEGG_G1_All$FDR<0.05,]
colnames(KEGG_G1_All)[2] <- "KEGGPathways"

KEGG_G2_All <- read.csv("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/t1dvsaab/alpha/KEGG.tsv", sep = "\t")
Term <- str_split_fixed(KEGG_G2_All$Term, ":", 2)
colnames(Term) <- c("Term1", "Term2")
KEGG_G2_All$Term <- Term[,2]
KEGG_G2_All$Term <- gsub("\\ -.*","",KEGG_G2_All$Term)
KEGG_G2_All <- KEGG_G2_All[KEGG_G2_All$FDR<0.05,]
colnames(KEGG_G2_All)[2] <- "KEGGPathways"


KEGG_G3_All <- read.csv("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/aabvsctrl/alpha/KEGG.tsv", sep = "\t")
Term <- str_split_fixed(KEGG_G3_All$Term, ":", 2)
colnames(Term) <- c("Term1", "Term2")
KEGG_G3_All$Term <- Term[,2]
KEGG_G3_All$Term <- gsub("\\ -.*","",KEGG_G3_All$Term)
KEGG_G3_All <- KEGG_G3_All[KEGG_G3_All$FDR<0.05,]
colnames(KEGG_G3_All)[2] <- "KEGGPathways"

## KEGG Combined plot
KEGG_G1_All <- KEGG_G1_All[c(1:10),]
KEGG_G1_All <- cbind("var"=rep("T1D-Ctrl",nrow(KEGG_G1_All)), KEGG_G1_All)

KEGG_G2_All <- KEGG_G2_All[c(1:10),]
KEGG_G2_All <- cbind("var"=rep("T1D-AAB",nrow(KEGG_G2_All)), KEGG_G2_All)

KEGG_G3_All <- KEGG_G3_All[c(1:10),]
KEGG_G3_All <- cbind("var"=rep("AAB-Ctrl",nrow(KEGG_G3_All)), KEGG_G3_All)

KEGG_All <- rbind(KEGG_G1_All, KEGG_G2_All, KEGG_G3_All)

kegg_subplot <- function(dat, keggtitle)
{
  ggplot(data = dat, aes(x = reorder(var, Count), y = reorder(KEGGPathways, -Count), size = Count, color = -log10(FDR)))+
    geom_point()+
    ggtitle(keggtitle)+
    scale_color_gradientn(colors = rainbow(5))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, color="black", size=10, face="bold"), 
          axis.text.x = element_text(color="black", size=10, face="bold"),
          axis.text.y = element_text(color="black", size=10, face="bold"),
          axis.title.x=element_blank(), axis.title.y=element_blank())
}
p8 <- kegg_subplot(KEGG_All, "Alpha cells- KEGG Pathways (FDR < 0.05)")
pdf("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/GO and KEGG Plots/KEGG_alpha.pdf", width = 8, height = 5)
p8
dev.off()
############################################################################################################################

## Ductal cells
KEGG_G1_All <- read.csv("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/t1dvsctrl/ductal/KEGG.tsv", sep = "\t")
Term <- str_split_fixed(KEGG_G1_All$Term, ":", 2)
colnames(Term) <- c("Term1", "Term2")
KEGG_G1_All$Term <- Term[,2]
KEGG_G1_All$Term <- gsub("\\ -.*","",KEGG_G1_All$Term)
KEGG_G1_All <- KEGG_G1_All[KEGG_G1_All$FDR<0.05,]
colnames(KEGG_G1_All)[2] <- "KEGGPathways"


KEGG_G2_All <- read.csv("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/t1dvsaab/ductal/KEGG.tsv", sep = "\t")
Term <- str_split_fixed(KEGG_G2_All$Term, ":", 2)
colnames(Term) <- c("Term1", "Term2")
KEGG_G2_All$Term <- Term[,2]
KEGG_G2_All$Term <- gsub("\\ -.*","",KEGG_G2_All$Term)
KEGG_G2_All <- KEGG_G2_All[KEGG_G2_All$FDR<0.05,]
colnames(KEGG_G2_All)[2] <- "KEGGPathways"


KEGG_G3_All <- read.csv("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/aabvsctrl/ductal/KEGG.tsv", sep = "\t")
Term <- str_split_fixed(KEGG_G3_All$Term, ":", 2)
colnames(Term) <- c("Term1", "Term2")
KEGG_G3_All$Term <- Term[,2]
KEGG_G3_All$Term <- gsub("\\ -.*","",KEGG_G3_All$Term)
KEGG_G3_All <- KEGG_G3_All[KEGG_G3_All$FDR<0.05,]
colnames(KEGG_G3_All)[2] <- "KEGGPathways"

## KEGG Combined plot
KEGG_G1_All <- KEGG_G1_All[c(1:10),]
KEGG_G1_All <- cbind("var"=rep("T1D-Ctrl",nrow(KEGG_G1_All)), KEGG_G1_All)

KEGG_G2_All <- KEGG_G2_All[c(1:10),]
KEGG_G2_All <- cbind("var"=rep("T1D-AAB",nrow(KEGG_G2_All)), KEGG_G2_All)

KEGG_G3_All <- KEGG_G3_All[c(1:10),]
KEGG_G3_All <- cbind("var"=rep("AAB-Ctrl",nrow(KEGG_G3_All)), KEGG_G3_All)

KEGG_All <- rbind(KEGG_G1_All, KEGG_G2_All, KEGG_G3_All)

kegg_subplot <- function(dat, keggtitle)
{
  ggplot(data = dat, aes(x = reorder(var, Count), y = reorder(KEGGPathways, -Count), size = Count, color = -log10(FDR)))+
    geom_point()+
    ggtitle(keggtitle)+
    scale_color_gradientn(colors = rainbow(5))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, color="black", size=10, face="bold"), 
          axis.text.x = element_text(color="black", size=10, face="bold"),
          axis.text.y = element_text(color="black", size=10, face="bold"),
          axis.title.x=element_blank(), axis.title.y=element_blank())
}
p8 <- kegg_subplot(KEGG_All, "Ductal cells- KEGG Pathways (FDR < 0.05)")
pdf("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/GO and KEGG Plots/KEGG_ductal.pdf", width = 10, height = 10)
p8
dev.off()
############################################################################################################################
## Immune cells

KEGG_G1_All <- read.csv("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/t1dvsctrl/immune/KEGG.tsv", sep = "\t")
Term <- str_split_fixed(KEGG_G1_All$Term, ":", 2)
colnames(Term) <- c("Term1", "Term2")
KEGG_G1_All$Term <- Term[,2]
KEGG_G1_All$Term <- gsub("\\ -.*","",KEGG_G1_All$Term)
KEGG_G1_All <- KEGG_G1_All[KEGG_G1_All$FDR<0.05,]
colnames(KEGG_G1_All)[2] <- "KEGGPathways"

KEGG_G2_All <- read.csv("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/t1dvsaab/immune/KEGG.tsv", sep = "\t")
Term <- str_split_fixed(KEGG_G2_All$Term, ":", 2)
colnames(Term) <- c("Term1", "Term2")
KEGG_G2_All$Term <- Term[,2]
KEGG_G2_All$Term <- gsub("\\ -.*","",KEGG_G2_All$Term)
KEGG_G2_All <- KEGG_G2_All[KEGG_G2_All$FDR<0.05,]
colnames(KEGG_G2_All)[2] <- "KEGGPathways"


KEGG_G3_All <- read.csv("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/aabvsctrl/immune/KEGG.tsv", sep = "\t")
Term <- str_split_fixed(KEGG_G3_All$Term, ":", 2)
colnames(Term) <- c("Term1", "Term2")
KEGG_G3_All$Term <- Term[,2]
KEGG_G3_All$Term <- gsub("\\ -.*","",KEGG_G3_All$Term)
KEGG_G3_All <- KEGG_G3_All[KEGG_G3_All$FDR<0.05,]
colnames(KEGG_G3_All)[2] <- "KEGGPathways"

## KEGG Combined plot
KEGG_G1_All <- KEGG_G1_All[c(1:10),]
KEGG_G1_All <- cbind("var"=rep("T1D-Ctrl",nrow(KEGG_G1_All)), KEGG_G1_All)

KEGG_G2_All <- KEGG_G2_All[c(1:10),]
KEGG_G2_All <- cbind("var"=rep("T1D-AAB",nrow(KEGG_G2_All)), KEGG_G2_All)

KEGG_G3_All <- KEGG_G3_All[c(1:10),]
KEGG_G3_All <- cbind("var"=rep("AAB-Ctrl",nrow(KEGG_G3_All)), KEGG_G3_All)

KEGG_All <- rbind(KEGG_G1_All, KEGG_G2_All, KEGG_G3_All)

kegg_subplot <- function(dat, keggtitle)
{
  ggplot(data = dat, aes(x = reorder(var, Count), y = reorder(KEGGPathways, -Count), size = Count, color = -log10(FDR)))+
    geom_point()+
    ggtitle(keggtitle)+
    scale_color_gradientn(colors = rainbow(5))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, color="black", size=10, face="bold"), 
          axis.text.x = element_text(color="black", size=10, face="bold"),
          axis.text.y = element_text(color="black", size=10, face="bold"),
          axis.title.x=element_blank(), axis.title.y=element_blank())
}
p8 <- kegg_subplot(KEGG_All, "Immune cells- KEGG Pathways (FDR < 0.05)")
pdf("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/GO and KEGG Plots/KEGG_Immune.pdf", width = 10, height = 10)
p8
dev.off()

############################################################################################################################

## Beta cells
KEGG_G1_All <- read.csv("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/t1dvsctrl/beta/KEGG.tsv", sep = "\t")
Term <- str_split_fixed(KEGG_G1_All$Term, ":", 2)
colnames(Term) <- c("Term1", "Term2")
KEGG_G1_All$Term <- Term[,2]
KEGG_G1_All$Term <- gsub("\\ -.*","",KEGG_G1_All$Term)
KEGG_G1_All <- KEGG_G1_All[KEGG_G1_All$FDR<0.05,]
colnames(KEGG_G1_All)[2] <- "KEGGPathways"

KEGG_G2_All <- read.csv("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/t1dvsaab/beta/KEGG.tsv", sep = "\t")
Term <- str_split_fixed(KEGG_G2_All$Term, ":", 2)
colnames(Term) <- c("Term1", "Term2")
KEGG_G2_All$Term <- Term[,2]
KEGG_G2_All$Term <- gsub("\\ -.*","",KEGG_G2_All$Term)
KEGG_G2_All <- KEGG_G2_All[KEGG_G2_All$FDR<0.05,]
colnames(KEGG_G2_All)[2] <- "KEGGPathways"


KEGG_G3_All <- read.csv("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/aabvsctrl/beta/KEGG.tsv", sep = "\t")
Term <- str_split_fixed(KEGG_G3_All$Term, ":", 2)
colnames(Term) <- c("Term1", "Term2")
KEGG_G3_All$Term <- Term[,2]
KEGG_G3_All$Term <- gsub("\\ -.*","",KEGG_G3_All$Term)
KEGG_G3_All <- KEGG_G3_All[KEGG_G3_All$FDR<0.05,]
colnames(KEGG_G3_All)[2] <- "KEGGPathways"

## KEGG Combined plot
KEGG_G1_All <- KEGG_G1_All[c(1:10),]
KEGG_G1_All <- cbind("var"=rep("T1D-Ctrl",nrow(KEGG_G1_All)), KEGG_G1_All)

KEGG_G2_All <- KEGG_G2_All[c(1:10),]
KEGG_G2_All <- cbind("var"=rep("T1D-AAB",nrow(KEGG_G2_All)), KEGG_G2_All)

KEGG_G3_All <- KEGG_G3_All[c(1:10),]
KEGG_G3_All <- cbind("var"=rep("AAB-Ctrl",nrow(KEGG_G3_All)), KEGG_G3_All)

KEGG_All <- rbind(KEGG_G1_All, KEGG_G2_All, KEGG_G3_All)

kegg_subplot <- function(dat, keggtitle)
{
  ggplot(data = dat, aes(x = reorder(var, Count), y = reorder(KEGGPathways, -Count), size = Count, color = -log10(FDR)))+
    geom_point()+
    ggtitle(keggtitle)+
    scale_color_gradientn(colors = rainbow(5))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, color="black", size=10, face="bold"), 
          axis.text.x = element_text(color="black", size=10, face="bold"),
          axis.text.y = element_text(color="black", size=10, face="bold"),
          axis.title.x=element_blank(), axis.title.y=element_blank())
}
p8 <- kegg_subplot(KEGG_All, "Beta cells- KEGG Pathways (FDR < 0.05)")
pdf("C:/Users/abhij/Desktop/Spring 2022/ML paper/Pathways/pathways/GO and KEGG Plots/KEGG_Beta.pdf", width = 10, height = 10)
p8
dev.off()
############################################################################################################################


