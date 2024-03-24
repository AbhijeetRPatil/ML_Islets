rm(list = ls())
library("ggplot2")
library("RobustRankAggreg")
library("dplyr")
library("org.Hs.eg.db")
library("clusterProfiler")
library("enrichplot")
library("GOSemSim")
library("ggplot2")
library("ggridges")

setwd("C:/Users/abhij/Desktop/RRA_New")

conv.gene.names <- function(df)
{
  tmp <- strsplit(df@result$core_enrichment, "/")
  new.tmp <- list()
  for(i in 1:length(tmp))
  {
    new.tmp[[i]] <- group1_tmp_list_fin$Genes[match(tmp[[i]], group1_tmp_list_fin$EntrezIDs$ENTREZID)]  
  }
  df@result$GeneNames <- vector(mode = "list",length=nrow(df@result)) 
  
  for(i in 1:length(new.tmp)){
    df@result$GeneNames[[i]] <- new.tmp[[i]]
  }
  df@result$GeneNames <- as.character(df@result$GeneNames)
  df@result$GeneNames <- gsub('\"', "", df@result$GeneNames, fixed = TRUE)
  df@result$GeneNames <- gsub('c(', "", df@result$GeneNames, fixed = TRUE)
  df@result$GeneNames <- gsub(')', "", df@result$GeneNames, fixed = TRUE)
  return(df)
}

## All t1d vs ctrl
xgb_t1dvsctrl <- read.csv("./rankedlists/RRA_Combined/T1DvsCTL_CellType_RRA.csv")
group1_tmp_list_fin <- data.frame(xgb_t1dvsctrl$X)
colnames(group1_tmp_list_fin) <- "Genes" 
hs <- org.Hs.eg.db
my.symbols <- group1_tmp_list_fin$Genes
group1_tmp_list_fin$EntrezIDs <- AnnotationDbi::select(hs, 
                                                       keys = my.symbols,
                                                       columns = c("ENTREZID", "SYMBOL"),
                                                       keytype = "SYMBOL")
# write.csv(group1_tmp_list_fin, "./rankedlists/xgb_t1dvsctrl_RRA_GL.csv")

group1_tmp_list_fin$Score <- -log10(xgb_t1dvsctrl$Score)
# group1_tmp_list_fin$Score[1] <- 500
x <- group1_tmp_list_fin$Score
names(x) <- group1_tmp_list_fin$EntrezIDs$ENTREZID

All_GO <- gseGO(geneList     = x,
                OrgDb        = org.Hs.eg.db,
                ont          = "BP",
                # minGSSize    = 100,
                # maxGSSize    = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "none",
                verbose      = FALSE)

All_GO <- conv.gene.names(All_GO)
write.csv(All_GO@result, "./rankedlists/RRA_Combined/T1DvsCTL_RRA_Combined.csv", row.names = F)

d <- godata('org.Hs.eg.db', ont="BP")
All_GO2 <- pairwise_termsim(All_GO, method="Wang", semData = d)
pdf(file = "./rankedlists/RRA_Combined/T1DvsCTL_All_GO_Cluster.pdf", width = 12, height = 12)
emapplot_cluster(All_GO2)
dev.off()

# # rm(list=ls())
# library("ggplot2")
# library("stringr")
# # setwd("C:/Users/abhij/Desktop/RRA_New")
# 
# #########################
# ## All cells
# GO_G1_All <- read.csv("./rankedlists/RRA_Combined/T1DvsCTL_RRA_Combined.csv")
# 
# GO_G1_All$Term <- GO_G1_All$Description
# GO_G1_All$Term <- gsub("\\ -.*","",GO_G1_All$Term)
# colnames(GO_G1_All)[2] <- "GOPathways"
# 
# ## GO Combined plot
# GO_G1_All <- cbind("var"=rep(" All Cells",nrow(GO_G1_All)), GO_G1_All)
# GO_G1_All <- GO_G1_All[1:20,]
# ############################################################################################################################
# 
# ############################################################################################################################
# 
# GO_All <- rbind(GO_G1_All)
# 
# GO_subplot <- function(dat, GOtitle)
# {
#   ggplot(data = dat, aes(x = var, y = GOPathways, size = setSize, color = -log10(p.adjust)))+
#     geom_point()+
#     ggtitle(GOtitle)+
#     scale_color_gradientn(colors = c("grey", "red") )+#rainbow(5) c("#001219", "#e9d8a6")
#     theme_bw()+
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           plot.title = element_text(hjust=0.5, color="black", size=9, face="bold"), 
#           axis.text.x = element_text(angle=45, hjust=1, color="black", size=9, face="bold"),
#           axis.text.y = element_text(color="black", size=9, face="bold"),
#           axis.title.x=element_blank(), axis.title.y=element_blank())
# }
# p8 <- GO_subplot(GO_All, "T1D vs. CTL GO:BP (FDR < 0.05) RRA Combined")
# pdf("./rankedlists/t1dvsctrl/GO_T1DvsCTL_Top20_FDR.pdf", width = 15, height = 20)
# p8
# dev.off()
