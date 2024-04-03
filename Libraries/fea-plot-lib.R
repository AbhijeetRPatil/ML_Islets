########################## KEGG plot function #########################################
## PVAL
KEGG_plot <- function(df, plot_title)
{
  ggplot(df, aes(x = reorder(KEGG_Pathways, Count),
                 y = Count,
                 fill = -log10(PValue))) +
    ggtitle(plot_title) +
    geom_bar(stat = "identity",  position = "dodge") + 
    coord_flip() + scale_fill_continuous(type = "viridis") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                plot.title = element_text(hjust = 0.5, color="black", size=15, face="bold"), 
                                                                                axis.text.x = element_text(color="black", size=10, face="bold"),
                                                                                axis.text.y = element_text(color="black", size=10, face="bold"),
                                                                                axis.title.x=element_blank(), axis.title.y=element_blank())
}
## FDR
KEGG_plot_FDR <- function(df, plot_title)
{
  ggplot(df, aes(x = reorder(KEGG_Pathways, Count),
                 y = Count,
                 fill = -log10(FDR))) +
    ggtitle(plot_title) +
    geom_bar(stat = "identity",  position = "dodge") + 
    coord_flip() + scale_fill_continuous(type = "viridis") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                plot.title = element_text(hjust = 0.5, color="black", size=15, face="bold"), 
                                                                                axis.text.x = element_text(color="black", size=10, face="bold"),
                                                                                axis.text.y = element_text(color="black", size=10, face="bold"),
                                                                                axis.title.x=element_blank(), axis.title.y=element_blank())
}

#KEGG_plot_FDR <- function(df, plot_title)
#{
#  ggplot(df) +
#    ggtitle(plot_title) +
#    geom_col(aes(Count, KEGG_Pathways), fill = -log10(FDR)) + #,  position = "dodge") +
#    coord_flip() + scale_fill_continuous(type = "viridis") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                                plot.title = element_text(hjust = 0.5, color="black", size=15, face="bold"),
#                                                                                axis.text.x = element_text(color="black", size=10, face="bold"),
#                                                                                axis.text.y = element_text(color="black", size=10, face="bold"),
#                                                                                axis.title.x=element_blank(), axis.title.y=element_blank())
#}



########################## REACTOME plot function #########################################
## PVAL
REACTOME_plot <- function(df, plot_title)
{
  ggplot(df, aes(x = reorder(REACTOME_Pathways, Count),
                 y = Count,
                 fill = -log10(PValue))) +
    ggtitle(plot_title) +
    geom_bar(stat = "identity",  position = "dodge") + 
    coord_flip() + scale_fill_continuous(type = "viridis") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                plot.title = element_text(hjust = 0.5, color="black", size=15, face="bold"), 
                                                                                axis.text.x = element_text(color="black", size=10, face="bold"),
                                                                                axis.text.y = element_text(color="black", size=10, face="bold"),
                                                                                axis.title.x=element_blank(), axis.title.y=element_blank())
}
## FDR
REACTOME_plot_FDR <- function(df, plot_title)
{
  ggplot(df, aes(x = reorder(REACTOME_Pathways, Count),
                 y = Count,
                 fill = -log10(FDR))) +
    ggtitle(plot_title) +
    geom_bar(stat = "identity",  position = "dodge") + 
    coord_flip() + scale_fill_continuous(type = "viridis") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                plot.title = element_text(hjust = 0.5, color="black", size=15, face="bold"), 
                                                                                axis.text.x = element_text(color="black", size=10, face="bold"),
                                                                                axis.text.y = element_text(color="black", size=10, face="bold"),
                                                                                axis.title.x=element_blank(), axis.title.y=element_blank())
}
########################## GO plot function #########################################
## PVAL
GO_plot <- function(df, plot_title)
{
  ggplot(df, aes(x = reorder(GO_Term, Count),
                 y = Count,
                 fill = -log10(PValue))) +
    ggtitle(plot_title) +
    geom_bar(stat = "identity",  position = "dodge") + 
    coord_flip() + scale_fill_continuous(type = "viridis") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                plot.title = element_text(hjust = 0.5, color="black", size=15, face="bold"), 
                                                                                axis.text.x = element_text(color="black", size=10, face="bold"),
                                                                                axis.text.y = element_text(color="black", size=10, face="bold"),
                                                                                axis.title.x=element_blank(), axis.title.y=element_blank())
}
## PVAL
GO_plot_FDR <- function(df, plot_title)
{
  ggplot(df, aes(x = reorder(GO_Term, Count),
                 y = Count,
                 fill = -log10(FDR))) +
    ggtitle(plot_title) +
    geom_bar(stat = "identity",  position = "dodge") + 
    coord_flip() + scale_fill_continuous(type = "viridis") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                plot.title = element_text(hjust = 0.5, color="black", size=15, face="bold"), 
                                                                                axis.text.x = element_text(color="black", size=10, face="bold"),
                                                                                axis.text.y = element_text(color="black", size=10, face="bold"),
                                                                                axis.title.x=element_blank(), axis.title.y=element_blank())
}
