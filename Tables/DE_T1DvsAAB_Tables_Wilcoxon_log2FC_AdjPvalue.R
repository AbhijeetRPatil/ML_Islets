rm(list=ls())
library('dplyr')
setwd("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/")

DE_T1DvsAAB <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/DE_T1DvsAAB.rds")
DE_T1DvsAAB <- DE_T1DvsAAB[abs(DE_T1DvsAAB$avg_log2FC) >= 0.5 & DE_T1DvsAAB$p_val_adj<0.05,]
# write.xlsx(DE_T1DvsAAB, 'Allcells_T1DvsAAB_Wilcox.xlsx')

DE_T1DvsAAB_acinar <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/DE_T1DvsAAB_acinar.rds")
DE_T1DvsAAB_acinar <- DE_T1DvsAAB_acinar[abs(DE_T1DvsAAB_acinar$avg_log2FC) >= 0.5 & DE_T1DvsAAB_acinar$p_val_adj<0.05,]
# writexl::write_xlsx(DE_T1DvsAAB_acinar, 'Acinar_T1DvsAAB_Wilcox.xlsx')

DE_T1DvsAAB_alpha <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/DE_T1DvsAAB_alpha.rds")
DE_T1DvsAAB_alpha <- DE_T1DvsAAB_alpha[abs(DE_T1DvsAAB_alpha$avg_log2FC) >= 0.5 & DE_T1DvsAAB_alpha$p_val_adj<0.05,]
# writexl::write_xlsx(DE_T1DvsAAB_alpha, 'Alpha_T1DvsAAB_Wilcox.xlsx')

DE_T1DvsAAB_beta <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/DE_T1DvsAAB_beta.rds")
DE_T1DvsAAB_beta <- DE_T1DvsAAB_beta[abs(DE_T1DvsAAB_beta$avg_log2FC) >= 0.5 & DE_T1DvsAAB_beta$p_val_adj<0.05,]
# writexl::write_xlsx(DE_T1DvsAAB_beta, 'Beta_T1DvsAAB_Wilcox.xlsx')

DE_T1DvsAAB_delta <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/DE_T1DvsAAB_delta.rds")
# writexl::write_xlsx(DE_T1DvsAAB_delta, 'Delta_T1DvsAAB_Wilcox.xlsx')
DE_T1DvsAAB_delta <- DE_T1DvsAAB_delta[abs(DE_T1DvsAAB_delta$avg_log2FC) >= 0.5 & DE_T1DvsAAB_delta$p_val_adj<0.05,]

DE_T1DvsAAB_ductal <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/DE_T1DvsAAB_ductal.rds")
DE_T1DvsAAB_ductal <- DE_T1DvsAAB_ductal[abs(DE_T1DvsAAB_ductal$avg_log2FC) >= 0.5 & DE_T1DvsAAB_ductal$p_val_adj<0.05,]
# writexl::write_xlsx(DE_T1DvsAAB_ductal, 'Ductal_T1DvsAAB_Wilcox.xlsx')

DE_T1DvsAAB_immune <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/DE_T1DvsAAB_immune.rds")
DE_T1DvsAAB_immune <- DE_T1DvsAAB_immune[abs(DE_T1DvsAAB_immune$avg_log2FC) >= 0.5 & DE_T1DvsAAB_immune$p_val_adj<0.05,]
# writexl::write_xlsx(DE_T1DvsAAB_immune, 'Immune_T1DvsAAB_Wilcox.xlsx')

DE_T1DvsAAB_endothelial <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/DE_T1DvsAAB_endothelial.rds")
DE_T1DvsAAB_endothelial <- DE_T1DvsAAB_endothelial[abs(DE_T1DvsAAB_endothelial$avg_log2FC) >= 0.5 & DE_T1DvsAAB_endothelial$p_val_adj<0.05,]
# writexl::write_xlsx(DE_T1DvsAAB_endothelial, 'Endothelial_T1DvsAAB_Wilcox.xlsx')

DE_T1DvsAAB_stellates <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/DE_T1DvsAAB_stellates.rds")
DE_T1DvsAAB_stellates <- DE_T1DvsAAB_stellates[abs(DE_T1DvsAAB_stellates$avg_log2FC) >= 0.5 & DE_T1DvsAAB_stellates$p_val_adj<0.05,]
# writexl::write_xlsx(DE_T1DvsAAB_stellates, 'Stellates_T1DvsAAB_Wilcox.xlsx')

library('xlsx')
write.xlsx(DE_T1DvsAAB, file="DE_T1DvsAAB_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="All Cells-DE T1DvsAAB", row.names=T)
write.xlsx(DE_T1DvsAAB_acinar, file="DE_T1DvsAAB_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="Acinar-DE T1DvsAAB", append=TRUE, row.names=T)
write.xlsx(DE_T1DvsAAB_alpha, file="DE_T1DvsAAB_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="Alpha-DE T1DvsAAB", append=TRUE, row.names=T)
write.xlsx(DE_T1DvsAAB_beta, file="DE_T1DvsAAB_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="Beta-DE T1DvsAAB", append=TRUE, row.names=T)
write.xlsx(DE_T1DvsAAB_delta, file="DE_T1DvsAAB_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="Delta-DE T1DvsAAB", append=TRUE, row.names=T)
write.xlsx(DE_T1DvsAAB_ductal, file="DE_T1DvsAAB_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="Ductal-DE T1DvsAAB", append=TRUE, row.names=T)
write.xlsx(DE_T1DvsAAB_endothelial, file="DE_T1DvsAAB_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="Endothelial-DE T1DvsAAB", append=TRUE, row.names=T)
write.xlsx(DE_T1DvsAAB_immune, file="DE_T1DvsAAB_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="Immune-DE T1DvsAAB", append=TRUE, row.names=T)
write.xlsx(DE_T1DvsAAB_stellates, file="DE_T1DvsAAB_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="Stellates-DE T1DvsAAB", append=TRUE, row.names=T)