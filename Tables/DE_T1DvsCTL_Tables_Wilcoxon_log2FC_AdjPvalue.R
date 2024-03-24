rm(list=ls())
library('dplyr')
setwd("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/")

DE_T1DvsCTL <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/DE_T1DvsCTL.rds")
DE_T1DvsCTL <- DE_T1DvsCTL[abs(DE_T1DvsCTL$avg_log2FC) >= 0.5 & DE_T1DvsCTL$p_val_adj<0.05,]
# write.xlsx(DE_T1DvsCTL, 'Allcells_T1DvsCTL_Wilcox.xlsx')

DE_T1DvsCTL_acinar <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/DE_T1DvsCTL_acinar.rds")
DE_T1DvsCTL_acinar <- DE_T1DvsCTL_acinar[abs(DE_T1DvsCTL_acinar$avg_log2FC) >= 0.5 & DE_T1DvsCTL_acinar$p_val_adj<0.05,]
# writexl::write_xlsx(DE_T1DvsCTL_acinar, 'Acinar_T1DvsCTL_Wilcox.xlsx')

DE_T1DvsCTL_alpha <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/DE_T1DvsCTL_alpha.rds")
DE_T1DvsCTL_alpha <- DE_T1DvsCTL_alpha[abs(DE_T1DvsCTL_alpha$avg_log2FC) >= 0.5 & DE_T1DvsCTL_alpha$p_val_adj<0.05,]
# writexl::write_xlsx(DE_T1DvsCTL_alpha, 'Alpha_T1DvsCTL_Wilcox.xlsx')

DE_T1DvsCTL_beta <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/DE_T1DvsCTL_beta.rds")
DE_T1DvsCTL_beta <- DE_T1DvsCTL_beta[abs(DE_T1DvsCTL_beta$avg_log2FC) >= 0.5 & DE_T1DvsCTL_beta$p_val_adj<0.05,]
# writexl::write_xlsx(DE_T1DvsCTL_beta, 'Beta_T1DvsCTL_Wilcox.xlsx')

DE_T1DvsCTL_delta <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/DE_T1DvsCTL_delta.rds")
# writexl::write_xlsx(DE_T1DvsCTL_delta, 'Delta_T1DvsCTL_Wilcox.xlsx')
DE_T1DvsCTL_delta <- DE_T1DvsCTL_delta[abs(DE_T1DvsCTL_delta$avg_log2FC) >= 0.5 & DE_T1DvsCTL_delta$p_val_adj<0.05,]

DE_T1DvsCTL_ductal <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/DE_T1DvsCTL_ductal.rds")
DE_T1DvsCTL_ductal <- DE_T1DvsCTL_ductal[abs(DE_T1DvsCTL_ductal$avg_log2FC) >= 0.5 & DE_T1DvsCTL_ductal$p_val_adj<0.05,]
# writexl::write_xlsx(DE_T1DvsCTL_ductal, 'Ductal_T1DvsCTL_Wilcox.xlsx')

DE_T1DvsCTL_immune <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/DE_T1DvsCTL_immune.rds")
DE_T1DvsCTL_immune <- DE_T1DvsCTL_immune[abs(DE_T1DvsCTL_immune$avg_log2FC) >= 0.5 & DE_T1DvsCTL_immune$p_val_adj<0.05,]
# writexl::write_xlsx(DE_T1DvsCTL_immune, 'Immune_T1DvsCTL_Wilcox.xlsx')

DE_T1DvsCTL_endothelial <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/DE_T1DvsCTL_endothelial.rds")
DE_T1DvsCTL_endothelial <- DE_T1DvsCTL_endothelial[abs(DE_T1DvsCTL_endothelial$avg_log2FC) >= 0.5 & DE_T1DvsCTL_endothelial$p_val_adj<0.05,]
# writexl::write_xlsx(DE_T1DvsCTL_endothelial, 'Endothelial_T1DvsCTL_Wilcox.xlsx')

DE_T1DvsCTL_stellates <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/DE_T1DvsCTL_stellates.rds")
DE_T1DvsCTL_stellates <- DE_T1DvsCTL_stellates[abs(DE_T1DvsCTL_stellates$avg_log2FC) >= 0.5 & DE_T1DvsCTL_stellates$p_val_adj<0.05,]
# writexl::write_xlsx(DE_T1DvsCTL_stellates, 'Stellates_T1DvsCTL_Wilcox.xlsx')

library('xlsx')
write.xlsx(DE_T1DvsCTL, file="DE_T1DvsCTL_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="All Cells-DE T1DvsCTL", row.names=T)
write.xlsx(DE_T1DvsCTL_acinar, file="DE_T1DvsCTL_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="Acinar-DE T1DvsCTL", append=TRUE, row.names=T)
write.xlsx(DE_T1DvsCTL_alpha, file="DE_T1DvsCTL_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="Alpha-DE T1DvsCTL", append=TRUE, row.names=T)
write.xlsx(DE_T1DvsCTL_beta, file="DE_T1DvsCTL_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="Beta-DE T1DvsCTL", append=TRUE, row.names=T)
write.xlsx(DE_T1DvsCTL_delta, file="DE_T1DvsCTL_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="Delta-DE T1DvsCTL", append=TRUE, row.names=T)
write.xlsx(DE_T1DvsCTL_ductal, file="DE_T1DvsCTL_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="Ductal-DE T1DvsCTL", append=TRUE, row.names=T)
write.xlsx(DE_T1DvsCTL_endothelial, file="DE_T1DvsCTL_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="Endothelial-DE T1DvsCTL", append=TRUE, row.names=T)
write.xlsx(DE_T1DvsCTL_immune, file="DE_T1DvsCTL_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="Immune-DE T1DvsCTL", append=TRUE, row.names=T)
write.xlsx(DE_T1DvsCTL_stellates, file="DE_T1DvsCTL_Wilcoxon_log2FC_AdjPvalue.xlsx", sheetName="Stellates-DE T1DvsCTL", append=TRUE, row.names=T)