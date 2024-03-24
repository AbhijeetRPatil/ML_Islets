rm(list=ls())
library('dplyr')
options(future.globals.maxSize = 1000000 * 1024^2)
setwd("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/")

DE_T1DvsAAB <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/PB_DE_T1DvsAAB_Allcells.rds")
DE_T1DvsAAB <- DE_T1DvsAAB[DE_T1DvsAAB$p_val_adj<0.05,]
DE_T1DvsAAB <- DE_T1DvsAAB[complete.cases(DE_T1DvsAAB),]

# write.xlsx(DE_T1DvsAAB, 'Allcells_T1DvsAAB_Wilcox.xlsx')

DE_T1DvsAAB_acinar <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/PB_DE_T1DvsAAB_Acinar.rds")
DE_T1DvsAAB_acinar <- DE_T1DvsAAB_acinar[DE_T1DvsAAB_acinar$p_val_adj<0.05,]
DE_T1DvsAAB_acinar <- DE_T1DvsAAB_acinar[complete.cases(DE_T1DvsAAB_acinar),]
# writexl::write_xlsx(DE_T1DvsAAB_acinar, 'Acinar_T1DvsAAB_Wilcox.xlsx')

DE_T1DvsAAB_alpha <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/PB_DE_T1DvsAAB_Alpha.rds")
DE_T1DvsAAB_alpha <- DE_T1DvsAAB_alpha[DE_T1DvsAAB_alpha$p_val_adj<0.05,]
DE_T1DvsAAB_alpha <- DE_T1DvsAAB_alpha[complete.cases(DE_T1DvsAAB_alpha),]
# writexl::write_xlsx(DE_T1DvsAAB_alpha, 'Alpha_T1DvsAAB_Wilcox.xlsx')

DE_T1DvsAAB_beta <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/PB_DE_T1DvsAAB_Beta.rds")
DE_T1DvsAAB_beta <- DE_T1DvsAAB_beta[DE_T1DvsAAB_beta$p_val_adj<0.05,]
DE_T1DvsAAB_beta <- DE_T1DvsAAB_beta[complete.cases(DE_T1DvsAAB_beta),]
# writexl::write_xlsx(DE_T1DvsAAB_beta, 'Beta_T1DvsAAB_Wilcox.xlsx')

DE_T1DvsAAB_delta <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/PB_DE_T1DvsAAB_Delta.rds")
DE_T1DvsAAB_delta <- DE_T1DvsAAB_delta[DE_T1DvsAAB_delta$p_val_adj<0.05,]
DE_T1DvsAAB_delta <- DE_T1DvsAAB_delta[complete.cases(DE_T1DvsAAB_delta),]
# writexl::write_xlsx(DE_T1DvsAAB_delta, 'Delta_T1DvsAAB_Wilcox.xlsx')

DE_T1DvsAAB_ductal <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/PB_DE_T1DvsAAB_Ductal.rds")
DE_T1DvsAAB_ductal <- DE_T1DvsAAB_ductal[DE_T1DvsAAB_ductal$p_val_adj<0.05,]
DE_T1DvsAAB_ductal <- DE_T1DvsAAB_ductal[complete.cases(DE_T1DvsAAB_ductal),]
# writexl::write_xlsx(DE_T1DvsAAB_ductal, 'Ductal_T1DvsAAB_Wilcox.xlsx')

DE_T1DvsAAB_immune <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/PB_DE_T1DvsAAB_Immune.rds")
DE_T1DvsAAB_immune <- DE_T1DvsAAB_immune[DE_T1DvsAAB_immune$p_val_adj<0.05,]
DE_T1DvsAAB_immune <- DE_T1DvsAAB_immune[complete.cases(DE_T1DvsAAB_immune),]
# writexl::write_xlsx(DE_T1DvsAAB_immune, 'Immune_T1DvsAAB_Wilcox.xlsx')

DE_T1DvsAAB_endothelial <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/PB_DE_T1DvsAAB_Endothelial.rds")
DE_T1DvsAAB_endothelial <- DE_T1DvsAAB_endothelial[DE_T1DvsAAB_endothelial$p_val_adj<0.05,]
DE_T1DvsAAB_endothelial <- DE_T1DvsAAB_endothelial[complete.cases(DE_T1DvsAAB_endothelial),]
# writexl::write_xlsx(DE_T1DvsAAB_endothelial, 'Endothelial_T1DvsAAB_Wilcox.xlsx')

DE_T1DvsAAB_stellates <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsAAB/PB_DE_T1DvsAAB_Stellates.rds")
DE_T1DvsAAB_stellates <- DE_T1DvsAAB_stellates[DE_T1DvsAAB_stellates$p_val_adj<0.05,]
DE_T1DvsAAB_stellates <- DE_T1DvsAAB_stellates[complete.cases(DE_T1DvsAAB_stellates),]
# writexl::write_xlsx(DE_T1DvsAAB_stellates, 'Stellates_T1DvsAAB_Wilcox.xlsx')

library('xlsx')
write.xlsx(DE_T1DvsAAB, file="DE_T1DvsAAB_DESeq2_AdjPvalue.xlsx", sheetName="All Cells-DE T1DvsAAB", row.names=T)
# rm(DE_T1DvsAAB)
# gc()
write.xlsx(DE_T1DvsAAB_acinar, file="DE_T1DvsAAB_DESeq2_AdjPvalue.xlsx", sheetName="Acinar-DE T1DvsAAB", append=TRUE, row.names=T)
# rm(DE_T1DvsAAB_acinar)
# gc()

write.xlsx(DE_T1DvsAAB_alpha, file="DE_T1DvsAAB_DESeq2_AdjPvalue.xlsx", sheetName="Alpha-DE T1DvsAAB", append=TRUE, row.names=T)
# rm(DE_T1DvsAAB_alpha)
# gc()

write.xlsx(DE_T1DvsAAB_beta, file="DE_T1DvsAAB_DESeq2_AdjPvalue.xlsx", sheetName="Beta-DE T1DvsAAB", append=TRUE, row.names=T)
# rm(DE_T1DvsAAB_beta)
# gc()

write.xlsx(DE_T1DvsAAB_delta, file="DE_T1DvsAAB_DESeq2_AdjPvalue.xlsx", sheetName="Delta-DE T1DvsAAB", append=TRUE, row.names=T)
# rm(DE_T1DvsAAB_delta)
# gc()

write.xlsx(DE_T1DvsAAB_ductal, file="DE_T1DvsAAB_DESeq2_AdjPvalue.xlsx", sheetName="Ductal-DE T1DvsAAB", append=TRUE, row.names=T)
# rm(DE_T1DvsAAB_ductal)
# gc()

write.xlsx(DE_T1DvsAAB_endothelial, file="DE_T1DvsAAB_DESeq2_AdjPvalue.xlsx", sheetName="Endothelial-DE T1DvsAAB", append=TRUE, row.names=T)
# rm(DE_T1DvsAAB_endothelial)
# gc()

write.xlsx(DE_T1DvsAAB_immune, file="DE_T1DvsAAB_DESeq2_AdjPvalue.xlsx", sheetName="Immune-DE T1DvsAAB", append=TRUE, row.names=T)
# rm(DE_T1DvsAAB_immune)
# gc()

write.xlsx(DE_T1DvsAAB_stellates, file="DE_T1DvsAAB_DESeq2_AdjPvalue.xlsx", sheetName="Stellates-DE T1DvsAAB", append=TRUE, row.names=T)
# rm(DE_T1DvsAAB_stellates)
# gc()
