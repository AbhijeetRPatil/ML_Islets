rm(list=ls())
library('dplyr')
options(future.globals.maxSize = 1000000 * 1024^2)
setwd("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/")

DE_T1DvsCTL <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/PB_DE_T1DvsCTL_Allcells.rds")
DE_T1DvsCTL <- DE_T1DvsCTL[DE_T1DvsCTL$p_val_adj<0.05,]
DE_T1DvsCTL <- DE_T1DvsCTL[complete.cases(DE_T1DvsCTL),]

# write.xlsx(DE_T1DvsCTL, 'Allcells_T1DvsCTL_Wilcox.xlsx')

DE_T1DvsCTL_acinar <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/PB_DE_T1DvsCTL_Acinar.rds")
DE_T1DvsCTL_acinar <- DE_T1DvsCTL_acinar[DE_T1DvsCTL_acinar$p_val_adj<0.05,]
DE_T1DvsCTL_acinar <- DE_T1DvsCTL_acinar[complete.cases(DE_T1DvsCTL_acinar),]
# writexl::write_xlsx(DE_T1DvsCTL_acinar, 'Acinar_T1DvsCTL_Wilcox.xlsx')

DE_T1DvsCTL_alpha <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/PB_DE_T1DvsCTL_Alpha.rds")
DE_T1DvsCTL_alpha <- DE_T1DvsCTL_alpha[DE_T1DvsCTL_alpha$p_val_adj<0.05,]
DE_T1DvsCTL_alpha <- DE_T1DvsCTL_alpha[complete.cases(DE_T1DvsCTL_alpha),]
# writexl::write_xlsx(DE_T1DvsCTL_alpha, 'Alpha_T1DvsCTL_Wilcox.xlsx')

DE_T1DvsCTL_beta <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/PB_DE_T1DvsCTL_Beta.rds")
DE_T1DvsCTL_beta <- DE_T1DvsCTL_beta[DE_T1DvsCTL_beta$p_val_adj<0.05,]
DE_T1DvsCTL_beta <- DE_T1DvsCTL_beta[complete.cases(DE_T1DvsCTL_beta),]
# writexl::write_xlsx(DE_T1DvsCTL_beta, 'Beta_T1DvsCTL_Wilcox.xlsx')

DE_T1DvsCTL_delta <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/PB_DE_T1DvsCTL_Delta.rds")
DE_T1DvsCTL_delta <- DE_T1DvsCTL_delta[DE_T1DvsCTL_delta$p_val_adj<0.05,]
DE_T1DvsCTL_delta <- DE_T1DvsCTL_delta[complete.cases(DE_T1DvsCTL_delta),]
# writexl::write_xlsx(DE_T1DvsCTL_delta, 'Delta_T1DvsCTL_Wilcox.xlsx')

DE_T1DvsCTL_ductal <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/PB_DE_T1DvsCTL_Ductal.rds")
DE_T1DvsCTL_ductal <- DE_T1DvsCTL_ductal[DE_T1DvsCTL_ductal$p_val_adj<0.05,]
DE_T1DvsCTL_ductal <- DE_T1DvsCTL_ductal[complete.cases(DE_T1DvsCTL_ductal),]
# writexl::write_xlsx(DE_T1DvsCTL_ductal, 'Ductal_T1DvsCTL_Wilcox.xlsx')

DE_T1DvsCTL_immune <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/PB_DE_T1DvsCTL_Immune.rds")
DE_T1DvsCTL_immune <- DE_T1DvsCTL_immune[DE_T1DvsCTL_immune$p_val_adj<0.05,]
DE_T1DvsCTL_immune <- DE_T1DvsCTL_immune[complete.cases(DE_T1DvsCTL_immune),]
# writexl::write_xlsx(DE_T1DvsCTL_immune, 'Immune_T1DvsCTL_Wilcox.xlsx')

DE_T1DvsCTL_endothelial <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/PB_DE_T1DvsCTL_Endothelial.rds")
DE_T1DvsCTL_endothelial <- DE_T1DvsCTL_endothelial[DE_T1DvsCTL_endothelial$p_val_adj<0.05,]
DE_T1DvsCTL_endothelial <- DE_T1DvsCTL_endothelial[complete.cases(DE_T1DvsCTL_endothelial),]
# writexl::write_xlsx(DE_T1DvsCTL_endothelial, 'Endothelial_T1DvsCTL_Wilcox.xlsx')

DE_T1DvsCTL_stellates <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_T1DvsCTL/PB_DE_T1DvsCTL_Stellates.rds")
DE_T1DvsCTL_stellates <- DE_T1DvsCTL_stellates[DE_T1DvsCTL_stellates$p_val_adj<0.05,]
DE_T1DvsCTL_stellates <- DE_T1DvsCTL_stellates[complete.cases(DE_T1DvsCTL_stellates),]
# writexl::write_xlsx(DE_T1DvsCTL_stellates, 'Stellates_T1DvsCTL_Wilcox.xlsx')

library('xlsx')
write.xlsx(DE_T1DvsCTL, file="DE_T1DvsCTL_DESeq2_AdjPvalue.xlsx", sheetName="All Cells-DE T1DvsCTL", row.names=T)
# rm(DE_T1DvsCTL)
# gc()
write.xlsx(DE_T1DvsCTL_acinar, file="DE_T1DvsCTL_DESeq2_AdjPvalue.xlsx", sheetName="Acinar-DE T1DvsCTL", append=TRUE, row.names=T)
# rm(DE_T1DvsCTL_acinar)
# gc()

write.xlsx(DE_T1DvsCTL_alpha, file="DE_T1DvsCTL_DESeq2_AdjPvalue.xlsx", sheetName="Alpha-DE T1DvsCTL", append=TRUE, row.names=T)
# rm(DE_T1DvsCTL_alpha)
# gc()

write.xlsx(DE_T1DvsCTL_beta, file="DE_T1DvsCTL_DESeq2_AdjPvalue.xlsx", sheetName="Beta-DE T1DvsCTL", append=TRUE, row.names=T)
# rm(DE_T1DvsCTL_beta)
# gc()

write.xlsx(DE_T1DvsCTL_delta, file="DE_T1DvsCTL_DESeq2_AdjPvalue.xlsx", sheetName="Delta-DE T1DvsCTL", append=TRUE, row.names=T)
# rm(DE_T1DvsCTL_delta)
# gc()

write.xlsx(DE_T1DvsCTL_ductal, file="DE_T1DvsCTL_DESeq2_AdjPvalue.xlsx", sheetName="Ductal-DE T1DvsCTL", append=TRUE, row.names=T)
# rm(DE_T1DvsCTL_ductal)
# gc()

write.xlsx(DE_T1DvsCTL_endothelial, file="DE_T1DvsCTL_DESeq2_AdjPvalue.xlsx", sheetName="Endothelial-DE T1DvsCTL", append=TRUE, row.names=T)
# rm(DE_T1DvsCTL_endothelial)
# gc()

write.xlsx(DE_T1DvsCTL_immune, file="DE_T1DvsCTL_DESeq2_AdjPvalue.xlsx", sheetName="Immune-DE T1DvsCTL", append=TRUE, row.names=T)
# rm(DE_T1DvsCTL_immune)
# gc()

write.xlsx(DE_T1DvsCTL_stellates, file="DE_T1DvsCTL_DESeq2_AdjPvalue.xlsx", sheetName="Stellates-DE T1DvsCTL", append=TRUE, row.names=T)
# rm(DE_T1DvsCTL_stellates)
# gc()
