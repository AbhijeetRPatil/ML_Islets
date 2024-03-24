rm(list=ls())
library('dplyr')
options(future.globals.maxSize = 1000000 * 1024^2)
setwd("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_AABvsCTL/")

DE_AABvsCTL <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_AABvsCTL/DE_AABvsCTL.rds")
DE_AABvsCTL <- DE_AABvsCTL[DE_AABvsCTL$p_val_adj<0.05,]
# write.xlsx(DE_AABvsCTL, 'Allcells_AABvsCTL_Wilcox.xlsx')

DE_AABvsCTL_acinar <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_AABvsCTL/DE_AABvsCTL_acinar.rds")
DE_AABvsCTL_acinar <- DE_AABvsCTL_acinar[DE_AABvsCTL_acinar$p_val_adj<0.05,]
# writexl::write_xlsx(DE_AABvsCTL_acinar, 'Acinar_AABvsCTL_Wilcox.xlsx')

DE_AABvsCTL_alpha <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_AABvsCTL/DE_AABvsCTL_alpha.rds")
DE_AABvsCTL_alpha <- DE_AABvsCTL_alpha[DE_AABvsCTL_alpha$p_val_adj<0.05,]
# writexl::write_xlsx(DE_AABvsCTL_alpha, 'Alpha_AABvsCTL_Wilcox.xlsx')

DE_AABvsCTL_beta <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_AABvsCTL/DE_AABvsCTL_beta.rds")
DE_AABvsCTL_beta <- DE_AABvsCTL_beta[DE_AABvsCTL_beta$p_val_adj<0.05,]
# writexl::write_xlsx(DE_AABvsCTL_beta, 'Beta_AABvsCTL_Wilcox.xlsx')

DE_AABvsCTL_delta <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_AABvsCTL/DE_AABvsCTL_delta.rds")
# writexl::write_xlsx(DE_AABvsCTL_delta, 'Delta_AABvsCTL_Wilcox.xlsx')
DE_AABvsCTL_delta <- DE_AABvsCTL_delta[DE_AABvsCTL_delta$p_val_adj<0.05,]

DE_AABvsCTL_ductal <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_AABvsCTL/DE_AABvsCTL_ductal.rds")
DE_AABvsCTL_ductal <- DE_AABvsCTL_ductal[DE_AABvsCTL_ductal$p_val_adj<0.05,]
# writexl::write_xlsx(DE_AABvsCTL_ductal, 'Ductal_AABvsCTL_Wilcox.xlsx')

DE_AABvsCTL_immune <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_AABvsCTL/DE_AABvsCTL_immune.rds")
DE_AABvsCTL_immune <- DE_AABvsCTL_immune[DE_AABvsCTL_immune$p_val_adj<0.05,]
# writexl::write_xlsx(DE_AABvsCTL_immune, 'Immune_AABvsCTL_Wilcox.xlsx')

DE_AABvsCTL_endothelial <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_AABvsCTL/DE_AABvsCTL_endothelial.rds")
DE_AABvsCTL_endothelial <- DE_AABvsCTL_endothelial[DE_AABvsCTL_endothelial$p_val_adj<0.05,]
# writexl::write_xlsx(DE_AABvsCTL_endothelial, 'Endothelial_AABvsCTL_Wilcox.xlsx')

DE_AABvsCTL_stellates <- readRDS("C:/Users/abhij/Desktop/Spring 2022/ML paper/Revision_Nov2023/DE_AABvsCTL/DE_AABvsCTL_stellates.rds")
DE_AABvsCTL_stellates <- DE_AABvsCTL_stellates[DE_AABvsCTL_stellates$p_val_adj<0.05,]
# writexl::write_xlsx(DE_AABvsCTL_stellates, 'Stellates_AABvsCTL_Wilcox.xlsx')

library('xlsx')
write.xlsx(DE_AABvsCTL, file="DE_AABvsCTL_Wilcoxon_AdjPvalue.xlsx", sheetName="All Cells-DE AABvsCTL", row.names=T)
rm(DE_AABvsCTL)
gc()
write.xlsx(DE_AABvsCTL_acinar, file="DE_AABvsCTL_Wilcoxon_AdjPvalue.xlsx", sheetName="Acinar-DE AABvsCTL", append=TRUE, row.names=T)
rm(DE_AABvsCTL_acinar)
gc()

write.xlsx(DE_AABvsCTL_alpha, file="DE_AABvsCTL_Wilcoxon_AdjPvalue.xlsx", sheetName="Alpha-DE AABvsCTL", append=TRUE, row.names=T)
rm(DE_AABvsCTL_alpha)
gc()

write.xlsx(DE_AABvsCTL_beta, file="DE_AABvsCTL_Wilcoxon_AdjPvalue.xlsx", sheetName="Beta-DE AABvsCTL", append=TRUE, row.names=T)
rm(DE_AABvsCTL_beta)
gc()

write.xlsx(DE_AABvsCTL_delta, file="DE_AABvsCTL_Wilcoxon_AdjPvalue.xlsx", sheetName="Delta-DE AABvsCTL", append=TRUE, row.names=T)
rm(DE_AABvsCTL_delta)
gc()

write.xlsx(DE_AABvsCTL_ductal, file="DE_AABvsCTL_Wilcoxon_AdjPvalue.xlsx", sheetName="Ductal-DE AABvsCTL", append=TRUE, row.names=T)
rm(DE_AABvsCTL_ductal)
gc()

write.xlsx(DE_AABvsCTL_endothelial, file="DE_AABvsCTL_Wilcoxon_AdjPvalue.xlsx", sheetName="Endothelial-DE AABvsCTL", append=TRUE, row.names=T)
rm(DE_AABvsCTL_endothelial)
gc()

write.xlsx(DE_AABvsCTL_immune, file="DE_AABvsCTL_Wilcoxon_AdjPvalue.xlsx", sheetName="Immune-DE AABvsCTL", append=TRUE, row.names=T)
rm(DE_AABvsCTL_immune)
gc()

write.xlsx(DE_AABvsCTL_stellates, file="DE_AABvsCTL_Wilcoxon_AdjPvalue.xlsx", sheetName="Stellates-DE AABvsCTL", append=TRUE, row.names=T)
rm(DE_AABvsCTL_stellates)
gc()
