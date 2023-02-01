
library(tidyverse)
library(MeltR)
library(viridis)
library(cowplot)

####CPEB3####

list.files("Tables/SI_Table_2_long_RNA_analysis")

df.CPEB3 = read.csv("Tables/SI_Table_2_long_RNA_analysis/js7008_CPEB3_formatted_data.csv")

fit.CPEB3 = meltR.A(df.CPEB3,
                    NucAcid = c("RNA", "GGGGGCCACAGCAGAAGCGUUCACGUCGCAGCCCCUGUCAGAUUCUGGUGAAUCUGCGAAUUCUGCUG"),
                    Mmodel = "Monomolecular.2State",
                    file_path = "Tables/SI_Table_2_long_RNA_analysis",
                    file_prefix = "CPEB3",
                    Save_results = "all")

BL.CPEB3 = BLTrimmer(fit.CPEB3,
                     file_path = "Tables/SI_Table_2_long_RNA_analysis",
                     file_prefix = "CPEB3",
                     Save_results = "all")

head(fit.CPEB3$Derivatives.data)


PA = ggplot(fit.CPEB3$Derivatives.data,
       aes(x = Temperature, y = dA.dT, color = log10(Ct))) +
  geom_point() +
  scale_color_viridis() +
  theme_classic() +
  xlab("Temperature (\u00B0C)") +
  ylab("dA/dT")

####Gua####

list.files("Tables/SI_Table_2_long_RNA_analysis")

df.Gua = read.csv("Tables/SI_Table_2_long_RNA_analysis/js7008_Gua_formatted_data.csv")

fit.Gua = meltR.A(df.Gua,
                    NucAcid = c("RNA", "GGGGGCCACAGCAGAAGCGUUCACGUCGCAGCCCCUGUCAGAUUCUGGUGAAUCUGCGAAUUCUGCUG"),
                    Mmodel = "Monomolecular.2State",
                    file_path = "Tables/SI_Table_2_long_RNA_analysis",
                    file_prefix = "Gua",
                    Save_results = "all")

BL.Gua = BLTrimmer(fit.Gua,
                     file_path = "Tables/SI_Table_2_long_RNA_analysis",
                     file_prefix = "Gua",
                     Save_results = "all")

head(fit.Gua$Derivatives.data)

PB = ggplot(fit.Gua$Derivatives.data,
       aes(x = Temperature, y = dA.dT, color = log10(Ct))) +
  geom_point() +
  scale_color_viridis() +
  theme_classic() +
  xlab("Temperature (\u00B0C)") +
  ylab("dA/dT")

####Final Plot####


P = plot_grid(PA, PB)

ggsave("Figures/SI_Figure_5_long_RNA_analysis/R_figure_1.svg", P, width = 5, height = 2, scale = 1.5, bg = "white")


####Format table####

table.CPEB3 = BL.CPEB3$Ensemble.energies
table.CPEB3$RNA = "CPEB3"
table.Gua = BL.Gua$Ensemble.energies
table.Gua$RNA = "Gua"

df.table = bind_rows(table.CPEB3, table.Gua)

df.table = df.table[-c(2,5),] %>% select(RNA, Method, dH, CI95.dH, dS, CI95.dS, dG, CI95.dG, Tm_at_0.1mM, CI95.Tm_at_0.1mM)

write.csv(df.table, "Tables/SI_Table_2_long_RNA_analysis/Table_S2.csv")
