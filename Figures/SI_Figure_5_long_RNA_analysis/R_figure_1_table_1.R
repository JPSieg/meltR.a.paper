
library(tidyverse)
library(MeltR)
library(viridis)
library(cowplot)

####CPEB3####

list.files("Revisions/R_figure_1_table_1")

df.CPEB3 = read.csv("Revisions/R_figure_1_table_1/js7008_CPEB3_formatted_data.csv")

fit.CPEB3 = meltR.A(df.CPEB3,
                    NucAcid = c("RNA", "GGGGGCCACAGCAGAAGCGUUCACGUCGCAGCCCCUGUCAGAUUCUGGUGAAUCUGCGAAUUCUGCUG"),
                    Mmodel = "Monomolecular.2State",
                    file_path = "Revisions/R_figure_1_table_1",
                    file_prefix = "CPEB3",
                    Save_results = "all")

BL.CPEB3 = BLTrimmer(fit.CPEB3,
                     file_path = "Revisions/R_figure_1_table_1",
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

list.files("Revisions/R_figure_1_table_1")

df.Gua = read.csv("Revisions/R_figure_1_table_1/js7008_Gua_formatted_data.csv")

fit.Gua = meltR.A(df.Gua,
                    NucAcid = c("RNA", "GGGGGCCACAGCAGAAGCGUUCACGUCGCAGCCCCUGUCAGAUUCUGGUGAAUCUGCGAAUUCUGCUG"),
                    Mmodel = "Monomolecular.2State",
                    file_path = "Revisions/R_figure_1_table_1",
                    file_prefix = "Gua",
                    Save_results = "all")

BL.Gua = BLTrimmer(fit.Gua,
                     file_path = "Revisions/R_figure_1_table_1",
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

ggsave("Revisions/R_figure_1_table_1/R_figure_1.svg", P, width = 5, height = 2, scale = 1.5, bg = "white")
