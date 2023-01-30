devtools::load_all()
library(MeltR)
library(tidyverse)

unique(df.absorbance$RNA)
list.files("Revisions")

?meltR.A

fit = meltR.A(df.absorbance %>% filter(RNA == "UAUAUAUA"),
        NucAcid = c("RNA", "UAUAUAUA"),
        Mmodel = "Homoduplex.2State",
        Save_results = "all",
        file_path = "Revisions/R_figure_1",
        file_prefix = "Low_melting_helix")

BL.fit = BLTrimmer(fit,
          Save_results = "all",
          file_path = "Revisions/R_figure_1",
          file_prefix = "Low_melting_helix")

df.1 = df.Meltwin %>% filter(Sequence == "UAUAUAUA") %>% select(!Sequence)
df.2 = BL.fit$Ensemble.energies

colnames(df.1)[c(3, 5, 7, 8)] = c("Error dH", "Error dS", "Error dG", "Tm_at_0.1mM")
colnames(df.2)[c(3, 5, 7)] = c("Error dH", "Error dS", "Error dG")
df.1$Program = "Meltwin"
df.2$Program = "MeltR"
df.2 = df.2[,-9]
df.1$`Error dH` = as.character(df.1$`Error dH`)
df.1$`Error dG` = as.character(df.1$`Error dG`)
df.1$`Error dS` = as.character(df.1$`Error dS`)
df = bind_rows(df.1, df.2) %>% select("Program", "Method", "dH", "Error dH", "dS", "Error dS", "dG", "Error dG", "Tm_at_0.1mM")

write.csv(df, "Revisions/R_figure_1/R_table_1.csv", row.names = F)
