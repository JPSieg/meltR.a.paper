devtools::load_all()
library(MeltR)
library(tidyverse)

unique(df.absorbance$RNA)
list.files("Revisions")

####Low melting####

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

####High melting####

?meltR.A

fit = meltR.A(df.absorbance %>% filter(RNA == "AGCCGGCU"),
              NucAcid = c("RNA", "AGCCGGCU"),
              Mmodel = "Homoduplex.2State",
              Save_results = "all",
              file_path = "Revisions/R_figure_1",
              file_prefix = "High_melting_helix")

BL.fit = BLTrimmer(fit,
                   Save_results = "all",
                   file_path = "Revisions/R_figure_1",
                   file_prefix = "High_melting_helix")
