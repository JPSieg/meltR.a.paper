library(MeltR)
library(tidyverse)
devtools::load_all()

df = df.absorbance %>% filter(Experiment == "CROWD DP5")

meltR.A.fit = meltR.A(df,
                      NucAcid = c("RNA", "ACCGGU"),
                      Mmodel = "Homoduplex.2State",
                      concT = 80,
                      fitTs = c(20, 65))

?BLTrimmer

list.files("Figures/Figure_3_Baseline_trimmer_internal")

Trimmed = BLTrimmer(meltR.A.fit,
                    n.ranges.float = 5,
                    range.step.fixed = 5,
                    n.combinations = 100,
                    Save_results = "all",
                    file_path = "Figures/Figure_3_Baseline_trimmer_internal")
