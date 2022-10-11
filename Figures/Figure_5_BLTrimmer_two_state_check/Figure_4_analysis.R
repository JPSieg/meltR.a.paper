library(tidyverse)
library(MeltR)

list.files("Raw_data")

?o3a.to.MeltR.csv

df = o3a.to.MeltR.csv("Raw_data", Remove_readings = c(182:184),
                      Pathlength = c(1, 1, 1, 0.1, 0.1, 0.1,
                                     0.5, 0.5, 0.5, 0.5,
                                     0.1, 0.1, 0.1))

ggplot(df, aes(x = Temperature, y = Absorbance, color = Sample)) +
  geom_point()

unique(df$Sample)

v.blanks = c("js6040-Cell 01.o3a", "js6040-Cell 04.o3a",
             "js6040-Cell 07.o3a", "js6040-Cell 11.o3a")

l.blanks = list(c("js6040-Cell 02.o3a", "js6040-Cell 01.o3a"),
                c("js6040-Cell 03.o3a", "js6040-Cell 01.o3a"),
                c("js6040-Cell 05.o3a", "js6040-Cell 04.o3a"),
                c("js6040-Cell 06.o3a", "js6040-Cell 04.o3a"),
                c("js6040-Cell 08.o3a", "js6040-Cell 07.o3a"),
                c("js6040-Cell 09.o3a", "js6040-Cell 07.o3a"),
                c("js6040-Cell 10.o3a", "js6040-Cell 07.o3a"),
                c("js6040-Cell 12.o3a", "js6040-Cell 11.o3a"),
                c("js6040-Cell 13.o3a", "js6040-Cell 11.o3a"))

list.df.bg = {}

for (i in 1:length(l.blanks)){
  df.A = df %>% filter(Sample == l.blanks[[i]][1])
  df.B = df %>% filter(Sample == l.blanks[[i]][2])
  df.A$Absorbance = df.A$Absorbance - df.B$Absorbance
  list.df.bg[[i]] = df.A
}

df = bind_rows(list.df.bg) %>% filter(Sample != "js6040-Cell 06.o3a")

ggplot(df, aes(x = Temperature, y = Absorbance, color = Sample)) +
  geom_point() +
  geom_vline(xintercept = c(5,75))

write.csv(df, "js6040_data.csv", row.names = FALSE)

?meltR.A

fit = meltR.A(df,
        NucAcid = c("Custom", 100300, 101100),
        Mmodel = "Heteroduplex.2State",
        fitTs = c(5, 75),
        Save_results = "all")

Trim = BLTrimmer(fit, n.combinations = 1000, Save_results = "all")

write.csv(Trim$Ensemble.energies, "js6040_result.csv", row.names = FALSE)
