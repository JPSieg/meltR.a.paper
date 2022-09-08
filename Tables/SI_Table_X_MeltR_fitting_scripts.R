library(MeltR)
library(tidyverse)
devtools::load_all()

list.fits = {}

####CROWD DP1####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "CROWD DP1")

?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(25, 75))

fit = meltR.A(df,
              NucAcid = c("RNA", "CGCGCG"),
              concT = 20,
              Mmodel = "Homoduplex.2State",
              fitTs = c(20, 80),
              Save_results = "all",
              file_path = "Tables/SI_Table_X_MeltR_fits/Fit_data",
              file_prefix = "CROWD DP1")

?BLTrimmer

BL.fit = BLTrimmer(fit,
                   quantile.threshold = 1,
                   n.combinations = 100,
                   Save_results = "all",
                   file_path = "Tables/SI_Table_X_MeltR_fits/Fit_data",
                   file_prefix = "CROWD DP1")

####CROWD DP5####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "CROWD DP5")

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(20, 65))

list(c(25, 75),
     c(20, 70),
     c(10, 55),
     c(20, 70),
     c(20, 65),
     c(20, 65),
     c(20, 65),
     c(20, 65),
     c(20, 65))

fit = meltR.A(df,
              NucAcid = c("RNA", "ACCGGU"),
              concT = 20,
              Mmodel = "Homoduplex.2State",
              Save_results = "all",
              fitTs = c(20, 65),
              file_path = "Tables/SI_Table_X_MeltR_fits/Fit_data",
              file_prefix = "CROWD DP5")

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   quantile.threshold = 1,
                   n.combinations = 100,
                   Save_results = "all",
                   file_path = "Tables/SI_Table_X_MeltR_fits/Fit_data",
                   file_prefix = "CROWD DP5")

####CROWD_DP9####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "CROWD_DP9")

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(20, 65))


fit = meltR.A(df,
              NucAcid = c("RNA", "CCAUGG"),
              concT = 20,
              Mmodel = "Homoduplex.2State",
              Save_results = "all",
              file_path = "Tables/SI_Table_X_MeltR_fits/Fit_data",
              file_prefix = "CROWD DP9")

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   quantile.threshold = 1,
                   n.combinations = 100,
                   Save_results = "all",
                   file_path = "Tables/SI_Table_X_MeltR_fits/Fit_data",
                   file_prefix = "CROWD DP9")

####CROWD_DR1####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "CROWD_DR1")

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(15, 65))


fit = meltR.A(df,
              NucAcid = c("RNA", "GAUAUAUC"),
              concT = 11,
              Mmodel = "Homoduplex.2State",
              Save_results = "all",
              file_path = "Tables/SI_Table_X_MeltR_fits/Fit_data",
              file_prefix = "CROWD DR1")

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   quantile.threshold = 1,
                   n.combinations = 100,
                   Save_results = "all",
                   file_path = "Tables/SI_Table_X_MeltR_fits/Fit_data",
                   file_prefix = "CROWD DR1")

####CROWD_DR17####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "CROWD_DR17")

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(5, 40))


fit = meltR.A(df,
              NucAcid = c("RNA", "UAUAUAUA"),
              concT = 11,
              Mmodel = "Homoduplex.2State",
              Save_results = "all",
              file_path = "Tables/SI_Table_X_MeltR_fits/Fit_data",
              file_prefix = "CROWD DR17")

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   range.step.float = 1,
                   quantile.threshold = 1,
                   n.combinations = 100,
                   Save_results = "all",
                   file_path = "Tables/SI_Table_X_MeltR_fits/Fit_data",
                   file_prefix = "CROWD DR17")
