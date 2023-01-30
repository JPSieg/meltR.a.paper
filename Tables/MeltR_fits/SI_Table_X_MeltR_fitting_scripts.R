library(MeltR)
library(tidyverse)
devtools::load_all()

list.df = {}

####CROWD DP1####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "CROWD DP1")

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(25, 75))

unique(df$Sample)

fit = meltR.A(df %>% filter(Sample != "11"),
              NucAcid = c("RNA", "CGCGCG"),
              wavelength = 280,
              concT = 75,
              Mmodel = "Homoduplex.2State",
              fitTs = c(20, 80),
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "CROWD DP1")

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "CROWD DP1")

list.df[[1]] = BL.fit$Ensemble.energies
list.df[[1]]$Helix =  "CGCGCG"

####CROWD DP5####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "CROWD DP5")

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(20, 65))

fit = meltR.A(df %>% filter(Sample != "24"),
              NucAcid = c("RNA", "ACCGGU"),
              concT = 85,
              Mmodel = "Homoduplex.2State",
              Save_results = "all",
              fitTs = c(20, 65),
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "CROWD DP5")

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "CROWD DP5")

list.df[[2]] = BL.fit$Ensemble.energies
list.df[[2]]$Helix =  "ACCGGU"

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
              concT = max(df$Temperature),
              fitTs = c(20,65),
              wavelength = 280,
              Mmodel = "Homoduplex.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "CROWD DP9")

 #?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "CROWD DP9")

list.df[[3]] = BL.fit$Ensemble.energies
list.df[[3]]$Helix =  "CCAUGG"

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
              concT = 55,
              Mmodel = "Homoduplex.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "CROWD DR1")

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "CROWD DR1")

list.df[[4]] = BL.fit$Ensemble.energies
list.df[[4]]$Helix =  "GAUAUAUC"

####CROWD_DR5####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "CROWD_DR5")

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(15, 65))


fit = meltR.A(df,
              NucAcid = c("RNA", "GCAAUUGC"),
              concT = 55,
              wavelength = 280,
              Mmodel = "Homoduplex.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "CROWD DR5")

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "CROWD DR5")

list.df[[5]] = BL.fit$Ensemble.energies
list.df[[5]]$Helix =  "GCAAUUGC"

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
              concT = 40,
              Mmodel = "Homoduplex.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "CROWD DR17")

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                    n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "CROWD DR17")

list.df[[6]] = BL.fit$Ensemble.energies
list.df[[6]]$Helix =  "UAUAUAUA"

####Helix A####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "Helix A")

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(5, 40))


fit = meltR.A(df,
              blank = 1,
              NucAcid = c("RNA", "CGAAAGGU", "ACCUUUCG"),
              concT = 85,
              Mmodel = "Heteroduplex.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "Helix A")

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "Helix A")

list.df[[7]] = BL.fit$Ensemble.energies
list.df[[7]]$Helix =  "CGAAAGGU/ACCUUUCG"

####Helix C####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "Helix C")

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(5, 40))

#View(df)

fit = meltR.A(df %>% filter(Sample != 3),
              blank = 2,
              NucAcid = c("RNA", "CUGAGUC", "GACUCAG"),
              concT = 85,
              Mmodel = "Heteroduplex.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "Helix C")

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "Helix C")

list.df[[8]] = BL.fit$Ensemble.energies
list.df[[8]]$Helix = "CUGAGUC/GACUCAG"

####Helix D####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "Helix D")

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(5, 40))

#View(df)

fit = meltR.A(df,
              blank = 3,
              NucAcid = c("RNA", "CGUUGC", "GCAACG"),
              concT = 85,
              Mmodel = "Heteroduplex.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "Helix D")

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "Helix D")

list.df[[9]] = BL.fit$Ensemble.energies
list.df[[9]]$Helix = "CGUUGC/GCAACG"

####Helix A labeled####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "Helix A labeled")

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(5, 40))

#View(df)

unique(df$Sample)
fit = meltR.A(df,
              NucAcid = c("RNA", "CGAAAGGU", "ACCUUUCG"),
              concT = 85,
              Mmodel = "Heteroduplex.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "Helix A labeled")

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "Helix A labeled")

list.df[[10]] = BL.fit$Ensemble.energies
list.df[[10]]$Helix = "FAMCGAAAGGU/ACCUUUCGBHQ1"

####Helix C labeled####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "Helix C labeled")

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(5, 40))

#View(df)

unique(df$Sample)
fit = meltR.A(df,
              NucAcid = c("RNA", "CUGAGUC", "GACUCAG"),
              concT = 85,
              Mmodel = "Heteroduplex.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "Helix C labeled")

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "Helix C labeled")

list.df[[11]] = BL.fit$Ensemble.energies
list.df[[11]]$Helix = "FAMCUGAGUC/GACUCAGBHQ1"

####Helix D labeled####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "Helix D labeled")

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(5, 40))

#View(df)

unique(df$Sample)
fit = meltR.A(df,
              NucAcid = c("RNA", "CGUUGC", "GCAACG"),
              concT = 85,
              Mmodel = "Heteroduplex.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "Helix C labeled")

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "Helix C labeled")

list.df[[12]] = BL.fit$Ensemble.energies
list.df[[12]]$Helix = "FAMCGUUGC/GCAACGBHQ1"

####Consolidate results####

df = bind_rows(list.df)

write.csv(df, "Tables/MeltR_fits/Fit_results.csv", row.names = FALSE)
