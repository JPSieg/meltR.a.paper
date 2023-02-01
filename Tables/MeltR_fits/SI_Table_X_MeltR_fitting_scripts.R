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

####CROWD_DR8####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "CROWD_DR8")

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(40, 90))


fit = meltR.A(df,
              NucAcid = c("RNA", "AGCCGGCU"),
              concT = 40,
              wavelength = 280,
              Mmodel = "Homoduplex.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "CROWD_DR8")

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "CROWD_DR8")

list.df[[7]] = BL.fit$Ensemble.energies
list.df[[7]]$Helix =  "AGCCGGCU"

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

list.df[[8]] = BL.fit$Ensemble.energies
list.df[[8]]$Helix =  "CGAAAGGU/ACCUUUCG"

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

list.df[[9]] = BL.fit$Ensemble.energies
list.df[[9]]$Helix = "CUGAGUC/GACUCAG"

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

list.df[[10]] = BL.fit$Ensemble.energies
list.df[[10]]$Helix = "CGUUGC/GCAACG"

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

list.df[[11]] = BL.fit$Ensemble.energies
list.df[[11]]$Helix = "FAMCGAAAGGU/ACCUUUCGBHQ1"

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

list.df[[12]] = BL.fit$Ensemble.energies
list.df[[12]]$Helix = "FAMCUGAGUC/GACUCAGBHQ1"

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

list.df[[13]] = BL.fit$Ensemble.energies
list.df[[13]]$Helix = "FAMCGUUGC/GCAACGBHQ1"

####HP4AP2####

unique(df.absorbance$Experiment)

df = df.absorbance %>%
  filter(Experiment == "HP4AP2")# %>%
  #filter(Sample %in% c(91,92,93,97,98,99)) #Remove samples with no lower baseline

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(40, 95))

#View(df)

unique(df$Sample)

fit = meltR.A(df,
              NucAcid = c("RNA", df$RNA[1]),
              concT = 85,
              wavelength = 280,
              Mmodel = "Monomolecular.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "HP4AP2",
              fitTs = c(30, 90))

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "HP4AP2")

list.df[[14]] = BL.fit$Ensemble.energies
list.df[[14]]$Helix = df$RNA[1]

####HP4AP3####

unique(df.absorbance$Experiment)

df = df.absorbance %>%
  filter(Experiment == "HP4AP3") %>%
  filter(Sample %in% c(84, 86, 87,90)) #Remove samples with no lower baseline

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(40, 95))

#View(df)

unique(df$Sample)

fit = meltR.A(df,
              NucAcid = c("RNA", df$RNA[1]),
              concT = 85,
              wavelength = 280,
              Mmodel = "Monomolecular.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "HP4AP3",
              fitTs = c(40, 90))

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "HP4AP3")

list.df[[15]] = BL.fit$Ensemble.energies
list.df[[15]]$Helix = df$RNA[1]

#Note I'm going to remove this helix from the dataset because I don't think 4, not great, melting curves are enough to make conclusions 

####HP4AP8####

unique(df.absorbance$Experiment)

df = df.absorbance %>%
  filter(Experiment == "HP4AP8") %>%
  filter(Sample %in% c(91,92,93,97,98,99)) #Remove samples with no lower baseline

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(30, 90))

#View(df)

unique(df$Sample)

fit = meltR.A(df,
              NucAcid = c("RNA", df$RNA[1]),
              concT = 85,
              Mmodel = "Monomolecular.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "HP4AP8",
              fitTs = c(30, 90))

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "HP4AP8")

list.df[[16]] = BL.fit$Ensemble.energies
list.df[[16]]$Helix = df$RNA[1]

####HP4AP14####

unique(df.absorbance$Experiment)

df = df.absorbance %>%
  filter(Experiment == "HP4AP14") %>%
  filter(Sample %in% c(100, 101, 102, 103, 104, 105)) #Remove samples with no lower baseline

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(30, 85))

#View(df)

unique(df$Sample)

fit = meltR.A(df,
              NucAcid = c("RNA", df$RNA[1]),
              concT = 85,
              wavelength = 280,
              Mmodel = "Monomolecular.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "HP4AP14",
              fitTs = c(30, 85))

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "HP4AP14")

list.df[[17]] = BL.fit$Ensemble.energies
list.df[[17]]$Helix = df$RNA[1]

####HP4AP18####

unique(df.absorbance$Experiment)

df = df.absorbance %>%
  filter(Experiment == "HP4AP18") #%>%
  #filter(Sample %in% c(100, 101, 102, 103, 104, 105)) #Remove samples with no lower baseline

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(30, 75))

#View(df)

unique(df$Sample)

fit = meltR.A(df,
              NucAcid = c("RNA", df$RNA[1]),
              concT = 70,
              wavelength = 280,
              Mmodel = "Monomolecular.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "HP4AP18",
              fitTs = c(30, 75))

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "HP4AP18")

list.df[[18]] = BL.fit$Ensemble.energies
list.df[[18]]$Helix = df$RNA[1]


####HP4AP25####

unique(df.absorbance$Experiment)

df = df.absorbance %>%
  filter(Experiment == "HP4AP25") #%>%
#filter(Sample %in% c(100, 101, 102, 103, 104, 105)) #Remove samples with no lower baseline

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(40, 90))

#View(df)

unique(df$Sample)

fit = meltR.A(df,
              NucAcid = c("RNA", df$RNA[1]),
              concT = 80,
              wavelength = 280,
              Mmodel = "Monomolecular.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "HP4AP25",
              fitTs = c(40, 90))

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "HP4AP25")

list.df[[19]] = BL.fit$Ensemble.energies
list.df[[19]]$Helix = df$RNA[1]

####HP3AQ1####

unique(df.absorbance$Experiment)

df = df.absorbance %>%
  filter(Experiment == "HP3AQ1") %>%
  filter(Sample %in% c(129, 132, 136, 137, 138, 140, 141)) #Remove samples with no lower baseline

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(40, 90))

#View(df)

unique(df$Sample)

fit = meltR.A(df,
              NucAcid = c("RNA", df$RNA[1]),
              concT = 85,
              wavelength = 280,
              Mmodel = "Monomolecular.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "HP3AQ1",
              fitTs = c(40, 90))

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "HP3AQ1")

list.df[[20]] = BL.fit$Ensemble.energies
list.df[[20]]$Helix = df$RNA[1]

####HP3BH1####

unique(df.absorbance$Experiment)

df = df.absorbance %>%
  filter(Experiment == "HP3BH1") %>%
  filter(Sample %in% c(144, 145, 147, 152, 153)) #Remove samples with no lower baseline

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(40, 90))

#View(df)

unique(df$Sample)

fit = meltR.A(df,
              NucAcid = c("RNA", df$RNA[1]),
              concT = 85,
              wavelength = 280,
              Mmodel = "Monomolecular.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "HP3BH1",
              fitTs = c(40, 90))

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "HP3BH1")

unique(df.absorbance$Experiment)

list.df[[21]] = BL.fit$Ensemble.energies
list.df[[21]]$Helix = df$RNA[1]

####HP3AQ2####

unique(df.absorbance$Experiment)

df = df.absorbance %>%
  filter(Experiment == "HP3AQ2") %>%
  filter(Sample %in% c(154, 155, 156, 158, 159, 161, 162)) #Remove samples with no lower baseline

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(40, 90))

#View(df)

unique(df$Sample)

fit = meltR.A(df,
              NucAcid = c("RNA", df$RNA[1]),
              concT = 90,
              wavelength = 280,
              Mmodel = "Monomolecular.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "HP3AQ2",
              fitTs = c(40, 90))

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "HP3AQ2")

list.df[[22]] = BL.fit$Ensemble.energies
list.df[[22]]$Helix = df$RNA[1]


####HP3S5####

unique(df.absorbance$Experiment)

df = df.absorbance %>%
  filter(Experiment == "HP3S5") %>%
  filter(Sample %in% c(163, 164, 165, 167, 168, 170, 171)) #Remove samples with no lower baseline

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(40, 90))

#View(df)

unique(df$Sample)

fit = meltR.A(df,
              NucAcid = c("RNA", df$RNA[1]),
              concT = 90,
              wavelength = 280,
              Mmodel = "Monomolecular.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "HP3S5",
              fitTs = c(40, 90))

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "HHP3S5")

list.df[[23]] = BL.fit$Ensemble.energies
list.df[[23]]$Helix = df$RNA[1]

####HP3AC2####

unique(df.absorbance$Experiment)

df = df.absorbance %>%
  filter(Experiment == "HP3AC2") %>%
  filter(Sample %in% c(172, 174, 175, 176, 177, 179, 180)) #Remove samples with no lower baseline

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(35, 90))

#View(df)

unique(df$Sample)

fit = meltR.A(df,
              NucAcid = c("RNA", df$RNA[1]),
              concT = 90,
              wavelength = 280,
              Mmodel = "Monomolecular.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "HP3AC2",
              fitTs = c(35, 90))

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "HP3AC2")

list.df[[24]] = BL.fit$Ensemble.energies
list.df[[24]]$Helix = df$RNA[1]

####HP3AC4####

unique(df.absorbance$Experiment)

df = df.absorbance %>%
  filter(Experiment == "HP3AC4") %>%
  filter(Sample %in% c(182, 184, 185, 177, 187, 188)) #Remove samples with no lower baseline

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(50, 95))

#View(df)

unique(df$Sample)

fit = meltR.A(df,
              NucAcid = c("RNA", df$RNA[1]),
              concT = 90,
              wavelength = 280,
              Mmodel = "Monomolecular.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "HP3AC4",
              fitTs = c(50, 95))

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "HP3AC4")


list.df[[25]] = BL.fit$Ensemble.energies
list.df[[25]]$Helix = df$RNA[1]

####HP5_4_2####

unique(df.absorbance$Experiment)

df = df.absorbance %>%
  filter(Experiment == "HP5_4_2") %>%
  filter(Sample %in% c(208, 212, 214, 216, 218, 220)) #Remove samples with no lower baseline

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(40, 90))

#View(df)

unique(df$Sample)

fit = meltR.A(df,
              NucAcid = c("RNA", df$RNA[1]),
              concT = 90,
              wavelength = 260,
              Mmodel = "Monomolecular.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "HP5_4_2",
              fitTs = c(50, 95))

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "HP5_4_2")


list.df[[26]] = BL.fit$Ensemble.energies
list.df[[26]]$Helix = df$RNA[1]

####HP5_4_4####

unique(df.absorbance$Experiment)

df = df.absorbance %>%
  filter(Experiment == "HP5_4_4") %>%
  filter(Sample %in% c(223, 224, 2230, 231, 235, 236)) #Remove samples with no lower baseline

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(30, 90))

#View(df)

unique(df$Sample)

fit = meltR.A(df,
              NucAcid = c("RNA", df$RNA[1]),
              concT = 90,
              wavelength = 260,
              Mmodel = "Monomolecular.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "HP5_4_4",
              fitTs = c(30, 90))

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "HP5_4_4")


list.df[[27]] = BL.fit$Ensemble.energies
list.df[[27]]$Helix = df$RNA[1]

####HP5_4_5####

unique(df.absorbance$Experiment)

df = df.absorbance %>%
  filter(Experiment == "HP5_4_5") %>%
  filter(Sample %in% c(238, 239, 247, 250, 251, 252)) #Remove samples with no lower baseline

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(30, 90))

#View(df)

unique(df$Sample)

fit = meltR.A(df,
              NucAcid = c("RNA", df$RNA[1]),
              concT = 90,
              wavelength = 260,
              Mmodel = "Monomolecular.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "HP5_4_5",
              fitTs = c(30, 90))

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "HP5_4_5")


list.df[[28]] = BL.fit$Ensemble.energies
list.df[[28]]$Helix = df$RNA[1]

####HP5_4_9####

unique(df.absorbance$Experiment)

df = df.absorbance %>%
  filter(Experiment == "HP5_4_9") %>%
  filter(Sample %in% c(254, 258, 259, 260, 262, 266, 267, 268)) #Remove samples with no lower baseline

#?meltR.A

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point() +
  geom_vline(xintercept = c(30, 90))

#View(df)

unique(df$Sample)

fit = meltR.A(df,
              NucAcid = c("RNA", df$RNA[1]),
              concT = 90,
              wavelength = 260,
              Mmodel = "Monomolecular.2State",
              Save_results = "all",
              file_path = "Tables/MeltR_fits/Fit_data",
              file_prefix = "HP5_4_9",
              fitTs = c(30, 90))

#?BLTrimmer

BL.fit = BLTrimmer(fit,
                   n.combinations = 1000,
                   Save_results = "all",
                   file_path = "Tables/MeltR_fits/Fit_data",
                   file_prefix = "HP5_4_9")


list.df[[29]] = BL.fit$Ensemble.energies
list.df[[29]]$Helix = df$RNA[1]

####Consolidate results####

df = bind_rows(list.df)

####Remove GCCUUCGGGC from the data set####

df = df %>% filter(Helix != "GCCUUCGGGC")

write.csv(df, "Tables/MeltR_fits/Fit_results.csv", row.names = FALSE)
