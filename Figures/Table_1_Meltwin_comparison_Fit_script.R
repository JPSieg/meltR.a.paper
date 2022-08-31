devtools::load_all()
library(MeltR)
library(tidyverse)

list.df.result = {}
list.df.indv.fit = {}

####"CROWD DP1####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "CROWD DP1")

head(df)

ggplot(df %>% filter(Sample == 12),
       aes(x = Temperature,
           y = Absorbance)) +
  geom_point() +
  geom_vline(xintercept = c(10,75))

unique(df$Sample)

T.range = list(c(27,83), #2
     c(27,75), #3
     c(15,75), #4
     c(25,75), #6
     c(25,75), #7
     c(20,75), #8
     c(15,85), #10
     c(20,70), #11
     c(15,70)) #12 

fit = meltR.A(df,
         NucAcid = c("RNA", "CGCGCG"),
         Mmodel = "Homoduplex.2State",
         blank = "none",
        concT = 10,
        fitTs = T.range,
        wavelength = 280,
        Save_results = "all",
        file_path = "Figures/Table_1_Meltwin_comparison/CROWD DP1")

?BLTrimmer
BLTrimmer(fit,
          Trim.method = "floating", 
          Assess.method = 2,
          n.combinations = 200,
          Save_results = "all")

BLTrimmer(fit,
          Assess.method = 2,
          n.combinations = 500,
          Save_results = "all")

BLTrimmer(fit,
          Assess.method = 3,
          n.combinations = 500,
          Save_results = "all")

BLTrimmer(fit,
          Trim.method = "fixed",
          Assess.method = 1,
          n.combinations = 500,
          Save_results = "all")

BLTrimmer(fit,
          Trim.method = "fixed",
          Assess.method = 2,
          n.combinations = 500,
          Save_results = "all")

BLTrimmer(fit,
          Trim.method = "fixed",
          Assess.method = 3,
          n.combinations = 500,
          Save_results = "all")


i = 1
list.df.result[[i]] = fit$Summary
list.df.result[[i]]$Experiment = df$Experiment[1]
list.df.indv.fit[[i]] = fit$Method.1.indvfits
list.df.indv.fit[[i]]$Experiment = df$Experiment[1]

####Helix A####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "Helix A")

head(df)

ggplot(df %>% filter(Sample == 10),
       aes(x = Temperature,
           y = Absorbance)) +
  geom_point() +
  geom_vline(xintercept = c(30,80))

unique(df$Sample)

T.range = list(c(15,70), #2
               c(15,80), #3
               c(15,80), #4
               c(15,85), #5
               c(25,77), #6
               c(25,80), #7
               c(25,80), #8
               c(30,80), #9
               c(30,80)) #10 

fit = meltR.A(df,
              NucAcid = c("RNA", "CGAAAGGU", "ACCUUUCG"),
              Mmodel = "Heteroduplex.2State",
              blank = 1,
              concT = 10,
              fitTs = T.range,
              wavelength = 260,
              Save_results = "all",
              file_path = "Figures/Table_1_Meltwin_comparison/Helix A")

i = i + 1
list.df.result[[i]] = fit$Summary
list.df.result[[i]]$Experiment = df$Experiment[1]
list.df.indv.fit[[i]] = fit$Method.1.indvfits
list.df.indv.fit[[i]]$Experiment = df$Experiment[1]

####Helix C####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "Helix C")

head(df)

unique(df$Sample)

ggplot(df %>% filter(Sample == 10),
       aes(x = Temperature,
           y = Absorbance)) +
  geom_point() +
  geom_vline(xintercept = c(25,80))

unique(df$Sample)

T.range = list(c(10,70), #1
               c(20,75), #3
               c(15,63), #4
               c(15,75), #5
               c(20,68), #6
               c(20,75), #7
               c(20,75), #8
               c(20,80), #9
               c(25,80)) #10 

fit = meltR.A(df,
              NucAcid = c("RNA", "CUGAGUC", "GACUCAG"),
              Mmodel = "Heteroduplex.2State",
              blank = 2,
              concT = 10,
              fitTs = T.range,
              wavelength = 260,
              Save_results = "all",
              file_path = "Figures/Table_1_Meltwin_comparison/Helix C")

i = i + 1
list.df.result[[i]] = fit$Summary
list.df.result[[i]]$Experiment = df$Experiment[1]
list.df.indv.fit[[i]] = fit$Method.1.indvfits
list.df.indv.fit[[i]]$Experiment = df$Experiment[1]

####Helix D####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "Helix D")

head(df)

unique(df$Sample)

ggplot(df %>% filter(Sample == 10),
       aes(x = Temperature,
           y = Absorbance)) +
  geom_point() +
  geom_vline(xintercept = c(24,70))

unique(df$Sample)

T.range = list(c(5,60), #1
               c(5,60), #2
               c(8,70), #4
               c(8,70), #5
               c(8,70), #6
               c(8,70), #7
               c(8,70), #8
               c(15,70), #9
               c(24,70)) #10 

fit = meltR.A(df,
              NucAcid = c("RNA", "CUGAGUC", "GACUCAG"),
              Mmodel = "Heteroduplex.2State",
              blank = 3,
              concT = 10,
              fitTs = T.range,
              wavelength = 260,
              Save_results = "all",
              file_path = "Figures/Table_1_Meltwin_comparison/Helix D")

i = i + 1
list.df.result[[i]] = fit$Summary
list.df.result[[i]]$Experiment = df$Experiment[1]
list.df.indv.fit[[i]] = fit$Method.1.indvfits
list.df.indv.fit[[i]]$Experiment = df$Experiment[1]

####Helix A labeled####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "Helix A labeled")

head(df)

unique(df$Sample)

ggplot(df %>% filter(Sample == 9),
       aes(x = Temperature,
           y = Absorbance)) +
  geom_point() +
  geom_vline(xintercept = c(40,95))

unique(df$Sample)

T.range = list(c(24,80), #2
               c(25,85), #3
               c(30,85), #4
               c(30,90), #5
               c(40,85), #6
               c(35,95), #7
               c(35,95), #8
               c(40,95)) #9 

head(df)

fit = meltR.A(df,
              NucAcid = c("RNA", "CGAAAGGU", "ACCUUUCG"),
              Mmodel = "Heteroduplex.2State",
              blank = "none",
              concT = 10,
              fitTs = T.range,
              wavelength = 260,
              Save_results = "all",
              file_path = "Figures/Table_1_Meltwin_comparison/Helix A labeled")

i = i + 1
list.df.result[[i]] = fit$Summary
list.df.result[[i]]$Experiment = df$Experiment[1]
list.df.indv.fit[[i]] = fit$Method.1.indvfits
list.df.indv.fit[[i]]$Experiment = df$Experiment[1]

####Helix C labeled####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "Helix C labeled")

head(df)

unique(df$Sample)

ggplot(df %>% filter(Sample == 9),
       aes(x = Temperature,
           y = Absorbance)) +
  geom_point() +
  geom_vline(xintercept = c(35,90))

unique(df$Sample)

T.range = list(c(25,85), #2
               c(25,85), #3
               c(30,85), #4
               c(30,90), #5
               c(35,90), #6
               c(35,81), #7
               c(35,90), #8
               c(35,90)) #9 

head(df)

fit = meltR.A(df,
              NucAcid = c("RNA", "CUGAGUC", "GACUCAG"),
              Mmodel = "Heteroduplex.2State",
              blank = "none",
              concT = 10,
              fitTs = T.range,
              wavelength = 260,
              Save_results = "all",
              file_path = "Figures/Table_1_Meltwin_comparison/Helix C labeled")

i = i + 1
list.df.result[[i]] = fit$Summary
list.df.result[[i]]$Experiment = df$Experiment[1]
list.df.indv.fit[[i]] = fit$Method.1.indvfits
list.df.indv.fit[[i]]$Experiment = df$Experiment[1]

####Helix D labeled####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "Helix D labeled")

head(df)

unique(df$Sample)

ggplot(df %>% filter(Sample == 9),
       aes(x = Temperature,
           y = Absorbance)) +
  geom_point() +
  geom_vline(xintercept = c(25,85))

unique(df$Sample)

T.range = list(c(20,75), #3
               c(20,75), #4
               c(20,80), #5
               c(20,80), #6
               c(25,80), #7
               c(30,80), #8
               c(25,85)) #9  

head(df)

fit = meltR.A(df,
              NucAcid = c("RNA", "CGUUGC", "GCAACG"),
              Mmodel = "Heteroduplex.2State",
              blank = "none",
              concT = 10,
              fitTs = T.range,
              wavelength = 260,
              Save_results = "all",
              file_path = "Figures/Table_1_Meltwin_comparison/Helix D labeled")

i = i + 1
list.df.result[[i]] = fit$Summary
list.df.result[[i]]$Experiment = df$Experiment[1]
list.df.indv.fit[[i]] = fit$Method.1.indvfits
list.df.indv.fit[[i]]$Experiment = df$Experiment[1]

####Consolidate MeltR data####

df.MeltR.summary = bind_rows(list.df.result)

usethis::use_data(df.MeltR.summary, overwrite = T)

df.MeltR.indv.fits = bind_rows(list.df.indv.fit)

usethis::use_data(df.MeltR.indv.fits, overwrite = T)
