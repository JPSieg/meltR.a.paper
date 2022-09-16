library(MeltR)
library(tidyverse)

####Read in data####

list.files("Figures/SI_Figure_X_BLtrimmer_modeled_data")

df = read.csv("Figures/SI_Figure_X_BLtrimmer_modeled_data/Xia_et_al_data.csv")

head(df)

list.df.model = {}
list.df.params = {}

?meltR.A.model

####model data####

for (i in 1:nrow(df)){
  
  if (df$F[i] == df$R[i]){
    RNA = c("RNA", df$F[i])
    Duplex = "Homoduplex.2State"
  }else{
    RNA = c("RNA", df$F[i], df$R[i])
    Duplex = "Heteroduplex.2State"
  }
  
  l.model = meltR.A.model(NucAcid = RNA,
                          Mmodel = Duplex,
                          Absorbance_error = 0.001)
  
  list.df.model[[i]] = l.model[[2]]
  list.df.model[[i]]$seqF =  df$F[i]
  list.df.model[[i]]$seqR =  df$R[i]
  list.df.params[[i]] = l.model[[1]][1,] %>% select(seqF, seqR, dG, dH, dS, Tm.at.0.1mM)
  
}

####Write data####

df.Modeled_absorbance = bind_rows(list.df.model)
df.Energies = bind_rows(list.df.params)

write.csv(df.Modeled_absorbance, "Figures/SI_Figure_X_BLtrimmer_modeled_data/Modeled_Xia_data.csv", row.names = FALSE)
write.csv(df.Energies, "Figures/SI_Figure_X_BLtrimmer_modeled_data/Modeled_Xia_energies.csv", row.names = FALSE)
