library(tidyverse)
library(MeltR)
library(doParallel)

####Read in data####

list.files("Figures/SI_Figure_X_BLtrimmer_modeled_data")

df.total = read.csv("Figures/SI_Figure_X_BLtrimmer_modeled_data/Modeled_Xia_data.csv")

Experiments = unique(df.total$seqF)

list.df = {}

pb = txtProgressBar(min = 1, max = length(Experiments), initial = 1, style = 3)

for (i in 1:length(Experiments)){
  df = df.total %>% filter(seqF == Experiments[i])
  #setTxtProgressBar(pb, i)
  if (df$seqF[i] == df$seqR[i]){
    RNA = c("RNA", df$seqF[i])
    Duplex = "Homoduplex.2State"
  }else{
    RNA = c("RNA", df$seqF[i], df$seqR[i])
    Duplex = "Heteroduplex.2State"
  }
  
  print(RNA)
  
  fit = meltR.A(df,
                NucAcid = RNA,
                Mmodel = Duplex,
                Silent = FALSE)
  
  BL.fit = BLTrimmer(fit,
                     Silent = FALSE)
  
  list.df[[i]] = BL.fit$Ensemble.energies
  list.df[[i]]$Helix = df$seqF[i]
}

df = bind_rows(list.df)

write.csv(df, "Figures/SI_Figure_X_BLtrimmer_modeled_data/Modeled_fit_results.csv", row.names = FALSE)
