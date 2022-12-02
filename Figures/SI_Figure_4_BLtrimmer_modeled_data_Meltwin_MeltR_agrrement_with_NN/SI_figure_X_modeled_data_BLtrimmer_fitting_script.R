library(tidyverse)
library(MeltR)
library(doParallel)

####Read in data####

list.files("Figures/")

df.total = read.csv("Figures/SI_Figure_4_BLtrimmer_modeled_data_Meltwin_MeltR_agrrement_with_NN/Modeled_Xia_data.csv")

Experiments = unique(df.total$seqF)

Experiments = Experiments[-2]

list.df = {}

pb = txtProgressBar(min = 1, max = length(Experiments), initial = 1, style = 3)

for (i in 1:length(Experiments)){
  df = df.total %>% filter(seqF == Experiments[i])
  #setTxtProgressBar(pb, i)
  if (df$seqF[i] == df$seqR[i]){
    RNA = c("RNA", df$seqF[i])
    Duplex = "Homoduplex.2State"
  }else{
    if (df$seqR[i] == ""){
      RNA = c("RNA", df$seqF[i])
      Duplex = "Monomolecular.2State"
    }else{
      RNA = c("RNA", df$seqF[i], df$seqR[i])
      Duplex = "Heteroduplex.2State"
    }
  }
  
  fit = meltR.A(df,
                NucAcid = RNA,
                Mmodel = Duplex,
                Silent = F)
  
  BL.fit = BLTrimmer(fit,
                     Silent = F)
  
  list.df[[i]] = BL.fit$Ensemble.energies
  list.df[[i]]$Helix = df$seqF[i]
}

df = bind_rows(list.df)

write.csv(df, "Figures/SI_Figure_4_BLtrimmer_modeled_data_Meltwin_MeltR_agrrement_with_NN/Modeled_fit_results.csv", row.names = FALSE)
