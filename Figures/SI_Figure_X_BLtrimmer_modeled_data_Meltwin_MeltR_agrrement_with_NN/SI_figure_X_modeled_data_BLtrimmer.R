library(MeltR)
library(tidyverse)

####Read in data####

list.files("Figures/")

df = read.csv("Figures/SI_Figure_X_BLtrimmer_modeled_data_Meltwin_MeltR_agrrement_with_NN/Xia_et_al_data.csv")

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
    if (df$R[i] == ""){
      RNA = c("RNA", df$F[i])
      Duplex = "Monomolecular.2State"
    }else{
      RNA = c("RNA", df$F[i], df$R[i])
      Duplex = "Heteroduplex.2State"
    }
  }
  
  if (Duplex == "Monomolecular.2State"){
    
    
    Extinct = MeltR::calc.extcoeff(RNA)
    Extinct = Extinct$Total
    
    Extinct.lower = Extinct*0.6
    
    Samples = 9
    
    Sample = 1:Samples
    
    start.Ct = 0.2/Extinct.lower
    end.Ct = 2/(Extinct*0.1)
    
    lnCt = seq(log(start.Ct), log(end.Ct), length.out = Samples)
    
    Ct = exp(lnCt)
    
    A.1cm = Ct*Extinct
    A.0.5cm = Ct*Extinct*0.5
    
    Pathlength = rep(NA, length(Ct))
    
    Pathlength[which(A.1cm <= 2)] = 1
    Pathlength[which(A.1cm >= 2)] = 0.5
    Pathlength[which(A.0.5cm >= 2)] = 0.1
    
    mED = c()
    mSS = c()
    bED = c()
    bSS = c()
    
    for (j in 1:Samples){
      mED[j] = rnorm(1, 0.001113342, 0.0007760626)
      mSS[j] = rnorm(1, 0.001113342, 0.0007760626)
      bED[j] = Ct[j]*Pathlength[j]*Extinct.lower - mED[j]*90
      bSS[j] = Ct[j]*Pathlength[j]*Extinct - mED[j]*90
    }
    
    seqF = RNA[2]
    seqR = NA
    dG = df$G[i]
    dH = df$H[i]
    dS = df$S[i]
    Tm.at.0.1mM = dH/(dS/1000) - 273.15
    
    df2 = data.frame(seqF, seqR, dG, dH, dS, Tm.at.0.1mM, Extinct, Extinct.lower, mED, mSS, bED, bSS, Sample, Ct,  Pathlength)
    
    output = list(df2)
    
    ####Model Absorbance data####
    
    list.df = {}
    Absorbance_error = 0.0005
    
    for (j in 1:nrow(df2)){
      Sample = j
      Temperature = seq(5, 95, 0.5)
      Pathlength = df2$Pathlength[j]
      f = 1/(exp((df2$dS[j]/1.9872) - (df2$dH[j]/(0.0019872*(273.15 + Temperature)))) + 1)
      Absorbance = (df2$mED[j]*Temperature + df2$bED[j])*(1 - f) +  (df2$mSS[j]*Temperature + df2$bSS[j])*f + rnorm(length(Temperature), 0, Absorbance_error)
      list.df[[j]] = data.frame(Sample, Temperature, Pathlength, Absorbance)
    }
    
    df3 = list.df[[1]]
    j = 2
    while(j < (1 + Sample)){
      df3 = rbind(df3, list.df[[j]])
      j = j + 1
    }
    
    output[[2]] = df3
    
    l.model  <- output
    
    list.df.model[[i]] = l.model[[2]]
    list.df.model[[i]]$seqF =  df$F[i]
    list.df.model[[i]]$seqR =  df$R[i]
    list.df.params[[i]] = l.model[[1]][1,] %>% select(seqF, seqR, dG, dH, dS, Tm.at.0.1mM)
    
  }else{
    l.model = meltR.A.model(NucAcid = RNA,
                            Mmodel = Duplex,
                            Absorbance_error = 0.0005)
    
    list.df.model[[i]] = l.model[[2]]
    list.df.model[[i]]$seqF =  df$F[i]
    list.df.model[[i]]$seqR =  df$R[i]
    list.df.params[[i]] = l.model[[1]][1,] %>% select(seqF, seqR, dG, dH, dS, Tm.at.0.1mM)
  }
  
}

####Write data####

df.Modeled_absorbance = bind_rows(list.df.model)
df.Energies = bind_rows(list.df.params)

write.csv(df.Modeled_absorbance, "Figures/SI_Figure_X_BLtrimmer_modeled_data_Meltwin_MeltR_agrrement_with_NN/Modeled_Xia_data.csv", row.names = FALSE)
write.csv(df.Energies, "Figures/SI_Figure_X_BLtrimmer_modeled_data_Meltwin_MeltR_agrrement_with_NN/Modeled_Xia_energies.csv", row.names = FALSE)
