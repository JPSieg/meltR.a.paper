library(tidyverse)
library(cowplot)
library(viridis)
library(MeltR)
devtools::load_all()

####Single state####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "CROWD DP1")

ggplot(df, aes(x = Temperature, y = Absorbance, color = Sample)) +
  geom_point() +
  geom_vline(xintercept = c(20,75))


#?meltR.A

fit = meltR.A(df,
              NucAcid = c("RNA", "CGCGCG"),
              Mmodel = "Homoduplex.2State",
              wavelength = 280,
              fitTs = c(20, 70))

v.log10n = seq(1, log10(5^length(unique(df$Sample))), length.out = 50)

v.n = floor(10^v.log10n)

v.n*0.5/1000

n.core = 3

sum(v.n[33:50])/n.core
sum(v.n[33:50])
sum(sum(v.n[c(40:46, 34,  39, 38)]), sum(v.n[c(33, 36,  48, 49, 35, 37)]), sum(v.n[c(47, 50)]))
#22243
sum(v.n[c(40:46, 34,  39, 38)])
sum(v.n[c(36,  48, 49, 35, 37)])
sum(v.n[c(47, 50)])

list.n = list(v.n[c(40:46, 34,  39, 38)],
              v.n[c(36,  48, 49, 35, 37)],
              v.n[c(47, 50)])

list.files("Figures/SI_Figure_X_Time_versus_n.comb")

BL = function(v){
  for ( i in 1:length(v)){
    Trim = BLTrimmer(fit, n.combinations =  v[[i]], Silent = T)
    df.fit = Trim$Ensemble.energies
    df.fit$S.time = Trim$System.time$all[1] + Trim$System.time$meltR.A[1]
    df.fit$n.comb = v[[i]]
    write.csv(df.fit,
              paste("Figures/SI_Figure_X_Time_versus_n.comb/Fit_results/",
                    v[[i]], "_BL_combinations_result.csv", sep = ""),
              row.names = FALSE)
  }
}
    
library(doParallel)

c1 = parallel::makeCluster(2)
doParallel::registerDoParallel(cores = n.core)
    
list.result = foreach::foreach(k = 1:4)  %dopar% {
  BL(list.n[[k]])
}

print("Done")

