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

v.n = v.n[1:39]

v.n*0.5/1000

n.core = 3
list.files("Figures/SI_Figure_3_Time_versus_n.comb")

BL = function(v){
    Trim = BLTrimmer(fit, n.combinations =  v, Silent = T)
    df.fit = Trim$Ensemble.energies
    df.fit$S.time = Trim$System.time$all[1] + Trim$System.time$meltR.A[1]
    df.fit$n.comb = v
    write.csv(df.fit,
              paste("Figures/SI_Figure_3_Time_versus_n.comb/Fit_results/",
                    v, "_BL_combinations_result.csv", sep = ""),
              row.names = FALSE)
}
    
library(doParallel)

c1 = parallel::makeCluster(2)
doParallel::registerDoParallel(cores = n.core)
    
list.result = foreach::foreach(k = 1:length(v.n))  %dopar% {
  BL(v.n[k])
}

print("Done")

