library(tidyverse)
library(MeltR)
devtools::load_all()

####Read in data####

list.files("Figures/SI_Figure_X_BLtrimmer_modeled_data")

df = read.csv("Figures/SI_Figure_X_BLtrimmer_modeled_data/Modeled_fit_results.csv")

df
