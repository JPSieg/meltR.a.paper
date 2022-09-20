library(tidyverse)
library(ggrepel)
library(cowplot)
devtools::load_all()

df.Meltwin

list.files("Tables/")

df.MeltR = read.csv("Tables/MeltR_fits/Fit_results.csv") 

unique(df.MeltR$Method)

df3 = df.MeltR %>% filter(Method == "3 Global fit") %>%
  arrange(Helix) 

Helix = df3$Helix

H3 = paste(df3$H, " (", df3$CI95.H, ")", sep = "") 
S3 = paste(df3$S, " (", df3$CI95.S, ")", sep = "") 
G3 = paste(df3$G, " (", df3$CI95.G, ")", sep = "") 
Tm3 = paste(df3$Tm_at_0.1mM, " (", df3$CI95.Tm_at_0.1mM, ")", sep = "") 

df.Meltwin = df.Meltwin %>%
  arrange(Sequence) %>%
  filter(Sequence != "UAUAUAUA")

df.Meltwin1 =  df.Meltwin %>% filter(Method == "1 individual fits")

H.MW1 = paste(df.Meltwin1$dH, " (\u00B1", df.Meltwin1$SE.dH, ")", sep = "")
S.MW1 = paste(df.Meltwin1$dS, " (\u00B1", df.Meltwin1$SE.dS, ")", sep = "")
G.MW1 = paste(df.Meltwin1$dG, " (\u00B1", df.Meltwin1$SE.dG, ")", sep = "")
Tm.MW1 = paste(df.Meltwin1$Tm)

pH = 100*abs(df3$H - df.Meltwin1$dH)/abs((df3$H + df.Meltwin1$dH)/2)
pS = 100*abs(df3$S - df.Meltwin1$dS)/abs((df3$S + df.Meltwin1$dS)/2)
pG = 100*abs(df3$G - df.Meltwin1$dG)/abs((df3$G + df.Meltwin1$dG)/2)
pTm = 100*abs(df3$Tm_at_0.1mM - df.Meltwin1$Tm)/abs((df3$Tm_at_0.1mM + df.Meltwin1$Tm)/2)

df.Table = data.frame(Helix, H3, H.MW1, S3, S.MW1, G3, G.MW1, Tm3, Tm.MW1, pH, pS, pG, pTm)

write.csv(df.Table, "Tables/SI_Table_X_BLtrimmer_results_comparing_MeltR_method_3_to_Meltwin_method_1/SI_Table_X_BLtrimmer_results_comparing_MeltR_method_3_to_Meltwin_method_1.csv", row.names = FALSE)
