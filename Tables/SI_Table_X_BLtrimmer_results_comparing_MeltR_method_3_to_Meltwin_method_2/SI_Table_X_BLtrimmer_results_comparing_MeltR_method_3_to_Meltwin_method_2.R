library(tidyverse)
library(ggrepel)
library(cowplot)
devtools::load_all()

df.Meltwin

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

df.Meltwin2 =  df.Meltwin %>% filter(Method == "2 Tm versus ln[Ct]")

H.MW2 = paste(df.Meltwin2$dH, " (\u00B1", df.Meltwin2$SE.dH, ")", sep = "")
S.MW2 = paste(df.Meltwin2$dS, " (\u00B1", df.Meltwin2$SE.dS, ")", sep = "")
G.MW2 = paste(df.Meltwin2$dG, " (\u00B1", df.Meltwin2$SE.dG, ")", sep = "")
Tm.MW2 = paste(df.Meltwin2$Tm)

pH = 100*abs(df3$H - df.Meltwin2$dH)/abs(mean(df3$H + df.Meltwin2$dH))
pS = 100*abs(df3$S - df.Meltwin2$dS)/abs(mean(df3$S + df.Meltwin2$dS))
pG = 100*abs(df3$G - df.Meltwin2$dG)/abs(mean(df3$G + df.Meltwin2$dG))
pTm = 100*abs(df3$Tm_at_0.1mM - df.Meltwin2$Tm)/abs(mean(df3$Tm_at_0.1mM + df.Meltwin2$Tm))

helix = unique(df.MeltR$Helix)

df.Table = data.frame(Helix, H3, H.MW2, S3, S.MW2, G3, G.MW2, Tm3, Tm.MW2, pH, pS, pG, pTm)

write.csv(df.Table, "Tables/SI_Table_X_BLtrimmer_results_comparing_MeltR_method_3_to_Meltwin_method_2/SI_Table_X_BLtrimmer_results_comparing_MeltR_method_3_to_Meltwin_method_2.csv", row.names = FALSE)
