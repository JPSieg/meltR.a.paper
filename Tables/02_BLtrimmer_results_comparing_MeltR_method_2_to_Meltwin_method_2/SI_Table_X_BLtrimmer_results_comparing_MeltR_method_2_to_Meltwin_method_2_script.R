library(tidyverse)
library(ggrepel)
library(cowplot)
devtools::load_all()

df.Meltwin

list.files("Tables/SI_Table_X_BLtrimmer_results_comparing_MeltR_method_2_to_Meltwin_method_2")

df.MeltR = read.csv("Tables/SI_Table_X_MeltR_fits/Fit_results.csv") 

unique(df.MeltR$Method)

df2 = df.MeltR %>% filter(Method == "2 Tm versus ln[Ct]") %>%
  arrange(Helix) 

Helix = df2$Helix

H2 = paste(df2$S, " (", df2$CI95.H, ")", sep = "") 
S2 = paste(df2$S, " (", df2$CI95.S, ")", sep = "") 
G2 = paste(df2$G, " (", df2$CI95.G, ")", sep = "") 
Tm2 = paste(df2$Tm_at_0.1mM, " (", df2$CI95.Tm_at_0.1mM, ")", sep = "") 

df.Meltwin = df.Meltwin %>%
  filter(Method == "2 Tm versus ln[Ct]") %>%
  arrange(Sequence) %>%
  filter(Sequence != "UAUAUAUA")

H.MW = paste(df.Meltwin$dH, " (\u00B1", df.Meltwin$SE.dH, ")", sep = "")
S.MW = paste(df.Meltwin$dS, " (\u00B1", df.Meltwin$SE.dS, ")", sep = "")
G.MW = paste(df.Meltwin$dG, " (\u00B1", df.Meltwin$SE.dG, ")", sep = "")
Tm.MW = paste(df.Meltwin$Tm)

pH = 100*abs(df2$H - df.Meltwin$dH)/abs((df2$H + df.Meltwin$dH)/2)
pS = 100*abs(df2$S - df.Meltwin$dS)/abs((df2$S + df.Meltwin$dS)/2)
pG = 100*abs(df2$G - df.Meltwin$dG)/abs((df2$G + df.Meltwin$dG)/2)
pTm = 100*abs(df2$Tm_at_0.1mM - df.Meltwin$Tm)/abs((df2$Tm_at_0.1mM + df.Meltwin$Tm)/2)

df.Table = data.frame(Helix, H2, H.MW, S2, S.MW, G2, G.MW, Tm2, Tm.MW, pH, pS, pG, pTm)

write.csv(df.Table, "Tables/SI_Table_X_BLtrimmer_results_comparing_MeltR_method_2_to_Meltwin_method_2/SI_Table_X_BLtrimmer_results_comparing_MeltR_method_2_to_Meltwin_method_2.csv", row.names = FALSE)

