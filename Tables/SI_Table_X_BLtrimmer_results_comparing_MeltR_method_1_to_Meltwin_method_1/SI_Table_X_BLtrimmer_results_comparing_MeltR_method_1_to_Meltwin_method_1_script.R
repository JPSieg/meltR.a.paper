library(tidyverse)
library(ggrepel)
library(cowplot)
devtools::load_all()

df.Meltwin

list.files("Tables/SI_Table_X_BLtrimmer_results_comparing_MeltR_method_1_to_Meltwin_method_1")

df.MeltR = read.csv("Tables/SI_Table_X_MeltR_fits/Fit_results.csv") 

unique(df.MeltR$Method)

df1 = df.MeltR %>% filter(Method == "1 individual fits") %>%
  arrange(Helix) 

Helix = df1$Helix

H1 = paste(df1$S, " (", df1$CI95.H, ")", sep = "") 
S1 = paste(df1$S, " (", df1$CI95.S, ")", sep = "") 
G1 = paste(df1$G, " (", df1$CI95.G, ")", sep = "") 
Tm1 = paste(df1$Tm_at_0.1mM, " (", df1$CI95.Tm_at_0.1mM, ")", sep = "") 

df.Meltwin = df.Meltwin %>%
  filter(Method == "1 individual fits") %>%
  arrange(Sequence) %>%
  filter(Sequence != "UAUAUAUA")

H.MW = paste(df.Meltwin$dH, " (\u00B1", df.Meltwin$SE.dH, ")", sep = "")
S.MW = paste(df.Meltwin$dS, " (\u00B1", df.Meltwin$SE.dS, ")", sep = "")
G.MW = paste(df.Meltwin$dG, " (\u00B1", df.Meltwin$SE.dG, ")", sep = "")
Tm.MW = paste(df.Meltwin$Tm)

pH = 100*abs(df1$H - df.Meltwin$dH)/abs((df1$H + df.Meltwin$dH)/2)
pS = 100*abs(df1$S - df.Meltwin$dS)/abs((df1$S + df.Meltwin$dS)/2)
pG = 100*abs(df1$G - df.Meltwin$dG)/abs((df1$G + df.Meltwin$dG)/2)
pTm = 100*abs(df1$Tm_at_0.1mM - df.Meltwin$Tm)/abs((df1$Tm_at_0.1mM + df.Meltwin$Tm)/2)

df.Table = data.frame(Helix, H1, H.MW, S1, S.MW, G1, G.MW, Tm1, Tm.MW, pH, pS, pG, pTm)

write.csv(df.Table, "Tables/SI_Table_X_BLtrimmer_results_comparing_MeltR_method_1_to_Meltwin_method_1/SI_Table_X_BLtrimmer_results_comparing_MeltR_method_1_to_Meltwin_method_1.csv", row.names = FALSE)

