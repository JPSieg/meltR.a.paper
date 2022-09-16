library(tidyverse)
library(ggrepel)
library(cowplot)
devtools::load_all()

df.Meltwin

list.files("Tables/SI_Table_X_BLtrimmer_results_comparing_MeltR_method_3_to_Meltwin_method_1_and_2")

df.MeltR = read.csv("Tables/SI_Table_X_MeltR_fits/Fit_results.csv") 

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

df.Meltwin2 =  df.Meltwin %>% filter(Method == "2 Tm versus ln[Ct]")

H.MW2 = paste(df.Meltwin2$dH, " (\u00B1", df.Meltwin2$SE.dH, ")", sep = "")
S.MW2 = paste(df.Meltwin2$dS, " (\u00B1", df.Meltwin2$SE.dS, ")", sep = "")
G.MW2 = paste(df.Meltwin2$dG, " (\u00B1", df.Meltwin2$SE.dG, ")", sep = "")
Tm.MW2 = paste(df.Meltwin2$Tm)

pH = c()
pS = c()
pG = c()
pTm = c()

helix = unique(df.MeltR$Helix)

for (i in 1:length(Helix)){
  df = df3 %>% filter(Helix == helix[i])
  dfM = df.Meltwin %>% filter(Sequence == helix[i])
  pH[i] = 100*max(abs(df$H[1] - dfM$dH))/abs(mean(c(df$H, dfM$dH)))
  pS[i] = 100*max(abs(df$S[1] - dfM$dS))/abs(mean(c(df$S, dfM$dS)))
  pG[i] = 100*max(abs(df$G[1] - dfM$dG))/abs(mean(c(df$G, dfM$dG)))
  pTm[i] = 100*max(abs(df$Tm_at_0.1mM[1] - dfM$Tm))/abs(mean(c(df$Tm_at_0.1mM, dfM$Tm)))
}

df.Table = data.frame(Helix, H3, H.MW1, H.MW2, S3, S.MW1, S.MW2, G3, G.MW1, G.MW2, Tm3, Tm.MW1, Tm.MW2, pH, pS, pG, pTm)

write.csv(df.Table, "Tables/SI_Table_X_BLtrimmer_results_comparing_MeltR_method_3_to_Meltwin_method_1_and_2/SI_Table_X_BLtrimmer_results_comparing_MeltR_method_3_to_Meltwin_method_1_and_2.csv", row.names = FALSE)
