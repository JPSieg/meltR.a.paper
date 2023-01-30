library(tidyverse)
devtools::load_all()

df.Meltwin

list.files("Tables")

df.MeltR = read.csv("Tables/MeltR_fits/Fit_results.csv") 

unique(df.MeltR$Method)

df3 = df.MeltR %>% filter(Method == "3 Global fit") %>%
  arrange(Helix) 

Helix = df3$Helix

H3 = paste(df3$dH, " (", df3$CI95.dH, ")", sep = "") 
S3 = paste(df3$dS, " (", df3$CI95.dS, ")", sep = "") 
G3 = paste(df3$dG, " (", df3$CI95.dG, ")", sep = "") 
Tm3 = paste(df3$Tm_at_0.1mM, " (", df3$CI95.Tm_at_0.1mM, ")", sep = "") 

df.Meltwin = df.Meltwin %>%
  arrange(Sequence)

df.Meltwin2 =  df.Meltwin %>% filter(Method == "2 Tm versus ln[Ct]")

H.MW2 = paste(df.Meltwin2$dH, " (\u00B1", df.Meltwin2$SE.dH, ")", sep = "")
S.MW2 = paste(df.Meltwin2$dS, " (\u00B1", df.Meltwin2$SE.dS, ")", sep = "")
G.MW2 = paste(df.Meltwin2$dG, " (\u00B1", df.Meltwin2$SE.dG, ")", sep = "")
Tm.MW2 = paste(df.Meltwin2$Tm)

pH = 100*abs(df3$dH - df.Meltwin2$dH)/abs((df3$dH + df.Meltwin2$dH)/2)
pS = 100*abs(df3$dS - df.Meltwin2$dS)/abs((df3$dS + df.Meltwin2$dS)/2)
pG = 100*abs(df3$dG - df.Meltwin2$dG)/abs((df3$dG + df.Meltwin2$dG)/2)
pTm = 100*abs(df3$Tm_at_0.1mM - df.Meltwin2$Tm)/abs((df3$Tm_at_0.1mM + df.Meltwin2$Tm)/2)

helix = unique(df.MeltR$Helix)

df.Table = data.frame(Helix, H3, H.MW2, S3, S.MW2, G3, G.MW2, Tm3, Tm.MW2, pH, pS, pG, pTm)

write.csv(df.Table, "Tables/04_BLtrimmer_results_comparing_MeltR_method_3_to_Meltwin_method_2/04_BLtrimmer_results_comparing_MeltR_method_3_to_Meltwin_method_2.csv", row.names = FALSE)
