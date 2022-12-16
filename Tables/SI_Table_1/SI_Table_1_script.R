library(tidyverse)
devtools::load_all()

list.files("Tables/")

df.MeltR = read.csv("Tables/MeltR_fits/Fit_results.csv") 

unique(df.MeltR$Method)

####MeltR method 1####

df1 = df.MeltR %>% filter(Method == "1 Individual fits") %>% arrange(Helix)

Helix = df1$Helix

H1 = paste(df1$dH, " (", df1$CI95.dH, ")", sep = "") 
S1 = paste(df1$dS, " (", df1$CI95.dS, ")", sep = "") 
G1 = paste(df1$dG, " (", df1$CI95.dG, ")", sep = "") 
Tm1 = paste(df1$Tm_at_0.1mM, " (", df1$CI95.Tm_at_0.1mM, ")", sep = "") 

df.1 = data.frame(H1, S1, G1, Tm1)

colnames(df.1) = paste("MeltR", colnames(df.1))

####MeltR method 2####

df2 = df.MeltR %>% filter(Method == "2 Tm versus ln[Ct]") %>% arrange(Helix)

H2 = paste(df2$dH, " (", df2$CI95.dH, ")", sep = "") 
S2 = paste(df2$dS, " (", df2$CI95.dS, ")", sep = "") 
G2 = paste(df2$dG, " (", df2$CI95.dG, ")", sep = "") 
Tm2 = paste(df2$Tm_at_0.1mM, " (", df2$CI95.Tm_at_0.1mM, ")", sep = "") 

df.2 = data.frame(H2, S2, G2, Tm2)

colnames(df.2) = paste("MeltR", colnames(df.2))

####MeltR Method 3####

df3 = df.MeltR %>% filter(Method == "3 Global fit") %>% arrange(Helix)

H3 = paste(df3$dH, " (", df3$CI95.dH, ")", sep = "") 
S3 = paste(df3$dS, " (", df3$CI95.dS, ")", sep = "") 
G3 = paste(df3$dG, " (", df3$CI95.dG, ")", sep = "") 
Tm3 = paste(df3$Tm_at_0.1mM, " (", df3$CI95.Tm_at_0.1mM, ")", sep = "") 

df.3 = data.frame(H3, S3, G3, Tm3)

colnames(df.3) = paste("MeltR", colnames(df.3))

####MeltWin Method 1####

df4 = df.Meltwin %>%
  filter(Method == "1 individual fits") %>%
  filter(Sequence != "UAUAUAUA") %>%
  arrange(Sequence)

H1 = paste(df4$dH, " (", df4$SE.dH, ")", sep = "") 
S1 = paste(df4$dS, " (", df4$SE.dS, ")", sep = "") 
G1 = paste(df4$dG, " (", df4$SE.dG, ")", sep = "") 
Tm1 = paste(df4$Tm, sep = "") 

df.4 = data.frame(H1, S1, G1, Tm1)

colnames(df.4) = paste("MeltWin", colnames(df.4))

####MeltWin Method 2####

df5 = df.Meltwin %>%
  filter(Method == "2 Tm versus ln[Ct]") %>%
  filter(Sequence != "UAUAUAUA") %>%
  arrange(Sequence)

H2 = paste(df5$dH, " (", df5$SE.dH, ")", sep = "") 
S2 = paste(df5$dS, " (", df5$SE.dS, ")", sep = "") 
G2 = paste(df5$dG, " (", df5$SE.dG, ")", sep = "") 
Tm2 = paste(df5$Tm, sep = "") 

df.5 = data.frame(H2, S2, G2, Tm2)

colnames(df.5) = paste("MeltWin", colnames(df.5))

####Consolidate data####

df.Table = bind_cols(df.1, df.2, df.3, df.4, df.5) %>%
  select("MeltR H1", "MeltR H2", "MeltR H3", "MeltWin H1", "MeltWin H2",
         "MeltR S1", "MeltR S2", "MeltR S3", "MeltWin S1", "MeltWin S2",
         "MeltR G1", "MeltR G2", "MeltR G3", "MeltWin G1", "MeltWin G2",
         "MeltR Tm1", "MeltR Tm2", "MeltR Tm3", "MeltWin Tm1", "MeltWin Tm2")

colnames(df.Table)

####Write the table####

write.csv(df.Table, "Tables/SI_Table_1/SI_Table_1.csv", row.names = FALSE)
