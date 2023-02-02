library(tidyverse)
devtools::load_all()

list.files("Tables/")

df.MeltR = read.csv("Tables/MeltR_fits/Fit_results.csv") 

df1 = df.MeltR %>% filter(Method == "1 Individual fits")

Helix = df1$Helix

H1 = paste(df1$dH, " (", df1$CI95.dH, ")", sep = "") 
S1 = paste(df1$dS, " (", df1$CI95.dS, ")", sep = "") 
G1 = paste(df1$dG, " (", df1$CI95.dG, ")", sep = "") 
Tm1 = paste(df1$Tm_at_0.1mM, " (", df1$CI95.Tm_at_0.1mM, ")", sep = "") 

df2 = df.MeltR %>% filter(Method == "2 Tm versus ln[Ct]")

H2 = paste(df2$dH, " (", df2$CI95.dH, ")", sep = "") 
S2 = paste(df2$dS, " (", df2$CI95.dS, ")", sep = "") 
G2 = paste(df2$dG, " (", df2$CI95.dG, ")", sep = "") 
Tm2 = paste(df2$Tm_at_0.1mM, " (", df2$CI95.Tm_at_0.1mM, ")", sep = "") 

df3 = df.MeltR %>% filter(Method == "3 Global fit")

H3 = paste(df3$dH, " (", df3$CI95.dH, ")", sep = "") 
S3 = paste(df3$dS, " (", df3$CI95.dS, ")", sep = "") 
G3 = paste(df3$dG, " (", df3$CI95.dG, ")", sep = "") 
Tm3 = paste(df3$Tm_at_0.1mM, " (", df3$CI95.Tm_at_0.1mM, ")", sep = "") 

df.Table = data.frame(Helix, H1, H2, H3, S1, S2, S3, G1, G2, G3, Tm1, Tm2, Tm3)

pH = c()
pS = c()
pG = c()
pTm = c()

helix = unique(df.MeltR$Helix)

for (i in 1:length(Helix)){
  df = df.MeltR %>% filter(Helix == helix[i])
  pH[i] = 100*abs(range(df %>% select(dH), na.rm = T)[1]-range(df %>% select(dH), na.rm = T)[2])/abs(mean(df$dH, na.rm = T))
  pS[i] = 100*abs(range(df %>% select(dS), na.rm = T)[1]-range(df %>% select(dS), na.rm = T)[2])/abs(mean(df$dS, na.rm = T))
  pG[i] = 100*abs(range(df %>% select(dG), na.rm = T)[1]-range(df %>% select(dG), na.rm = T)[2])/abs(mean(df$dG, na.rm = T))
  pTm[i] = 100*abs(range(df %>% select(Tm_at_0.1mM), na.rm = T)[1]-range(df %>% select(Tm_at_0.1mM), na.rm = T)[2])/abs(mean(df$Tm_at_0.1mM, na.rm = T))
}

df.Table$pH = pH
df.Table$pS = pS
df.Table$pG = pG
df.Table$pTm = pTm

head(df.Table)

write.csv(df.Table, "Tables/05_BLtrimmer_results_comparing_MeltR_methods_1_2_3/05_BLtrimmer_results_comparing_MeltR_methods_1_2_3.csv", row.names = FALSE)
