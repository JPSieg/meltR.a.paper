library(tidyverse)
library(ggrepel)
library(cowplot)
devtools::load_all()

list.files("Tables/SI_Table_X_BLtrimmer_results_comparing_MeltR_methods_1_2_3")

df.MeltR = read.csv("Tables/SI_Table_X_MeltR_fits/Fit_results.csv") 

unique(df.MeltR$Method)


df1 = df.MeltR %>% filter(Method == "1 individual fits")

Helix = df1$Helix

H1 = paste(df1$H, " (", df1$CI95.H, ")", sep = "") 
S1 = paste(df1$S, " (", df1$CI95.S, ")", sep = "") 
G1 = paste(df1$G, " (", df1$CI95.G, ")", sep = "") 
Tm1 = paste(df1$Tm_at_0.1mM, " (", df1$CI95.Tm_at_0.1mM, ")", sep = "") 

df2 = df.MeltR %>% filter(Method == "2 Tm versus ln[Ct]")

H2 = paste(df2$H, " (", df2$CI95.H, ")", sep = "") 
S2 = paste(df2$S, " (", df2$CI95.S, ")", sep = "") 
G2 = paste(df2$G, " (", df2$CI95.G, ")", sep = "") 
Tm2 = paste(df2$Tm_at_0.1mM, " (", df2$CI95.Tm_at_0.1mM, ")", sep = "") 

df3 = df.MeltR %>% filter(Method == "3 Global fit")

H3 = paste(df3$H, " (", df3$CI95.H, ")", sep = "") 
S3 = paste(df3$S, " (", df3$CI95.S, ")", sep = "") 
G3 = paste(df3$G, " (", df3$CI95.G, ")", sep = "") 
Tm3 = paste(df3$Tm_at_0.1mM, " (", df3$CI95.Tm_at_0.1mM, ")", sep = "") 

df.Table = data.frame(Helix, H1, H2, H3, S1, S2, S3, G1, G2, G3, Tm1, Tm2, Tm3)

pH = c()
pS = c()
pG = c()
pTm = c()

helix = unique(df.MeltR$Helix)

for (i in 1:length(Helix)){
  df = df.MeltR %>% filter(Helix == helix[i])
  pH[i] = 100*abs(range(df %>% select(H))[1]-range(df %>% select(H))[2])/abs(mean((df %>% select(H))[[1]]))
  pS[i] = 100*abs(range(df %>% select(S))[1]-range(df %>% select(S))[2])/abs(mean((df %>% select(S))[[1]]))
  pG[i] = 100*abs(range(df %>% select(G))[1]-range(df %>% select(G))[2])/abs(mean((df %>% select(G))[[1]]))
  pTm[i] = 100*abs(range(df %>% select(Tm_at_0.1mM))[1]-range(df %>% select(Tm_at_0.1mM))[2])/abs(mean((df %>% select(Tm_at_0.1mM))[[1]]))
}

df.Table$pH = pH
df.Table$pS = pS
df.Table$pG = pG
df.Table$pTm = pTm

head(df.Table)

write.csv(df.Table, "Tables/SI_Table_X_BLtrimmer_results_comparing_MeltR_methods_1_2_3/SI_Table_X_BLtrimmer_results_comparing_MeltR_methods_1_2_3.csv", row.names = FALSE)

