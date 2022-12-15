library(ggbeeswarm)
library(tidyverse)
library(ggpubr)
devtools::load_all()
library(MeltR)

####Check error from NN model####

df.MeltR = read.csv("Tables/MeltR_fits/Fit_results.csv")

#?Helix.energy

Helix.energy("CGAAAGGU","ACCUUUCG")
Helix.energy("CGUUGC", "GCAACG")
Helix.energy("CUGAGUC", "GACUCAG")
Helix.energy("CGAAAGGU","ACCUUUCG", F.Q = TRUE)
Helix.energy("CGUUGC", "GCAACG", F.Q = TRUE)
Helix.energy("CUGAGUC", "GACUCAG", F.Q = TRUE)
Helix.energy("ACCGGU", "ACCGGU",
             AA.UU = -0.88,
             AU.AU = -1.13,
             UA.UA = -1.36,
             CU.AG = -2.03,
             CA.UG = -1.91,
             GU.AC = -2.36,
             GA.UC = -2.36,
             CG.CG = -2.19,
             GG.CC = -3.32,
             GC.GC = -3.44,
             Initiation = 4.63,
             Term.AU = 0.55,)
Helix.energy("ACCGGU", "ACCGGU",
             AA.UU = -6.3,
             AU.AU = -6.7,
             UA.UA = -12.7,
             CU.AG = -8.8,
             CA.UG = -13.0,
             GU.AC = -12.8,
             GA.UC = -14.5,
             CG.CG = -9.0,
             GG.CC = -14.6,
             GC.GC = -18.9,
             Initiation = 5.3,
             Term.AU = 2.8)

Helix.energy("CGAAAGGU", "ACCUUUCG",
             AA.UU = -6.82,
             AU.AU = -9.38,
             UA.UA = -7.69,
             CU.AG = -10.48,
             CA.UG = -10.44,
             GU.AC = -11.40,
             GA.UC = -12.44,
             CG.CG = -10.64,
             GG.CC = -13.39,
             GC.GC = -14.88,
             Initiation = 3.61,
             Term.AU = 3.72,
             FAMC.GBHQ1 = -13.07,
             F.Q = TRUE)

df.known = read.csv("Figures/Figure_3_Check_meltR.A/NN_numbers.csv")

colnames(df.MeltR)

df.MeltR$Known.dG = NA
df.MeltR$Known.dH = NA
df.MeltR$Known.dS = NA

for (i in 1:length(unique(df.MeltR$Helix))){
  helix = unique(df.MeltR$Helix)[i]
  df.MeltR$Known.dG[which(df.MeltR$Helix == helix)] = df.known$dG[which(df.known$Helix == helix)]
  df.MeltR$Known.dH[which(df.MeltR$Helix == helix)] = df.known$dH[which(df.known$Helix == helix)]
  df.MeltR$Known.dS[which(df.MeltR$Helix == helix)] = df.known$dS[which(df.known$Helix == helix)]
}

df.MeltR$error.dG = (df.MeltR$dG - df.MeltR$Known.dG)/df.MeltR$Known.dG
df.MeltR$error.dH = (df.MeltR$dH - df.MeltR$Known.dH)/df.MeltR$Known.dH
df.MeltR$error.dS = (df.MeltR$dS - 1000*df.MeltR$Known.dS)/(1000*df.MeltR$Known.dS)


df1 = df.MeltR %>%
  select(Helix, Method, error.dG, error.dH, error.dS)
df1$Program = "MeltR"

df.Meltwin$Known.dG = NA
df.Meltwin$Known.dH = NA
df.Meltwin$Known.dS = NA

colnames(df.Meltwin)[1] = "Helix"

df.Meltwin = df.Meltwin %>% filter(Helix != "UAUAUAUA")

for (i in 1:length(unique(df.Meltwin$Helix))){
  helix = unique(df.Meltwin$Helix)[i]
  print(i)
  df.Meltwin$Known.dG[which(df.Meltwin$Helix == helix)] = df.known$dG[which(df.known$Helix == helix)]
  df.Meltwin$Known.dH[which(df.Meltwin$Helix == helix)] = df.known$dH[which(df.known$Helix == helix)]
  df.Meltwin$Known.dS[which(df.Meltwin$Helix == helix)] = df.known$dS[which(df.known$Helix == helix)]
}

df.Meltwin$error.dG = (df.Meltwin$dG - df.Meltwin$Known.dG)/df.Meltwin$Known.dG
df.Meltwin$error.dH = (df.Meltwin$dH - df.Meltwin$Known.dH)/df.Meltwin$Known.dH
df.Meltwin$error.dS = (df.Meltwin$dS - 1000*df.Meltwin$Known.dS)/(1000*df.Meltwin$Known.dS)


df2 = df.Meltwin %>%
  select(Helix, Method, error.dG, error.dH, error.dS)
df2$Program = "MeltWin"

df = bind_rows(df1, df2)

df$Program = paste(df$Program, df$Method, sep = "\n")

unique(df$Program)

df$Program = factor(df$Program,
                    levels = c("MeltR\n1 individual fits", "MeltWin\n1 individual fits",
                               "MeltR\n2 Tm versus ln[Ct]", "MeltWin\n2 Tm versus ln[Ct]",
                               "MeltR\n3 Global fit"))
l.comp = list(c("MeltR\n1 individual fits", "MeltWin\n1 individual fits"),
              c("MeltR\n2 Tm versus ln[Ct]", "MeltWin\n2 Tm versus ln[Ct]"),
              c("MeltR\n3 Global fit", "MeltWin\n1 individual fits"),
              c("MeltR\n3 Global fit", "MeltWin\n2 Tm versus ln[Ct]"))

length(unique(df$Helix))

ggplot(df, aes(x = Program, y = 100*error.dG)) +
  geom_hline(yintercept =  0) +
  geom_line(mapping = aes(group = Helix)) + 
  stat_compare_means(label.x = 2, label.y = 10) +
  #stat_compare_means(comparisons = l.comp) +
  geom_beeswarm() +
  theme_classic() +
  #scale_y_continuous(breaks = c(-10, -5, 0, 5,  10), limits = c(-10, 10)) +
  theme(axis.text.x = element_text(color = "black",
                                    angle = 45,
                                   hjust = 1),
        axis.text.y = element_text(color = "black")) +
  xlab("") +
  ylab("%error \u0394G\u00B037 (kcal/mol)")

list.files("Figures")

ggsave("Figures/SI_Figure_4_BLtrimmer_modeled_data_Meltwin_MeltR_agrrement_with_NN/SI_Figure_X_Meltwin_MeltR_agrrement_with_NN.svg",
       scale = 2,
       height = 2,
       width = 4)
