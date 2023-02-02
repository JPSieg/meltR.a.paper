library(ggbeeswarm)
library(tidyverse)
library(ggpubr)
devtools::load_all()
library(MeltR)

####Check error from NN model####

df.MeltR = read.csv("Tables/MeltR_fits/Fit_results.csv")


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

df.Meltwin = df.Meltwin

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
                    levels = c("MeltR\n1 Individual fits", "MeltWin\n1 individual fits",
                               "MeltR\n3 Global fit",
                               "MeltR\n2 Tm versus ln[Ct]", "MeltWin\n2 Tm versus ln[Ct]"))
l.comp = list(c("MeltR\n1 Individual fits", "MeltWin\n1 individual fits"),
              c("MeltR\n2 Tm versus ln[Ct]", "MeltWin\n2 Tm versus ln[Ct]"),
              c("MeltR\n3 Global fit", "MeltWin\n1 individual fits"),
              c("MeltR\n3 Global fit", "MeltWin\n2 Tm versus ln[Ct]"))

length(unique(df$Helix))

Mmodel = c()

for (i in 1:nrow(df)){
  Mmodel[i] = df.known$Mmodel[which(df.known$Helix == df$Helix[i])]
}

df$Mmodel = Mmodel

View(df)

ggplot(df, aes(x = Program, y = 100*error.dG)) +
  geom_hline(yintercept =  0) +
  geom_line(mapping = aes(group = Helix, color = Mmodel)) + 
  stat_compare_means(label.x = 2, label.y = 10) +
  scale_color_manual(values = viridis::viridis(4)) +
  #scale_y_continuous(limits = c(-10, 10), breaks = c(-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10)) +
  #stat_compare_means(comparisons = l.comp) +
  geom_beeswarm(mapping = aes(color = Mmodel, shape = Mmodel), size = 3) +
  theme_classic() +
  scale_y_continuous(breaks = c(-50, -25, -10, -5, 0, 5,  10), limits = c(-50, 10)) +
  theme(axis.text.x = element_text(color = "black",
                                    angle = 45,
                                   hjust = 1),
        axis.text.y = element_text(color = "black")) +
  xlab("") +
  ylab("%error \u0394G\u00B037")

list.files("Figures")

ggsave("Figures/SI_Figure_4_BLtrimmer_modeled_data_Meltwin_MeltR_agrrement_with_NN/SI_Figure_X_Meltwin_MeltR_agrrement_with_NN.svg",
       scale = 2.5,
       height = 1.8,
       width = 4)

write.csv(df, "Figures/SI_Figure_4_BLtrimmer_modeled_data_Meltwin_MeltR_agrrement_with_NN/NN_Agreement.csv", row.names = F)
