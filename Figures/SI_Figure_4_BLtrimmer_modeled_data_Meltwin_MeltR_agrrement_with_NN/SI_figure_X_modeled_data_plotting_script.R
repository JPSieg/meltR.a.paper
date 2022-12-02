library(tidyverse)
library(MeltR)
library(cowplot)
devtools::load_all()

####Read in data####

list.files("Figures/")

df = read.csv("Figures/SI_Figure_4_BLtrimmer_modeled_data_Meltwin_MeltR_agrrement_with_NN/Modeled_fit_results.csv")

####Convert confidence intervals to nummeric####

Hmin = c()
Hmax = c()
Smin = c()
Smax = c()
Gmin = c()
Gmax = c()
Tmin = c()
Tmax = c()

for (i in 1:nrow(df)){
  Hmin[i] = as.numeric(strsplit(df$CI95.dH[i], split = " to ")[[1]][1])
  Hmax[i] = as.numeric(strsplit(df$CI95.dH[i], split = " to ")[[1]][2])
  Smin[i] = as.numeric(strsplit(df$CI95.dS[i], split = " to ")[[1]][1])
  Smax[i] = as.numeric(strsplit(df$CI95.dS[i], split = " to ")[[1]][2])
  Gmin[i] = as.numeric(strsplit(df$CI95.dG[i], split = " to ")[[1]][1])
  Gmax[i] = as.numeric(strsplit(df$CI95.dG[i], split = " to ")[[1]][2])
  Tmin[i] = as.numeric(strsplit(df$CI95.Tm_at_0.1mM[i], split = " to ")[[1]][1])
  Tmax[i] = as.numeric(strsplit(df$CI95.Tm_at_0.1mM[i], split = " to ")[[1]][2])
}

df$Hmin = Hmin
df$Hmax = Hmax
df$Smin = Smin
df$Smax = Smax
df$Gmin = Gmin
df$Gmax = Gmax
df$Tmin = Tmin
df$Tmax = Tmax

####Add in known parameters####

list.files("Figures/SI_Figure_4_BLtrimmer_modeled_data_Meltwin_MeltR_agrrement_with_NN")

df.known = read.csv("Figures/SI_Figure_4_BLtrimmer_modeled_data_Meltwin_MeltR_agrrement_with_NN/Modeled_Xia_energies.csv")

helices = unique(df$Helix)

df$Known.H = NA
df$Known.S = NA
df$Known.G = NA
df$known.Tm = NA

for (i in helices){
  df$Known.H[which(df$Helix == i)] = df.known$dH[which(df.known$seqF == i)]
  df$Known.S[which(df$Helix == i)] = df.known$dS[which(df.known$seqF == i)]
  df$Known.G[which(df$Helix == i)] = df.known$dG[which(df.known$seqF == i)]
  df$known.Tm[which(df$Helix == i)] = df.known$Tm.at.0.1mM[which(df.known$seqF == i)]
}

####dH plot####

colnames(df)

PA = ggplot(df, aes(x = Known.H,
                    y = dH, ymin = Hmin, ymax = Hmax,
                    label = Helix,
                    color = Method)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_pointrange() +
  theme_classic() +
  coord_fixed() +
  xlim(-85, -30) +
  ylim(-85, -30) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  ylab("MeltR calculated \u0394H\u00B0 (kcal/mol)") +
  xlab("Known \u0394H\u00B0 (kcal/mol)") +
  scale_color_manual(values = viridis::viridis(3, end = 0.8),
                     name = "Meltwin\nmethod")

PA  

####dS plot####

PB = ggplot(df, aes(x = 1000*Known.S,
                    y = dS, ymin = Smin, ymax = Smax,
                    label = Helix,
                    color = Method)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_pointrange() +
  theme_classic() +
  coord_fixed() +
  xlim(-230, -80) +
  ylim(-230, -80)+
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  ylab("MeltR calculated \u0394S\u00B0 (cal/mol/K)") +
  xlab("Known \u0394S\u00B0 (cal/mol/K)")+
  scale_color_manual(values = viridis::viridis(3, end = 0.8),
                     name = "Meltwin\nmethod")

PB

PC = ggplot(df, aes(x = Known.G,
                    y = dG, ymin = Gmin, ymax = Gmax,
                    label = Helix,
                    color = Method)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_pointrange() +
  theme_classic() +
  coord_fixed() +
  #xlim(-11, -5) +
  #ylim(-11, -5)+
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  ylab("MeltR calculated \u0394G\u00B037 (kcal/mol)") +
  xlab("Known \u0394G\u00B037 (kcal/mol)")+
  scale_color_manual(values = viridis::viridis(3, end = 0.8),
                     name = "Meltwin\nmethod")

PC

PD = ggplot(df, aes(x = known.Tm,
                    y = Tm_at_0.1mM, ymin = Tmin, ymax = Tmax,
                    label = Helix,
                    color = Method)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_pointrange() +
  theme_classic() +
  coord_fixed() +
  #xlim(30, 70) +
  #ylim(30, 70)+
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.75, 0.25))  +
  xlab("Known Tm at 0.1 mM (\u00B0C)") +
  ylab("MeltR calculated Tm at 0.1 mM (\u00B0C)")+
  scale_color_manual(values = viridis::viridis(3, end = 0.8),
                     name = "MeltR\nmethod")
PD

P = plot_grid(PA, PB, PC, PD,
              labels = c("A", "B", "C", "D"), label_size = 16)  

list.files("Figures/")

ggsave("Figures/SI_Figure_4_BLtrimmer_modeled_data_Meltwin_MeltR_agrrement_with_NN/SI_Figure_X_BLtrimmer_modeled_data.svg", P,
       width = 5.3, height = 5.3, units = "in", scale = 1.5, bg = "white")


####Rewrite results file####

write.csv(df, "Figures/SI_Figure_4_BLtrimmer_modeled_data_Meltwin_MeltR_agrrement_with_NN/Modeled_fit_results.csv", row.names = FALSE)

df$error = abs(df$dH -df$Known.H)

#View(df)

  