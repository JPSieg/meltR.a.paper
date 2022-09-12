library(tidyverse)
library(ggrepel)
library(cowplot)

df.Meltwin

list.files("Tables/SI_Table_X_MeltR_fits")

df.MeltR = read.csv("Tables/SI_Table_X_MeltR_fits/Fit_results.csv") %>% filter(Method == "3 Global fit")

head(df.MeltR)

Hmin = c()
Hmax = c()
Smin = c()
Smax = c()
Gmin = c()
Gmax = c()
Tmin = c()
Tmax = c()

for (i in 1:nrow(df.MeltR)){
  Hmin[i] = as.numeric(strsplit(df.MeltR$CI95.H[i], split = " to ")[[1]][1])
  Hmax[i] = as.numeric(strsplit(df.MeltR$CI95.H[i], split = " to ")[[1]][2])
  Smin[i] = as.numeric(strsplit(df.MeltR$CI95.S[i], split = " to ")[[1]][1])
  Smax[i] = as.numeric(strsplit(df.MeltR$CI95.S[i], split = " to ")[[1]][2])
  Gmin[i] = as.numeric(strsplit(df.MeltR$CI95.G[i], split = " to ")[[1]][1])
  Gmax[i] = as.numeric(strsplit(df.MeltR$CI95.G[i], split = " to ")[[1]][2])
  Tmin[i] = as.numeric(strsplit(df.MeltR$CI95.Tm_at_0.1mM[i], split = " to ")[[1]][1])
  Tmax[i] = as.numeric(strsplit(df.MeltR$CI95.Tm_at_0.1mM[i], split = " to ")[[1]][2])
}

df.MeltR$Hmin = Hmin
df.MeltR$Hmax = Hmax
df.MeltR$Smin = Smin
df.MeltR$Smax = Smax
df.MeltR$Gmin = Gmin
df.MeltR$Gmax = Gmax
df.MeltR$Tmin = Tmin
df.MeltR$Tmax = Tmax

df.M1 = cbind(df.MeltR %>%
                arrange(Helix) %>%
                select(Helix, H, S, G, Hmin, Hmax, Smin, Smax, Gmin, Gmax, Tm_at_0.1mM, Tmin, Tmax),
              df.Meltwin %>%
                filter(Sequence != "UAUAUAUA") %>%
                filter(Method == "1 individual fits") %>%
                arrange(Sequence) %>%
                select(dH, SE.dH, dS, SE.dS, dG, SE.dG, Tm))

df.M1$Method = "1 individual fits"

df.M2 = cbind(df.MeltR %>%
                arrange(Helix) %>%
                select(Helix, H, S, G, Hmin, Hmax, Smin, Smax, Gmin, Gmax, Tm_at_0.1mM, Tmin, Tmax),
              df.Meltwin %>%
                filter(Sequence != "UAUAUAUA") %>%
                filter(Method == "2 Tm versus ln[Ct]") %>%
                arrange(Sequence) %>%
                select(dH, SE.dH, dS, SE.dS, dG, SE.dG, Tm))

df.M2$Method = "2 Tm versus ln[Ct]"

df = bind_rows(df.M1, df.M2)

####dH plot####

df$Labels = df$Helix

df$Labels[which(df$Method == "2 Tm versus ln[Ct]")] = NA

View(df)

PA = ggplot(df, aes(x = dH, xmin = dH - SE.dH, xmax = dH + SE.dH,
               y = H, ymin = Hmin, ymax = Hmax,
               label = Helix,
               color = Method)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_pointrange() +
  geom_errorbarh() +
  #geom_text_repel(size = 2) +
  theme_classic() +
  coord_fixed() +
  xlim(-85, -45) +
  ylim(-85, -45) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  ylab("MeltR Global fit \u0394H\u00B0 (kcal/mol)") +
  xlab("Meltwin \u0394H\u00B0 (kcal/mol)") +
  scale_color_manual(values = viridis::viridis(2, end = 0.8),
                     name = "Meltwin\nmethod")
  
  
PB = ggplot(df, aes(x = dS, xmin = dS - SE.dS, xmax = dS + SE.dS,
               y = S, ymin = Smin, ymax = Smax,
               label = Helix,
               color = Method)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_pointrange() +
  geom_errorbarh() +
  #geom_text_repel(size = 2) +
  theme_classic() +
  coord_fixed() +
  xlim(-250, -110) +
  ylim(-250, -110)+
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  ylab("MeltR Global fit \u0394S\u00B0 (cal/mol/K)") +
  xlab("Meltwin \u0394S\u00B0 (cal/mol/K)")+
  scale_color_manual(values = viridis::viridis(2, end = 0.8),
                     name = "Meltwin\nmethod")

PC = ggplot(df, aes(x = dG, xmin = dG - SE.dG, xmax = dG + SE.dG,
               y = G, ymin = Gmin, ymax = Gmax,
               label = Helix,
               color = Method)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_pointrange() +
  geom_errorbarh() +
  theme_classic() +
  coord_fixed() +
  xlim(-11, -5) +
  ylim(-11, -5)+
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  ylab("MeltR Global fit \u0394G\u00B037 (kcal/mol)") +
  xlab("Meltwin \u0394G\u00B037 (kcal/mol)")+
  scale_color_manual(values = viridis::viridis(2, end = 0.8),
                     name = "Meltwin\nmethod")

PD = ggplot(df, aes(x = Tm,
               y = Tm_at_0.1mM, ymin = Tmin, ymax = Tmax,
               label = Helix,
               color = Method)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_pointrange() +
  theme_classic() +
  coord_fixed() +
  xlim(30, 70) +
  ylim(30, 70)+
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.75, 0.15))  +
  xlab("Meltwin Tm at 0.1 mM (\u00B0C)") +
  ylab("MeltR Global fit Tm at 0.1 mM (\u00B0C)")+
  scale_color_manual(values = viridis::viridis(2, end = 0.8),
                     name = "Meltwin\nmethod")

P = plot_grid(PA, PB, PC, PD)  

list.files("Figures/Figure_4_Check_meltR.A/")

ggsave("Figures/Figure_4_Check_meltR.A/Figure_4_Meltwin_comparison.svg", P,
       width = 5.3, height = 5.3, units = "in", scale = 1.5, bg = "white")
