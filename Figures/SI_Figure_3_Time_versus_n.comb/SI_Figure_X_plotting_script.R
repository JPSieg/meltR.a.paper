library(tidyverse)
library(cowplot)
library(viridis)
library(MeltR)
library(ggpubr)
devtools::load_all()

####Read in data####

list.files("Figures/SI_Figure_3_Time_versus_n.comb/Fit_results")

list.df = lapply(paste("Figures/SI_Figure_3_Time_versus_n.comb/Fit_results",
                       list.files("Figures/SI_Figure_3_Time_versus_n.comb/Fit_results"),
                       sep = "/"), read.csv)

df = bind_rows(list.df)

####Make CIs characters####

df$CI95.dH = as.character(df$CI95.dH)
df$CI95.dS = as.character(df$CI95.dS)
df$CI95.dG = as.character(df$CI95.dG)
df$CI95.Tm_at_0.1mM = as.character(df$CI95.Tm_at_0.1mM)



####Pull out error bars####

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

####Plot Time v n.comb####

head(df)

A = ggplot(df %>% filter(Method == "3 Global fit"),
       aes(x = n.comb, y = S.time/60)) +
  geom_smooth(method = "lm") +
  stat_regline_equation(aes(label = ..eq.label..)) +
  geom_point() +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  xlab("Number of combinations") +
  ylab("System time (min)") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10")

####parameter convergence####

B = ggplot(df, aes(x = n.comb, y = dH, color = Method, shape = Method,
               ymin = Hmin, ymax = Hmax)) +
  geom_vline(xintercept = 1000) +
  geom_pointrange() +
  geom_line() +
  theme_classic() +
  scale_x_continuous(trans = "log10") +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.8, 0.15)) +
  scale_color_manual(values = viridis(5)[2:4]) +
  xlab("Number of combinations") +
  ylab("\u0394H\u00B0 (kcal/mol)")

C = ggplot(df, aes(x = n.comb, y = dS, color = Method, shape = Method,
               ymin = Smin, ymax = Smax)) +
  geom_vline(xintercept = 1000) +
  geom_pointrange() +
  geom_line() +
  theme_classic() +
  scale_x_continuous(trans = "log10") +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_color_manual(values = viridis(5)[2:4]) +
  xlab("Number of combinations") +
  ylab("\u0394S\u00B0 (cal/mol/K)")

D = ggplot(df, aes(x = n.comb, y = dG, color = Method, shape = Method,
               ymin = Gmin, ymax = Gmax)) +
  geom_vline(xintercept = 1000) +
  geom_pointrange() +
  geom_line() +
  theme_classic() +
  scale_x_continuous(trans = "log10") +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_color_manual(values = viridis(5)[2:4]) +
  xlab("Number of combinations") +
  ylab("\u0394G\u00B037 (kcal/mol)")

E = ggplot(df, aes(x = n.comb, y = Tm_at_0.1mM, color = Method, shape = Method,
               ymin = Tmin, ymax = Tmax)) +
  geom_vline(xintercept = 1000) +
  geom_pointrange() +
  geom_line() +
  theme_classic() +
  scale_x_continuous(trans = "log10") +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  scale_color_manual(values = viridis(5)[2:4]) +
  xlab("Number of combinations") +
  ylab("Tm at 0.1 mM (\u00B0C)")

####Make final plot####

P = plot_grid(A,
          plot_grid(B, C, D, E,
                    labels = c("B", "C", "D", "E"),
                    label_size = 16),
          ncol = 1, rel_heights = c(1, 2.5),
          labels = c("A", ""),
          label_size = 16)

ggsave("Figures/SI_Figure_3_Time_versus_n.comb/SI_Figure_X_Time_versus_ncomb.svg", P,
       width = 4, height = 5, units = "in", scale = 2.5, bg = "white")
