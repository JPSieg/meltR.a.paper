library(tidyverse)
library(ggrepel)
library(cowplot)
devtools::load_all()

df.Meltwin

File = "Tables/MeltR_fits/Fit_results.csv"

list.files("Tables/MeltR_fits/Fit_results")

df.MeltR = read.csv(File) %>% filter(Method == "3 Global fit")

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
  Hmin[i] = as.numeric(strsplit(df.MeltR$CI95.dH[i], split = " to ")[[1]][1])
  Hmax[i] = as.numeric(strsplit(df.MeltR$CI95.dH[i], split = " to ")[[1]][2])
  Smin[i] = as.numeric(strsplit(df.MeltR$CI95.dS[i], split = " to ")[[1]][1])
  Smax[i] = as.numeric(strsplit(df.MeltR$CI95.dS[i], split = " to ")[[1]][2])
  Gmin[i] = as.numeric(strsplit(df.MeltR$CI95.dG[i], split = " to ")[[1]][1])
  Gmax[i] = as.numeric(strsplit(df.MeltR$CI95.dG[i], split = " to ")[[1]][2])
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

colnames(df.MeltR)[c(2, 4, 6)] = c("H", "S", "G")

df.M1 = cbind(df.MeltR %>%
                arrange(Helix) %>%
                select(Helix, H, S, G, Hmin, Hmax, Smin, Smax, Gmin, Gmax, Tm_at_0.1mM, Tmin, Tmax),
              df.Meltwin %>%
                filter(Method == "1 individual fits") %>%
                arrange(Sequence) %>%
                select(dH, SE.dH, dS, SE.dS, dG, SE.dG, Tm))

df.M1$Method = "1 Individual fits"

df.M2 = cbind(df.MeltR %>%
                arrange(Helix) %>%
                select(Helix, H, S, G, Hmin, Hmax, Smin, Smax, Gmin, Gmax, Tm_at_0.1mM, Tmin, Tmax),
              df.Meltwin %>%
                filter(Method == "2 Tm versus ln[Ct]") %>%
                arrange(Sequence) %>%
                select(dH, SE.dH, dS, SE.dS, dG, SE.dG, Tm))

df.M2$Method = "2 Tm versus ln[Ct]"

df = bind_rows(df.M1, df.M2)

####dH plot####

df$Labels = df$Helix

df$Labels[which(df$Method == "2 Tm versus ln[Ct]")] = NA

colnames(df)

PA = ggplot(df, aes(y = dH, ymin = dH - SE.dH, ymax = dH + SE.dH,
               x = H, xmin = Hmin, xmax = Hmax,
               label = Helix,
               color = Method)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_pointrange() +
  geom_errorbarh() +
  #geom_text_repel(size = 2) +
  theme_classic() +
  coord_fixed() +
  xlim(-85, -20) +
  ylim(-85, -20) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.75, 0.2)) +
  xlab("MeltR Global fit \u0394H\u00B0 (kcal/mol)") +
  ylab("MeltWin \u0394H\u00B0 (kcal/mol)") +
  scale_color_manual(values = viridis::viridis(2, end = 0.8),
                     name = "MeltWin\nmethod")

PB = ggplot(df, aes(y = dS, ymin = dS - SE.dS, ymax = dS + SE.dS,
               x = S, xmin = Smin, xmax = Smax,
               label = Helix,
               color = Method)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_pointrange() +
  geom_errorbarh() +
  #geom_text_repel(size = 2) +
  theme_classic() +
  coord_fixed() +
  xlim(-250, -75) +
  ylim(-250, -75)+
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  xlab("MeltR Global fit \u0394S\u00B0 (cal/mol*K)") +
  ylab("MeltWin \u0394S\u00B0 (cal/mol/K)")+
  scale_color_manual(values = viridis::viridis(2, end = 0.8),
                     name = "MeltWin\nmethod")

PC = ggplot(df, aes(x = dG, xmin = dG - SE.dG, xmax = dG + SE.dG,
               y = G, ymin = Gmin, ymax = Gmax,
               label = Helix,
               color = Method)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_pointrange() +
  geom_errorbarh() +
  theme_classic() +
  coord_fixed() +
  xlim(-11, 0) +
  ylim(-11, 0)+
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  xlab("MeltR Global fit \u0394G\u00B037 (kcal/mol)") +
  ylab("MeltWin \u0394G\u00B037 (kcal/mol)")+
  scale_color_manual(values = viridis::viridis(2, end = 0.8),
                     name = "MeltWin\nmethod")

PD = ggplot(df, aes(y = Tm,
               x = Tm_at_0.1mM, xmin = Tmin, xmax = Tmax,
               label = Helix,
               color = Method)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_pointrange() +
  theme_classic() +
  coord_fixed() +
  xlim(30, 80) +
  ylim(30, 80)+
  theme(axis.text = element_text(color = "black"),
        legend.position = "none")  +
  ylab("MeltWin Tm at 0.1 mM (\u00B0C)") +
  xlab("MeltR Global fit Tm at 0.1 mM (\u00B0C)")+
  scale_color_manual(values = viridis::viridis(2, end = 0.8),
                     name = "MeltWin\nmethod")

####Percent error heat map Bimolecular####

df.MeltR = read.csv(File) 

Comparisons = list(c(1, 1),
                   c(1, 2),
                   c(2, 1),
                   c(2, 2),
                   c(3, 1),
                   c(3, 2))

df.MeltR$MModel = NA
df.Meltwin$MModel = NA
df.MeltR$MModel[-which(df.MeltR$Helix %in% c("GGCUAACGGCC", "GCCGUGAGGC", "GCCUUCGGGC", "GCGGCAACGC", "GCUGAAAGGC", "GCUGAGAGGC", "GCCUAACGGC", "GCCUUUUAGGC", "GGCCACAGGCC", "GCCGUUUCGGC", "GGCGAGACGCC","GGCAAAAUGCC"))] = "Bimolecular"
df.MeltR$MModel[which(df.MeltR$Helix %in% c("GGCUAACGGCC", "GCCGUGAGGC", "GCCUUCGGGC", "GCGGCAACGC", "GCUGAAAGGC", "GCUGAGAGGC", "GCCUAACGGC", "GCCUUUUAGGC", "GGCCACAGGCC", "GCCGUUUCGGC", "GGCGAGACGCC","GGCAAAAUGCC"))] = "Monomolecular"
df.Meltwin$MModel[-which(df.Meltwin$Sequence %in% c("GGCUAACGGCC", "GCCGUGAGGC", "GCCUUCGGGC", "GCGGCAACGC", "GCUGAAAGGC", "GCUGAGAGGC", "GCCUAACGGC", "GCCUUUUAGGC", "GGCCACAGGCC", "GCCGUUUCGGC", "GGCGAGACGCC","GGCAAAAUGCC"))] = "Bimolecular"
df.Meltwin$MModel[which(df.Meltwin$Sequence %in% c("GGCUAACGGCC", "GCCGUGAGGC", "GCCUUCGGGC", "GCGGCAACGC", "GCUGAAAGGC", "GCUGAGAGGC", "GCCUAACGGC", "GCCUUUUAGGC", "GGCCACAGGCC", "GCCGUUUCGGC", "GGCGAGACGCC","GGCAAAAUGCC"))] = "Monomolecular"

list.df.comparison = {}

for (i in 1:length(Comparisons)){
  colnames(df.MeltR)
  unique(df.MeltR$Method)
  
  if (Comparisons[[i]][1] == 1){method = "1 Individual fits"}
  if (Comparisons[[i]][1] == 2){method = "2 Tm versus ln[Ct]"}
  if (Comparisons[[i]][1] == 3){method = "3 Global fit"}
  
  df.method.MeltR = df.MeltR %>% filter(Method == method) %>% filter(MModel == "Bimolecular") %>% arrange(Helix) 
  
  colnames(df.Meltwin)
  unique(df.Meltwin$Method)
  
  if (Comparisons[[i]][2] == 1){method = "1 individual fits"}
  if (Comparisons[[i]][2] == 2){method = "2 Tm versus ln[Ct]"}
  
  df.method.Meltwin = df.Meltwin %>% filter(Method == method) %>% filter(MModel == "Bimolecular")  %>% arrange(Sequence)
  
  H = mean(100*abs(df.method.MeltR$dH - df.method.Meltwin$dH)/abs((df.method.MeltR$dH + df.method.Meltwin$dH)/2), na.rm = T)
  S = mean(100*abs(df.method.MeltR$dS - df.method.Meltwin$dS)/abs((df.method.MeltR$dS + df.method.Meltwin$dS)/2), na.rm = T)
  G = mean(100*abs(df.method.MeltR$dG - df.method.Meltwin$dG)/abs((df.method.MeltR$dG + df.method.Meltwin$dG)/2), na.rm = T)
  Tm = mean(100*abs(df.method.MeltR$Tm_at_0.1mM - df.method.Meltwin$Tm)/abs((df.method.MeltR$Tm_at_0.1mM + df.method.Meltwin$Tm)/2), na.rm = T)
  
  df.H = data.frame("Parameter" =  "dH",
                    "MeltR.method" = Comparisons[[i]][1],
                    "MeltWin.method" = Comparisons[[i]][2],
                    "Error" = H)
  df.S = data.frame("Parameter" =  "dS",
                    "MeltR.method" = Comparisons[[i]][1],
                    "MeltWin.method" = Comparisons[[i]][2],
                    "Error" = S)
  df.G = data.frame("Parameter" =  "dG",
                    "MeltR.method" = Comparisons[[i]][1],
                    "MeltWin.method" = Comparisons[[i]][2],
                    "Error" = G)
  df.Tm = data.frame("Parameter" =  "Tm",
                    "MeltR.method" = Comparisons[[i]][1],
                    "MeltWin.method" = Comparisons[[i]][2],
                    "Error" = Tm)
  
  list.df.comparison[[i]] = bind_rows(df.H, df.S, df.G, df.Tm)
}

df.compare = bind_rows(list.df.comparison)


list.files("Figures/Figure_3_Check_meltR.A/")

write.csv(df.compare, "Figures/Figure_3_Check_meltR.A/MeltR_MeltWin_percent_error_by_method_bimolecular.csv")

head(df.compare)

#View(df.compare)

df.compare$Parameter = factor(df.compare$Parameter,
                              levels = c("dH", "dS", "dG", "Tm"),
                              labels = c("\u0394H\u00B0", "\u0394S\u00B0", "\u0394G\u00B037", "Tm"))

df.compare$Label = paste(round(df.compare$Error, digits = 1), "%", sep = "")


P1 = ggplot(df.compare, aes(x = MeltR.method, y = MeltWin.method,
                       fill = Error, label = Label)) +
  geom_tile() +
  geom_text() +
  facet_wrap(~Parameter) +
  scale_y_continuous(breaks = c(1, 2)) +
  scale_fill_viridis_c(minor_breaks = 10, limits = c(0, 5.1),
                       name = "%error") +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  xlab("MeltR method") +
  ylab("MeltWin method")

####Percent error heat map Monomolecular####

list.df.comparison = {}

for (i in 1:length(Comparisons)){
  colnames(df.MeltR)
  unique(df.MeltR$Method)
  
  if (Comparisons[[i]][1] == 1){method = "1 Individual fits"}
  if (Comparisons[[i]][1] == 2){method = "2 Tm versus ln[Ct]"}
  if (Comparisons[[i]][1] == 3){method = "3 Global fit"}
  
  df.method.MeltR = df.MeltR %>% filter(Method == method) %>% filter(!MModel == "Bimolecular") %>% arrange(Helix) 
  
  colnames(df.Meltwin)
  unique(df.Meltwin$Method)
  
  if (Comparisons[[i]][2] == 1){method = "1 individual fits"}
  if (Comparisons[[i]][2] == 2){method = "2 Tm versus ln[Ct]"}
  
  df.method.Meltwin = df.Meltwin %>% filter(Method == method) %>% filter(!MModel == "Bimolecular")  %>% arrange(Sequence)
  
  H = mean(100*abs(df.method.MeltR$dH - df.method.Meltwin$dH)/abs((df.method.MeltR$dH + df.method.Meltwin$dH)/2), na.rm = T)
  S = mean(100*abs(df.method.MeltR$dS - df.method.Meltwin$dS)/abs((df.method.MeltR$dS + df.method.Meltwin$dS)/2), na.rm = T)
  G = mean(100*abs(df.method.MeltR$dG - df.method.Meltwin$dG)/abs((df.method.MeltR$dG + df.method.Meltwin$dG)/2), na.rm = T)
  Tm = mean(100*abs(df.method.MeltR$Tm_at_0.1mM - df.method.Meltwin$Tm)/abs((df.method.MeltR$Tm_at_0.1mM + df.method.Meltwin$Tm)/2), na.rm = T)
  
  df.H = data.frame("Parameter" =  "dH",
                    "MeltR.method" = Comparisons[[i]][1],
                    "MeltWin.method" = Comparisons[[i]][2],
                    "Error" = H)
  df.S = data.frame("Parameter" =  "dS",
                    "MeltR.method" = Comparisons[[i]][1],
                    "MeltWin.method" = Comparisons[[i]][2],
                    "Error" = S)
  df.G = data.frame("Parameter" =  "dG",
                    "MeltR.method" = Comparisons[[i]][1],
                    "MeltWin.method" = Comparisons[[i]][2],
                    "Error" = G)
  df.Tm = data.frame("Parameter" =  "Tm",
                     "MeltR.method" = Comparisons[[i]][1],
                     "MeltWin.method" = Comparisons[[i]][2],
                     "Error" = Tm)
  
  list.df.comparison[[i]] = bind_rows(df.H, df.S, df.G, df.Tm)
}

df.compare = bind_rows(list.df.comparison)


list.files("Figures/Figure_3_Check_meltR.A/")

write.csv(df.compare, "Figures/Figure_3_Check_meltR.A/MeltR_MeltWin_percent_error_by_method_monomolecular.csv")

head(df.compare)

#View(df.compare)

df.compare$Parameter = factor(df.compare$Parameter,
                              levels = c("dH", "dS", "dG", "Tm"),
                              labels = c("\u0394H\u00B0", "\u0394S\u00B0", "\u0394G\u00B037", "Tm"))

df.compare$Label = paste(round(df.compare$Error, digits = 1), "%", sep = "")


P2 = ggplot(df.compare, aes(x = MeltR.method, y = MeltWin.method,
                            fill = Error, label = Label)) +
  geom_tile() +
  geom_text() +
  facet_wrap(~Parameter) +
  scale_y_continuous(breaks = c(1, 2)) +
  scale_fill_viridis_c(minor_breaks = 10, limits = c(0, 5.1),
                       name = "%error") +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  xlab("MeltR method") +
  ylab("MeltWin method")

####Consolidate plots####

P.bottom = plot_grid(PA, PB, PC, PD,
              labels = c("B", "C", "D", "F"), label_size = 16)  

P.top = plot_grid(P1, P2)


P = plot_grid(P.top, P.bottom, ncol = 1,
              labels = c("A", ""), label_size = 16,
              rel_heights = c(1, 3))

list.files("Figures/Figure_3_Check_meltR.A/")

ggsave("Figures/Figure_3_Check_meltR.A/Figure_3_Meltwin_comparison.svg", P,
       width = 5.5, height = 6, units = "in", scale = 1.5, bg = "white")

