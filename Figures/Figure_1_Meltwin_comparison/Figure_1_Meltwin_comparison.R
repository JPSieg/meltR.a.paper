library(tidyverse)
library(ggpubr)
library(cowplot)

####Data for A B and C####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "CROWD DP1")

head(df)

ggplot(df %>% filter(Sample == 12),
       aes(x = Temperature,
           y = Absorbance)) +
  geom_point() +
  geom_vline(xintercept = c(10,75))

unique(df$Sample)

T.range = list(c(27,83), #2
               c(27,75), #3
               c(15,75), #4
               c(25,75), #6
               c(25,75), #7
               c(20,75), #8
               c(15,85), #10
               c(20,70), #11
               c(15,70)) #12 

fit = meltR.A(df,
              NucAcid = c("RNA", "CGCGCG"),
              Mmodel = "Homoduplex.2State",
              blank = "none",
              concT = 10,
              fitTs = T.range,
              wavelength = 280,
              Save_results = "all",
              file_path = "Figures/Table_1_Meltwin_comparison/CROWD DP1")

####A####

ggplot(fit$Method.1.data, aes(x = Temperature)) +
  geom_point(mapping = aes(y = Absorbance, color = Ct))



####B####

####C####

####D####

df.D = cbind(df.MeltR.indv.fits %>% arrange(Experiment) %>% select(Experiment, H, S, G, Tm),
df.Meltwin.indv.fits %>% arrange(Experiment) %>% select(dH.kcal.mol, dS.kcal.molK, dG.kcal.mol, Tm.degC))

P1 = ggplot(df.D, aes(x = H, y = dH.kcal.mol, color = Experiment)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.8) +
  stat_regline_equation(aes(label = ..rr.label.., color = NULL)) +
  theme_classic() +
  ggtitle("dH (kcal/mol)") +
  xlab("MeltR") +
  ylab("Meltwin") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

P2 = ggplot(df.D, aes(x = 1000*S, y = dS.kcal.molK, color = Experiment)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.8) +
  stat_regline_equation(aes(label = ..rr.label.., color = NULL)) +
  theme_classic() +
  ggtitle("dS (cal/mol/K)") +
  xlab("MeltR") +
  ylab("Meltwin") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
 
P3 = ggplot(df.D, aes(x = G, y = dG.kcal.mol, color = Experiment)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.8) +
  stat_regline_equation(aes(label = ..rr.label.., color = NULL)) +
  theme_classic() +
  ggtitle("dG (kcal/mol)") +
  xlab("MeltR") +
  ylab("Meltwin") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

P4 = ggplot(df.D, aes(x = Tm, y = Tm.degC, color = Experiment)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.8) +
  stat_regline_equation(aes(label = ..rr.label.., color = NULL)) +
  theme_classic() +
  ggtitle("Tm (degC)") +
  xlab("MeltR") +
  ylab("Meltwin") +
  theme(legend.position = c(0.8, 0.35),
        legend.background = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

PD = plot_grid(P1, P2, P3, P4)

PD

####E####

df.E = cbind(df.MeltR.summary %>% arrange(Experiment) %>% filter(Method != "3 Global fit") %>% select(Experiment, Method, H, S, G, Tm_at_0.1mM),
             df.Meltwin.summary %>% arrange(Experiment) %>% select(dH.kcal.mol, dS.kcal.molK, dG.kcal.mol, Tm.degC))
df.E = df.E %>% filter(Method != "2 Tm versus ln[Ct]")

P1 = ggplot(df.E, aes(x = H, y = dH.kcal.mol, color = Experiment)) +
  stat_regline_equation(aes(label = ..rr.label.., color = NULL)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  ggtitle("dH (kcal/mol)") +
  xlab("MeltR") +
  ylab("Meltwin") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

P2 = ggplot(df.E, aes(x = S, y = dS.kcal.molK, color = Experiment)) +
  stat_regline_equation(aes(label = ..rr.label.., color = NULL)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  ggtitle("dS (cal/mol/K)") +
  xlab("MeltR") +
  ylab("Meltwin") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

P3 = ggplot(df.E, aes(x = G, y = dG.kcal.mol, color = Experiment)) +
  stat_regline_equation(aes(label = ..rr.label.., color = NULL)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  ggtitle("dG (kcal/mol)") +
  xlab("MeltR") +
  ylab("Meltwin") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

P4 = ggplot(df.E, aes(x = Tm_at_0.1mM, y = Tm.degC, color = Experiment)) +
  stat_regline_equation(aes(label = ..rr.label.., color = NULL)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  ggtitle("Tm (degC)") +
  xlab("MeltR") +
  ylab("Meltwin") +
  theme(legend.position = c(0.8, 0.35),
        legend.background = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

PE = plot_grid(P1, P2, P3, P4)

####F####


df.F = cbind(df.MeltR.summary %>% arrange(Experiment) %>% filter(Method != "3 Global fit") %>% select(Experiment, Method, H, S, G, Tm_at_0.1mM),
             df.Meltwin.summary %>% arrange(Experiment) %>% select(dH.kcal.mol, dS.kcal.molK, dG.kcal.mol, Tm.degC))
df.F = df.F %>% filter(Method == "2 Tm versus ln[Ct]")

P1 = ggplot(df.F, aes(x = H, y = dH.kcal.mol, color = Experiment)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.8) +
  stat_regline_equation(aes(label = ..rr.label.., color = NULL)) +
  theme_classic() +
  ggtitle("dH (kcal/mol)") +
  xlab("MeltR") +
  ylab("Meltwin") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

P2 = ggplot(df.F, aes(x = S, y = dS.kcal.molK, color = Experiment)) +
  geom_abline(intercept = 0, slope = 1) +
  stat_regline_equation(aes(label = ..rr.label.., color = NULL)) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  ggtitle("dS (cal/mol/K)") +
  xlab("MeltR") +
  ylab("Meltwin") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

P3 = ggplot(df.F, aes(x = G, y = dG.kcal.mol, color = Experiment)) +
  stat_regline_equation(aes(label = ..rr.label.., color = NULL)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  ggtitle("dG (kcal/mol)") +
  xlab("MeltR") +
  ylab("Meltwin") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

P4 = ggplot(df.F, aes(x = Tm_at_0.1mM, y = Tm.degC, color = Experiment)) +
  stat_regline_equation(aes(label = ..rr.label.., color = NULL)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  ggtitle("Tm (degC)") +
  xlab("MeltR") +
  ylab("Meltwin") +
  theme(legend.position = c(0.8, 0.35),
        legend.background = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

PF = plot_grid(P1, P2, P3, P4)

####G####

df.G1 = cbind(df.MeltR.summary %>% arrange(Experiment) %>% filter(Method != "3 Global fit") %>% select(Experiment, Method, H, S, G, Tm_at_0.1mM),
             df.Meltwin.summary %>% arrange(Experiment) %>% select(dH.kcal.mol, dS.kcal.molK, dG.kcal.mol, Tm.degC))
df.G1 = df.G1 %>% filter(Method != "2 Tm versus ln[Ct]")
df.G2 = cbind(df.MeltR.summary %>% arrange(Experiment) %>% filter(Method != "3 Global fit") %>% select(Experiment, Method, H, S, G, Tm_at_0.1mM),
              df.Meltwin.summary %>% arrange(Experiment) %>% select(dH.kcal.mol, dS.kcal.molK, dG.kcal.mol, Tm.degC))
df.G2 = df.F %>% filter(Method == "2 Tm versus ln[Ct]")
df.G = rbind(df.G1, df.G2)

P1 = ggplot(df.G, aes(x = H, y = dH.kcal.mol, color = Method)) +
  stat_regline_equation(aes(label = ..rr.label..)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  ggtitle("dH (kcal/mol)") +
  xlab("MeltR") +
  ylab("Meltwin") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

P2 = ggplot(df.G, aes(x = S, y = dS.kcal.molK, color = Method)) +
  stat_regline_equation(aes(label = ..rr.label..)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  ggtitle("dS (cal/mol/K)") +
  xlab("MeltR") +
  ylab("Meltwin") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

P3 = ggplot(df.G, aes(x = G, y = dG.kcal.mol, color = Method)) +
  stat_regline_equation(aes(label = ..rr.label..)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  ggtitle("dG (kcal/mol)") +
  xlab("MeltR") +
  ylab("Meltwin") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

P4 = ggplot(df.G, aes(x = Tm_at_0.1mM, y = Tm.degC, color = Method)) +
  stat_regline_equation(aes(label = ..rr.label..)) +
  geom_abline(intercept = 0, slope = 1) +
  stat_regline_equation(aes(label = ..rr.label..)) +
  geom_point(alpha = 0.8) +
  theme_classic() +
  ggtitle("Tm (degC)") +
  xlab("MeltR") +
  ylab("Meltwin") +
  theme(legend.position = c(0.8, 0.35),
        legend.background = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

PG = plot_grid(P1, P2, P3, P4)

