library(tidyverse)
library(cowplot)
library(viridis)
library(MeltR)
devtools::load_all()

####Multi state analysis####

list.files("Figures/Figure_4_BLTrimmer_two_state_check")

df = read.csv("Figures/Figure_4_BLTrimmer_two_state_check/js6040_multistate_data.csv")

ggplot(df, aes(x = Temperature, y = Absorbance, color = Sample)) +
  geom_point() +
  geom_vline(xintercept = c(5,75))

?meltR.A

fit = meltR.A(df,
        NucAcid = c("Custom", 100300, 101100),
        Mmodel = "Heteroduplex.2State",
        fitTs = c(5, 75),
        Save_results = "all",
        file_path = "Figures/Figure_4_BLTrimmer_two_state_check",
        file_prefix = "Multistate")

Trim = BLTrimmer(meltR.A.fit = fit, n.combinations = 1000, Save_results = "all",
                 file_path = "Figures/Figure_4_BLTrimmer_two_state_check",
                 file_prefix = "Multistate",
                 memory.light = F)

Trim$System.time

head(Trim$Baseline.data)

df.MS = Trim$Baseline.data
df.MS$Helix = "Multi-state" 

####Single state####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "CROWD DP1")

ggplot(df, aes(x = Temperature, y = Absorbance, color = Sample)) +
  geom_point() +
  geom_vline(xintercept = c(20,75))


#?meltR.A

fit = meltR.A(df,
              NucAcid = c("RNA", "CGCGCG"),
              Mmodel = "Homoduplex.2State",
              wavelength = 280,
              fitTs = c(20, 70),
              Save_results = "all",
              file_path = "Figures/Figure_4_BLTrimmer_two_state_check",
              file_prefix = "Two-state")

Trim = BLTrimmer(fit, n.combinations = 1000, Save_results = "all",
                 file_path = "Figures/Figure_4_BLTrimmer_two_state_check",
                 file_prefix = "Two-state",
                 memory.light = F)

Trim$System.time

head(Trim$Baseline.data)

df.2S = Trim$Baseline.data
df.2S$Helix = "Two-state" 

####Consolidate and format####

df = bind_rows(df.2S,
               df.MS)

unique(df$Helix)

df$Helix = factor(df$Helix,
                  levels = c("Two-state", "Multi-state"),
                  labels = c("Two-state", "Non two-state"))

colnames(df)

####Make B####

unique(df$Helix)

A = ggplot()+
  geom_histogram(data = df %>% filter(Helix == "Two-state"),
                 mapping = aes(x = 100*frac.dH1.dH2.error, fill = Helix),
                 bins = 50, alpha = 1, fill = viridis(5)[1]) +
  geom_histogram(data = df %>% filter(Helix != "Two-state"),
                 mapping = aes(x = 100*frac.dH1.dH2.error, fill = Helix),
                 bins = 50, alpha = 0.8, fill = viridis(5)[3]) +
  xlab("%Difference in \u0394H\u00B0 between\n Method 1 and Method 2") +
  ylab("Frequency") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 35), breaks = c(0, 5, 10, 15, 20, 25, 30,5)) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.35, 0.75),
        legend.background = element_blank()) 

####Make C####

B = ggplot(df) +
  geom_point(aes(x = 100*frac.dH1.dH2.error, y = dH1), color = viridis(5)[3], alpha = 0.1) +
  geom_point(aes(x = 100*frac.dH1.dH2.error, y = dH2), color = viridis(5)[2], alpha = 0.1) +
  facet_wrap(~Helix, scales = "free") +
  xlab("%Difference in \u0394H\u00B0 between\n Method 1 and Method 2") +
  ylab("\u0394H\u00B0 (kcal/mol") +
  theme_classic() +
  scale_fill_manual(values = viridis(5), name = "") +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.5, 0.5)) 

####Save plots####

plot_grid(A, ncol = 1, labels = c("C"), label_size = 14)

list.files("Figures/Figure_4_BLTrimmer_two_state_check")

ggsave("Figures/Figure_4_BLTrimmer_two_state_check/Figure_4_Two_state_test.svg", scale = 1.2,
       width = 5, height = 2.5)
