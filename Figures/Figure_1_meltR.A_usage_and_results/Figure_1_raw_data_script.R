
library(tidyverse)
library(cowplot)
library(viridis)
library(MeltR)
devtools::load_all()

####Single state####

unique(df.absorbance$Experiment)

df = df.absorbance %>% filter(Experiment == "CROWD DP1")

fit = meltR.A(df,
              NucAcid = c("RNA", "CGCGCG"),
              Mmodel = "Homoduplex.2State",
              wavelength = 280,
              fitTs = c(20, 70))

unique(df$Sample)

df$Sample = factor(df$Sample,
       levels = c("2",  "3",  "4",  "6",  "7",  "8",  "10", "11", "12"),
       labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9"))

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  geom_point(alpha = 0.1) +
  facet_wrap(~Sample) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_y_continuous(limits = c(0, 2), breaks = c(0, 0.5, 1, 1.5, 2)) +
  xlab("Temperature (\u00B0C)") +
  ylab("Absorbance (At 280 nm)") +
  theme(axis.text = element_text(color = "black"))

list.files("Figures/Figure_1_meltR.A_usage_and_results")

ggsave("Figures/Figure_1_meltR.A_usage_and_results/Figure_1_raw_data.svg",
       width = 1.731,
       height = 1.04,
       units = "in",
       scale = 3.5)
