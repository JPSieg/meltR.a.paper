library(MeltR)
library(tidyverse)
library(cowplot)
devtools::load_all()

fit = meltR.A(df.absorbance %>% filter(Experiment == "Helix A"),
        blank = 1,
        NucAcid = c("RNA", "CGAAAGGU", "ACCUUUCG"),
        Mmodel = "Heteroduplex.2State",
        fitTs = c(25, 80))

Trim = BLTrimmer(fit)

dH.range = c(-80, -60)
Tm.range = c(40, 55)

####Define functions####

Mmodel_names <- c("Monomolecular.2State",
                  "Monomolecular.3State",
                  "Heteroduplex.2State",
                  "Homoduplex.2State")
Mmodels <- list(function(K){ (K/(1+K)) },
                function(K1, K2, Ct){ 1/(1 + K1 + (K1*K2)) },
                function(K, Ct){ ( (2/(K*Ct)) + 2 - sqrt(((((2/(K*Ct)) + 2)^2) - 4)))/2 },
                function(K, Ct){ ((1/(2*K*Ct)) + 2 - sqrt(((((1/(2*K*Ct)) + 2)^2) - 4)))/2 })
names(Mmodels) <- Mmodel_names

G_VH = function(H, S, Temperature){exp((S/0.0019872) - ((1/((Temperature + 273.15)*0.0019872))*H))}

f = function(H, S, Temperature, Ct){
  K <- G_VH(H = H, S = S, Temperature = Temperature)
  model <- Mmodels$Heteroduplex.2State(K = K, Ct = Ct)
  return(model)
}

TmModel = function(H, S, lnCt){
  ((0.0019872/H)*lnCt) + ((S - 0.0019872*log(4))/H)
}

calcS = function(H, Tm, Ct){ (H/(273.15 + Tm)) + (0.0019872*log(4/Ct)) }

####Fit a real but ideal melting curve####

df = df.absorbance %>% filter(Experiment == "Helix A")

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point()

df = df %>% filter(Sample == 4)

fit = meltR.A(df,
              NucAcid = c("RNA", "CGAAAGGU", "ACCUUUCG"),
              Mmodel = "Heteroduplex.2State",
              methods = c(TRUE, FALSE, FALSE))


df$f = f(fit$Summary$dH[1], fit$Summary$dS[1]/1000, df$Temperature, fit$Method.1.indvfits$Ct)

plot(df$Temperature, df$f)

low = df$Temperature[which.min(abs(df$f - 0.9))]
high = df$Temperature[which.min(abs(df$f - 0.1))]
abline(v = low)
abline(h = 0.9)
abline(h = 0.1)
abline(v = high)

low
high

v.lows = seq(low, (low - 25), by = -1)
v.highs = seq(high, (high + 25), by = +1)

df.ranges = data.frame(v.lows,
                       v.highs)
df.ranges$Width = df.ranges$v.highs - df.ranges$v.lows
df.ranges$Baseline.length = df.ranges$v.highs - high

mDS = c()
bDS = c()
mSS = c()
bSS = c()
dH = c()
dS = c()
dG = c()
Tm = c()

for (i in 1:nrow(df.ranges)){
  tryCatch({
    fit = meltR.A(df,
                  NucAcid = c("RNA", "CGAAAGGU", "ACCUUUCG"),
                  Mmodel = "Heteroduplex.2State",
                  methods = c(TRUE, FALSE, FALSE),
                  fitTs = c(df.ranges$v.lows[i], df.ranges$v.highs[i]),
                  Silent = TRUE)
    mDS[i] = coef(fit$Method.1.fit[[1]])[3]
    bDS[i] = coef(fit$Method.1.fit[[1]])[4]
    mSS[i] = coef(fit$Method.1.fit[[1]])[5]
    bSS[i] = coef(fit$Method.1.fit[[1]])[6]
    dH[i] = fit$Method.1.indvfits$dH[1]
    dS[i] = fit$Method.1.indvfits$dS[1]
    dG[i] = fit$Method.1.indvfits$dG[1]
    Tm[i] = fit$Method.1.indvfits$Tm[1]
    
    
  }, error = function(e){})
  
}

df.ranges$mDS = mDS
df.ranges$bDS = bDS
df.ranges$mSS = mSS
df.ranges$bSS = bSS
df.ranges$dH = dH
df.ranges$dS = dS
df.ranges$dG = dG
df.ranges$Tm = Tm

calcTM = function(H, S, Ct){(H/((S/1000) - 0.0019872*log(4/Ct))) - 273.15}

dH.known = -68.40
dS.known = -187.89
dG.known = -10.13
Tm.known = calcTM(dH.known, dS.known, fit$Method.1.indvfits$Ct[1])

top = ggplot() +
  geom_abline(data = df.ranges, mapping = aes(slope = mDS, intercept = bDS, color = Width)) +
  geom_abline(data = df.ranges, mapping = aes(slope = mSS, intercept = bSS, color = Width)) +
  geom_point(data = df, mapping = aes(x = Temperature, y = Absorbance),
             color = "black", alpha = 0.2) +
  theme_classic() +
  scale_color_viridis_b() +
  theme(legend.position = "none")  +
  xlab("Temperature (\u00B0C)") +
  theme(axis.text = element_text(color = "black"))

dH = ggplot(df.ranges) +
  geom_hline(yintercept = dH.known) +
  geom_point(mapping = aes(x = Baseline.length, y = dH, color = Baseline.length)) +
  theme_classic() +
  scale_y_continuous(limits = dH.range) +
  scale_color_viridis_b() +
  theme(legend.position = "none")  +
  xlab("Baseline length (\u00B0C)") +
  ylab("\u0394H\u00B0 (kcal/mol)") +
  theme(axis.text = element_text(color = "black"))

dS.range = c(-187.89-20, -187.89+20)

dS = ggplot(df.ranges) +
  geom_hline(yintercept = dS.known) +
  geom_point(mapping = aes(x = Baseline.length, y = dS, color = Baseline.length)) +
  theme_classic() +
  scale_y_continuous(limits = dS.range) +
  scale_color_viridis_b() +
  theme(legend.position = "none")  +
  xlab("Baseline length (\u00B0C)") +
  ylab("\u0394S\u00B0 (cal/mol/K)") +
  theme(axis.text = element_text(color = "black"))

dG.range = c(-10.13-2, -10.13+2)

dG = ggplot(df.ranges) +
  geom_hline(yintercept = dG.known) +
  geom_point(mapping = aes(x = Baseline.length, y = dG, color = Baseline.length)) +
  theme_classic() +
  scale_y_continuous(limits = dG.range) +
  scale_color_viridis_b() +
  theme(legend.position = "none")  +
  xlab("Baseline length (\u00B0C)") +
  ylab("\u0394S\u00B0 (cal/mol/K)") +
  theme(axis.text = element_text(color = "black"))

Tm = ggplot(df.ranges) +
  geom_hline(yintercept = Tm.known) +
  geom_point(mapping = aes(x = Baseline.length, y = Tm, color = Baseline.length)) +
  theme_classic() +
  scale_y_continuous(limits = Tm.range) +
  scale_color_viridis_b() +
  theme(legend.position = "none") +
  xlab("Baseline length (\u00B0C)") +
  ylab("Tm (\u00B0C)") +
  theme(axis.text = element_text(color = "black"))

Ideal = plot_grid(top, dH, dS, dG, Tm,
                      align = "v",
                      ncol = 1,
                  rel_heights = c(2, 1, 1, 1, 1))

####Make modeled data####


df = df.absorbance %>% filter(Experiment == "Helix A")

head(df)

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point()

df = df %>% filter(Sample == 4)

fit = meltR.A(df,
        NucAcid = c("RNA", "CGAAAGGU", "ACCUUUCG"),
        Mmodel = "Heteroduplex.2State",
        methods = c(TRUE, FALSE, FALSE))


fit$Summary$dS[1]
fit$Summary$dH[1]

df$f = f(dH.known, dS.known/1000, df$Temperature, fit$Method.1.indvfits$Ct)

plot(df$Temperature, df$f)

Tm = dH.known/((dS.known/1000) - 0.00198720425864083*log(4/fit$Method.1.indvfits$Ct))

mDS = df.ranges$mDS[15]
bDS = df.ranges$bDS[15]
mSS = df.ranges$mSS[15]
bSS = df.ranges$bSS[15]
Ct = fit$Method.1.indvfits$Ct[1]

Tm = fit$Summary$H[1]/((fit$Summary$S[1]/1000) - 0.00198720425864083*log(4/Ct))
Tm-273.15

df.model = df

df.model$Absorbance = (10^5)*((mDS*df.model$Temperature + bDS)*Ct*df.model$f + (mSS*df.model$Temperature + bSS)*Ct*(1 - df.model$f)) + rnorm(nrow(df.model), 0, 0.00005)

plot(df.model$Temperature, df.model$Absorbance)

####Fit modeled data####

plot(df.model$Temperature, df.model$f)

low = df.model$Temperature[which.min(abs(df.model$f - 0.9))]
high = df.model$Temperature[which.min(abs(df.model$f - 0.1))]
abline(v = low)
abline(h = 0.9)
abline(h = 0.1)
abline(v = high)

low
high

v.lows = seq(low, (low - 25), by = -1)
v.highs = seq(high, (high + 25), by = +1)

df.model.ranges = data.frame(v.lows,
                       v.highs)
df.model.ranges$Width = df.model.ranges$v.highs - df.model.ranges$v.lows
df.model.ranges$Baseline.length = df.model.ranges$v.highs - high

mDS = c()
bDS = c()
mSS = c()
bSS = c()
dH = c()
dS = c()
dG = c()
Tm = c()

for (i in 1:nrow(df.ranges)){
  tryCatch({
    fit = meltR.A(df.model,
                  NucAcid = c("RNA", "CGAAAGGU", "ACCUUUCG"),
                  Mmodel = "Heteroduplex.2State",
                  methods = c(TRUE, FALSE, FALSE),
                  fitTs = c(df.ranges$v.lows[i], df.ranges$v.highs[i]),
                  Silent = TRUE)
    mDS[i] = coef(fit$Method.1.fit[[1]])[3]
    bDS[i] = coef(fit$Method.1.fit[[1]])[4]
    mSS[i] = coef(fit$Method.1.fit[[1]])[5]
    bSS[i] = coef(fit$Method.1.fit[[1]])[6]
    dH[i] = fit$Method.1.indvfits$dH[1]
    dS[i] = fit$Method.1.indvfits$dS[1]
    dG[i] = fit$Method.1.indvfits$dG[1]
    Tm[i] = fit$Method.1.indvfits$Tm[1]
    
    
  }, error = function(e){})
  
}

df.ranges$mDS = mDS
df.ranges$bDS = bDS
df.ranges$mSS = mSS
df.ranges$bSS = bSS
df.ranges$dH = dH
df.ranges$dS = dS
df.ranges$dG = dG
df.ranges$Tm = Tm

calcTM = function(H, S, Ct){(H/((S/1000) - 0.0019872*log(4/Ct))) - 273.15}

dH.known = -68.40
dS.known = -187.89
dG.known = -10.13
Tm.known = calcTM(dH.known, dS.known, fit$Method.1.indvfits$Ct[1])

top = ggplot() +
  geom_abline(data = df.ranges, mapping = aes(slope = mDS, intercept = bDS, color = Width)) +
  geom_abline(data = df.ranges, mapping = aes(slope = mSS, intercept = bSS, color = Width)) +
  geom_point(data = df.model, mapping = aes(x = Temperature, y = Absorbance),
             color = "black", alpha = 0.2) +
  theme_classic() +
  scale_color_viridis_b() +
  theme(legend.position = "none")  +
  xlab("Temperature (\u00B0C)") +
  theme(axis.text = element_text(color = "black"))

dH = ggplot(df.ranges) +
  geom_hline(yintercept = dH.known) +
  geom_point(mapping = aes(x = Baseline.length, y = dH, color = Baseline.length)) +
  theme_classic() +
  scale_y_continuous(limits = dH.range) +
  scale_color_viridis_b() +
  theme(legend.position = "none")  +
  xlab("Baseline length (\u00B0C)") +
  ylab("\u0394H\u00B0 (kcal/mol)") +
  theme(axis.text = element_text(color = "black"))

dS.range = c(-187.89-20, -187.89+20)

dS = ggplot(df.ranges) +
  geom_hline(yintercept = dS.known) +
  geom_point(mapping = aes(x = Baseline.length, y = dS, color = Baseline.length)) +
  theme_classic() +
  scale_y_continuous(limits = dS.range) +
  scale_color_viridis_b() +
  theme(legend.position = "none")  +
  xlab("Baseline length (\u00B0C)") +
  ylab("\u0394S\u00B0 (cal/mol/K)") +
  theme(axis.text = element_text(color = "black"))

dG.range = c(-10.13-2, -10.13+2)

dG = ggplot(df.ranges) +
  geom_hline(yintercept = dG.known) +
  geom_point(mapping = aes(x = Baseline.length, y = dG, color = Baseline.length)) +
  theme_classic() +
  scale_y_continuous(limits = dG.range) +
  scale_color_viridis_b() +
  theme(legend.position = "none")  +
  xlab("Baseline length (\u00B0C)") +
  ylab("\u0394S\u00B0 (cal/mol/K)") +
  theme(axis.text = element_text(color = "black"))

Tm = ggplot(df.ranges) +
  geom_hline(yintercept = Tm.known) +
  geom_point(mapping = aes(x = Baseline.length, y = Tm, color = Baseline.length)) +
  theme_classic() +
  scale_y_continuous(limits = Tm.range) +
  scale_color_viridis_b() +
  theme(legend.position = "none") +
  xlab("Baseline length (\u00B0C)") +
  ylab("Tm (\u00B0C)") +
  theme(axis.text = element_text(color = "black"))

Model = plot_grid(top, dH, dS, dG, Tm,
                  align = "v",
                  ncol = 1,
                  rel_heights = c(2, 1, 1, 1, 1))


####Fit a real but non-ideal melting curve####

df = df.absorbance %>% filter(Experiment == "Helix A")

ggplot(df, aes(x = Temperature, y = Absorbance)) +
  facet_wrap(~Sample, scales = "free") +
  geom_point()

df = df %>% filter(Sample == 5)


fit = meltR.A(df,
              NucAcid = c("RNA", "CGAAAGGU", "ACCUUUCG"),
              Mmodel = "Heteroduplex.2State",
              methods = c(TRUE, FALSE, FALSE))


df$f = f(fit$Summary$dH[1], fit$Summary$dS[1]/1000, df$Temperature, fit$Method.1.indvfits$Ct)

plot(df$Temperature, df$f)

low = df$Temperature[which.min(abs(df$f - 0.9))]
high = df$Temperature[which.min(abs(df$f - 0.1))]
abline(v = low)
abline(h = 0.9)
abline(h = 0.1)
abline(v = high)

low
high

v.lows = seq(low, (low - 25), by = -1)
v.highs = seq(high, (high + 25), by = +1)

df.ranges = data.frame(v.lows,
                       v.highs)
df.ranges$Width = df.ranges$v.highs - df.ranges$v.lows
df.ranges$Baseline.length = df.ranges$v.highs - high

mDS = c()
bDS = c()
mSS = c()
bSS = c()
dH = c()
dS = c()
dG = c()
Tm = c()

for (i in 1:nrow(df.ranges)){
  tryCatch({
    fit = meltR.A(df,
                  NucAcid = c("RNA", "CGAAAGGU", "ACCUUUCG"),
                  Mmodel = "Heteroduplex.2State",
                  methods = c(TRUE, FALSE, FALSE),
                  fitTs = c(df.ranges$v.lows[i], df.ranges$v.highs[i]),
                  Silent = TRUE)
    mDS[i] = coef(fit$Method.1.fit[[1]])[3]
    bDS[i] = coef(fit$Method.1.fit[[1]])[4]
    mSS[i] = coef(fit$Method.1.fit[[1]])[5]
    bSS[i] = coef(fit$Method.1.fit[[1]])[6]
    dH[i] = fit$Method.1.indvfits$dH[1]
    dS[i] = fit$Method.1.indvfits$dS[1]
    dG[i] = fit$Method.1.indvfits$dG[1]
    Tm[i] = fit$Method.1.indvfits$Tm[1]
    
    
  }, error = function(e){})
  
}

df.ranges$mDS = mDS
df.ranges$bDS = bDS
df.ranges$mSS = mSS
df.ranges$bSS = bSS
df.ranges$dH = dH
df.ranges$dS = dS
df.ranges$dG = dG
df.ranges$Tm = Tm

calcTM = function(H, S, Ct){(H/((S/1000) - 0.0019872*log(4/Ct))) - 273.15}

dH.known = -68.40
dS.known = -187.89
dG.known = -10.13
Tm.known = calcTM(dH.known, dS.known, fit$Method.1.indvfits$Ct[1])

top = ggplot() +
  geom_abline(data = df.ranges, mapping = aes(slope = mDS, intercept = bDS, color = Width)) +
  geom_abline(data = df.ranges, mapping = aes(slope = mSS, intercept = bSS, color = Width)) +
  geom_point(data = df, mapping = aes(x = Temperature, y = Absorbance),
             color = "black", alpha = 0.2) +
  theme_classic() +
  scale_color_viridis_b() +
  theme(legend.position = "none")  +
  xlab("Temperature (\u00B0C)") +
  theme(axis.text = element_text(color = "black"))

dH = ggplot(df.ranges) +
  geom_hline(yintercept = dH.known) +
  geom_point(mapping = aes(x = Baseline.length, y = dH, color = Baseline.length)) +
  theme_classic() +
  scale_y_continuous(limits = dH.range) +
  scale_color_viridis_b() +
  theme(legend.position = "none")  +
  xlab("Baseline length (\u00B0C)") +
  ylab("\u0394H\u00B0 (kcal/mol)") +
  theme(axis.text = element_text(color = "black"))

dS.range = c(-187.89-20, -187.89+20)

dS = ggplot(df.ranges) +
  geom_hline(yintercept = dS.known) +
  geom_point(mapping = aes(x = Baseline.length, y = dS, color = Baseline.length)) +
  theme_classic() +
  scale_y_continuous(limits = dS.range) +
  scale_color_viridis_b() +
  theme(legend.position = "none")  +
  xlab("Baseline length (\u00B0C)") +
  ylab("\u0394S\u00B0 (cal/mol/K)") +
  theme(axis.text = element_text(color = "black"))

dG.range = c(-10.13-2, -10.13+2)

dG = ggplot(df.ranges) +
  geom_hline(yintercept = dG.known) +
  geom_point(mapping = aes(x = Baseline.length, y = dG, color = Baseline.length)) +
  theme_classic() +
  scale_y_continuous(limits = dG.range) +
  scale_color_viridis_b() +
  theme(legend.position = "none")  +
  xlab("Baseline length (\u00B0C)") +
  ylab("\u0394S\u00B0 (cal/mol/K)") +
  theme(axis.text = element_text(color = "black"))

Tm = ggplot(df.ranges) +
  geom_hline(yintercept = Tm.known) +
  geom_point(mapping = aes(x = Baseline.length, y = Tm, color = Baseline.length)) +
  theme_classic() +
  scale_y_continuous(limits = Tm.range) +
  scale_color_viridis_b() +
  theme(legend.position = "none") +
  xlab("Baseline length (\u00B0C)") +
  ylab("Tm (\u00B0C)") +
  theme(axis.text = element_text(color = "black"))

Non.ideal = plot_grid(top, dH, dS, dG, Tm,
                  align = "v",
                  ncol = 1,
                  rel_heights = c(2, 1, 1, 1, 1))


####Consolidate plot####

P = plot_grid(Model, Ideal, Non.ideal, nrow = 1,
              labels = c("A", "B", "C"), label_size = 16)

list.files("Figures/SI_Figure_X_BLtrimming_theory")

ggsave("Figures/SI_Figure_X_BLtrimming_theory/SI_Figure_X_BLtrimming_theory.svg", P,
       width = 5.3, height = 4.5, units = "in", scale = 2)
