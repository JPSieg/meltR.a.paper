library(tidyverse)

list.files("data-raw")

df.index = read.csv("data-raw/Znosko_data_index.csv")

head(df.index)

####Read in Znosko data####

list.df = {}

for (i in 1:nrow(df.index)){
  #print(i)
  path = paste("data-raw", df.index$File[i], sep = "/")
  
  Con = file(path)
  Lines = readLines(Con)
  close(Con)
  
  start = which(Lines == paste("CORRECTED ABSORBANCES @ ", df.index$Wavelength[i],".0 nm", sep = ""))
  Lines = Lines[(start+5):(length(Lines))]
  
  Cell = as.integer(strsplit(df.index$Cell[i], split = " ")[[1]][2]) + 1 
  
  Temperature = c()
  Absorbance = c()
  
  for (j in 1:length(Lines)){
    Temperature[j] = as.numeric(strsplit(Lines[j], split = ",,")[[1]][1])
    Absorbance[j] = as.numeric(strsplit(Lines[j], split = ",,")[[1]][Cell])
  }
  
  df = data.frame(Temperature, Absorbance)
  
  df$Experiment = df.index$Experiment[i]
  df$RNA = df.index$RNA[i]
  df$Cell = df.index$Cell[i]
  df$Buffer = df.index$Buffer[i]
  df$File = df.index$File[i]
  df$Blank = df.index$Blank[i]
  df$Pathlength = df.index$Pathlength[i]
  df$Wavelength = df.index$Wavelength[i]
  
  list.df[[i]] = df
  
}

####Subtract Znosko background####

df = bind_rows(list.df)

for (i in 1:length(list.df)){
  df.bg = df %>%
    filter(Experiment == list.df[[i]]$Experiment[1]) %>%
    filter(File == list.df[[i]]$File[1]) %>%
    filter(Cell == list.df[[i]]$Blank[1])
  list.df[[i]]$Absorbance = list.df[[i]]$Absorbance - (df.bg$Absorbance*list.df[[i]]$Pathlength)
  list.df[[i]]$Sample = i
}

df.Znosko = bind_rows(list.df) %>% filter(Absorbance != 0)

####Read in Bevilacqua data####

list.files("data-raw")

df.index = read.csv("data-raw/Bavilacqua_data_index.csv")

list.df = {}

for (i in 1:nrow(df.index)){
  df = read.csv(paste("data-raw/", df.index$File[i], sep = ""))
  df$Experiment = df.index$Experiment[i]
  df$RNA = df.index$RNA[i]
  df$Buffer = df.index$Buffer[i]
  df$File = df.index$File[i]
  df$Blank = df.index$Blank[i]
  df$Wavelength = df.index$Wavelength[i]
  list.df[[i]] = df
}

df.Bevilacqua = bind_rows(list.df)

table(df.Bevilacqua$Experiment)

colnames(df.Znosko)
colnames(df.Bevilacqua)

df.absorbance = bind_rows(df.Znosko %>% select("Experiment", "RNA", "Buffer", "File",
                     "Blank", "Wavelength", "Sample", "Temperature", "Absorbance",  "Pathlength"),
          df.Bevilacqua %>% select("Experiment", "RNA", "Buffer", "File",
                               "Blank", "Wavelength", "Sample", "Temperature", "Absorbance",  "Pathlength"))

usethis::use_data(df.absorbance, overwrite = T)

####Read in indv fit data from Meltwin####

list.files("data-raw")

df.index = read.csv("data-raw/Meltwin_indv_fits_data_index.csv")

list.df = {}

for (i in 1:nrow(df.index)){
  df = read.csv(paste("data-raw", df.index$File[i], sep = "/"))
  df$Experiment = df.index$Experiment[i]
  df$RNA = df.index$RNA[i]
  df$Buffer = df.index$Buffer[i]
  df$File = df.index$File[i]
  df$Molecular_model = df.index$Molecular_model[i]
  df$Wavelength = df.index$Wavelength[i]
  df$Source = df.index$Source[i]
  list.df[[i]] = df
}

df.Meltwin.indv.fits = bind_rows(list.df)

usethis::use_data(df.Meltwin.indv.fits, overwrite = T)

####Read in data from Meltwin####

list.files("data-raw")

df.index = read.csv("data-raw/Meltwin_summary_data_index.csv")

list.df = {}

for (i in 1:nrow(df.index)){
  df = read.csv(paste("data-raw", df.index$File[i], sep = "/"))
  df$Experiment = df.index$Experiment[i]
  df$RNA = df.index$RNA[i]
  df$Buffer = df.index$Buffer[i]
  df$File = df.index$File[i]
  df$Molecular_model = df.index$Molecular_model[i]
  df$Wavelength = df.index$Wavelength[i]
  df$Source = df.index$Source[i]
  list.df[[i]] = df
}

df.Meltwin.summary = bind_rows(list.df)

usethis::use_data(df.Meltwin.summary, overwrite = T)

####Read in Adams data####

df.index = read.csv("data-raw/Adams_data_index.csv")

df.index$File = gsub("MirandaAData", "", gsub("\\", "/", df.index$File, fixed = TRUE), fixed = TRUE)

df.index$File  = paste("data-raw", df.index$File, sep = "")

list.df = {}

for (i in 1:nrow(df.index)){
  df = read.csv(df.index$File[i]) %>% select(Sample, Pathlength, Temperature, Absorbance)
  df$Experiment = df.index$Experiment[i]
  df$RNA = df.index$RNA[i]
  df$Buffer = df.index$Buffer[i]
  df$File = df.index$File[i]
  df$Blank = df.index$Blank[i]
  df$Wavelength = df.index$Wavelength[i]
  df$Sample = paste("Sample", df$Sample)
  Samples = paste("Sample", unique(df$Sample))
  blank = unique(df$Blank)
  list.df.samples = {}
  for (j in Samples){
    df.blank = df %>% filter(Sample == blank)
    df.A = df %>% filter(Sample == j)
    df.A$Absorbance = df.A$Absorbance - df.A$Pathlength*df.blank$Absorbance/df.blank$Pathlength
    
  }
}
