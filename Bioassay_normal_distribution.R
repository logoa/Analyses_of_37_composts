# Check normal distribution for all data sets
# Batches and concentrations separatly
# Author: Anja Logo
# Date: 16.08.12

# Here were all the composts included for all the diseases systems!
# CHALLENGE: Has to be repeated if we decide to exclude data (F.e. K6, K10, K27, K36)
# Descion was based on mainly on qqplot, whereas the shapiro.test was taken as a control

#---------------------------------Preparations-------------------------------
# Load packages
library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(factoextra)
library(multcompView)
library(corrplot)
library(rstatix)
library(Hmisc)
library(knitr)
library(reshape2)
library(kableExtra)
library(ggpubr)
library(dunn.test)
library(readxl)

setwd("C:/Users/anja.logo/ownCloud/Kompostmikrobio/Statistics/Bioassays/")
source("../20230227_Bioassays_funct.R", local = knitr::knit_global()) # Loading Functions

# Load data
# Gu-cress
df_cp <-read.csv("Compost_ds_gp/20230728_df_cp.csv", header =TRUE, sep=";")
df_cp$X <-NULL
factors <- colnames(df_cp[,1:9])
df_cp[factors] <-lapply(df_cp[factors], factor)
df_cp$treatment <- as.character(df_cp$treatment)
df_cp$treatment <- factor(df_cp$treatment, levels=unique(df_cp$treatment))
df_cp$batch <- as.integer(df_cp$batch)

# Gu-cucumber rep 1
df_cu_gu <- read.csv("Compost_ds_gp/20230227_df_cu_gu.csv", header = TRUE, sep = ";")
df_cu_gu$X <-NULL
factors <- colnames(df_cu_gu[,c(2:7)])
df_cu_gu[factors] <-lapply(df_cu_gu[factors], factor)
df_cu_gu$treatment <- as.character(df_cu_gu$treatment)
df_cu_gu$treatment <- factor(df_cu_gu$treatment, levels=unique(df_cu_gu$treatment)) # Do not delete, sets the order of the levels
df_cu_gu$treatment <-relevel(df_cu_gu$treatment, "Std5")
df_cu_gu$treatment <-relevel(df_cu_gu$treatment, "Std4")
df_cu_gu$treatment <-relevel(df_cu_gu$treatment, "Std3")
df_cu_gu$treatment <-relevel(df_cu_gu$treatment, "Std2")
df_cu_gu$batch <- as.integer(df_cu_gu$batch)

# Gu-cucumber rep 2
df_cu_gu_r_all <-read.csv("Compost_ds_gp/20230227_df_cu_gu_R.csv", header = TRUE, sep = ";") # repetition batch III
df_cu_gu_r_all$X <-NULL
factors <- colnames(df_cu_gu_r_all[,c(2:7)])
df_cu_gu_r_all[factors] <-lapply(df_cu_gu_r_all[factors], factor)
df_cu_gu_r_all$treatment <- as.character(df_cu_gu_r_all$treatment)
df_cu_gu_r_all$treatment <- factor(df_cu_gu_r_all$treatment, levels=unique(df_cu_gu_r_all$treatment)) # Do not delete, sets the order of the levels
df_cu_gu_r_all$treatment <-relevel(df_cu_gu_r_all$treatment, "Std5")
df_cu_gu_r_all$treatment <-relevel(df_cu_gu_r_all$treatment, "Std4")
df_cu_gu_r_all$treatment <-relevel(df_cu_gu_r_all$treatment, "Std3")
df_cu_gu_r_all$treatment <-relevel(df_cu_gu_r_all$treatment, "Std2")
df_cu_gu_r_all$batch <- as.integer(df_cu_gu_r_all$batch)

# Rs-cucumber
df_cu_rs <-read.csv("Compost_ds_gp/20230227_df_cu_rs.csv", header = TRUE, sep = ";")
df_cu_rs$X <-NULL
df_cu_rs$treatment <- as.character(df_cu_rs$treatment)
df_cu_rs$treatment <- factor(df_cu_rs$treatment, levels=unique(df_cu_rs$treatment))
df_cu_rs$treatment <-relevel(df_cu_rs$treatment, "Std5")
df_cu_rs$treatment <-relevel(df_cu_rs$treatment, "Std4")
df_cu_rs$treatment <-relevel(df_cu_rs$treatment, "Std3")
df_cu_rs$treatment <-relevel(df_cu_rs$treatment, "Std2")

rm(factors)

df.normal <- read_excel("Bioassay_normal_distribution.xlsx")

#---------------------------------Batches overall---------
df.normal.all <- df.normal[df.normal$batch == "all",]

# Function to test normality batch overall

test.normal.dis = function(data, factor){
  data$x = data$treatment
  data$y = factor
  ANOVA = aov (y ~ x, data = data)
  par(mfrow=c(2,2))
  plot(ANOVA)
  s.test =shapiro.test(resid(ANOVA))
  print(s.test$p.value)
}

# Growth promotion
# in larger data sets shapiro.test is not reliable anymore!

data =df_cp[df_cp$conc== 0,]
s1 = test.normal.dis(data, data$gp)
data =df_cu_gu[df_cu_gu$conc== 0,]
s2 = test.normal.dis(data, data$gp)
data =df_cu_rs[df_cu_rs$conc ==0,]
s3 = test.normal.dis(data, (data$gp))

df.normal.all[df.normal.all$variable =="gp",]$shapiro = c(s1,s2,s3)
df.normal.all[df.normal.all$variable =="gp",]$qqplot = c(1,1,0)
df.normal.all[df.normal.all$variable =="gp",]$transformation = rep("none",3)
df.normal.all[df.normal.all$variable =="gp",]$test = c("ANOVA","ANOVA", "KRUSKAL")
df.normal.all[df.normal.all$variable =="gp",]$variance = rep(1,3)
rm(s1,s2,s3, data)

#  disease suppression
# gu cress
data =df_cp %>% 
  filter(
      (batch ==1 & conc == 1.35)|
      (batch ==2 & conc == 0.45)|
      (batch ==3 & conc == 1.35)|
      (batch ==4 & conc == 1.35))

s1 = test.normal.dis(data, log(data$biomass+1))
s2 = test.normal.dis(data, sqrt(data$rel))

# gu cucumber
data =df_cu_gu %>% 
  filter(
    (batch ==1 & conc == 1.35)|
      (batch ==2 & conc == 0.45)|
      (batch ==3 & conc == 1.35)|
      (batch ==4 & conc == 1.35))
s3 = test.normal.dis(data, data$biomass)
s4 = test.normal.dis(data, data$rel)
s5 = test.normal.dis(data, data$surv)

# rs cucumber
data =df_cu_rs %>% 
  filter(
    (batch ==1 & conc == 1.4)|
      (batch ==2 & conc == 0.8)|
      (batch ==3 & conc == 0.8)|
      (batch ==4 & conc == 0.8))
s6 = test.normal.dis(data, data$biomass)
s7 = test.normal.dis(data, data$rel)
s8 = test.normal.dis(data, data$surv)

df.normal.all[df.normal.all$conc =="mixed",]$shapiro = c(s1,s2,s3,s4,s5,s6,s7,s8)
df.normal.all[df.normal.all$conc =="mixed",]$qqplot = c(1,1,1,1,0,0,0,0)
df.normal.all[df.normal.all$conc =="mixed",]$variance = c(1,1,1,1,0,0,0,0)
df.normal.all[df.normal.all$conc =="mixed",]$transformation = c("log+1", "sqrt", rep("none",6))
df.normal.all[df.normal.all$conc =="mixed" & df.normal.all$qqplot == 1,]$test = "ANOVA"
df.normal.all[df.normal.all$conc =="mixed" & df.normal.all$qqplot == 0,]$test = "KRUSKAL"

rm(s1, s2, s3, s4, s5, s6, s7, s8)
#---------------------------------G.ultimum-cress----------------------------------------

# Loop to generate ANOVA and check assumptions of ANOVA

df.normal.cr <-df.normal[df.normal$batch != "all" & df.normal$plant.pathogen =="globi-cress",]

# Loop shapiro-test and qqplot
for (i in 1:nrow(df.normal.cr)) {
  data = df_cp[df_cp$batch ==  df.normal.cr$batch[i] & df_cp$conc == df.normal.cr$conc[i],]
  data$x = data$treatment
  data$y = data[, df.normal.cr$variable[i]] 
  ANOVA = aov (y ~ x, data = data)
  par(mfrow=c(2,2))
  plot(ANOVA)
  mtext(paste("batch", df.normal.cr$batch[i], "conc", df.normal.cr$conc[i], df.normal.cr$variable[i]), side =3, line= -21, outer =TRUE)
  s.test =shapiro.test(resid(ANOVA))
  df.normal.cr[i,"shapiro"] <- s.test$p.value
  rm(data, s.test, ANOVA)
}

df.normal.cr$qqplot   <-c(0,1,1,1,1,
                          0,1,1,1,0,
                          1,1,0,1,1,
                          1,1,1,1,1,
                          0,1,1,0,1,
                          1)
df.normal.cr$variance <-c(1,1,1,1,1,
                          1,1,1,1,1,
                          1,1,1,1,1,
                          1,1,1,1,1,
                          1,1,1,1,1,
                          1)
for (i in 1:nrow(df.normal.cr)) {
  if (df.normal.cr[i,"qqplot"] == 1) {df.normal.cr[i, "transformation"] = "none"
  }
}

df.not.normal <-df.normal.cr[df.normal.cr$qqplot ==0,]

# log tranformation, sqrt does not improve for normal distribution!
for (i in 1:nrow(df.not.normal)) {
  data = df_cp[df_cp$batch ==  df.not.normal$batch[i] & df_cp$conc == df.not.normal$conc[i],]
  data$x = data$treatment
  data$y = data[, df.not.normal$variable[i]] 
  ANOVA = aov (log(y+ 1) ~ x, data = data)
  par(mfrow=c(2,2))
  plot(ANOVA)
  mtext(paste("batch", df.not.normal$batch[i], "conc", df.not.normal$conc[i], df.not.normal$variable[i]), side =3, line= -21, outer =TRUE)
  s.test =shapiro.test(resid(ANOVA))
  df.not.normal[i,"shapiro"] <- s.test$p.value
  rm(data, s.test, ANOVA)
}
df.not.normal$transformation <- c("none", "log+1", "log+1",
                                  "none", "log+1", "log+1")
df.not.normal$qqplot <- c(0,1,1,
                          0,1,1)
df.normal.cr[df.normal.cr$qqplot ==0,] <-df.not.normal

# Select the accuarate test
for (i in 1:nrow(df.normal.cr)) {
  if (df.normal.cr[i,"qqplot"] == 1) {df.normal.cr[i,"test"] <- "ANOVA"
  } else df.normal.cr[i,"test"] <- "KRUSKAL" 
}
rm(i, df.not.normal)

#---------------------------------G.ultimum-cucumber----------------------------------------

df.normal.gu.cuc <-df.normal[df.normal$batch != "all" & df.normal$plant.pathogen =="globi-cuc",]

for (i in 1:nrow(df.normal.gu.cuc)) {
  data = df_cu_gu[df_cu_gu$batch ==  df.normal.gu.cuc$batch[i] & df_cu_gu$conc == df.normal.gu.cuc$conc[i],]
  data$x = data$treatment
  data$y = data[, df.normal.gu.cuc$variable[i]] 
  ANOVA = aov (y ~ x, data = data)
  par(mfrow=c(2,2))
  plot(ANOVA)
  mtext(paste("batch", df.normal.gu.cuc$batch[i], "conc", df.normal.gu.cuc$conc[i], df.normal.gu.cuc$variable[i]), side =3, line= -21, outer =TRUE)
  s.test =shapiro.test(resid(ANOVA))
  df.normal.gu.cuc[i,"shapiro"] <- round(s.test$p.value,3)
  rm(data, s.test, ANOVA)
}
df.normal.gu.cuc$qqplot   <-c(1,1,1,1,1,
                              0,1,1,1,0,
                              1,0,1,1,1,
                              1,0,1,1,1,
                              1,1,0,1,1,
                              1, rep(0,4),
                              1,rep(0,3),1,
                              1,1)
df.normal.gu.cuc$variance <-c(1,1,1,1,1,
                              1,1,1,1,1,
                              1,1,1,1,1,
                              1,1,1,1,1,
                              1,1,1,1,1,
                              1,rep(1,11))

for (i in 1:nrow(df.normal.gu.cuc)) {
  if (df.normal.gu.cuc[i,"qqplot"] == 1) {df.normal.gu.cuc[i, "transformation"] = "none"
  }
}

df.not.normal <-df.normal.gu.cuc[df.normal.gu.cuc$qqplot ==0,]

# sqrt() or log tranformation do note improve the distribtution
for (i in 1:nrow(df.not.normal)) {
  data = df_cu_gu[df_cu_gu$batch ==  df.not.normal$batch[i] & df_cu_gu$conc == df.not.normal$conc[i],]
  data$x = data$treatment
  data$y = data[, df.not.normal$variable[i]] 
  ANOVA = aov (log(y+1) ~ x, data = data)
  par(mfrow=c(2,2))
  plot(ANOVA)
  mtext(paste("batch", df.not.normal$batch[i], "conc", df.not.normal$conc[i], df.not.normal$variable[i]), side =3, line= -21, outer =TRUE)
  s.test =shapiro.test(resid(ANOVA))
  df.not.normal[i,"shapiro"] <- s.test$p.value
  rm(data, s.test, ANOVA)
}

df.not.normal$qqplot = c(0,1,rep(0, nrow(df.not.normal)-2))

df.not.normal[df.not.normal$qqplot== 0,]$transformation <- "none"
df.not.normal[df.not.normal$qqplot!= 0,]$transformation <- "log+1"

df.normal.gu.cuc[df.normal.gu.cuc$qqplot ==0,] <-df.not.normal

df.normal.gu.cuc[df.normal.gu.cuc$qqplot ==0,]$test <- "KRUSKAL"
df.normal.gu.cuc[df.normal.gu.cuc$qqplot ==1,]$test <- "ANOVA"

rm(i, df.not.normal)

#---------------------------------R.solani-cucumber----------------------------------------

df.normal.rs.cuc <-df.normal[df.normal$batch != "all" & df.normal$plant.pathogen =="rsolani-cuc",]

for (i in 1:nrow(df.normal.rs.cuc)) {
  data = df_cu_rs[df_cu_rs$batch ==  df.normal.rs.cuc$batch[i] & df_cu_rs$conc == df.normal.rs.cuc$conc[i],]
  data$x = data$treatment
  data$y = data[, df.normal.rs.cuc$variable[i]] 
  ANOVA = aov (y ~ x, data = data)
  par(mfrow=c(2,2))
  plot(ANOVA)
  mtext(paste("batch", df.normal.rs.cuc$batch[i], "conc", df.normal.rs.cuc$conc[i], df.normal.rs.cuc$variable[i]), side =3, line= -21, outer =TRUE)
  s.test =shapiro.test(resid(ANOVA))
  df.normal.rs.cuc[i,"shapiro"] <- round(s.test$p.value,3)
  rm(data, s.test, ANOVA)
}
df.normal.rs.cuc$qqplot  <-c(1,0,0,0,1,
                              1,1,0,1,0,
                              0,1,1,1,1,
                              0,1,0,0,1,
                              0,1,0,0,1,
                              0,0,0)
df.normal.rs.cuc$variance <-c(1,1,1,1,1,
                              1,1,1,1,1,
                              1,0,1,1,1,
                              1,1,1,1,0,
                              1,1,1,1,1,
                              1,1,0)

for (i in 1:nrow(df.normal.rs.cuc)) {
  if (df.normal.rs.cuc[i,"qqplot"] == 1) {df.normal.rs.cuc[i, "transformation"] = "none"
  }
}
df.not.normal <-df.normal.rs.cuc[df.normal.rs.cuc$qqplot ==0,]

# log tranformation 
for (i in 1:nrow(df.not.normal)) {
  data = df_cu_rs[df_cu_rs$batch ==  df.not.normal$batch[i] & df_cu_rs$conc == df.not.normal$conc[i],]
  data$x = data$treatment
  data$y = data[, df.not.normal$variable[i]] 
  ANOVA = aov (log(y+1) ~ x, data = data)
  par(mfrow=c(2,2))
  plot(ANOVA)
  mtext(paste("batch", df.not.normal$batch[i], "conc", df.not.normal$conc[i], df.not.normal$variable[i]), side =3, line= -21, outer =TRUE)
  s.test =shapiro.test(resid(ANOVA))
  df.not.normal[i,"shapiro"] <- s.test$p.value
  rm(data, s.test, ANOVA)
}
df.not.normal$qqplot <- c(0,0,0,0,1,
                          0,0,1,0,0,
                          1,0,1,0,1)
df.not.normal[df.not.normal$qqplot ==1,]$transformation <- "log+1"

df.normal.rs.cuc[df.normal.rs.cuc$qqplot == 0, ] <- df.not.normal # substitut

df.not.normal = df.normal.rs.cuc[df.normal.rs.cuc$qqplot == 0,] # reisolate again the ones that were not normal ditributed

# sqrt transformation
for (i in 1:nrow(df.not.normal)) {
  data = df_cu_rs[df_cu_rs$batch ==  df.not.normal$batch[i] & df_cu_rs$conc == df.not.normal$conc[i],]
  data$x = data$treatment
  data$y = data[, df.not.normal$variable[i]] 
  ANOVA = aov (sqrt(y) ~ x, data = data)
  par(mfrow=c(2,2))
  plot(ANOVA)
  mtext(paste("batch", df.not.normal$batch[i], "conc", df.not.normal$conc[i], df.not.normal$variable[i]), side =3, line= -21, outer =TRUE)
  s.test =shapiro.test(resid(ANOVA))
  df.not.normal[i,"shapiro"] <- s.test$p.value
  rm(data, s.test, ANOVA)
}

df.not.normal$qqplot <-c(1,1,1,0,0,
                         0,0,0,0,0)
df.not.normal[df.not.normal$qqplot ==1,]$transformation <- "sqrt"
df.not.normal[df.not.normal$qqplot ==0,]$transformation <- "none"

df.normal.rs.cuc[df.normal.rs.cuc$qqplot ==0,] <- df.not.normal
df.normal.rs.cuc[df.normal.rs.cuc$qqplot ==0,]$test <- "KRUSKAL"
df.normal.rs.cuc[df.normal.rs.cuc$qqplot ==1,]$test <- "ANOVA"

df.not.normal = df.normal.rs.cuc[df.normal.rs.cuc$qqplot ==0,]

# back to non transformed
for (i in 1:nrow(df.not.normal)) {
  data = df_cu_rs[df_cu_rs$batch ==  df.not.normal$batch[i] & df_cu_rs$conc == df.not.normal$conc[i],]
  data$x = data$treatment
  data$y = data[, df.not.normal$variable[i]] 
  ANOVA = aov (y ~ x, data = data)
  par(mfrow=c(2,2))
  plot(ANOVA)
  mtext(paste("batch", df.not.normal$batch[i], "conc", df.not.normal$conc[i], df.not.normal$variable[i]), side =3, line= -21, outer =TRUE)
  s.test =shapiro.test(resid(ANOVA))
  df.not.normal[i,"shapiro"] <- s.test$p.value
  rm(data, s.test, ANOVA)
}
df.normal.rs.cuc[df.normal.rs.cuc$qqplot ==0,] <- df.not.normal
rm(df.not.normal, i)


#---------------------------------Create csv file-----------
data =rbind(df.normal.all, df.normal.cr, df.normal.gu.cuc, df.normal.rs.cuc)
data$shapiro = round(data$shapiro, 3)

#write.csv(data, file ="Bioassay_normal_distribution_filled.csv")

