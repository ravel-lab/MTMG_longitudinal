library("nlme")
require("multcomp")
setwd("~/bioinformatics_analysis/PGATES/MT/MG_comparison/")
rela_data <- read.csv("species_relative_exp.csv",sep=",")

rela_data_trim <- rela_data[rela_data$Read_keep_5mil != 0, ]

AR1<-lme(Ratio ~Species + LogAbund + LogAbund*Species, random=~1|SID/UID,
       data=rela_data_trim, na.action = (na.omit), method = "REML")
summary(AR1)
anova(AR1)
TukeyHSD(AR1)
