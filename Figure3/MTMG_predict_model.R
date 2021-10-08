setwd("~/bioinformatics_analysis/MT_MG_paper/Figures/Predict_updown/")
predict_data <- read.csv("all_species_expr_predict.csv")
predict_data$Compar <- as.factor(predict_data$Compar)
predict_data$Species <- as.factor(predict_data$Species)

AR2<-lme(LogChange_abs ~ LogPriorExpr + Species + LogPriorExpr*Species, random=~1|SID/Compar,
         data=predict_data, na.action = (na.omit), method = "REML")

summary(AR2)
anova(AR2)
post_hoc <- summary(glht(AR2, linfct=mcp(Species="Tukey")), test = adjusted(type = "hochberg"))
post_hoc_cld <- cld(post_hoc,level=0.05)
post_hoc_cld