###libraries####
library(ggplot2)
library(cowplot)
library(broom)
# library(lattice)
library(car)
library(reshape2)
# library(reshape)
# library(TSA)
library(plyr)
library(dplyr)
library(tidyverse)
# library(tidyr)
# library(rmarkdown)
library(visreg)
# library(MASS)
library(modEvA)
# library(BiodiversityR)
# library(gridExtra)
library(AICcmodavg)
# library(nlme)
library(mgcv)
library(lme4)
# library(lsmeans)
library(multcomp)
library(MuMIn)

# library(VTrack)
# library(igraph)
# library(MARSS)
library(splitstackshape) ###package to use cSplit and do text to column like excel
library(chron)
library(lubridate)
# library(rgdal)#package for geospatial data
# library(RInSp)#package for intraspecific niche variation
# library(boot)#package for boostrapping operations
library(vegan)

############################################################




##########Set Data#######################################################################

setwd("E:/.....")

###
#Data 1 - Removal by spp | treatment | reef
###
df <- read.csv("Lirman.data.csv", stringsAsFactors = FALSE)
df$Reef <- factor(df$Reef)
df <- dplyr::rename(df, one.week = X..CORALS.REMOVED.1.WEEK,
                    six.months = X..REMOVED.6.MO,
                    prop.bites.1w = X..OF.REMAINING.CORALS.WITH.BITES.1.WEEK)



###################################################################





####################Models for paper############################


###
#First, data wrangling to add time of removal as treatment
###

#data with prop.bites removed and removal period as treatment
df.1.1 <- df %>% 
  dplyr::select(.,-prop.bites.1w) %>% 
  gather(., removal.period, proportion, 5:6)



####Model 1 - Removal probability------------------------

###
#Option 1 - puck data
###

df.1.1b <- df.1.1 %>%  
  filter(., Treatment == "puck")
df.1.1b$SPP <- factor(df.1.1b$SPP)
df.1.1b$Reef <- factor(df.1.1b$Reef)
df.1.1b$removal.period <- factor(df.1.1b$removal.period)

ggplot(df.1.1b, aes(Reef, proportion, colour = SPP, fill = SPP)) +
  geom_bar(stat = "identity") +
  facet_wrap(~removal.period)

glm.removal <- glm(proportion ~ SPP*removal.period + Reef*removal.period, weigh = N.corals, family = binomial (link = "logit"), data = df.1.1b)
glm.removal.2 <- glm(proportion ~ SPP*removal.period + Reef, weigh = N.corals, family = binomial (link = "logit"), data = df.1.1b)
glm.removal.3 <- glm(proportion ~ SPP + Reef*removal.period, weigh = N.corals, family = binomial (link = "logit"), data = df.1.1b)
glm.removal.4 <- glm(proportion ~ SPP + Reef + removal.period, weigh = N.corals, family = binomial (link = "logit"), data = df.1.1b)
glm.removal.5 <- glm(proportion ~ SPP + Reef, weigh = N.corals, family = binomial (link = "logit"), data = df.1.1b)

q1.modTL.list<-list(glm.removal, glm.removal.2, glm.removal.3, glm.removal.4, glm.removal.5)
q1.modTL.sel<-model.sel(q1.modTL.list)
capture.output(q1.modTL.sel, file="LirmanPaper.Output/q1.model.sel.txt")

###
#According AIC from procedure above (model.sel), glm.removal.4 best model
glm.summ.1 <- summary(glm.removal.4)
capture.output(glm.summ.1, file="LirmanPaper.Output/glm.summ.1.txt")

anova.glm1 <- anova(glm.removal.4, test = "Chisq")
capture.output(anova.glm1, file="LirmanPaper.Output/anova.glm1.txt")


visreg(glm.removal.4, "SPP", by = "removal.period", scale = "response")
visreg(glm.removal.4, "SPP", by = "Reef", scale = "response")

Dsquared(glm.removal.4)



visreg(glm.removal.4, "SPP", by = "Reef", scale = "response")

#####################################################################################



####################Working with best model - glm.removal.4########################
library(car)
library(relaimpo) #package to estimate relative importance of regressors in linear models
library(multcomp)

####Diagnostics----------------------------------

###
#residual plots
###

#model 1
residualPlots(glm.removal.4)


###
#Tukey Tests
###

#model 1
spp_coef<-glht(glm.removal.4, mcp(SPP = "Tukey"))$linfct
tidy(cld(glht(glm.removal.4, mcp(SPP = "Tukey"))))

reef_coef<-glht(glm.removal.4, mcp(Reef = "Tukey"))$linfct
tidy(cld(glht(glm.removal.4, mcp(Reef = "Tukey"))))

period_coef<-glht(glm.removal.4, mcp(removal.period = "Tukey"))$linfct
tidy(cld(glht(glm.removal.4, mcp(removal.period = "Tukey"))))

tukey_pairwise_summ.1 <- summary(glht(glm.removal.4, linfct = rbind(spp_coef, reef_coef, period_coef)))

capture.output(tukey_pairwise_summ.1, file="LirmanPaper.Output/tukey_pairwise_summ.1b.txt")



###
#Model Plots
###

#Model 1

new.df1 <- data.frame(expand.grid(c("Mcav","Ofav", "Pcliv", "Pstri"), 
                                  c("1","2", "3"),
                                  c("one.week", "six.months")))

colnames(new.df1) <- c("SPP", "Reef", "removal.period")

dftempo1 <- data.frame(predict(glm.removal.4, newdata = new.df1, se.fit = TRUE, type = "response"))
quantile_norm <- qnorm(p = 0.975)
dftempo1$fit_lower <- dftempo1$fit - quantile_norm * dftempo1$se.fit
dftempo1$fit_upper <- dftempo1$fit + quantile_norm * dftempo1$se.fit
dftempo1$SPP <- new.df1$SPP
dftempo1$Reef <- new.df1$Reef
dftempo1$removal.period <- new.df1$removal.period


#Plotting the model estimates

ggplot(dftempo1, aes(Reef, fit, colour = SPP, fill = SPP)) + 
  geom_crossbar(aes(ymin = fit_lower, ymax = fit_upper), width = 0.8, size = 1, alpha = 0.4, position = position_dodge(width = 0.98)) + 
  theme_bw()+
  facet_wrap(~removal.period)+
  scale_fill_viridis_d(option = "inferno", name = "Spp") + 
  scale_colour_viridis_d(option = "inferno", name = "Spp") +
  labs(x = "Reef", y = "Probability of Removal") + 
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold"), 
        plot.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"))
ggsave("LirmanPaper.Output/Q1.TL.Mod.png", width = 25, height = 15, units = "cm")


