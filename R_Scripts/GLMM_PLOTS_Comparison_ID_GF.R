
library(tidyverse)


####################################################################
# LOADING DATA FOR species ID's (ID_simple)
####################################################################

fitness_final_ID <- read_csv("data_models_phenol_overlap.csv")

fitness_final_ID$Line <- NA

for (i in 1:nrow(fitness_final_ID)){
  if(fitness_final_ID$Plot[i] %in% c(1,2,3)){fitness_final_ID$Line[i] <- 1}
  else if(fitness_final_ID$Plot[i] %in% c(4,5,6)){fitness_final_ID$Line[i] <- 2}
  else{fitness_final_ID$Line[i] <- 3}
} 

fitness_final_ID$Plot <- as.factor(fitness_final_ID$Plot)
fitness_final_ID$Line <- as.factor(fitness_final_ID$Line)
fitness_final_ID$Subplot <- as.factor(fitness_final_ID$Subplot)
fitness_final_ID$ID <- as.factor(fitness_final_ID$ID)
fitness_final_ID$Plant_Simple <- as.factor(fitness_final_ID$Plant_Simple)

str(fitness_final_ID)

fitness_final_ID <- fitness_final_ID %>% dplyr::filter(Seeds_GF>0) %>%
  dplyr::arrange(Plot,Subplot,Plant_Simple,ID)

#############################################
# DATA PLANT SPECIES
##############################################

fitness_LEMA_ID <- fitness_final_ID %>% filter(Plant_Simple=="LEMA")
fitness_CHFU_ID <- fitness_final_ID %>% filter(Plant_Simple=="CHFU")
fitness_PUPA_ID <- fitness_final_ID %>% filter(Plant_Simple=="PUPA")
fitness_ME_ID <- fitness_final_ID %>% filter(Plant_Simple=="ME")
fitness_CHMI_ID <- fitness_final_ID %>% filter(Plant_Simple=="CHMI")

####################################################################
# LOADING DATA FOR species ID's (ID_simple)
####################################################################

fitness_final_GF <- read_csv("data_models_phenol_overlap_GF.csv")

fitness_final_GF$Line <- NA

for (i in 1:nrow(fitness_final_GF)){
  if(fitness_final_GF$Plot[i] %in% c(1,2,3)){fitness_final_GF$Line[i] <- 1}
  else if(fitness_final_GF$Plot[i] %in% c(4,5,6)){fitness_final_GF$Line[i] <- 2}
  else{fitness_final_GF$Line[i] <- 3}
} 

fitness_final_GF$Plot <- as.factor(fitness_final_GF$Plot)
fitness_final_GF$Line <- as.factor(fitness_final_GF$Line)
fitness_final_GF$Subplot <- as.factor(fitness_final_GF$Subplot)
fitness_final_GF$ID <- as.factor(fitness_final_GF$ID)
fitness_final_GF$Plant_Simple <- as.factor(fitness_final_GF$Plant_Simple)

str(fitness_final_GF)

fitness_final_GF <- fitness_final_GF %>% dplyr::filter(Seeds_GF>0) %>%
  dplyr::arrange(Plot,Subplot,Plant_Simple,ID)

#############################################
# DATA PLANT SPECIES (GF)
##############################################

fitness_LEMA_GF <- fitness_final_GF %>% filter(Plant_Simple=="LEMA")
fitness_CHFU_GF <- fitness_final_GF %>% filter(Plant_Simple=="CHFU")
fitness_PUPA_GF <- fitness_final_GF %>% filter(Plant_Simple=="PUPA")
fitness_ME_GF <- fitness_final_GF %>% filter(Plant_Simple=="ME")
fitness_CHMI_GF <- fitness_final_GF %>% filter(Plant_Simple=="CHMI")

###############################################
# MODELS
###############################################

library(MASS)
library(fitdistrplus)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(bbmle) ## for AICtab
library(usdm)
library(DHARMa)

########################################
# NEGATIVE BINOMIAL
########################################

# In glmmTMB
# nbinom1 (also called quasi-poisson) variance = µ * phi
# where µ is the mean and phi is the over-dispersion parameter
# Negative binomial distribution: linear parameterization
# 
# nbinom2 (the default negative binomial in most packages) variance = µ(1+µ/k)
# also written µ + (µ^2)/k Negative binomial distribution: quadratic parameterization


####################
# LEMA
######################

m2.nbinom_LEMA_ID <- glmmTMB(Seeds_GF ~  scale(homo_motif) * scale(DegreeIn) * scale(hete_motif) + (1|ID),
                           ziformula = ~ 1,
                           family = nbinom2(),
                           data = fitness_LEMA_ID)

summary(m2.nbinom_LEMA_ID)

LEMA_ID_comp <- MuMIn::dredge(m2.nbinom_LEMA_ID)


# get residuals ID
simulationOutput <- simulateResiduals(fittedModel = m2.nbinom_LEMA_ID, n = 500)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)
testUniformity(simulationOutput = simulationOutput)
plot(simulationOutput)

plotResiduals(simulationOutput, fitness_LEMA_ID$homo_motif)
plotResiduals(simulationOutput, fitness_LEMA_ID$hete_motif)
plotResiduals(simulationOutput, fitness_LEMA_ID$ID)
levels(fitness_LEMA_ID$ID)
plotResiduals(simulationOutput, fitness_LEMA_ID$DegreeIn)


fitness_LEMA_ID %>% group_by(DegreeIn) %>% count()
fitness_LEMA_ID %>% filter(DegreeIn>5)

m2.nbinom_LEMA_GF <- glmmTMB(Seeds_GF ~  scale(homo_motif) * scale(DegreeIn) * scale(hete_motif) +(1|ID),
                             ziformula= ~ 1,
                             family = nbinom2(),
                             data = fitness_LEMA_GF)

summary(m2.nbinom_LEMA_GF)

LEMA_GF_comp <- MuMIn::dredge(m2.nbinom_LEMA_GF)

# get residuals GF
simulationOutput <- simulateResiduals(fittedModel = m2.nbinom_LEMA_GF, n = 1500)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)
testUniformity(simulationOutput = simulationOutput)
plot(simulationOutput)

plotResiduals(simulationOutput, fitness_LEMA_GF$homo_motif)
plotResiduals(simulationOutput, fitness_LEMA_GF$hete_motif)
plotResiduals(simulationOutput, fitness_LEMA_GF$ID)
levels(fitness_LEMA_GF$ID)
plotResiduals(simulationOutput, fitness_LEMA_GF$DegreeIn)


##########################################
# CHFU
##########################################
m2.nbinom_CHFU_ID <- glmmTMB(Seeds_GF ~  scale(homo_motif) * scale(DegreeIn) * scale(hete_motif) + (1|ID),
                             ziformula = ~ 1,
                             family = nbinom2(),
                             data = fitness_CHFU_ID)

summary(m2.nbinom_CHFU_ID)

CHFU_ID_comp <- MuMIn::dredge(m2.nbinom_CHFU_ID)


# get residuals ID
simulationOutput <- simulateResiduals(fittedModel = m2.nbinom_CHFU_ID, n = 500)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)
testUniformity(simulationOutput = simulationOutput)
plot(simulationOutput)

plotResiduals(simulationOutput, fitness_CHFU_ID$homo_motif)
plotResiduals(simulationOutput, fitness_CHFU_ID$hete_motif)
plotResiduals(simulationOutput, fitness_CHFU_ID$ID)
levels(fitness_CHFU_ID$ID)
plotResiduals(simulationOutput, fitness_CHFU_ID$DegreeIn)


fitness_CHFU_ID %>% group_by(DegreeIn) %>% count()
fitness_CHFU_ID %>% filter(DegreeIn>5)

m2.nbinom_CHFU_GF <- glmmTMB(Seeds_GF ~  scale(homo_motif) * scale(DegreeIn) * scale(hete_motif) +(1|ID),
                             ziformula= ~ 1,
                             family = nbinom2(),
                             data = fitness_CHFU_GF)

summary(m2.nbinom_CHFU_GF)

CHFU_GF_comp <- MuMIn::dredge(m2.nbinom_CHFU_GF)

# get residuals GF
simulationOutput <- simulateResiduals(fittedModel = m2.nbinom_CHFU_GF, n = 1500)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)
testUniformity(simulationOutput = simulationOutput)
plot(simulationOutput)

plotResiduals(simulationOutput, fitness_CHFU_GF$homo_motif)
plotResiduals(simulationOutput, fitness_CHFU_GF$hete_motif)
plotResiduals(simulationOutput, fitness_CHFU_GF$ID)
levels(fitness_CHFU_GF$ID)
plotResiduals(simulationOutput, fitness_CHFU_GF$DegreeIn)


##########################################
# PUPA
##########################################

m2.nbinom_PUPA_ID <- glmmTMB(Seeds_GF ~  scale(homo_motif) * scale(DegreeIn) + (1|ID),
                             ziformula = ~ 1,
                             family = nbinom2(),
                             data = fitness_PUPA_ID)

summary(m2.nbinom_PUPA_ID)

PUPA_ID_comp <- MuMIn::dredge(m2.nbinom_PUPA_ID)


# get residuals ID
simulationOutput <- simulateResiduals(fittedModel = m2.nbinom_PUPA_ID, n = 500)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)
testUniformity(simulationOutput = simulationOutput)
plot(simulationOutput)

plotResiduals(simulationOutput, fitness_PUPA_ID$homo_motif)
plotResiduals(simulationOutput, fitness_PUPA_ID$hete_motif)
plotResiduals(simulationOutput, fitness_PUPA_ID$ID)
levels(fitness_PUPA_ID$ID)
plotResiduals(simulationOutput, fitness_PUPA_ID$DegreeIn)


fitness_PUPA_ID %>% group_by(DegreeIn) %>% count()
fitness_PUPA_ID %>% filter(DegreeIn>5)

m2.nbinom_PUPA_GF <- glmmTMB(Seeds_GF ~  scale(homo_motif) * scale(DegreeIn)  +(1|ID),
                             ziformula= ~ 1,
                             family = nbinom2(),
                             data = fitness_PUPA_GF)

summary(m2.nbinom_PUPA_GF)

PUPA_GF_comp <- MuMIn::dredge(m2.nbinom_PUPA_GF)

# get residuals GF
simulationOutput <- simulateResiduals(fittedModel = m2.nbinom_PUPA_GF, n = 1500)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)
testUniformity(simulationOutput = simulationOutput)
plot(simulationOutput)

plotResiduals(simulationOutput, fitness_PUPA_GF$homo_motif)
plotResiduals(simulationOutput, fitness_PUPA_GF$ID)
levels(fitness_PUPA_GF$ID)
plotResiduals(simulationOutput, fitness_PUPA_GF$DegreeIn)

##########################################
# ME
##########################################

m2.nbinom_ME_ID <- glmmTMB(Seeds_GF ~  scale(homo_motif) * scale(DegreeIn) * scale(hete_motif) + (1|ID),
                             ziformula = ~ 1,
                             family = nbinom2(),
                             data = fitness_ME_ID)

summary(m2.nbinom_ME_ID)

ME_ID_comp <- MuMIn::dredge(m2.nbinom_ME_ID)


# get residuals ID
simulationOutput <- simulateResiduals(fittedModel = m2.nbinom_ME_ID, n = 500)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)
testUniformity(simulationOutput = simulationOutput)
plot(simulationOutput)

plotResiduals(simulationOutput, fitness_ME_ID$homo_motif)
plotResiduals(simulationOutput, fitness_ME_ID$hete_motif)
plotResiduals(simulationOutput, fitness_ME_ID$ID)
levels(fitness_ME_ID$ID)
plotResiduals(simulationOutput, fitness_ME_ID$DegreeIn)


fitness_ME_ID %>% group_by(DegreeIn) %>% count()
fitness_ME_ID %>% filter(DegreeIn>5)

m2.nbinom_ME_GF <- glmmTMB(Seeds_GF ~  scale(homo_motif) * scale(DegreeIn) * scale(hete_motif) +(1|ID),
                             ziformula= ~ 1,
                             family = nbinom2(),
                             data = fitness_ME_GF)

summary(m2.nbinom_ME_GF)

ME_GF_comp <- MuMIn::dredge(m2.nbinom_ME_GF)

# get residuals GF
simulationOutput <- simulateResiduals(fittedModel = m2.nbinom_ME_GF, n = 1500)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)
testUniformity(simulationOutput = simulationOutput)
plot(simulationOutput)

plotResiduals(simulationOutput, fitness_ME_GF$homo_motif)
plotResiduals(simulationOutput, fitness_ME_GF$hete_motif)
plotResiduals(simulationOutput, fitness_ME_GF$ID)
levels(fitness_ME_GF$ID)
plotResiduals(simulationOutput, fitness_ME_GF$DegreeIn)


##########################################
# CHMI
##########################################

# Every pupa has hete_motif = 0

vif(as.data.frame(dplyr::select(fitness_CHMI,StrengthIn,homo_motif,hete_motif,Hub,Katz)))
vif(as.data.frame(dplyr::select(fitness_CHMI,DegreeIn, homo_motif,hete_motif)))
vif(as.data.frame(dplyr::select(fitness_CHMI,StrengthIn, homo_motif,hete_motif)))

fitness_CHMI %>% group_by(Plot) %>% count()


m2.nbinom_CHFI_ZI <- glmmTMB(Seeds_GF ~ (DegreeIn) + (1|Plot/ID),
                     ziformula= ~ 1,
                     family= nbinom2(),
                     data = fitness_CHMI)
summary(m2.nbinom_CHFI_ZI)


m2.nbinom_CHFI <- glmmTMB(Seeds_GF ~ (DegreeIn) + (1|Plot/ID),
                          ziformula= ~ 0,
                          family= nbinom2(),
                          data = fitness_CHMI)

AIC(m2.nbinom_CHFI_ZI,m2.nbinom_CHFI)


# get residuals
simulationOutput <- simulateResiduals(fittedModel = m2.nbinom_CHFI, n = 2500)
testDispersion(simulationOutput)
plot(simulationOutput)
testZeroInflation(simulationOutput)
