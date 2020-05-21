# glm models - DGC

library(tidyverse)
library(MASS)
library(fitdistrplus)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(bbmle) ## for AICtab
library(DHARMa)

# read data ---------------------------------------------------------------

fitness_orig <- read.csv(file = "data_models_phenol_overlap_GF.csv",
                         header = TRUE,
                         stringsAsFactors = FALSE)

# pre-analysis ------------------------------------------------------------

# remove fitness = 0

fitness.data <- subset(fitness_orig,Seeds_GF > 1)

# class of every column
# fitness.data %>% map_chr(class)

fitness_LEMA <- subset(fitness.data,Plant_Simple == "LEMA")
fitness_PUPA <- subset(fitness.data,Plant_Simple == "PUPA")
fitness_CHFU <- subset(fitness.data,Plant_Simple == "CHFU")

# scale sometimes gives problems
summary(scale(fitness_LEMA$homo_motif))
summary(scale(fitness_LEMA$hete_motif))
summary(scale(fitness_LEMA$DegreeIn))

summary(scale(fitness_PUPA$homo_motif))
summary(scale(fitness_PUPA$hete_motif)) # fail 
summary(scale(fitness_PUPA$DegreeIn))

summary(scale(fitness_CHFU$homo_motif))
summary(scale(fitness_CHFU$hete_motif))
summary(scale(fitness_CHFU$DegreeIn))


# models ------------------------------------------------------------------

# negbin - GF no space --------------------------------------------------

GF_LEMA_01 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                           scale(hete_motif) + 
                           scale(DegreeIn) + 
                           (1|ID),
                         family = nbinom2(),
                         data = fitness_LEMA)

GF_PUPA_01 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                           # scale(hete_motif) + 
                           scale(DegreeIn) + 
                           (1|ID),
                         family = nbinom2(),
                         data = fitness_PUPA)

GF_CHFU_01 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                           scale(hete_motif) + 
                           scale(DegreeIn) + 
                           (1|ID),
                         family = nbinom2(),
                         data = fitness_CHFU)

# negbin - GF ----------------------------------------------------------

GF_LEMA_02 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                           scale(hete_motif) + 
                           scale(DegreeIn) + 
                           (1|ID) + (1|Plot),
                      family = nbinom2(),
                      data = fitness_LEMA)

GF_PUPA_02 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                           # scale(hete_motif) + 
                           scale(DegreeIn) + 
                           (1|ID) + (1|Plot),
                      family = nbinom2(),
                      data = fitness_PUPA)

GF_CHFU_02 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                           scale(hete_motif) + 
                           scale(DegreeIn) + 
                           (1|ID) + (1|Plot),
                      family = nbinom2(),
                      data = fitness_CHFU)

# check residuals ---------------------------------------------------------

# get residuals

res_LEMA_01 <- simulateResiduals(fittedModel = GF_LEMA_01, n = 1500)
res_LEMA_02 <- simulateResiduals(fittedModel = GF_LEMA_02, n = 1500)

res_PUPA_01 <- simulateResiduals(fittedModel = GF_PUPA_01, n = 1500)
res_PUPA_02 <- simulateResiduals(fittedModel = GF_PUPA_02, n = 1500)

res_CHFU_01 <- simulateResiduals(fittedModel = GF_CHFU_01, n = 1500)
res_CHFU_02 <- simulateResiduals(fittedModel = GF_CHFU_02, n = 1500)

# Checking Residuals 
testDispersion(res_LEMA_01)
testDispersion(res_LEMA_02)

testDispersion(res_PUPA_01)
testDispersion(res_PUPA_02)

testDispersion(res_CHFU_01)
testDispersion(res_CHFU_02)

# plot residuals
# generally ok, adding plot does not help
# (models 02 are equal or worse)

plot(res_LEMA_01)
plot(res_LEMA_02)

plot(res_PUPA_01)
plot(res_PUPA_02)

plot(res_CHFU_01)
plot(res_CHFU_02)

# more specific plots
plotResiduals(res_LEMA_01, fitness_LEMA$homo_motif)
plotResiduals(res_LEMA_01, fitness_LEMA$hete_motif)
plotResiduals(res_LEMA_01, fitness_LEMA$DegreeIn)
plotResiduals(res_LEMA_01, fitness_LEMA$ID)

# 
plotResiduals(res_PUPA_01, fitness_PUPA$homo_motif)
# plotResiduals(res_PUPA_01, fitness_PUPA$hete_motif)
plotResiduals(res_PUPA_01, fitness_PUPA$DegreeIn)
plotResiduals(res_PUPA_01, fitness_PUPA$ID)

#
plotResiduals(res_CHFU_01, fitness_CHFU$homo_motif)
plotResiduals(res_CHFU_01, fitness_CHFU$hete_motif)
plotResiduals(res_CHFU_01, fitness_CHFU$DegreeIn)
plotResiduals(res_CHFU_01, fitness_CHFU$ID)

# check package GLMMadaptive ----------------------------------------------
# did this because it uses different integration method, 
# apparently more robust to non-gaussian residuals
# in this case, results are similar to glmmTMB

library(GLMMadaptive)

GF_LEMA_01_ad <- mixed_model(Seeds_GF ~ scale(homo_motif) + 
                           scale(hete_motif) + 
                           scale(DegreeIn), random = ~ 1|ID,
                         family = negative.binomial(),
                         data = fitness_LEMA)

GF_CHFU_01_ad <- mixed_model(Seeds_GF ~ scale(homo_motif) + 
                               scale(hete_motif) + 
                               scale(DegreeIn), random = ~ 1|ID,
                             family = negative.binomial(),
                             data = fitness_CHFU)



