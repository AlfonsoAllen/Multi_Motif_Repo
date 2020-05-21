
library(tidyverse)


# read data ---------------------------------------------------------------

fitness_final_aux <- read.csv(file = "data_models_phenol_overlap.csv",
                         header = TRUE,
                         stringsAsFactors = FALSE)

# Add G_F

G_F_list <- read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv") %>%
  dplyr::select(G_F,ID_Simple) %>% rename(ID=ID_Simple) %>% unique()

G_F_list <- bind_rows(G_F_list,tibble(G_F="None",ID="None"))

# Fix "Odontomyia_sp."

G_F_list$G_F[G_F_list$ID=="Odontomyia_sp."] <- "Small_flies"
G_F_list <- unique(G_F_list)

# Sanity check
G_F_list %>% group_by(ID) %>% count() %>% filter(n>1)

fitness_orig <- fitness_final_aux %>% dplyr::left_join(G_F_list,by = "ID")

#write_csv(fitness_final,"data_models_phenol_overlap_Random_GF.csv")

# Turn ID, GF and Plot into factors
fitness_orig$Plot <- as.factor(fitness_orig$Plot)
fitness_orig$ID <- as.factor(fitness_orig$ID)
fitness_orig$G_F <- as.factor(fitness_orig$G_F)

library(tidyverse)
library(MASS)
library(fitdistrplus)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(bbmle) ## for AICtab
library(DHARMa)


# pre-analysis ------------------------------------------------------------

# remove fitness = 0

fitness.data <- subset(fitness_orig,Seeds_GF > 0)

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

GF_all_sp <- glmmTMB(Seeds_GF ~ scale(homo_motif) + #Convergence problem
                       scale(hete_motif) + 
                       scale(DegreeIn) + 
                       (1 | Plant_Simple) + 
                       (1 | G_F)+ 
                       (1 | Plot),
                     family = nbinom2(),
                     data = fitness.data)


summary(GF_all_sp)

GF_all_sp_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                       scale(hete_motif) + 
                       scale(DegreeIn) + 
                       (1 | Plant_Simple) + 
                       (1 | G_F)+ 
                       (1 | Plot),
                      ziformula = ~1,
                     family = nbinom2(),
                     data = fitness.data)


summary(GF_all_sp_ZI)

GF_all_sp_GF <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                       scale(hete_motif) + 
                       scale(DegreeIn) + 
                       (1 | Plant_Simple) + 
                       (1 | G_F),
                     family = nbinom2(),
                     data = fitness.data)


summary(GF_all_sp_GF)

GF_all_sp_Plot <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                            scale(hete_motif) + 
                            scale(DegreeIn) + 
                            (1 | Plant_Simple) + 
                            (1 | Plot),
                          family = nbinom2(),
                          data = fitness.data)


summary(GF_all_sp_Plot)

GF_all_sp_Plant <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                            scale(hete_motif) + 
                            scale(DegreeIn) + 
                            (1 | Plant_Simple),
                          family = nbinom2(),
                          data = fitness.data)


summary(GF_all_sp_Plant)



# check residuals ---------------------------------------------------------

# get residuals

res_all <- simulateResiduals(fittedModel = GF_all_sp, n = 1500)
res_all_ZI <- simulateResiduals(fittedModel = GF_all_sp_ZI, n = 1500)
res_all_plot <- simulateResiduals(fittedModel = GF_all_sp_Plot, n = 1500)
res_all_GF <- simulateResiduals(fittedModel = GF_all_sp_GF, n = 1500)
res_all_plant <- simulateResiduals(fittedModel = GF_all_sp_Plant, n = 1500)

# Checking Residuals 
testDispersion(res_all)
testDispersion(res_all_ZI)
testDispersion(res_all_plot)
testDispersion(res_all_GF)
testDispersion(res_all_plant)

# plot residuals
# results show dispersion and heter.

plot(res_all) #KS + heter.
plot(res_all_ZI)#KS + heter.
plot(res_all_plot)#KS + heter.
plot(res_all_GF) #KS + heter.
plot(res_all_plant)#KS + heter.



# check differences between models with and without ZI, when seed > 0 -----
AIC(GF_all_sp,GF_all_sp_ZI,GF_all_sp_Plot,GF_all_sp_GF,GF_all_sp_Plant) 





#########################################
# INDIVIDUAL MODEL FOR EACH PLANT SPECIES
#########################################



# negbin - GF no space --------------------------------------------------

GF_LEMA_01 <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                        scale(hete_motif) + 
                        scale(DegreeIn) + 
                        (1|G_F),
                      family = nbinom2(),
                      data = fitness_LEMA)

GF_LEMA_02 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + # Convergence problems
                        scale(hete_motif) + 
                        scale(DegreeIn) + 
                        (1 + scale(homo_motif) + scale(hete_motif)|G_F),
                      family = nbinom2(),
                      data = fitness_LEMA)

GF_PUPA_01 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                        # scale(hete_motif) + 
                        scale(DegreeIn) + 
                        (1|G_F),
                      family = nbinom2(),
                      data = fitness_PUPA)

GF_CHFU_01 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                        scale(hete_motif) + 
                        scale(DegreeIn) + 
                        (1|G_F),
                      family = nbinom2(),
                      data = fitness_CHFU)

# negbin - ZI GF no space --------------------------------------------------

GF_LEMA_01_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) +  #Convergence problems without degree
                           scale(hete_motif) + 
                           #scale(DegreeIn) + 
                           (1|G_F),
                         ziformula = ~1,
                         family = nbinom2(),
                         data = fitness_LEMA)

GF_PUPA_01_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) + #Convergence problems with ID
                           # scale(hete_motif) + 
                           scale(DegreeIn) + 
                           (1|G_F),
                         ziformula = ~1,
                         family = nbinom2(),
                         data = fitness_PUPA)

GF_CHFU_01_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) +  # Convergence problems with G_F
                           scale(hete_motif) + 
                           scale(DegreeIn) + 
                           (1|G_F),
                         ziformula = ~1,
                         family = nbinom2(),
                         data = fitness_CHFU)



# negbin - GF ----------------------------------------------------------

GF_LEMA_02 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                        scale(hete_motif) + 
                        scale(DegreeIn) + 
                        (1|G_F) + (1|Plot),
                      #ziformula = ~1,
                      family = nbinom2(),
                      data = fitness_LEMA)

GF_PUPA_02 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + # Convergence problem without ZI 
                        # scale(hete_motif) + 
                        scale(DegreeIn) + 
                        (1|G_F) + (1|Plot),
                      #ziformula = ~1,
                      family = nbinom2(),
                      data = fitness_PUPA)

GF_CHFU_02 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + # Convergence problem with ZI
                        scale(hete_motif) + 
                        scale(DegreeIn) + 
                        (1|G_F) + (1|Plot),
                      #ziformula = ~1,
                      family = nbinom2(),
                      data = fitness_CHFU)


# negbin - space no GF ----------------------------------------------------------

GF_LEMA_03 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                        scale(hete_motif) + 
                        scale(DegreeIn) + 
                        (1|Plot),
                      #ziformula = ~1,
                      family = nbinom2(),
                      data = fitness_LEMA)

GF_PUPA_03 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + # Convergence problem without ZI 
                        # scale(hete_motif) + 
                        scale(DegreeIn) + 
                        (1|Plot),
                      #ziformula = ~1,
                      family = nbinom2(),
                      data = fitness_PUPA)

GF_CHFU_03 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + # Convergence problem with ZI
                        scale(hete_motif) + 
                        scale(DegreeIn) + 
                        (1|Plot),
                      #ziformula = ~1,
                      family = nbinom2(),
                      data = fitness_CHFU)

# check residuals ---------------------------------------------------------

# get residuals

res_LEMA_01 <- simulateResiduals(fittedModel = GF_LEMA_01, n = 1500)
res_LEMA_02 <- simulateResiduals(fittedModel = GF_LEMA_01_ZI, n = 1500)
res_LEMA_both <- simulateResiduals(fittedModel = GF_LEMA_02, n = 1500)
res_LEMA_03 <- simulateResiduals(fittedModel = GF_LEMA_03, n = 1500)

res_PUPA_01 <- simulateResiduals(fittedModel = GF_PUPA_01, n = 1500)
res_PUPA_02 <- simulateResiduals(fittedModel = GF_PUPA_01_ZI, n = 1500)
res_PUPA_both <- simulateResiduals(fittedModel = GF_PUPA_02, n = 1500)
res_PUPA_03 <- simulateResiduals(fittedModel = GF_PUPA_03, n = 1500)


res_CHFU_01 <- simulateResiduals(fittedModel = GF_CHFU_01, n = 1500)
res_CHFU_02 <- simulateResiduals(fittedModel = GF_CHFU_01_ZI, n = 1500)
res_CHFU_both <- simulateResiduals(fittedModel = GF_CHFU_02, n = 1500)
res_CHFU_03 <- simulateResiduals(fittedModel = GF_CHFU_03, n = 1500)

# Checking Residuals 
testDispersion(res_LEMA_01)
testDispersion(res_LEMA_02)
testDispersion(res_LEMA_both)
testDispersion(res_LEMA_03)

testDispersion(res_PUPA_01)
testDispersion(res_PUPA_02)
testDispersion(res_PUPA_both)
testDispersion(res_PUPA_03)

testDispersion(res_CHFU_01)
testDispersion(res_CHFU_02)
testDispersion(res_CHFU_both)
testDispersion(res_CHFU_03)
# plot residuals
# homo+hete by ID and G_F as random factor worsen redisuals
# PUPA develops heter.
# results for LEMA improve adding Plot as cross effect

# Without degree hete. gets wild

plot(res_LEMA_01) #KS + heter.
plot(res_LEMA_02)#KS + heter.
plot(res_LEMA_both)#KS
plot(res_LEMA_03) #KS

plot(res_PUPA_01) #heter.
plot(res_PUPA_02) #heter.
plot(res_PUPA_both) #heter.
plot(res_PUPA_03) #heter.

plot(res_CHFU_01) #KS
plot(res_CHFU_02) #KS + heter.
plot(res_CHFU_both) #KS + High heter.
plot(res_CHFU_03)#KS + High heter.

# more specific plots
plotResiduals(res_LEMA_01, fitness_LEMA$homo_motif)
plotResiduals(res_LEMA_01, fitness_LEMA$hete_motif)
plotResiduals(res_LEMA_01, fitness_LEMA$DegreeIn)
plotResiduals(res_LEMA_01, fitness_LEMA$G_F)

# 
plotResiduals(res_PUPA_01, fitness_PUPA$homo_motif)
# plotResiduals(res_PUPA_01, fitness_PUPA$hete_motif)
plotResiduals(res_PUPA_01, fitness_PUPA$DegreeIn)
plotResiduals(res_PUPA_01, fitness_PUPA$G_F)

#
plotResiduals(res_CHFU_01, fitness_CHFU$homo_motif)
plotResiduals(res_CHFU_01, fitness_CHFU$hete_motif)
plotResiduals(res_CHFU_01, fitness_CHFU$DegreeIn)
plotResiduals(res_CHFU_01, fitness_CHFU$G_F)


# check differences between models with and without ZI, when seed > 0 -----
AIC(GF_LEMA_01,GF_LEMA_01_ZI, GF_LEMA_02, GF_LEMA_03) 
AIC(GF_CHFU_01,GF_CHFU_01_ZI, GF_CHFU_02, GF_CHFU_03)
AIC(GF_PUPA_01,GF_PUPA_01_ZI, GF_PUPA_02, GF_PUPA_03)

