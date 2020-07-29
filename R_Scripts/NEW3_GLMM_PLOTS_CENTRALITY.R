
library(tidyverse)


# read data ---------------------------------------------------------------

fitness_final_aux0 <- read.csv(file = "data_models_phenol_overlap.csv",
                         header = TRUE,
                         stringsAsFactors = FALSE)

PageRank_results <- read.csv(file = "PageRank_results.csv",
                              header = TRUE,
                              stringsAsFactors = FALSE)

PageRank_results1 <- PageRank_results %>% separate(Label,c("Subplot","Plant0")," ") %>% dplyr::select(-Plant0) 

# NOTE THAT PREVIOUS CENTRALITY MEASURES WERE WRONG BECAUSE THEY WEREN'T NORMALIZED

fitness_final_aux <- fitness_final_aux0 %>% 
  left_join(PageRank_results1,by=c("Plot","Subplot","Plant_Simple"))

fitness_final_aux[is.na(fitness_final_aux)] <- 0
fitness_final_aux$Ratio[fitness_final_aux$ID=="None"] <- 1
fitness_final_aux$pecentage_same_plants[fitness_final_aux$ID=="None"] <- 1

fitness_final_aux %>% filter(ID!="None") %>% count()
fitness_final_aux %>% filter(ID!="None") %>% count(wt=visits_GF)
fitness_final_aux %>% filter(Real_PR_Multi>0) %>% count()

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

fitness.data %>% count(wt=visits_GF)

# class of every column
# fitness.data %>% map_chr(class)

fitness_LEMA <- subset(fitness.data,Plant_Simple == "LEMA")
fitness_PUPA <- subset(fitness.data,Plant_Simple == "PUPA")
fitness_CHFU <- subset(fitness.data,Plant_Simple == "CHFU")


# Data of LEMA plants with R > 1

# x <- fitness_LEMA %>% dplyr::select(Plot,Subplot,G_F,Plant_Simple,Ratio) %>% filter(Ratio>1)
# write_csv(x,"data_LEMA_PLANTS_R.csv")

# scale sometimes gives problems
summary(scale(fitness_LEMA$homo_motif))
summary(scale(fitness_LEMA$hete_motif))
summary(scale(fitness_LEMA$DegreeIn))
summary(scale(fitness_LEMA$Real_PR_Layer))
summary(scale(fitness_LEMA$Delta))
summary(scale(fitness_LEMA$Ratio))
summary(scale(fitness_LEMA$pecentage_same_plants))

summary(scale(fitness_PUPA$homo_motif))
summary(scale(fitness_PUPA$hete_motif)) # fail
summary(scale(fitness_PUPA$DegreeIn))
summary(scale(fitness_PUPA$Real_PR_Layer))
summary(scale(fitness_PUPA$Delta))
summary(scale(fitness_PUPA$Ratio))
summary(scale(fitness_PUPA$pecentage_same_plants)) # fail

summary(scale(fitness_CHFU$homo_motif))
summary(scale(fitness_CHFU$hete_motif))
summary(scale(fitness_CHFU$DegreeIn))
summary(scale(fitness_CHFU$Real_PR_Layer))
summary(scale(fitness_CHFU$Ratio))
summary(scale(fitness_CHFU$Delta))
summary(scale(fitness_CHFU$pecentage_same_plants))


# models ------------------------------------------------------------------

GF_all_sp <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                       scale(hete_motif) + 
                       scale(StrengthIn) + scale(Delta) + 
                       (1 | Plant_Simple) + 
                       (1 | G_F)+ 
                       (1 | Plot),
                     family = nbinom2(),
                     data = fitness.data)


summary(GF_all_sp)

GF_all_sp_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                       scale(hete_motif) + 
                        scale(StrengthIn) + scale(Delta) +
                       (1 | Plant_Simple) + 
                       (1 | G_F)+ 
                       (1 | Plot),
                      ziformula = ~1,
                     family = nbinom2(),
                     data = fitness.data)


summary(GF_all_sp_ZI)

GF_all_sp_GF <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                       scale(hete_motif) + 
                      scale(StrengthIn) + scale(Delta) + 
                       (1 | Plant_Simple) + 
                       (1 | G_F),
                     family = nbinom2(),
                     data = fitness.data)


summary(GF_all_sp_GF)

GF_all_sp_Plot <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                            scale(hete_motif) + 
                            scale(StrengthIn) + scale(Delta) + 
                            (1 | Plant_Simple) + 
                            (1 | Plot),
                          family = nbinom2(),
                          data = fitness.data)


summary(GF_all_sp_Plot)

GF_all_sp_Plant <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                            scale(hete_motif) + 
                             scale(StrengthIn) + scale(Delta) + 
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


library(usdm)

vif(as.data.frame(dplyr::select(fitness_LEMA,StrengthIn,homo_motif,hete_motif,Ratio)))
vif(as.data.frame(dplyr::select(fitness_PUPA,StrengthIn,homo_motif,Ratio)))
vif(as.data.frame(dplyr::select(fitness_CHFU,StrengthIn,homo_motif,hete_motif,Ratio)))

vif(as.data.frame(dplyr::select(fitness_LEMA,StrengthIn,Real_PR_Layer,DegreeIn,Ratio)))
vif(as.data.frame(dplyr::select(fitness_LEMA,DegreeIn,homo_motif,hete_motif,individuals)))
vif(as.data.frame(dplyr::select(fitness_LEMA,PageRank,homo_motif,hete_motif,individuals)))


# negbin - GF no space --------------------------------------------------

GF_LEMA_01 <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                        scale(hete_motif) + 
                        scale(StrengthIn) + scale(Delta) +
                        (1|G_F),
                      family = nbinom2(),
                      data = fitness_LEMA)

GF_LEMA_02 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + # Convergence problems
                        scale(hete_motif) + 
                        scale(StrengthIn) + scale(Delta) +
                        (1 + scale(homo_motif) + scale(hete_motif)|G_F),
                      family = nbinom2(),
                      data = fitness_LEMA)

GF_PUPA_01 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                        # scale(hete_motif) + 
                        scale(StrengthIn) + scale(Delta) + 
                        (1|G_F),
                      family = nbinom1(),
                      data = fitness_PUPA)

GF_CHFU_01 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                        scale(hete_motif) + 
                        scale(StrengthIn) + scale(Delta) +
                        (1|G_F),
                      family = nbinom1(),
                      data = fitness_CHFU)

# negbin - ZI GF no space --------------------------------------------------

GF_LEMA_01_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) +  #Convergence problems without degree
                           scale(hete_motif) + 
                           scale(StrengthIn) + scale(Delta) + 
                           (1|G_F),
                         ziformula = ~1,
                         family = nbinom2(),
                         data = fitness_LEMA)

GF_PUPA_01_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) + #Convergence problems with ID
                           # scale(hete_motif) + 
                           scale(StrengthIn) + scale(Delta) +
                           (1|G_F),
                         ziformula = ~1,
                         family = nbinom1(),
                         data = fitness_PUPA)

GF_CHFU_01_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) +  
                           scale(hete_motif) + 
                           scale(StrengthIn) + scale(Delta) + 
                           (1|G_F),
                         ziformula = ~1,
                         family = nbinom1(),
                         data = fitness_CHFU)



# negbin - GF ----------------------------------------------------------

GF_LEMA_02 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                        scale(hete_motif) + 
                        scale(StrengthIn) + scale(Delta) + 
                        (1|G_F) + (1|Plot),
                      #ziformula = ~1,
                      family = nbinom2(),
                      data = fitness_LEMA)

GF_PUPA_02 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + # Convergence problem without ZI 
                        # scale(hete_motif) + 
                        scale(StrengthIn) + scale(Delta) + 
                        (1|G_F) + (1|Plot),
                      #ziformula = ~1,
                      family = nbinom1(),
                      data = fitness_PUPA)

GF_CHFU_02 <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                        scale(hete_motif) + 
                        scale(StrengthIn) + scale(Delta) + 
                        (1|G_F) + (1|Plot),
                      #ziformula = ~1,
                      family = nbinom1(),
                      data = fitness_CHFU)


# negbin - space no GF ----------------------------------------------------------

GF_LEMA_03 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                        scale(hete_motif) + 
                        scale(StrengthIn) + scale(Delta) + 
                        (1|Plot),
                      #ziformula = ~1,
                      family = nbinom2(),
                      data = fitness_LEMA)

GF_PUPA_03 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                        # scale(hete_motif) + 
                        scale(StrengthIn) + scale(Delta) + 
                        (1|Plot),
                      #ziformula = ~1,
                      family = nbinom1(),
                      data = fitness_PUPA)

GF_CHFU_03 <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                        scale(hete_motif) + 
                        scale(StrengthIn) + scale(Delta) +
                        (1|Plot),
                      #ziformula = ~1,
                      family = nbinom1(),
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

testZeroInflation(res_LEMA_01)
testZeroInflation(res_LEMA_02)
testZeroInflation(res_LEMA_both)
testZeroInflation(res_LEMA_03)

testZeroInflation(res_PUPA_01)
testZeroInflation(res_PUPA_02)
testZeroInflation(res_PUPA_both)
testZeroInflation(res_PUPA_03)

testZeroInflation(res_CHFU_01)
testZeroInflation(res_CHFU_02)
testZeroInflation(res_CHFU_both)
testZeroInflation(res_CHFU_03)

# plot residuals
# homo+hete by ID and G_F as random factor worsen redisuals
# PUPA develops heter.
# results for LEMA improve adding Plot as cross effect

# Without degree hete. gets wild

plot(res_LEMA_01) #heter.
plot(res_LEMA_02)#KS + heter.
plot(res_LEMA_both)#KS+ heter
plot(res_LEMA_03) #KS+ heter

plot(res_PUPA_01) 
plot(res_PUPA_02)
plot(res_PUPA_both)
plot(res_PUPA_03)

plot(res_CHFU_01) #KS 
plot(res_CHFU_02) #KS + heter.
plot(res_CHFU_both) #KS + High heter.
plot(res_CHFU_03)#KS + High heter.

fitness.data %>% group_by(G_F) %>% count()

# more specific plots
plotResiduals(res_LEMA_01, scale(fitness_LEMA$homo_motif))
plotResiduals(res_LEMA_01, scale(fitness_LEMA$hete_motif))
plotResiduals(res_LEMA_01, scale(fitness_LEMA$StrengthIn))
plotResiduals(res_LEMA_01, scale(fitness_LEMA$Ratio))
plotResiduals(res_LEMA_01, scale(fitness_LEMA$Delta))

# 
plotResiduals(res_PUPA_02, scale(fitness_PUPA$homo_motif))
# plotResiduals(res_PUPA_01, fitness_PUPA$hete_motif)
plotResiduals(res_PUPA_02, scale(fitness_PUPA$StrengthIn))
plotResiduals(res_PUPA_02, scale(fitness_PUPA$Ratio))
plotResiduals(res_PUPA_02, scale(fitness_PUPA$Delta))
#
plotResiduals(res_CHFU_02, scale(fitness_CHFU$homo_motif))
plotResiduals(res_CHFU_02, scale(fitness_CHFU$hete_motif))
plotResiduals(res_CHFU_02, scale(fitness_CHFU$StrengthIn))
plotResiduals(res_CHFU_02, scale(fitness_CHFU$Ratio))
plotResiduals(res_CHFU_02, scale(fitness_CHFU$Delta))



# check differences between models with and without ZI, when seed > 0 -----
AIC(GF_LEMA_01,GF_LEMA_01_ZI, GF_LEMA_02, GF_LEMA_03) 
AIC(GF_CHFU_01,GF_CHFU_01_ZI, GF_CHFU_02, GF_CHFU_03)
AIC(GF_PUPA_01,GF_PUPA_01_ZI, GF_PUPA_02, GF_PUPA_03)


library(performance)
r2(GF_LEMA_01)
r2(GF_PUPA_01)
r2(GF_CHFU_01)

summary(GF_LEMA_01)
summary(GF_PUPA_01_ZI)
summary(GF_CHFU_01_ZI)

ggplot(fitness_LEMA)+
  geom_point(aes(x=scale(Delta),y=Seeds_GF))

fitness_LEMA2 <- fitness_LEMA
fitness_LEMA$Ratio <- scale(fitness_LEMA$Ratio)
library(ggeffects)
ggpredict(GF_LEMA_01, c("") ) %>% plot(rawdata = TRUE, jitter = .01)


glmmTMB(Seeds_GF ~ scale(homo_motif) + 
          # scale(hete_motif) + 
          scale(StrengthIn) + scale(Delta) + 
          (1|Plot),
        #ziformula = ~1,
        family = nbinom1(),
        data = fitness_PUPA)


fitness_CHFU %>% group_by(G_F) %>% count()
fitness_LEMA %>% group_by(G_F) %>% count()
fitness_PUPA %>% group_by(G_F) %>% count()








############################################
# WITHOUT SCALING + CHFU binomial2
############################################

# negbin - GF no space --------------------------------------------------

GF_LEMA_01 <- glmmTMB(Seeds_GF ~ (homo_motif) +
                        (hete_motif) + 
                        (StrengthIn) + (Delta) +
                        (1|G_F),
                      family = nbinom2(),
                      data = fitness_LEMA)

GF_LEMA_02 <- glmmTMB(Seeds_GF ~ (homo_motif) + # Convergence problems
                        (hete_motif) + 
                        (StrengthIn) + (Delta) +
                        (1 + (homo_motif) + (hete_motif)|G_F),
                      family = nbinom2(),
                      data = fitness_LEMA)

GF_PUPA_01 <- glmmTMB(Seeds_GF ~ (homo_motif) + 
                        # scale(hete_motif) + 
                        (StrengthIn) + (Delta) + 
                        (1|G_F),
                      family = nbinom1(),
                      data = fitness_PUPA)

GF_CHFU_01 <- glmmTMB(Seeds_GF ~ (homo_motif) + 
                        (hete_motif) + 
                        (StrengthIn) + (Delta) +
                        (1|G_F),
                      family = nbinom2(),
                      data = fitness_CHFU)

# negbin - ZI GF no space --------------------------------------------------

GF_LEMA_01_ZI <- glmmTMB(Seeds_GF ~ (homo_motif) +  
                           (hete_motif) + 
                           (StrengthIn) + (Delta) + 
                           (1|G_F),
                         ziformula = ~1,
                         family = nbinom2(),
                         data = fitness_LEMA)

GF_PUPA_01_ZI <- glmmTMB(Seeds_GF ~ (homo_motif) + #Convergence problems with ID
                           # scale(hete_motif) + 
                           (StrengthIn) +(Delta) +
                           (1|G_F),
                         ziformula = ~1,
                         family = nbinom1(),
                         data = fitness_PUPA)

GF_CHFU_01_ZI <- glmmTMB(Seeds_GF ~ (homo_motif) +  
                           (hete_motif) + 
                           (StrengthIn) + (Delta) + 
                           (1|G_F),
                         ziformula = ~1,
                         family = nbinom2(),
                         data = fitness_CHFU)



# negbin - GF ----------------------------------------------------------

GF_LEMA_02 <- glmmTMB(Seeds_GF ~ (homo_motif) + 
                        (hete_motif) + 
                        (StrengthIn) + (Delta) + 
                        (1|G_F) + (1|Plot),
                      #ziformula = ~1,
                      family = nbinom2(),
                      data = fitness_LEMA)

GF_PUPA_02 <- glmmTMB(Seeds_GF ~ (homo_motif) + 
                        # scale(hete_motif) + 
                        (StrengthIn) + (Delta) + 
                        (1|G_F) + (1|Plot),
                      #ziformula = ~1,
                      family = nbinom1(),
                      data = fitness_PUPA)

GF_CHFU_02 <- glmmTMB(Seeds_GF ~ (homo_motif) +
                        (hete_motif) + 
                        (StrengthIn) + (Delta) + 
                        (1|G_F) + (1|Plot),
                      #ziformula = ~1,
                      family = nbinom2(),
                      data = fitness_CHFU)


# negbin - space no GF ----------------------------------------------------------

GF_LEMA_03 <- glmmTMB(Seeds_GF ~ (homo_motif) + 
                        (hete_motif) + 
                        (StrengthIn) + (Delta) + 
                        (1|Plot),
                      #ziformula = ~1,
                      family = nbinom2(),
                      data = fitness_LEMA)

GF_PUPA_03 <- glmmTMB(Seeds_GF ~ (homo_motif) + 
                        # scale(hete_motif) + 
                        (StrengthIn) + (Delta) + 
                        (1|Plot),
                      #ziformula = ~1,
                      family = nbinom1(),
                      data = fitness_PUPA)

GF_CHFU_03 <- glmmTMB(Seeds_GF ~ (homo_motif) +
                        (hete_motif) + 
                        (StrengthIn) + (Delta) +
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

testZeroInflation(res_LEMA_01)
testZeroInflation(res_LEMA_02)
testZeroInflation(res_LEMA_both)
testZeroInflation(res_LEMA_03)

testZeroInflation(res_PUPA_01)
testZeroInflation(res_PUPA_02)
testZeroInflation(res_PUPA_both)
testZeroInflation(res_PUPA_03)

testZeroInflation(res_CHFU_01)
testZeroInflation(res_CHFU_02)
testZeroInflation(res_CHFU_both)
testZeroInflation(res_CHFU_03)

# plot residuals
# homo+hete by ID and G_F as random factor worsen redisuals
# PUPA develops heter.
# results for LEMA improve adding Plot as cross effect

# Without degree hete. gets wild

plot(res_LEMA_01) #heter.
plot(res_LEMA_02)#KS + heter.
plot(res_LEMA_both)#KS+ heter
plot(res_LEMA_03) #KS+ heter

plot(res_PUPA_01) 
plot(res_PUPA_02)
plot(res_PUPA_both)
plot(res_PUPA_03)

plot(res_CHFU_01) #KS 
plot(res_CHFU_02) #KS + heter.
plot(res_CHFU_both) #KS + High heter.
plot(res_CHFU_03)#KS + High heter.

fitness.data %>% group_by(G_F) %>% count()

# more specific plots
plotResiduals(res_LEMA_01, (fitness_LEMA$homo_motif))
plotResiduals(res_LEMA_01, (fitness_LEMA$hete_motif))
plotResiduals(res_LEMA_01, (fitness_LEMA$StrengthIn))
plotResiduals(res_LEMA_01, (fitness_LEMA$Ratio))
plotResiduals(res_LEMA_01, (fitness_LEMA$Delta))

# 
plotResiduals(res_PUPA_02, scale(fitness_PUPA$homo_motif))
# plotResiduals(res_PUPA_01, fitness_PUPA$hete_motif)
plotResiduals(res_PUPA_02, scale(fitness_PUPA$StrengthIn))
plotResiduals(res_PUPA_02, scale(fitness_PUPA$Ratio))
plotResiduals(res_PUPA_02, scale(fitness_PUPA$Delta))
#
plotResiduals(res_CHFU_02, (fitness_CHFU$homo_motif))
plotResiduals(res_CHFU_02, (fitness_CHFU$hete_motif))
plotResiduals(res_CHFU_02, (fitness_CHFU$StrengthIn))
plotResiduals(res_CHFU_02, (fitness_CHFU$Ratio))
plotResiduals(res_CHFU_02, (fitness_CHFU$Delta))



# check differences between models with and without ZI, when seed > 0 -----
AIC(GF_LEMA_01,GF_LEMA_01_ZI, GF_LEMA_02, GF_LEMA_03) 
AIC(GF_CHFU_01,GF_CHFU_01_ZI, GF_CHFU_02, GF_CHFU_03)
AIC(GF_PUPA_01,GF_PUPA_01_ZI, GF_PUPA_02, GF_PUPA_03)


library(performance)
r2(GF_LEMA_01)
r2(GF_PUPA_01)
r2(GF_CHFU_01)

summary(GF_LEMA_01)
summary(GF_PUPA_01_ZI)
summary(GF_CHFU_01_ZI)


############################################
# WITHOUT SCALING + CHFU binomial2 +log SEEDS + gaussian
############################################

# negbin - GF no space --------------------------------------------------

GF_LEMA_01 <- glmmTMB(log(Seeds_GF) ~ (homo_motif) +
                        (hete_motif) + 
                        (StrengthIn) + (Delta) +
                        (1|G_F),
                      family = gaussian(),
                      data = fitness_LEMA)

GF_LEMA_02 <- glmmTMB(log(Seeds_GF) ~ (homo_motif) + # Convergence problems
                        (hete_motif) + 
                        (StrengthIn) + (Delta) +
                        (1 + (homo_motif) + (hete_motif)|G_F),
                      family = gaussian(),
                      data = fitness_LEMA)

GF_PUPA_01 <- glmmTMB(log(Seeds_GF) ~ (homo_motif) + 
                        # scale(hete_motif) + 
                        (StrengthIn) + (Delta) + 
                        (1|G_F),
                      family = gaussian(),
                      data = fitness_PUPA)

GF_CHFU_01 <- glmmTMB(log(Seeds_GF) ~ (homo_motif) + 
                        (hete_motif) + 
                        (StrengthIn) + (Delta) +
                        (1|G_F),
                      family = gaussian(),
                      data = fitness_CHFU)

# negbin - ZI GF no space --------------------------------------------------

GF_LEMA_01_ZI <- glmmTMB(log(Seeds_GF) ~ (homo_motif) +  
                           (hete_motif) + 
                           (StrengthIn) + (Delta) + 
                           (1|G_F),
                         ziformula = ~1,
                         family = gaussian(),
                         data = fitness_LEMA)

GF_PUPA_01_ZI <- glmmTMB(log(Seeds_GF) ~ (homo_motif) + #Convergence problems with ID
                           # scale(hete_motif) + 
                           (StrengthIn) +(Delta) +
                           (1|G_F),
                         ziformula = ~1,
                         family = gaussian(),
                         data = fitness_PUPA)

GF_CHFU_01_ZI <- glmmTMB(log(Seeds_GF) ~ (homo_motif) +  
                           (hete_motif) + 
                           (StrengthIn) + (Delta) + 
                           (1|G_F),
                         ziformula = ~1,
                         family = gaussian(),
                         data = fitness_CHFU)



# negbin - GF ----------------------------------------------------------

GF_LEMA_02 <- glmmTMB(log(Seeds_GF) ~ (homo_motif) + 
                        (hete_motif) + 
                        (StrengthIn) + (Delta) + 
                        (1|G_F) + (1|Plot),
                      #ziformula = ~1,
                      family = gaussian(),
                      data = fitness_LEMA)

GF_PUPA_02 <- glmmTMB(log(Seeds_GF) ~ (homo_motif) + 
                        # scale(hete_motif) + 
                        (StrengthIn) + (Delta) + 
                        (1|G_F) + (1|Plot),
                      #ziformula = ~1,
                      family = gaussian(),
                      data = fitness_PUPA)

GF_CHFU_02 <- glmmTMB(log(Seeds_GF) ~ (homo_motif) +
                        (hete_motif) + 
                        (StrengthIn) + (Delta) + 
                        (1|G_F) + (1|Plot),
                      #ziformula = ~1,
                      family = gaussian(),
                      data = fitness_CHFU)


# negbin - space no GF ----------------------------------------------------------

GF_LEMA_03 <- glmmTMB(log(Seeds_GF) ~ (homo_motif) + 
                        (hete_motif) + 
                        (StrengthIn) + (Delta) + 
                        (1|Plot),
                      #ziformula = ~1,
                      family = gaussian(),
                      data = fitness_LEMA)

GF_PUPA_03 <- glmmTMB(log(Seeds_GF) ~ (homo_motif) + 
                        # scale(hete_motif) + 
                        (StrengthIn) + (Delta) + 
                        (1|Plot),
                      #ziformula = ~1,
                      family = gaussian(),
                      data = fitness_PUPA)

GF_CHFU_03 <- glmmTMB(log(Seeds_GF) ~ (homo_motif) +
                        (hete_motif) + 
                        (StrengthIn) + (Delta) +
                        (1|Plot),
                      #ziformula = ~1,
                      family = gaussian(),
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

testZeroInflation(res_LEMA_01)
testZeroInflation(res_LEMA_02)
testZeroInflation(res_LEMA_both)
testZeroInflation(res_LEMA_03)

testZeroInflation(res_PUPA_01)
testZeroInflation(res_PUPA_02)
testZeroInflation(res_PUPA_both)
testZeroInflation(res_PUPA_03)

testZeroInflation(res_CHFU_01)
testZeroInflation(res_CHFU_02)
testZeroInflation(res_CHFU_both)
testZeroInflation(res_CHFU_03)

# plot residuals
# homo+hete by ID and G_F as random factor worsen redisuals
# PUPA develops heter.
# results for LEMA improve adding Plot as cross effect

# Without degree hete. gets wild

plot(res_LEMA_01) #heter.
plot(res_LEMA_02)#KS + heter.
plot(res_LEMA_both)#KS+ heter
plot(res_LEMA_03) #KS+ heter

plot(res_PUPA_01) 
plot(res_PUPA_02)
plot(res_PUPA_both)
plot(res_PUPA_03)

plot(res_CHFU_01) #KS 
plot(res_CHFU_02) #KS + heter.
plot(res_CHFU_both) #KS + High heter.
plot(res_CHFU_03)#KS + High heter.

fitness.data %>% group_by(G_F) %>% count()

# more specific plots
plotResiduals(res_LEMA_01, (fitness_LEMA$homo_motif))
plotResiduals(res_LEMA_01, (fitness_LEMA$hete_motif))
plotResiduals(res_LEMA_01, (fitness_LEMA$StrengthIn))
plotResiduals(res_LEMA_01, (fitness_LEMA$Ratio))
plotResiduals(res_LEMA_01, (fitness_LEMA$Delta))

# 
plotResiduals(res_PUPA_02, scale(fitness_PUPA$homo_motif))
# plotResiduals(res_PUPA_01, fitness_PUPA$hete_motif)
plotResiduals(res_PUPA_02, scale(fitness_PUPA$StrengthIn))
plotResiduals(res_PUPA_02, scale(fitness_PUPA$Ratio))
plotResiduals(res_PUPA_02, scale(fitness_PUPA$Delta))
#
plotResiduals(res_CHFU_02, (fitness_CHFU$homo_motif))
plotResiduals(res_CHFU_02, (fitness_CHFU$hete_motif))
plotResiduals(res_CHFU_02, (fitness_CHFU$StrengthIn))
plotResiduals(res_CHFU_02, (fitness_CHFU$Ratio))
plotResiduals(res_CHFU_02, (fitness_CHFU$Delta))



# check differences between models with and without ZI, when seed > 0 -----
AIC(GF_LEMA_01,GF_LEMA_01_ZI, GF_LEMA_02, GF_LEMA_03) 
AIC(GF_CHFU_01,GF_CHFU_01_ZI, GF_CHFU_02, GF_CHFU_03)
AIC(GF_PUPA_01,GF_PUPA_01_ZI, GF_PUPA_02, GF_PUPA_03)


library(performance)
r2(GF_LEMA_01)
r2(GF_PUPA_01)
r2(GF_CHFU_01)

summary(GF_LEMA_01)
summary(GF_PUPA_01_ZI)
summary(GF_CHFU_01_ZI)

ggplot(fitness_LEMA)+
  geom_point(aes(x=scale(Delta),y=Seeds_GF))

fitness_LEMA2 <- fitness_LEMA
fitness_LEMA$Ratio <- scale(fitness_LEMA$Ratio)
library(ggeffects)
ggpredict(GF_LEMA_01, c("") ) %>% plot(rawdata = TRUE, jitter = .01)


glmmTMB(Seeds_GF ~ scale(homo_motif) + 
          # scale(hete_motif) + 
          scale(StrengthIn) + scale(Delta) + 
          (1|Plot),
        #ziformula = ~1,
        family = nbinom1(),
        data = fitness_PUPA)


fitness_CHFU %>% group_by(G_F) %>% count()
fitness_LEMA %>% group_by(G_F) %>% count()
fitness_PUPA %>% group_by(G_F) %>% count()



































































##########################################
#GLMs
##########################################

fitness_LEMA <- fitness_LEMA %>% mutate(scl_homo=scale(homo_motif),
                                        scl_hete=scale(hete_motif),
                                        scl_strength=scale(StrengthIn), 
                                        scl_ratio=scale(Delta))

LEMA_glm1 <- glm.nb(Seeds_GF ~ scale(homo_motif) + scale(hete_motif) + scale(StrengthIn) +  
                         scale(Delta), data = fitness_LEMA)

LEMA_glm2 <- glmmTMB(Seeds_GF ~ scale(homo_motif) + scale(hete_motif)+
                    scale(StrengthIn) + scale(Delta),
                  #ziformula = ~1,
                  family = nbinom2(),
                  data = fitness_LEMA)

LEMA_glm3 <- glmmTMB(Seeds_GF ~ scl_homo +scl_hete+
                       scl_strength + scl_ratio,
                     #ziformula = ~1,
                     family = nbinom2(),
                     data = fitness_LEMA)


# Quasipoisson is not supported by DARHMA!!!!!! <- USAMOS NEG. BIN.
PUPA_glm1 <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                       scale(StrengthIn) + scale(Delta),
                     #ziformula = ~1,
                     #family = quasipoisson,
                     data = fitness_PUPA)

# Quasipoisson is not supported by DARHMA!!!!!! <- USAMOS NEG. BIN.
PUPA_glm11 <- glm(Seeds_GF ~ scale(homo_motif) +
                      scale(StrengthIn) + scale(Delta),
                    #ziformula = ~1,
                    family = quasipoisson,
                    data = fitness_PUPA)


PUPA_glm2 <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                      scale(StrengthIn) + scale(Delta),
                    #ziformula = ~1,
                    family = nbinom1(),
                    data = fitness_PUPA)

CHFU_glm1 <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                      scale(hete_motif) + 
                      scale(StrengthIn) + scale(Delta),
                    #ziformula = ~1,
                    #family = quasipoisson,
                    data = fitness_CHFU)

CHFU_glm11 <- glm(Seeds_GF ~ scale(homo_motif) +
                      scale(hete_motif) + 
                      scale(StrengthIn) + scale(Delta),
                    #ziformula = ~1,
                    family = quasipoisson,
                    data = fitness_CHFU)

CHFU_glm2 <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                      scale(hete_motif) + 
                      scale(StrengthIn) + scale(Delta),
                    #ziformula = ~1,
                    family = nbinom1(),
                    data = fitness_CHFU)


summary(LEMA_glm1)
summary(LEMA_glm2)

summary(PUPA_glm1)
summary(PUPA_glm2)


summary(CHFU_glm1)
summary(CHFU_glm2)

# glmmTMB residuals show a better behavior

res_LEMA_04 <- simulateResiduals(fittedModel = LEMA_glm1, n = 1500)
res_PUPA_04 <- simulateResiduals(fittedModel = PUPA_glm1, n = 1500)
res_CHFU_04 <- simulateResiduals(fittedModel = CHFU_glm1, n = 1500)

res_LEMA_05 <- simulateResiduals(fittedModel = LEMA_glm2, n = 1500)
res_PUPA_05 <- simulateResiduals(fittedModel = PUPA_glm2, n = 1500)
res_CHFU_05 <- simulateResiduals(fittedModel = CHFU_glm2, n = 1500)

res_LEMA_06 <- simulateResiduals(fittedModel = LEMA_glm3, n = 1500)

# Residuals 

testDispersion(res_LEMA_04)
testDispersion(res_LEMA_05)
testDispersion(res_LEMA_06)
testDispersion(res_PUPA_04)
testDispersion(res_PUPA_05)
testDispersion(res_CHFU_04)
testDispersion(res_CHFU_05)

testZeroInflation(res_LEMA_04)
testZeroInflation(res_LEMA_05)
testZeroInflation(res_LEMA_06)
testZeroInflation(res_PUPA_04)
testZeroInflation(res_PUPA_05)
testZeroInflation(res_CHFU_04)
testZeroInflation(res_CHFU_05)

plot(res_LEMA_04)
plot(res_LEMA_05)
plot(res_LEMA_06)
plot(res_PUPA_04) #Neg. bin induces hete.
plot(res_PUPA_05)
plot(res_CHFU_04)#Neg. bin induces hete.
plot(res_CHFU_05)

# more specific plots
plotResiduals(res_LEMA_04, scale(fitness_LEMA$homo_motif))
plotResiduals(res_LEMA_04, scale(fitness_LEMA$hete_motif))
plotResiduals(res_LEMA_04, scale(fitness_LEMA$StrengthIn))
plotResiduals(res_LEMA_04, scale(fitness_LEMA$Ratio))

res_LEMA_04_P <- recalculateResiduals(res_LEMA_04,group = fitness_LEMA$G_F)
plot(res_LEMA_04_P)

plotResiduals(res_LEMA_05, scale(fitness_LEMA$homo_motif))
plotResiduals(res_LEMA_05, scale(fitness_LEMA$hete_motif))
plotResiduals(res_LEMA_05, scale(fitness_LEMA$StrengthIn))
plotResiduals(res_LEMA_05, scale(fitness_LEMA$Ratio))

library(car)

Anova(LEMA_glm2)
Anova(PUPA_glm2)
Anova(CHFU_glm2)


library(performance)
performance::r2(LEMA_glm1) #0.1838923
performance::r2(LEMA_glm2) #0.065
performance::r2(LEMA_glm3) #0.065
performance::r2(PUPA_glm1) #0.2111155
performance::r2(PUPA_glm11) # 0.9999705
performance::r2(PUPA_glm2) #0.792
performance::r2(CHFU_glm1) #0.2130732
performance::r2(CHFU_glm11) #1
performance::r2(CHFU_glm2) #0.908

# Conditional plots

library(visreg)
visreg(LEMA_glm1)#,trans=exp)
visreg(LEMA_glm2)
visreg(LEMA_glm2)

visreg(LEMA_glm3,"scl_homo",xlab="Ho-triplet",ylab="Seeds per individual")
visreg(LEMA_glm3,"scl_ratio",xlab="scl_ratio",ylab="Seeds per individual")

visreg(PUPA_glm1)
visreg(PUPA_glm11)
visreg(PUPA_glm2)

visreg(CHFU_glm1)
visreg(CHFU_glm11)
visreg(CHFU_glm2)
visreg(CHFU_glm2)


library(lattice)
m <- matrix(c(0.1,0.2,0.3,0.4), 2, 2)
col.l <- colorRampPalette(c('blue', 'green'))(30)

plot.new()
par(mfrow=c(1,3))
visreg(LEMA_glm2,"homo_motif",xlab="Ho-triplet",ylab="Seeds per individual",trans=exp)
mtext(side=3, "LEMA", outer=F)
visreg(PUPA_glm2,"homo_motif",xlab="Ho-triplet",ylab="Seeds per individual",trans=exp)
mtext(side=3, "PUPA", outer=F)
visreg(CHFU_glm2,"homo_motif",xlab="Ho-triplet",ylab="Seeds per individual",trans=exp)
mtext(side=3, "CHFU", outer=F)
dev.off()

plot.new()
par(mfrow=c(1,2))
visreg(LEMA_glm2,"hete_motif",xlab="He-triplet",ylab="Seeds per individual",trans=exp)
mtext(side=3, "LEMA", outer=F)
visreg(CHFU_glm2,"hete_motif",xlab="He-triplet",ylab="Seeds per individual",trans=exp)
mtext(side=3, "CHFU", outer=F)
dev.off()

plot.new()
par(mfrow=c(1,3))
visreg(LEMA_glm2,"StrengthIn",xlab="In-Strength",ylab="Seeds per individual",trans=exp)
mtext(side=3, "LEMA", outer=F)
visreg(PUPA_glm2,"StrengthIn",xlab="In-Strength",ylab="Seeds per individual",trans=exp)
mtext(side=3, "PUPA", outer=F)
visreg(CHFU_glm2,"StrengthIn",xlab="In-Strength",ylab="Seeds per individual",trans=exp)
mtext(side=3, "CHFU", outer=F)
dev.off()

plot.new()
par(mfrow=c(1,3))
visreg(LEMA_glm2,"Ratio",xlab="Ratio",ylab="Seeds per individual",trans=exp)
mtext(side=3, "LEMA", outer=F)
visreg(PUPA_glm2,"Ratio",xlab="Ratio",ylab="Seeds per individual",trans=exp)
mtext(side=3, "PUPA", outer=F)
visreg(CHFU_glm2,"Ratio",xlab="Ratio",ylab="Seeds per individual",trans=exp)
mtext(side=3, "CHFU", outer=F)
dev.off()

a=0.035
b=0.05
c = 0.05
# make labels and margins smaller
plot.new()
#par(cex=0.7, mai=c(0.1,0.1,0.2,0.1))

# define area for the boxplot
par(fig=c(0+a,0.33-a,0.25+b,0.50-c), new=TRUE)
visreg(LEMA_glm2,"Ratio",xlab="Ratio",ylab="Seeds per individual",trans=exp)
mtext(side=3, "LEMA", outer=F)
# define area for the boxplot
par(fig=c(0.33+a,0.66-a,0.25+b,0.50-c), new=TRUE)
visreg(PUPA_glm2,"Ratio",xlab="Ratio",ylab="Seeds per individual",trans=exp)
mtext(side=3, "PUPA", outer=F)
# define area for the boxplot
par(fig=c(0.66+a,1-a,0.25+b,0.5-c), new=TRUE)
visreg(CHFU_glm2,"Ratio",xlab="Ratio",ylab="Seeds per individual",trans=exp)
mtext(side=3, "CHFU", outer=F)

# define area for the boxplot
par(fig=c(0.17+a,0.5-a,0.0+b,0.25-c), new=TRUE)
visreg(LEMA_glm2,"hete_motif",xlab="He-triplet",ylab="Seeds per individual",trans=exp)
mtext(side=3, "LEMA", outer=F)
# define area for the boxplot
par(fig=c(0.5+a,.83-a,0.0+b,0.25-c), new=TRUE)
visreg(CHFU_glm2,"hete_motif",xlab="He-triplet",ylab="Seeds per individual",trans=exp)
mtext(side=3, "CHFU", outer=F)




####################################################3
#####################################################

# En LEMA (que gasta una binomial negativa "de toda la vida") los resultados bayesianos 
# y glmmTMB son muy similares (en términos de coeficientes y parámetros), y sus R2 
# están dentro del mismo orden de magnitud (0.06, 0.11).
# Por otro lado, con CHFU, la cuasi-poisson de glmmTMB (una binomial negativa de 
# "parametrización lineal") y la binomial negativa bayesiana ajustan coeficientes y 
# parámetros similares. Sin embargo, R2 de la quasipoisson es 1, mientras que la bayesiana
# es 0.2. Si ajusto CHFU con glmmTMB y una binomial negativa obtengo R2=0.2. La gracia es 
# que Dharma muestra que aparecen problemas con los residuos (edited) 
# 
# NOTA: calculo R2 con el paquete "performance". Al usarlo con la quasi-poisson 
# no protesta ni indica que esté sucediendo algo chungo, y los resultados para la 
# cuasi-poisson de MASS y glmmTMB son parecidos: 0.99 y 0.9 (edited) 

library(brms)

LEMA_by <- brm(Seeds_GF ~ scale(homo_motif) +
                 scale(hete_motif) + 
                 scale(StrengthIn) + scale(Delta), data = fitness_LEMA,
               family = negbinomial(link = "log"))
summary(LEMA_by)
r2(LEMA_by)
r2(LEMA_glm2)
r2(LEMA_glm1)

CHFU_by <- brm(Seeds_GF ~ scale(homo_motif) +
                 scale(hete_motif) + 
                 scale(StrengthIn) + scale(Delta), data = fitness_CHFU,
               family = negbinomial(link = "log"))
summary(CHFU_by)
r2(CHFU_by)
r2(CHFU_glm1)
r2(CHFU_glm11)



res_CHFU_06 <- simulateResiduals(fittedModel = CHFU_by, n = 1500)
plot(res_CHFU_06)

posterior <- as.array(CHFU_by)
dim(posterior)
dimnames(posterior)



library("bayesplot")
library("ggplot2")
library("rstanarm")

color_scheme_set("blue")
mcmc_trace(posterior, pars = c("b_scalehomo_motif", "b_scalehete_motif","b_scaleStrengthIn","b_scaleRatio"))