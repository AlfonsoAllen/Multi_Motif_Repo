
# Fit GLMMs for multi-motif study in Caracoles (2020)
# Focal individuals with insect visits have their own value of seeds/fruit
library(tidyverse)
source("R_Scripts/functions.R")

# Load data for models
fitness.data.GF <- load_data_models_2020_without_agg()
fitness_orig <- load_data_models_2020() 
fitness.data <- subset(fitness_orig,Seeds_GF > 0)

# Differences in visits received and focal plants examined due to seed=0
visits_plant_orig <- fitness_orig %>% ungroup() %>% group_by(Plant_Simple) %>%
  count(wt=visits_GF)
visits_plant_data <- fitness.data %>% ungroup() %>% group_by(Plant_Simple) %>%
  count(wt=visits_GF) %>% rename(nn=n) %>% 
  left_join(visits_plant_orig, by ="Plant_Simple") %>% mutate(dif=n-nn)

fitness.data %>% ungroup() %>% filter(DegreeIn==0) %>% group_by(Plant_Simple) %>%
  count() %>% arrange(desc(n))
fitness.data %>% ungroup() %>% filter(DegreeIn!=0) %>% group_by(Plant_Simple) %>%
  count() %>% arrange(desc(n))


# We correct the values of PageRank obtained for multilayers without isolated nodes
corrections_PR <- corrections_pagerank()

for (Plot_i in 1:9){
  
  #Correct data with seeds that may be 0---
  
  fitness_orig$Real_PR_Multi[fitness_orig$Plot==as.character(Plot_i) &
                               fitness_orig$DegreeIn!=0 ] <-  pull(corrections_PR[Plot_i,3])*
    fitness_orig$Real_PR_Multi[fitness_orig$Plot==as.character(Plot_i) &
                                 fitness_orig$DegreeIn!=0 ]
  
  fitness_orig$Real_PR_Layer[fitness_orig$Plot==as.character(Plot_i) &
                               fitness_orig$DegreeIn!=0 ] <-  pull(corrections_PR[Plot_i,3])*
    fitness_orig$Real_PR_Layer[fitness_orig$Plot==as.character(Plot_i) &
                                 fitness_orig$DegreeIn!=0 ]
  
  fitness_orig$Real_PR_Multi[fitness_orig$Plot==as.character(Plot_i) &
                               fitness_orig$DegreeIn==0 ] <- pull(corrections_PR[Plot_i,2])
  
  #Correct data with seeds > 0---
  
  fitness.data$Real_PR_Multi[fitness.data$Plot==as.character(Plot_i) &
                               fitness.data$DegreeIn!=0 ] <-  pull(corrections_PR[Plot_i,3])*
    fitness.data$Real_PR_Multi[fitness.data$Plot==as.character(Plot_i) &
                                 fitness.data$DegreeIn!=0 ]
  
  fitness.data$Real_PR_Layer[fitness.data$Plot==as.character(Plot_i) &
                               fitness.data$DegreeIn!=0 ] <-  pull(corrections_PR[Plot_i,3])*
    fitness.data$Real_PR_Layer[fitness.data$Plot==as.character(Plot_i) &
                                 fitness.data$DegreeIn!=0 ]
  
  fitness.data$Real_PR_Multi[fitness.data$Plot==as.character(Plot_i) &
                               fitness.data$DegreeIn==0 ] <- pull(corrections_PR[Plot_i,2])
}



# Add community
modules_final = NULL

for (Plot_i in 1:9){
  
  module_i <- read_csv(paste0("Processed_data/Modularity_Pheno_Overlap/2020_NN_Modularity_Plot",Plot_i,".csv")  )
  
  module_i <- module_i %>% filter(type == "plant") %>%
    separate(species,c("Subplot","Plant_Simple")," ") %>%
    dplyr::select(Plot,Subplot,Plant_Simple,module)
  
  modules_final=bind_rows(modules_final,module_i)

  
}

modules_final$Plot <- as.factor(modules_final$Plot)


fitness_orig <- fitness_orig %>% left_join(modules_final,
                                           by=c("Plot","Subplot","Plant_Simple"))

fitness_orig$module[is.na(fitness_orig$module)] <- "isolated"
fitness_orig$module <- as.factor(fitness_orig$module)

fitness.data <- fitness.data %>% left_join(modules_final,
                                               by=c("Plot","Subplot","Plant_Simple"))

fitness.data$module[is.na(fitness.data$module)] <- "isolated"
fitness.data$module <- as.factor(fitness.data$module)

# Load libraries for analysis-----


library(MASS)
library(fitdistrplus)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(bbmle) ## for AICtab
library(DHARMa)
library(performance)
library(visreg)
library(RColorBrewer)
library(nlme)
library(usdm)

###############################
# MODEL ASSUMPTIONS: There is correlation between visits and seeds

##########################
# SEED PRODUCTION ~ VISITS
############################
library(ggpmisc)
my.formula <- y ~ x
ggplot(fitness.data, aes(x = visits_GF,y=log(Seeds_GF))) +
  geom_point(alpha=0.3)+geom_smooth(aes(color=Plant_Simple),method = "lm", se = F)+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, size = 3) +
  facet_wrap(~Plant_Simple)

ggplot(fitness.data, aes(x = visits_GF,y=log(Seeds_GF))) +
  geom_point(alpha=0.3)+geom_smooth(aes(color=Plant_Simple),method = "lm", se = F)+
  stat_fit_glance(method = 'lm', method.args = list(formula = my.formula),
                  aes(label = paste("P-value = ",
                                    signif(..p.value.., digits = 4), sep = "")),
                  parse = TRUE, size = 1.8,label.y.npc = 'bottom')+
  facet_wrap(~Plant_Simple)+
  theme(strip.text.y = element_text(size = 6, colour = "red", angle = 90))

# For LEMA, CHFU and PUPA correlations are positive and significant -> OK!

cor(log(fitness.data$visits_GF[fitness.data$Plant_Simple=="LEMA"]),
    fitness.data$homo_motif[fitness.data$Plant_Simple=="LEMA"],method = "spearman")
cor(log(fitness.data$visits_GF[fitness.data$Plant_Simple=="CHFU"]),
    fitness.data$hete_motif[fitness.data$Plant_Simple=="CHFU"],method = "spearman")
cor(log(fitness.data$visits_GF[fitness.data$Plant_Simple=="PUPA"]),
    fitness.data$StrengthIn[fitness.data$Plant_Simple=="PUPA"],method = "spearman")

GF_MIX_NB_Visits <- glmmTMB(Seeds_GF ~ scale(visits_GF) +
                              (scale(visits_GF)|Plot) +(scale(visits_GF)|Plant_Simple) ,
                            family = nbinom2(),
                            data = fitness.data)

GF_MIX_LIN_Visits <- lmer(log(Seeds_GF) ~ (visits_GF) +
                            (1|Plot) +(1|Plant_Simple) ,
                          #family = gaussian(),
                          data = fitness.data)

summary(GF_MIX_NB_Visits)
summary(GF_MIX_LIN_Visits)# Intercept model shows a significant relation between seeds and
# visits

###############################
# CORRELATIONS BETWEEN EACH CENTRALITY METRIC AND VISITS

cor(fitness.data$visits_GF,fitness.data$homo_motif,method = "spearman")
cor(fitness.data$visits_GF,fitness.data$hete_motif,method = "spearman")
cor(fitness.data$visits_GF,fitness.data$StrengthIn,method = "spearman")
cor(fitness.data$visits_GF,fitness.data$Ratio,method = "spearman") #This looks odd -> Check
#values without isolated plant nodes

# Correlations only for visited focals

cor(fitness.data$visits_GF[fitness.data$DegreeIn!=0],
    fitness.data$homo_motif[fitness.data$DegreeIn!=0],method = "spearman")
cor.test(fitness.data$visits_GF[fitness.data$DegreeIn!=0],
         fitness.data$homo_motif[fitness.data$DegreeIn!=0],method = "spearman")
cor(fitness.data$visits_GF[fitness.data$DegreeIn!=0],
    fitness.data$hete_motif[fitness.data$DegreeIn!=0],method = "spearman")
cor.test(fitness.data$visits_GF[fitness.data$DegreeIn!=0],
    fitness.data$hete_motif[fitness.data$DegreeIn!=0],method = "spearman")

cor(fitness.data$visits_GF[fitness.data$DegreeIn!=0],
    fitness.data$StrengthIn[fitness.data$DegreeIn!=0],method = "spearman")
cor.test(fitness.data$visits_GF[fitness.data$DegreeIn!=0],
    fitness.data$StrengthIn[fitness.data$DegreeIn!=0],method = "spearman")

cor(fitness.data$visits_GF[fitness.data$DegreeIn!=0],
    fitness.data$Real_PR_Multi[fitness.data$DegreeIn!=0],method = "spearman") #OK!
cor.test(fitness.data$visits_GF[fitness.data$DegreeIn!=0],
         fitness.data$Real_PR_Multi[fitness.data$DegreeIn!=0],method = "spearman")

cor(fitness.data$visits_GF[fitness.data$DegreeIn!=0],
    fitness.data$Ratio[fitness.data$DegreeIn!=0],method = "spearman") #OK!
cor.test(fitness.data$visits_GF[fitness.data$DegreeIn!=0],
    fitness.data$Ratio[fitness.data$DegreeIn!=0],method = "spearman")

plot(fitness.data$visits_GF[fitness.data$DegreeIn!=0],
     fitness.data$homo_motif[fitness.data$DegreeIn!=0])
plot(fitness.data$visits_GF[fitness.data$DegreeIn!=0],
     fitness.data$hete_motif[fitness.data$DegreeIn!=0])
plot(fitness.data$visits_GF[fitness.data$DegreeIn!=0],
     fitness.data$StrengthIn[fitness.data$DegreeIn!=0])
plot(fitness.data$visits_GF[fitness.data$DegreeIn!=0],
     fitness.data$Ratio[fitness.data$DegreeIn!=0])
plot(fitness.data$visits_GF,(fitness.data$Real_PR_Multi))


###############################
# ALL PLANT SPECIES MODEL: 2 Centrality Index + Homo_motif + Hetero_motifs
# (WITH and WITHOUT ZERO INFLATION FACTOR)
###############################

hist(log(fitness.data$Seeds_GF)) # Response variable looks aprox. normal.


fitness.data %>% group_by(Plant_Simple) %>% count()
fitness.data %>% group_by(Plant_Simple) %>% count()

#CHECK--------
vif(as.data.frame(dplyr::select(fitness.data %>% ungroup(),
                                homo_motif,hete_motif,StrengthIn, Ratio)))
cor(as.data.frame(dplyr::select(fitness.data %>% ungroup(),
                                homo_motif,hete_motif,StrengthIn, Ratio)))

GF_MIX_NB_intercept_Plot_Plant <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                       scale(hete_motif) + 
                                       scale(StrengthIn) + scale(Ratio) +
                                         Plot+(1|Plant_Simple),
                                     #ziformula = ~1,
                                     family = nbinom2(),
                                     data = fitness_orig)


GF_MIX_LIN_intercept_Plot_Plant <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                       scale(hete_motif) + 
                                       scale(StrengthIn) + scale(Ratio) +
                                         Plot+(1|Plant_Simple),
                                     #ziformula = ~1,
                                     family = gaussian(),
                                     data = fitness.data)

GF_MIX_NB_intercept_Plot_Plant_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                            scale(hete_motif) + 
                                            scale(StrengthIn) + scale(Ratio) +
                                              Plot+(1|Plant_Simple),
                                          ziformula = ~1,
                                          family = nbinom2(),
                                          data = fitness_orig)


# Summaries

summary(GF_MIX_NB_intercept_Plot_Plant)
summary(GF_MIX_NB_intercept_Plot_Plant_ZI)
summary(GF_MIX_LIN_intercept_Plot_Plant)

# Simulating residuals with Dahrma

res_GF_MIX_NB_intercept_Plot_Plant <- simulateResiduals(fittedModel = GF_MIX_NB_intercept_Plot_Plant, n = 500)
res_GF_MIX_LIN_intercept_Plot_Plant <- simulateResiduals(fittedModel = GF_MIX_LIN_intercept_Plot_Plant, n = 500)
res_GF_MIX_NB_intercept_Plot_Plant_ZI <- simulateResiduals(fittedModel = GF_MIX_NB_intercept_Plot_Plant_ZI, n = 500)

# Checking Residuals 
testZeroInflation(res_GF_MIX_NB_intercept_Plot_Plant)
testDispersion(res_GF_MIX_NB_intercept_Plot_Plant)

testZeroInflation(res_GF_MIX_NB_intercept_Plot_Plant_ZI)
testDispersion(res_GF_MIX_NB_intercept_Plot_Plant_ZI)

plot(res_GF_MIX_NB_intercept_Plot_Plant) #KS + 2 Quant desv
plot(res_GF_MIX_NB_intercept_Plot_Plant_ZI) #KS +2 Quant desv
plot(res_GF_MIX_LIN_intercept_Plot_Plant) #KS + 2 Quant desv


plotResiduals(res_GF_MIX_LIN_intercept_Plot_Plant, fitness.data$Plot)
plotResiduals(res_GF_MIX_LIN_intercept_Plot_Plant, fitness.data$Plant_Simple)
plotResiduals(res_GF_MIX_LIN_intercept_Plot_Plant, fitness.data$homo_motif) # It's not OK
plotResiduals(res_GF_MIX_LIN_intercept_Plot_Plant, fitness.data$hete_motif) #2D



################
# LEMA
################

fitness_orig_LEMA <- fitness_orig %>% filter(Plant_Simple=="LEMA")
fitness.data_LEMA <- fitness_orig %>% filter(Plant_Simple=="LEMA")

fitness.data_LEMA$Plot %>% unique()
fitness.data_LEMA$hete_motif %>% unique()

###############
# EXPLORING 

hist(fitness_orig_LEMA$Seeds_GF,50)
hist(log(fitness.data_LEMA$Seeds_GF),50)

fitness.data_LEMA$l_seeds <-  log(fitness.data_LEMA$Seeds_GF)

ggplot(fitness.data_LEMA)+geom_histogram(aes(Seeds_GF))+facet_wrap(~Plot)
ggplot(fitness.data_LEMA)+geom_histogram(aes(l_seeds))+facet_wrap(~Plot)

# Unusual frequent values in plots 4 and 5 are caused by inferred seeds (58)

ggplot(fitness.data_LEMA %>%
         filter(!(Seeds_GF==58 & Plot %in% c("4","5"))))+geom_histogram(aes(l_seeds))+facet_wrap(~Plot)

# We remove them
# 
# fitness_orig_LEMA <- fitness_orig_LEMA %>%
#   filter(!(Seeds_GF==58 & Plot %in% c("4","5")))
# fitness.data_LEMA <- fitness.data_LEMA %>%
#   filter(!(Seeds_GF==58 & Plot %in% c("4","5")))


library(lattice)
Z_LEMA <- cbind(log(fitness.data_LEMA$Seeds_GF),scale(fitness.data_LEMA$Fruit_GF),
           scale(fitness.data_LEMA$visits_GF),scale(fitness.data_LEMA$homo_motif),
           scale(fitness.data_LEMA$hete_motif),scale(fitness.data_LEMA$StrengthIn),
           scale(fitness.data_LEMA$DegreeIn),scale(fitness.data_LEMA$Ratio))

colnames(Z_LEMA) <- c("seeds", "fruits","visits", "homo",
                 "hete", "strengthIn", "DegreeIn", "Ratio")

dotplot(as.matrix(Z_LEMA), groups = FALSE,
        strip = strip.custom(bg = 'white',
                             par.strip.text = list(cex = 0.8)),
        scales = list(x = list(relation = "free"),
                      y = list(relation = "free"),
                      draw = T),
        col = 1, cex  = 0.5, pch = 16,
        xlab = "Value of the variable",
        ylab = "Order of the data from text file")

# There maybe outlayers according to the Cleveland plots

fitness.data_LEMA$s_hete <-  scale(fitness.data_LEMA$hete_motif)
fitness.data_LEMA$s_ratio <-  scale(fitness.data_LEMA$Ratio)
fitness.data_LEMA$s_fruits <-  scale(fitness.data_LEMA$Fruit_GF)
fitness.data_LEMA %>% filter(s_hete>5)
fitness.data.GF %>% filter(Plot=="3",
                         Subplot %in% c("B5","E2"),
                         Plant_Simple=="LEMA")

fitness.data_LEMA %>% filter(s_ratio>4)
fitness.data.GF %>% filter(Plot=="9",
                         Subplot %in% c("B5"),
                         Plant_Simple=="LEMA") #Seeds extrapolated

fitness.data_LEMA %>% filter(s_fruits>3)
fitness.data.GF %>% filter(Plot=="1",
                         Subplot %in% c("B2","B4","B6","C2"),
                         Plant_Simple=="LEMA")

# We calculate Cook distances to check it 

LM_LEMA_LIN_intercept_Plot_Plant <- lmer(log(Seeds_GF) ~ scale(homo_motif) +
                                              scale(hete_motif) + 
                                              scale(StrengthIn) + 
                                           scale(Ratio)+(1|Plot),
                                            data = fitness.data_LEMA)

library(influence.ME)
infl <- influence(LM_LEMA_LIN_intercept_Plot_Plant, obs = TRUE)
# Calculate Cook's distance:

cooks.distance(infl)

#Plot Cook's distance:
plot(infl, which = "cook")

#####################
#MODELS 
#No pooling----
GF_LEMA_NB_intercept_Plot_Plant <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                            scale(hete_motif) +
                                            scale(StrengthIn) + scale(Ratio) +
                                             Plot,
                                          #ziformula = ~1,
                                          family = nbinom2(),
                                          data = fitness_orig_LEMA)


GF_LEMA_LIN_intercept_Plot_Plant <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                             scale(hete_motif) +
                                             scale(StrengthIn) + scale(Ratio) +
                                              Plot,
                                           #ziformula = ~1,
                                           family = gaussian(),
                                           data = fitness.data_LEMA)

GF_LEMA_LIN_intercept_Plot_Plant_lm <- lm(log(Seeds_GF) ~ scale(homo_motif) +
                                              scale(hete_motif) +
                                              scale(StrengthIn) + scale(Ratio) +
                                              Plot,
                                            data = fitness.data_LEMA)

GF_LEMA_NB_intercept_Plot_Plant_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                               scale(hete_motif) + 
                                               scale(StrengthIn) + scale(Ratio) +
                                                Plot,
                                             ziformula = ~1,
                                             family = nbinom2(),
                                             data = fitness_orig_LEMA)

# Random factor Plot---

GF_LEMA_NB_intercept_Plot_Plant_RD <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                             scale(hete_motif) +
                                             scale(StrengthIn) + scale(Ratio) +
                                             (1|Plot),
                                           #ziformula = ~1,
                                           family = nbinom2(),
                                           data = fitness_orig_LEMA)


GF_LEMA_LIN_intercept_Plot_Plant_RD <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                              scale(hete_motif) +
                                              scale(StrengthIn) + scale(Ratio) +
                                                (1|Plot),
                                            #ziformula = ~1,
                                            family = gaussian(),
                                            data = fitness.data_LEMA)

GF_LEMA_LIN_intercept_Plot_Plant_lm_RD <- lmer(log(Seeds_GF) ~ scale(homo_motif) +
                                            scale(hete_motif) +
                                            scale(StrengthIn) + scale(Ratio) +
                                              (1|Plot),
                                          data = fitness.data_LEMA)

GF_LEMA_NB_intercept_Plot_Plant_ZI_RD <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                                scale(hete_motif) + 
                                                scale(StrengthIn) + scale(Ratio) +
                                                  (1|Plot),
                                              ziformula = ~1,
                                              family = nbinom2(),
                                              data = fitness_orig_LEMA)


# Random factor Plot/module---

GF_LEMA_NB_intercept_Plot_Plant_RD_mod <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                                scale(hete_motif) +
                                                scale(StrengthIn) + scale(Ratio) +
                                                (1|Plot/module),
                                              #ziformula = ~1,
                                              family = nbinom2(),
                                              data = fitness_orig_LEMA)


GF_LEMA_LIN_intercept_Plot_Plant_RD_mod <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                                 scale(hete_motif) +
                                                 scale(StrengthIn) + scale(Ratio) +
                                                 (1|Plot/module),
                                               #ziformula = ~1,
                                               family = gaussian(),
                                               data = fitness.data_LEMA)

GF_LEMA_LIN_intercept_Plot_Plant_lm_RD_mod <- lmer(log(Seeds_GF) ~ scale(homo_motif) +
                                                 scale(hete_motif) +
                                                 scale(StrengthIn) + scale(Ratio) +
                                                 (1|Plot/module),
                                               data = fitness.data_LEMA)

GF_LEMA_NB_intercept_Plot_Plant_ZI_RD_mod <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                                   scale(hete_motif) + 
                                                   scale(StrengthIn) + scale(Ratio) +
                                                   (1|Plot/module),
                                                 ziformula = ~1,
                                                 family = nbinom2(),
                                                 data = fitness_orig_LEMA)



summary(GF_LEMA_NB_intercept_Plot_Plant)
summary(GF_LEMA_NB_intercept_Plot_Plant_ZI)
summary(GF_LEMA_LIN_intercept_Plot_Plant)
summary(GF_LEMA_NB_intercept_Plot_Plant_RD)
summary(GF_LEMA_NB_intercept_Plot_Plant_ZI_RD)
summary(GF_LEMA_LIN_intercept_Plot_Plant_RD)
summary(GF_LEMA_NB_intercept_Plot_Plant_RD_mod)
summary(GF_LEMA_NB_intercept_Plot_Plant_ZI_RD_mod)
summary(GF_LEMA_LIN_intercept_Plot_Plant_RD_mod)

# Results seem similar

# Simulating residuals with Dahrma

res_GF_LEMA_NB_intercept_Plot_Plant <- simulateResiduals(fittedModel = GF_LEMA_NB_intercept_Plot_Plant, n = 500)
res_GF_LEMA_LIN_intercept_Plot_Plant <- simulateResiduals(fittedModel = GF_LEMA_LIN_intercept_Plot_Plant, n = 500)
res_GF_LEMA_NB_intercept_Plot_Plant_ZI<- simulateResiduals(fittedModel = GF_LEMA_NB_intercept_Plot_Plant_ZI, n = 500)
res_GF_LEMA_NB_intercept_Plot_Plant_RD <- simulateResiduals(fittedModel = GF_LEMA_NB_intercept_Plot_Plant_RD, n = 500)
res_GF_LEMA_LIN_intercept_Plot_Plant_RD <- simulateResiduals(fittedModel = GF_LEMA_LIN_intercept_Plot_Plant_RD, n = 500)
res_GF_LEMA_NB_intercept_Plot_Plant_ZI_RD <- simulateResiduals(fittedModel = GF_LEMA_NB_intercept_Plot_Plant_ZI_RD, n = 500)
res_GF_LEMA_NB_intercept_Plot_Plant_RD_mod <- simulateResiduals(fittedModel = GF_LEMA_NB_intercept_Plot_Plant_RD_mod, n = 500)
res_GF_LEMA_LIN_intercept_Plot_Plant_RD_mod <- simulateResiduals(fittedModel = GF_LEMA_LIN_intercept_Plot_Plant_RD_mod, n = 500)
res_GF_LEMA_NB_intercept_Plot_Plant_ZI_RD_mod <- simulateResiduals(fittedModel = GF_LEMA_NB_intercept_Plot_Plant_ZI_RD_mod, n = 500)

# Checking Residuals 
testZeroInflation(res_GF_LEMA_NB_intercept_Plot_Plant)
testDispersion(res_GF_LEMA_NB_intercept_Plot_Plant)

testZeroInflation(res_GF_LEMA_NB_intercept_Plot_Plant_ZI)
testDispersion(res_GF_LEMA_NB_intercept_Plot_Plant_ZI)

testZeroInflation(res_GF_LEMA_NB_intercept_Plot_Plant_RD)
testDispersion(res_GF_LEMA_NB_intercept_Plot_Plant_RD)

testZeroInflation(res_GF_LEMA_NB_intercept_Plot_Plant_ZI_RD)
testDispersion(res_GF_LEMA_NB_intercept_Plot_Plant_ZI_RD)

plot(res_GF_LEMA_NB_intercept_Plot_Plant) #2D
plot(res_GF_LEMA_NB_intercept_Plot_Plant_ZI) #2D
plot(res_GF_LEMA_LIN_intercept_Plot_Plant) #1D

plot(res_GF_LEMA_NB_intercept_Plot_Plant_RD) #3D + outer Newton did not converge fully
plot(res_GF_LEMA_NB_intercept_Plot_Plant_ZI_RD) #3D
plot(res_GF_LEMA_LIN_intercept_Plot_Plant_RD) #KS + 3D

plot(res_GF_LEMA_NB_intercept_Plot_Plant_RD_mod) #3D 
plot(res_GF_LEMA_NB_intercept_Plot_Plant_ZI_RD_mod) #3D
plot(res_GF_LEMA_LIN_intercept_Plot_Plant_RD_mod) #1D

# Residuals for the linear no pooling model look better than those for other models

plotResiduals(res_GF_LEMA_LIN_intercept_Plot_Plant, fitness.data_LEMA$homo_motif)
plotResiduals(res_GF_LEMA_LIN_intercept_Plot_Plant, fitness.data_LEMA$hete_motif) 
plotResiduals(res_GF_LEMA_LIN_intercept_Plot_Plant, fitness.data_LEMA$StrengthIn)
plotResiduals(res_GF_LEMA_LIN_intercept_Plot_Plant, fitness.data_LEMA$Ratio)
plotResiduals(res_GF_LEMA_LIN_intercept_Plot_Plant, fitness.data_LEMA$Plot)
plotResiduals(res_GF_LEMA_LIN_intercept_Plot_Plant, fitness.data_LEMA$Delta)

plotResiduals(res_GF_LEMA_NB_intercept_Plot_Plant, fitness.data_LEMA$homo_motif) #1D
plotResiduals(res_GF_LEMA_NB_intercept_Plot_Plant, fitness.data_LEMA$hete_motif) 
plotResiduals(res_GF_LEMA_NB_intercept_Plot_Plant, fitness.data_LEMA$StrengthIn)
plotResiduals(res_GF_LEMA_NB_intercept_Plot_Plant, fitness.data_LEMA$Ratio)
plotResiduals(res_GF_LEMA_NB_intercept_Plot_Plant, fitness.data_LEMA$Plot)
plotResiduals(res_GF_LEMA_NB_intercept_Plot_Plant, fitness.data_LEMA$Delta)

# Testing invariance of residuals by Plot
fitness.data_LEMA$resid_LIN <- resid(GF_LEMA_LIN_intercept_Plot_Plant_RD, type = "pearson")

hist(fitness.data_LEMA$resid_LIN)

ggResidpanel::resid_panel(GF_LEMA_LIN_intercept_Plot_Plant_lm_RD, plots = "all")
ggResidpanel::resid_panel(GF_LEMA_LIN_intercept_Plot_Plant_lm, plots = "all")

Test_LEMA <- lm(resid_LIN ~ Plot, data = fitness.data_LEMA)
drop1(Test_LEMA, test = "F") # Plot is not significant -> OK!



################
# CHFU

fitness_orig_CHFU <- fitness_orig %>% filter(Plant_Simple=="CHFU")
fitness.data_CHFU <- fitness.data %>% filter(Plant_Simple=="CHFU")
fitness.data_CHFU$Plot %>% unique()

# Readjust factor levels
fitness_orig_CHFU$Plot <- as.factor(as.numeric(fitness_orig_CHFU$Plot))
fitness.data_CHFU$Plot <- as.factor(as.numeric(fitness.data_CHFU$Plot))

###############
# EXPLORING 

hist(fitness_orig_CHFU$Seeds_GF,50)
hist(log(fitness.data_CHFU$Seeds_GF),50)

fitness.data_CHFU$l_seeds <-  log(fitness.data_CHFU$Seeds_GF)

ggplot(fitness.data_CHFU %>% filter(Seeds_GF<2000))+
  geom_histogram(aes(Seeds_GF))+facet_wrap(~Plot)
ggplot(fitness.data_CHFU%>% filter(Seeds_GF<2000))+
  geom_histogram(aes(l_seeds))+facet_wrap(~Plot)

# Unusual frequent values (32) are mainly produced by non-visited focals 

library(lattice)
Z_CHFU <- cbind(log(fitness.data_CHFU$Seeds_GF),scale(fitness.data_CHFU$Fruit_GF),
                scale(fitness.data_CHFU$visits_GF),scale(fitness.data_CHFU$homo_motif),
                scale(fitness.data_CHFU$hete_motif),scale(fitness.data_CHFU$StrengthIn),
                scale(fitness.data_CHFU$DegreeIn),scale(fitness.data_CHFU$Ratio))

colnames(Z_CHFU) <- c("seeds", "fruits","visits", "homo",
                      "hete", "strengthIn", "DegreeIn", "Ratio")

dotplot(as.matrix(Z_CHFU), groups = FALSE,
        strip = strip.custom(bg = 'white',
                             par.strip.text = list(cex = 0.8)),
        scales = list(x = list(relation = "free"),
                      y = list(relation = "free"),
                      draw = T),
        col = 1, cex  = 0.5, pch = 16,
        xlab = "Value of the variable",
        ylab = "Order of the data from text file")

# There maybe outlayers according to the Cleveland plots

fitness.data_CHFU$s_hete <-  scale(fitness.data_CHFU$hete_motif)
fitness.data_CHFU$s_ratio <-  scale(fitness.data_CHFU$Ratio)
fitness.data_CHFU$s_strenght <-  scale(fitness.data_CHFU$StrengthIn)
fitness.data_CHFU$s_fruits <-  scale(fitness.data_CHFU$Fruit_GF)
fitness.data_CHFU %>% filter(s_hete>5)
fitness.data.GF %>% filter(Plot=="9",
                         Subplot %in% c("B1"),
                         Plant_Simple=="CHFU")

fitness.data_CHFU %>% filter(s_strenght>5)
fitness.data.GF %>% filter(Plot %in% c("9"),
                         Subplot %in% c("F6"),
                         Plant_Simple=="CHFU")

fitness.data_CHFU %>% filter(s_ratio > 4)
fitness.data.GF %>% filter(Plot %in% c("3","9"),
                         Subplot %in% c("B6","C6","D2"),
                         Plant_Simple=="CHFU") #Seeds extrapolated

fitness.data_CHFU %>% filter(s_fruits>5)
fitness.data.GF %>% filter(Plot=="1",
                         Subplot %in% c("B1"),
                         Plant_Simple=="CHFU")

# We calculate Cook distances to check it 

LM_CHFU_LIN_intercept_Plot_Plant <- lmer(log(Seeds_GF) ~ scale(homo_motif) +
                                           scale(hete_motif) + 
                                           scale(StrengthIn) + 
                                           scale(Ratio)+(1|Plot),
                                         data = fitness.data_CHFU)

library(influence.ME)
infl <- influence(LM_CHFU_LIN_intercept_Plot_Plant, obs = TRUE)

# Calculate Cook's distance:
cooks.distance(infl)
cook_dist_CHFU <- cooks.distance(infl)
#Plot Cook's distance:
plot(infl, which = "cook")

cook_dist_CHFU[which(cook_dist_CHFU>0.15)]
fitness.data_CHFU[which(cook_dist_CHFU>0.15),]

#Observation 160 reaches 0.18: 9	F6	CHFU (This is the strength outlayer)

# That observation was not inferred
# We keep it inside
#####################
#MODELS 

# No pooling---
GF_CHFU_NB_intercept_Plot_Plant <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                             scale(hete_motif) + 
                                             scale(StrengthIn) + scale(Ratio)+
                                             Plot,
                                           #ziformula = ~1,
                                           family = nbinom2(),
                                           data = fitness_orig_CHFU)


GF_CHFU_LIN_intercept_Plot_Plant <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                              scale(hete_motif) + 
                                              scale(StrengthIn) + scale(Ratio)+
                                              Plot,
                                            #ziformula = ~1,
                                            family = gaussian(),
                                            data = fitness.data_CHFU)

GF_CHFU_NB_intercept_Plot_Plant_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                                scale(hete_motif) + 
                                                scale(StrengthIn) + scale(Ratio)+
                                                Plot,
                                              ziformula = ~1,
                                              family = nbinom2(),
                                              data = fitness_orig_CHFU)

# Random factor---
GF_CHFU_NB_intercept_Plot_Plant_RD <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                             scale(hete_motif) + 
                                             scale(StrengthIn) + scale(Ratio)+
                                             (1|Plot),
                                           #ziformula = ~1,
                                           family = nbinom2(),
                                           data = fitness_orig_CHFU)


GF_CHFU_LIN_intercept_Plot_Plant_RD <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                              scale(hete_motif) + 
                                              scale(StrengthIn) + scale(Ratio)+
                                                (1|Plot),
                                            #ziformula = ~1,
                                            family = gaussian(),
                                            data = fitness.data_CHFU)

GF_CHFU_NB_intercept_Plot_Plant_ZI_RD <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                                scale(hete_motif) + 
                                                scale(StrengthIn) + scale(Ratio)+
                                                  (1|Plot),
                                              ziformula = ~1,
                                              family = nbinom2(),
                                              data = fitness_orig_CHFU)

# Random factor module---
GF_CHFU_NB_intercept_Plot_Plant_RD_mod <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                                scale(hete_motif) + 
                                                scale(StrengthIn) + scale(Ratio)+
                                                (1|Plot/module),
                                              #ziformula = ~1,
                                              family = nbinom2(),
                                              data = fitness_orig_CHFU)


GF_CHFU_LIN_intercept_Plot_Plant_RD_mod <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                                 scale(hete_motif) + 
                                                 scale(StrengthIn) + scale(Ratio)+
                                                 (1|Plot/module),
                                               #ziformula = ~1,
                                               family = gaussian(),
                                               data = fitness.data_CHFU)

GF_CHFU_NB_intercept_Plot_Plant_ZI_RD_mod <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                                   scale(hete_motif) + 
                                                   scale(StrengthIn) + scale(Ratio)+
                                                   (1|Plot/module),
                                                 ziformula = ~1,
                                                 family = nbinom2(),
                                                 data = fitness_orig_CHFU)
summary(GF_CHFU_NB_intercept_Plot_Plant)
summary(GF_CHFU_NB_intercept_Plot_Plant_ZI)
summary(GF_CHFU_LIN_intercept_Plot_Plant)

summary(GF_CHFU_NB_intercept_Plot_Plant_RD)
summary(GF_CHFU_NB_intercept_Plot_Plant_ZI_RD)
summary(GF_CHFU_LIN_intercept_Plot_Plant_RD)

summary(GF_CHFU_NB_intercept_Plot_Plant_RD_mod)
summary(GF_CHFU_NB_intercept_Plot_Plant_ZI_RD_mod)
summary(GF_CHFU_LIN_intercept_Plot_Plant_RD_mod)


# Simulating residuals with Dahrma

res_GF_CHFU_NB_intercept_Plot_Plant <- simulateResiduals(fittedModel = GF_CHFU_NB_intercept_Plot_Plant, n = 500)
res_GF_CHFU_LIN_intercept_Plot_Plant <- simulateResiduals(fittedModel = GF_CHFU_LIN_intercept_Plot_Plant, n = 500)
res_GF_CHFU_NB_intercept_Plot_Plant_ZI <- simulateResiduals(fittedModel = GF_CHFU_NB_intercept_Plot_Plant_ZI, n = 500)

res_GF_CHFU_NB_intercept_Plot_Plant_RD <- simulateResiduals(fittedModel = GF_CHFU_NB_intercept_Plot_Plant_RD, n = 500)
res_GF_CHFU_LIN_intercept_Plot_Plant_RD <- simulateResiduals(fittedModel = GF_CHFU_LIN_intercept_Plot_Plant_RD, n = 500)
res_GF_CHFU_NB_intercept_Plot_Plant_ZI_RD <- simulateResiduals(fittedModel = GF_CHFU_NB_intercept_Plot_Plant_ZI_RD, n = 500)

res_GF_CHFU_NB_intercept_Plot_Plant_RD_mod <- simulateResiduals(fittedModel = GF_CHFU_NB_intercept_Plot_Plant_RD_mod, n = 500)
res_GF_CHFU_LIN_intercept_Plot_Plant_RD_mod <- simulateResiduals(fittedModel = GF_CHFU_LIN_intercept_Plot_Plant_RD_mod, n = 500)
res_GF_CHFU_NB_intercept_Plot_Plant_ZI_RD_mod <- simulateResiduals(fittedModel = GF_CHFU_NB_intercept_Plot_Plant_ZI_RD_mod, n = 500)


# Checking Residuals 
testZeroInflation(res_GF_CHFU_NB_intercept_Plot_Plant)
testDispersion(res_GF_CHFU_NB_intercept_Plot_Plant)

testZeroInflation(res_GF_CHFU_NB_intercept_Plot_Plant_ZI)
testDispersion(res_GF_CHFU_NB_intercept_Plot_Plant_ZI) #BAD

testZeroInflation(res_GF_CHFU_NB_intercept_Plot_Plant_RD)
testDispersion(res_GF_CHFU_NB_intercept_Plot_Plant_RD)

testZeroInflation(res_GF_CHFU_NB_intercept_Plot_Plant_ZI_RD)
testDispersion(res_GF_CHFU_NB_intercept_Plot_Plant_ZI_RD) #BAD

plot(res_GF_CHFU_NB_intercept_Plot_Plant) #Dev
plot(res_GF_CHFU_NB_intercept_Plot_Plant_ZI) #Dev+1D
plot(res_GF_CHFU_LIN_intercept_Plot_Plant) #OK

plot(res_GF_CHFU_NB_intercept_Plot_Plant_RD) #OK + pattern?
plot(res_GF_CHFU_NB_intercept_Plot_Plant_ZI_RD) #OK + pattern
plot(res_GF_CHFU_LIN_intercept_Plot_Plant_RD) #1D + pattern

plot(res_GF_CHFU_NB_intercept_Plot_Plant_RD_mod) #OK + pattern?
plot(res_GF_CHFU_NB_intercept_Plot_Plant_ZI_RD_mod) #OK + pattern
plot(res_GF_CHFU_LIN_intercept_Plot_Plant_RD_mod) #OK + pattern

# Again the residuals of the no pooling model look better

plotResiduals(res_GF_CHFU_LIN_intercept_Plot_Plant, fitness.data_CHFU$homo_motif)
plotResiduals(res_GF_CHFU_LIN_intercept_Plot_Plant, fitness.data_CHFU$hete_motif)
plotResiduals(res_GF_CHFU_LIN_intercept_Plot_Plant, fitness.data_CHFU$StrengthIn)
plotResiduals(res_GF_CHFU_LIN_intercept_Plot_Plant, fitness.data_CHFU$Ratio)
plotResiduals(res_GF_CHFU_LIN_intercept_Plot_Plant, fitness.data_CHFU$Plot)
plotResiduals(res_GF_CHFU_LIN_intercept_Plot_Plant, fitness.data_CHFU$Delta)


fitness.data_CHFU$resid_LIN <- resid(GF_CHFU_LIN_intercept_Plot_Plant_RD, type = "response")

hist(fitness.data_CHFU$resid_LIN)

Test_CHFU <- lm(resid_LIN ~ Plot, data = fitness.data_CHFU)
drop1(Test_CHFU, test = "F") # Plot is not significant -> OK!

# fitness_orig_CHFU$resid_NB <- resid(GF_CHFU_NB_intercept_Plot_Plant, type = "response")
# 
# hist(fitness_orig_CHFU$resid_NB)
# 
# Test_CHFU_NB <- lm(resid_NB ~ Plot, data = fitness_orig_CHFU)
# drop1(Test_CHFU_NB, test = "F") # Plot is not significant -> OK!

################
# PUPA

fitness_orig_PUPA <- fitness_orig %>% filter(Plant_Simple=="PUPA")
fitness.data_PUPA <- fitness.data %>% filter(Plant_Simple=="PUPA")

fitness.data_PUPA$Plot %>% unique() 

###############
# EXPLORING 

hist(fitness_orig_PUPA$Seeds_GF,50)
hist(log(fitness.data_PUPA$Seeds_GF,50))

fitness.data_PUPA$l_seeds <-  log(fitness.data_PUPA$Seeds_GF)

ggplot(fitness.data_PUPA %>% filter(Seeds_GF<1000) )+
  geom_histogram(aes(Seeds_GF))+facet_wrap(~Plot)
ggplot(fitness.data_PUPA%>% filter(Seeds_GF<1000))+
  geom_histogram(aes(l_seeds))+facet_wrap(~Plot)

fitness.data_PUPA %>% ungroup() %>% filter(Seeds_GF<1000) %>% dplyr::select(Seeds_GF,Plot) %>%
   group_by(Plot) %>% summarise_all(mean)
# The  mean values of the distribution are in the interval [100,400]
# Maybe we should not transform the data because the diff. among plots will become
# smaller
# Update: the residuals depend on plot when (1|Plot) is NOT included

library(lattice)
Z_PUPA <- cbind(log(fitness.data_PUPA$Seeds_GF),scale(fitness.data_PUPA$Fruit_GF),
                scale(fitness.data_PUPA$visits_GF),scale(fitness.data_PUPA$homo_motif),
                scale(fitness.data_PUPA$hete_motif),scale(fitness.data_PUPA$StrengthIn),
                scale(fitness.data_PUPA$DegreeIn),scale(fitness.data_PUPA$Ratio))

colnames(Z_PUPA) <- c("seeds", "fruits","visits", "homo",
                      "hete", "strengthIn", "DegreeIn", "Ratio")

dotplot(as.matrix(Z_PUPA), groups = FALSE,
        strip = strip.custom(bg = 'white',
                             par.strip.text = list(cex = 0.8)),
        scales = list(x = list(relation = "free"),
                      y = list(relation = "free"),
                      draw = T),
        col = 1, cex  = 0.5, pch = 16,
        xlab = "Value of the variable",
        ylab = "Order of the data from text file")

# There maybe outlayers according to the Cleveland plots

fitness.data_PUPA$s_hete <-  scale(fitness.data_PUPA$hete_motif)
fitness.data_PUPA$s_ratio <-  scale(fitness.data_PUPA$Ratio)
fitness.data_PUPA$s_strenght <-  scale(fitness.data_PUPA$StrengthIn)
fitness.data_PUPA$s_degree <-  scale(fitness.data_PUPA$DegreeIn)
fitness.data_PUPA$s_fruits <-  scale(fitness.data_PUPA$Fruit_GF)
fitness.data_PUPA$l_seeds <-  log(fitness.data_PUPA$Seeds_GF)
fitness.data_PUPA %>% filter(s_hete>5)
fitness.data.GF %>% filter(Plot=="2",
                         Subplot %in% c("C6"),
                         Plant_Simple=="PUPA")

fitness.data_PUPA %>% filter(s_strenght>5)
fitness.data.GF %>% filter(Plot %in% c("2"),
                         Subplot %in% c("C6"),
                         Plant_Simple=="PUPA")

fitness.data_PUPA %>% filter(s_degree>5)
fitness.data.GF %>% filter(Plot %in% c("2"),
                         Subplot %in% c("C6"),
                         Plant_Simple=="PUPA")

fitness.data_PUPA %>% filter(l_seeds < 2.2)
fitness.data.GF %>% filter(Plot %in% c("2","8"),
                         Subplot %in% c("C2","B4","F5"),
                         Plant_Simple=="PUPA")

# We calculate Cook distances to check it 

LM_PUPA_LIN_intercept_Plot_Plant <- lmer(log(Seeds_GF) ~ scale(homo_motif) +
                                           scale(hete_motif) + 
                                           scale(StrengthIn) + 
                                           scale(Ratio)+(1|Plot),
                                         data = fitness.data_PUPA)

library(influence.ME)
infl <- influence(LM_PUPA_LIN_intercept_Plot_Plant, obs = TRUE)

# Calculate Cook's distance:
cooks.distance(infl)
cook_results <- cooks.distance(infl)
#Plot Cook's distance:
plot(infl, which = "cook")
cook_results[which(cook_results>.15)]

PUPA_Cook <- fitness.data_PUPA[which(cook_results>.15),] #These entries come from real obs.


# index_outlayer_orig <- which((fitness_orig_PUPA$Plot=="8" &
#                                 fitness_orig_PUPA$Subplot=="F5"))
# index_outlayer_data <- which((fitness.data_PUPA$Plot=="8" &
#                                 fitness.data_PUPA$Subplot=="F5"))
# fitness_orig_PUPA <- fitness_orig_PUPA[-index_outlayer_orig,]
# fitness.data_PUPA <- fitness.data_PUPA[-index_outlayer_data,]


# MODELS
#No pooling------
GF_PUPA_NB_intercept_Plot_Plant <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                             scale(hete_motif) + 
                                             scale(StrengthIn) + scale(Ratio)+
                                             Plot,
                                           #ziformula = ~1,
                                           family = nbinom2(),
                                           data = fitness_orig_PUPA)


GF_PUPA_LIN_intercept_Plot_Plant <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                              scale(hete_motif) + 
                                              scale(StrengthIn) + scale(Ratio)+
                                              Plot,
                                            #ziformula = ~1,
                                            family = gaussian(),
                                            data = fitness.data_PUPA)

GF_PUPA_LIN_intercept_Plot_Plant_lmer <- lm(log(Seeds_GF) ~ scale(homo_motif) +
                                              scale(hete_motif) + 
                                              scale(StrengthIn) + scale(Ratio)+
                                                Plot,
                                            data = fitness.data_PUPA)


GF_PUPA_NB_intercept_Plot_Plant_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                                scale(hete_motif) + 
                                                scale(StrengthIn) + scale(Ratio)+
                                                Plot,
                                              ziformula = ~1,
                                              family = nbinom2(),
                                              data = fitness_orig_PUPA)

#Random factor Plot------
GF_PUPA_NB_intercept_Plot_Plant_RD <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                             scale(hete_motif) + 
                                             scale(StrengthIn) + scale(Ratio)+
                                             (1|Plot),
                                           #ziformula = ~1,
                                           family = nbinom2(),
                                           data = fitness_orig_PUPA)


GF_PUPA_LIN_intercept_Plot_Plant_RD <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                              scale(hete_motif) + 
                                              scale(StrengthIn) + scale(Ratio)+
                                                (1|Plot),
                                            #ziformula = ~1,
                                            family = gaussian(),
                                            data = fitness.data_PUPA)

GF_PUPA_LIN_intercept_Plot_Plant_lmer_RD <- lmer(log(Seeds_GF) ~ scale(homo_motif) +
                                              scale(hete_motif) + 
                                              scale(StrengthIn) + scale(Ratio)+
                                                (1|Plot),
                                            data = fitness.data_PUPA)


GF_PUPA_NB_intercept_Plot_Plant_ZI_RD <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                                scale(hete_motif) + 
                                                scale(StrengthIn) + scale(Ratio)+
                                                  (1|Plot),
                                              ziformula = ~1,
                                              family = nbinom2(),
                                              data = fitness_orig_PUPA)

#Random factor Plot/module------
GF_PUPA_NB_intercept_Plot_Plant_RD_mod <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                                scale(hete_motif) + 
                                                scale(StrengthIn) + scale(Ratio)+
                                                (1|Plot/module),
                                              #ziformula = ~1,
                                              family = nbinom2(),
                                              data = fitness_orig_PUPA)


GF_PUPA_LIN_intercept_Plot_Plant_RD_mod <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                                 scale(hete_motif) + 
                                                 scale(StrengthIn) + scale(Ratio)+
                                                 (1|Plot/module),
                                               #ziformula = ~1,
                                               family = gaussian(),
                                               data = fitness.data_PUPA)

GF_PUPA_LIN_intercept_Plot_Plant_lmer_RD_mod <- lmer(log(Seeds_GF) ~ scale(homo_motif) +
                                                   scale(hete_motif) + 
                                                   scale(StrengthIn) + scale(Ratio)+
                                                   (1|Plot/module),
                                                 data = fitness.data_PUPA)


GF_PUPA_NB_intercept_Plot_Plant_ZI_RD_mod <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                                   scale(hete_motif) + 
                                                   scale(StrengthIn) + scale(Ratio)+
                                                   (1|Plot/module),
                                                 ziformula = ~1,
                                                 family = nbinom2(),
                                                 data = fitness_orig_PUPA)

summary(GF_PUPA_NB_intercept_Plot_Plant)
summary(GF_PUPA_NB_intercept_Plot_Plant_ZI)
summary(GF_PUPA_LIN_intercept_Plot_Plant)
summary(GF_PUPA_LIN_intercept_Plot_Plant_lmer)

summary(GF_PUPA_NB_intercept_Plot_Plant_RD)
summary(GF_PUPA_NB_intercept_Plot_Plant_ZI_RD)
summary(GF_PUPA_LIN_intercept_Plot_Plant_RD)
summary(GF_PUPA_LIN_intercept_Plot_Plant_lmer_RD)

summary(GF_PUPA_NB_intercept_Plot_Plant_RD_mod)
summary(GF_PUPA_NB_intercept_Plot_Plant_ZI_RD_mod)
summary(GF_PUPA_LIN_intercept_Plot_Plant_RD_mod)
summary(GF_PUPA_LIN_intercept_Plot_Plant_lmer_RD_mod)

# Simulating residuals with Dahrma : BEST BEHAVIOR LN_ZI

res_GF_PUPA_NB_intercept_Plot_Plant <- simulateResiduals(fittedModel = GF_PUPA_NB_intercept_Plot_Plant,
                                                         n = 1000)
res_GF_PUPA_LIN_intercept_Plot_Plant <- simulateResiduals(fittedModel = GF_PUPA_LIN_intercept_Plot_Plant, 
                                                          n = 1000)

res_GF_PUPA_NB_intercept_Plot_Plant_ZI<- simulateResiduals(fittedModel = GF_PUPA_NB_intercept_Plot_Plant_ZI, 
                                                           n = 1000)

res_GF_PUPA_NB_intercept_Plot_Plant_RD <- simulateResiduals(fittedModel = GF_PUPA_NB_intercept_Plot_Plant_RD,
                                                         n = 1000)

res_GF_PUPA_LIN_intercept_Plot_Plant_RD <- simulateResiduals(fittedModel = GF_PUPA_LIN_intercept_Plot_Plant_RD, 
                                                          n = 1000)

res_GF_PUPA_NB_intercept_Plot_Plant_ZI_RD<- simulateResiduals(fittedModel = GF_PUPA_NB_intercept_Plot_Plant_ZI_RD, 
                                                           n = 1000)


res_GF_PUPA_NB_intercept_Plot_Plant_RD_mod <- simulateResiduals(fittedModel = GF_PUPA_NB_intercept_Plot_Plant_RD_mod,
                                                            n = 1000)

res_GF_PUPA_LIN_intercept_Plot_Plant_RD_mod <- simulateResiduals(fittedModel = GF_PUPA_LIN_intercept_Plot_Plant_RD_mod, 
                                                             n = 1000)

res_GF_PUPA_NB_intercept_Plot_Plant_ZI_RD_mod <- simulateResiduals(fittedModel = GF_PUPA_NB_intercept_Plot_Plant_ZI_RD_mod, 
                                                              n = 1000)

# Checking Residuals 
testZeroInflation(res_GF_PUPA_NB_intercept_Plot_Plant) #BAD
testDispersion(res_GF_PUPA_NB_intercept_Plot_Plant)

testZeroInflation(res_GF_PUPA_NB_intercept_Plot_Plant_ZI)
testDispersion(res_GF_PUPA_NB_intercept_Plot_Plant_ZI)

testZeroInflation(res_GF_PUPA_NB_intercept_Plot_Plant_RD) #BAD
testDispersion(res_GF_PUPA_NB_intercept_Plot_Plant_RD)

testZeroInflation(res_GF_PUPA_NB_intercept_Plot_Plant_ZI_RD)
testDispersion(res_GF_PUPA_NB_intercept_Plot_Plant_ZI_RD)

plot(res_GF_PUPA_NB_intercept_Plot_Plant) #KS +1 Quant desv
plot(res_GF_PUPA_NB_intercept_Plot_Plant_ZI) # OK
plot(res_GF_PUPA_LIN_intercept_Plot_Plant) #OK

plot(res_GF_PUPA_NB_intercept_Plot_Plant_RD) #KS +1 Quant desv + pattern
plot(res_GF_PUPA_NB_intercept_Plot_Plant_ZI_RD) # 1D + pattern
plot(res_GF_PUPA_LIN_intercept_Plot_Plant_RD) #1D + pattern

plot(res_GF_PUPA_NB_intercept_Plot_Plant_RD_mod) #KS +1 Quant desv + pattern
plot(res_GF_PUPA_NB_intercept_Plot_Plant_ZI_RD_mod) # 1D + pattern
plot(res_GF_PUPA_LIN_intercept_Plot_Plant_RD_mod) #KS + 1D + pattern

# Residuals look better in the no poling model

plotResiduals(res_GF_PUPA_LIN_intercept_Plot_Plant, fitness.data_PUPA$homo_motif)
plotResiduals(res_GF_PUPA_LIN_intercept_Plot_Plant, fitness.data_PUPA$hete_motif) 
plotResiduals(res_GF_PUPA_LIN_intercept_Plot_Plant, fitness.data_PUPA$StrengthIn) 
plotResiduals(res_GF_PUPA_LIN_intercept_Plot_Plant, fitness.data_PUPA$Ratio)
plotResiduals(res_GF_PUPA_LIN_intercept_Plot_Plant, fitness.data_PUPA$Plot)
plotResiduals(res_GF_PUPA_LIN_intercept_Plot_Plant, fitness.data_PUPA$Delta)#1D


fitness.data_PUPA$resid_LIN <- resid(GF_PUPA_LIN_intercept_Plot_Plant_RD, type = "response")

hist(fitness.data_PUPA$resid_LIN)

Test_PUPA <- lm(resid_LIN ~ Plot, data = fitness.data_PUPA)
drop1(Test_PUPA, test = "F") # fPlot is not significant -> OK!


##################################3
# FINAL MODELS

LEMA_LIN_intercept_Plot_Plant <- lm(log(Seeds_GF) ~ scale(homo_motif) +
                                      scale(hete_motif) + 
                                      scale(StrengthIn) + scale(Ratio)+
                                      Plot,
                                    data = fitness.data_LEMA)
CHFU_LIN_intercept_Plot_Plant <- lm(log(Seeds_GF) ~ scale(homo_motif) +
                                      scale(hete_motif) + 
                                      scale(StrengthIn) + scale(Ratio)+
                                      Plot,
                                    data = fitness.data_CHFU)
PUPA_LIN_intercept_Plot_Plant <- lm(log(Seeds_GF) ~ scale(homo_motif) + 
                                      scale(hete_motif) + 
                                      scale(StrengthIn) + scale(Ratio)+
                                      Plot,
                                    data = fitness.data_PUPA)

PUPA_LIN_intercept_Plot_Plant_int <- lm(log(Seeds_GF) ~ scale(homo_motif) *
                                      scale(hete_motif) + 
                                      scale(StrengthIn) + scale(Ratio)+
                                      Plot,
                                    data = fitness.data_PUPA)

summary(LEMA_LIN_intercept_Plot_Plant)
summary(CHFU_LIN_intercept_Plot_Plant)
summary(PUPA_LIN_intercept_Plot_Plant)
summary(PUPA_LIN_intercept_Plot_Plant_int)


plot(effects::allEffects(PUPA_LIN_intercept_Plot_Plant_int), multiline=TRUE, ci.style="bars")

drop1(LEMA_LIN_intercept_Plot_Plant, test = "F")
drop1(CHFU_LIN_intercept_Plot_Plant, test = "F")
drop1(PUPA_LIN_intercept_Plot_Plant, test = "F")
drop1(PUPA_LIN_intercept_Plot_Plant_int, test = "F")

car::vif(LEMA_LIN_intercept_Plot_Plant) # All the ( GVIF^(1/(2*Df)) )^2 < 5 
# (similar to an ordinary VIF smaller than 5 for one-coefficient variables)
# values for scale(homo_motif) are close to multicoll. threshold
car::vif(CHFU_LIN_intercept_Plot_Plant) # OK
car::vif(PUPA_LIN_intercept_Plot_Plant) # OK, although values for scale(homo_motif) & scale(hete_motif)
# are close to multicoll. threshold

############################
jtools::summ(LEMA_LIN_intercept_Plot_Plant,confint = TRUE,digits = 3)
jtools::summ(CHFU_LIN_intercept_Plot_Plant,confint = TRUE,digits = 3)
jtools::summ(PUPA_LIN_intercept_Plot_Plant,confint = TRUE,digits = 3)
jtools::summ(PUPA_LIN_intercept_Plot_Plant_int,confint = TRUE,digits = 3)

############################
#Visualization of slopes by using visreg

p1 <- visreg(LEMA_LIN_intercept_Plot_Plant,"homo_motif",xlab="Homo triplet",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="L. maroccanus")
p2 <- visreg(CHFU_LIN_intercept_Plot_Plant,"homo_motif",xlab="Homo triplet",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="C. fuscatum")
p3 <- visreg(PUPA_LIN_intercept_Plot_Plant,"homo_motif",xlab="Homo triplet",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="P. paludosa")


p5 <- visreg(LEMA_LIN_intercept_Plot_Plant,"hete_motif",xlab="Hetero triplet",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="L. maroccanus")
p6 <- visreg(CHFU_LIN_intercept_Plot_Plant,"hete_motif",xlab="Hetero triplet",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="C. fuscatum")
p7 <- visreg(PUPA_LIN_intercept_Plot_Plant,"hete_motif",xlab="Hetero triplet",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="P. paludosa")


p9 <- visreg(LEMA_LIN_intercept_Plot_Plant,"StrengthIn",xlab="Within layer\n centrality",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="L. maroccanus")
p10 <- visreg(CHFU_LIN_intercept_Plot_Plant,"StrengthIn",xlab="Within layer\n centrality",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="C. fuscatum")
p11 <- visreg(PUPA_LIN_intercept_Plot_Plant,"StrengthIn",xlab="Within layer\n centrality",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="P. paludosa")


p13 <- visreg(LEMA_LIN_intercept_Plot_Plant,"Ratio",xlab="Among layer\n centrality ratio",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="L. maroccanus")
p14 <- visreg(CHFU_LIN_intercept_Plot_Plant,"Ratio",xlab="Among layer\n centrality ratio",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="C. fuscatum")
p15 <- visreg(PUPA_LIN_intercept_Plot_Plant,"Ratio",xlab="Among layer\n centrality ratio",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="P. paludosa")


p17 <- visreg(LEMA_LIN_intercept_Plot_Plant,"Plot",xlab="Plot",ylab=NULL,gg = TRUE,  partial=T, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="L. maroccanus")

p18 <- visreg(CHFU_LIN_intercept_Plot_Plant,"Plot",xlab="Plot",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="C. fuscatum")
p19 <- visreg(PUPA_LIN_intercept_Plot_Plant,"Plot",xlab="Plot",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="P. paludosa")

library(gridExtra)
library(gtable)
library(grid)


grid.arrange(
  grobs = list(p2, p1, p3, p6,p5,p7,p10,p9,
               p11, p14, p13,p15),
  widths = c(1, 1, 1),
  layout_matrix = rbind(c(1, 2, 3),
                        c(4,5, 6),
                        c(7, 8, 9),
                        c(10, 11, 12)),
  left = textGrob("log(Seeds per individual)", rot = 90, vjust = 1)
)
#Save in 800x1000

sjPlot::plot_model(LEMA_LIN_intercept_Plot_Plant, type = "est")
sjPlot::plot_model(CHFU_LIN_intercept_Plot_Plant, type = "est")
sjPlot::plot_model(PUPA_LIN_intercept_Plot_Plant, type = "est")


#######
visreg2d(PUPA_LIN_intercept_Plot_Plant_int,"homo_motif","hete_motif",xlab="Homo triplet",ylab="Hetero triplet",zlab ="log(Seeds)")
