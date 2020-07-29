
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

########################
# SPATIAL LOCATION
##########################

set.seed(123)

x_r <- runif(nrow(fitness_orig), min=0, max=0.001)
y_r <- runif(nrow(fitness_orig), min=0, max=0.001)

coordinates <- read_csv2("Raw_Data/Caracoles_allplotsposition.csv") %>%
  rename(x = x_coor2, y = y_coor2) 


coordinates$Plot <- as.factor(coordinates$Plot)

fitness_orig <- fitness_orig %>% left_join(coordinates,by=c("Plot","Subplot")) %>% mutate(x=x+x_r,y=y+y_r)

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
#fitness.data <- subset(fitness.data,ID!="None")

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
summary(scale(fitness_LEMA$Ratio))
summary(scale(fitness_LEMA$pecentage_same_plants))

summary(scale(fitness_PUPA$homo_motif))
summary(scale(fitness_PUPA$hete_motif)) # fail
summary(scale(fitness_PUPA$DegreeIn))
summary(scale(fitness_PUPA$Real_PR_Layer))
summary(scale(fitness_PUPA$Ratio))
summary(scale(fitness_PUPA$pecentage_same_plants)) # fail

summary(scale(fitness_CHFU$homo_motif))
summary(scale(fitness_CHFU$hete_motif))
summary(scale(fitness_CHFU$DegreeIn))
summary(scale(fitness_CHFU$Real_PR_Layer))
summary(scale(fitness_CHFU$Ratio))
summary(scale(fitness_CHFU$pecentage_same_plants))

###################################################
# ALL SPECIES MODELS
###################################################

#################################################################
# NEGATIVE BINOMIAL--------RD SLOPE---------NO-ZI

GF_MIX_NB_slope_GF <- glmmTMB((Seeds_GF) ~ scale(homo_motif) +
                                scale(hete_motif) + 
                                scale(StrengthIn) + scale(Ratio) +
                                (0+scale(homo_motif)|G_F) + (1+scale(homo_motif)|Plant_Simple),
                              family = nbinom2(),
                              data = fitness.data)

GF_MIX_NB_slope_GF_Plot <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                     scale(hete_motif) + 
                                     scale(StrengthIn) + scale(Ratio) + 
                                     (0+scale(homo_motif)|G_F) + (0+scale(homo_motif)|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                   #ziformula = ~1,
                                   family = nbinom2(),
                                   data = fitness.data)

GF_MIX_NB_slope_Plot <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                  scale(hete_motif) + 
                                  scale(StrengthIn) + scale(Ratio) + 
                                  (0+scale(homo_motif)|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                #ziformula = ~1,
                                family = nbinom2(),
                                data = fitness.data)

################################################
# NEGATIVE BINOMIAL--------RD SLOPE---------ZI

GF_MIX_NB_slope_GF_ZI <- glmmTMB((Seeds_GF) ~ scale(homo_motif) +
                                   scale(hete_motif) + 
                                   scale(StrengthIn) + scale(Ratio) +
                                   (0+scale(homo_motif)|G_F)+ (1+scale(homo_motif)|Plant_Simple),
                                 ziformula = ~1,
                                 family = nbinom2(),
                                 data = fitness.data)

GF_MIX_NB_slope_GF_Plot_ZI <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                        scale(hete_motif) + 
                                        scale(StrengthIn) + scale(Ratio) + 
                                        (0+scale(homo_motif)|G_F) + (1+scale(homo_motif)|Plot)+ (0+scale(homo_motif)|Plant_Simple),
                                      ziformula = ~1,
                                      family = nbinom2(),
                                      data = fitness.data)

GF_MIX_NB_slope_Plot_ZI <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                     scale(hete_motif) + 
                                     scale(StrengthIn) + scale(Ratio) + 
                                     (0+scale(homo_motif)|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                   ziformula = ~1,
                                   family = nbinom2(),
                                   data = fitness.data)

#########################################3
# # NEGATIVE BINOMIAL--------RD INTERCEPT---------NO-ZI

GF_MIX_NB_intercept_GF <- glmmTMB((Seeds_GF) ~ scale(homo_motif) +
                                    scale(hete_motif) + 
                                    scale(StrengthIn) + scale(Ratio) +
                                    (1|G_F)+(1+scale(homo_motif)|Plant_Simple),
                                  family = nbinom2(),
                                  data = fitness.data)

GF_MIX_NB_intercept_GF_Plot <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                         scale(hete_motif) + 
                                         scale(StrengthIn) + scale(Ratio) + 
                                         (1|G_F) + (1|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                       #ziformula = ~1,
                                       family = nbinom2(),
                                       data = fitness.data)

GF_MIX_NB_intercept_Plot <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                      scale(hete_motif) + 
                                      scale(StrengthIn) + scale(Ratio) + 
                                      (1|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                    #ziformula = ~1,
                                    family = nbinom2(),
                                    data = fitness.data)

#########################################################
# NEGATIVE BINOMIAL--------RD INTERCEPT---------ZI

GF_MIX_NB_intercept_GF_ZI <- glmmTMB((Seeds_GF) ~ scale(homo_motif) +
                                       scale(hete_motif) + 
                                       scale(StrengthIn) + scale(Ratio) +
                                       (1|G_F)+ (1+scale(homo_motif)|Plant_Simple),
                                     ziformula = ~1,
                                     family = nbinom2(),
                                     data = fitness.data)

GF_MIX_NB_intercept_GF_Plot_ZI <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                            scale(hete_motif) + 
                                            scale(StrengthIn) + scale(Ratio) + 
                                            (1|G_F) + (1|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                          ziformula = ~1,
                                          family = nbinom2(),
                                          data = fitness.data)

GF_MIX_NB_intercept_Plot_ZI <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                         scale(hete_motif) + 
                                         scale(StrengthIn) + scale(Ratio) + 
                                         (1|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                       ziformula = ~1,
                                       family = nbinom2(),
                                       data = fitness.data)

#########################################3
# # NEGATIVE BINOMIAL--------RD inter_slop---------NO-ZI

GF_MIX_NB_inter_slop_GF <- glmmTMB((Seeds_GF) ~ scale(homo_motif) +
                                     scale(hete_motif) + 
                                     scale(StrengthIn) + scale(Ratio) +
                                     (1+scale(homo_motif)|G_F)+ (1+scale(homo_motif)|Plant_Simple),
                                   family = nbinom2(),
                                   data = fitness.data)

GF_MIX_NB_inter_slop_GF_Plot <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                          scale(hete_motif) + 
                                          scale(StrengthIn) + scale(Ratio) + 
                                          (1+scale(homo_motif)|G_F) + (1+scale(homo_motif)|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                        #ziformula = ~1,
                                        family = nbinom2(),
                                        data = fitness.data)

GF_MIX_NB_inter_slop_Plot <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                       scale(hete_motif) + 
                                       scale(StrengthIn) + scale(Ratio) + 
                                       (1+scale(homo_motif)|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                     #ziformula = ~1,
                                     family = nbinom2(),
                                     data = fitness.data)

#########################################################
# NEGATIVE BINOMIAL--------RD inter_slop---------ZI

GF_MIX_NB_inter_slop_GF_ZI <- glmmTMB((Seeds_GF) ~ scale(homo_motif) +
                                        scale(hete_motif) + 
                                        scale(StrengthIn) + scale(Ratio) +
                                        (1+scale(homo_motif)|G_F)+ (1+scale(homo_motif)|Plant_Simple),
                                      ziformula = ~1,
                                      family = nbinom2(),
                                      data = fitness.data)

GF_MIX_NB_inter_slop_GF_Plot_ZI <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                             scale(hete_motif) + 
                                             scale(StrengthIn) + scale(Ratio) + 
                                             (1+scale(homo_motif)|G_F) + (1+scale(homo_motif)|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                           ziformula = ~1,
                                           family = nbinom2(),
                                           data = fitness.data)

GF_MIX_NB_inter_slop_Plot_ZI <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                          scale(hete_motif) + 
                                          scale(StrengthIn) + scale(Ratio) + 
                                          (1+scale(homo_motif)|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                        ziformula = ~1,
                                        family = nbinom2(),
                                        data = fitness.data)


############################
# LINEAL MODELS LOG(SEEDS)--------RD SLOPE---------NO-ZI

GF_MIX_LIN_slope_GF <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                 scale(hete_motif) + 
                                 scale(StrengthIn) + scale(Ratio) +
                                 (0+scale(homo_motif)|G_F)+ (1+scale(homo_motif)|Plant_Simple),
                               family = gaussian(),
                               data = fitness.data)

GF_MIX_LIN_slope_GF_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                      scale(hete_motif) + 
                                      scale(StrengthIn) + scale(Ratio) + 
                                      (0+scale(homo_motif)|G_F) + (0+scale(homo_motif)|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                    #ziformula = ~1,
                                    family = gaussian(),
                                    data = fitness.data)

GF_MIX_LIN_slope_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                   scale(hete_motif) + 
                                   scale(StrengthIn) + scale(Ratio) + 
                                   (0+scale(homo_motif)|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                 #ziformula = ~1,
                                 family = gaussian(),
                                 data = fitness.data)

##################################
# LINEAL--------RD SLOPE---------ZI

GF_MIX_LIN_slope_GF_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                    scale(hete_motif) + 
                                    scale(StrengthIn) + scale(Ratio) +
                                    (0+scale(homo_motif)|G_F)+ (1+scale(homo_motif)|Plant_Simple),
                                  ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness.data)

GF_MIX_LIN_slope_GF_Plot_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                         scale(hete_motif) + 
                                         scale(StrengthIn) + scale(Ratio) + 
                                         (0+scale(homo_motif)|G_F) + (0+scale(homo_motif)|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                       ziformula = ~1,
                                       family = gaussian(),
                                       data = fitness.data)


GF_MIX_LIN_slope_Plot_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                      scale(hete_motif) + 
                                      scale(StrengthIn) + scale(Ratio) + 
                                      (0+scale(homo_motif)|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                    ziformula = ~1,
                                    family = gaussian(),
                                    data = fitness.data)


#########################################3
# LINEAL--------RD INTERCEPT---------NO-ZI

GF_MIX_LIN_intercept_GF <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                     scale(hete_motif) + 
                                     scale(StrengthIn) + scale(Ratio) +
                                     (1|G_F)+ (1+scale(homo_motif)|Plant_Simple),
                                   family = gaussian(),
                                   data = fitness.data)


GF_MIX_LIN_intercept_GF_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                          scale(hete_motif) + 
                                          scale(StrengthIn) + scale(Ratio) + 
                                          (1|G_F) + (1|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                        #ziformula = ~1,
                                        family = gaussian(),
                                        data = fitness.data)

GF_MIX_LIN_intercept_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                       scale(hete_motif) + 
                                       scale(StrengthIn) + scale(Ratio) + 
                                       (1|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                     #ziformula = ~1,
                                     family = gaussian(),
                                     data = fitness.data)


# LINEAL--------RD INTER---------ZI

GF_MIX_LIN_intercept_GF_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                        scale(hete_motif) + 
                                        scale(StrengthIn) + scale(Ratio) +
                                        (1|G_F)+ (1+scale(homo_motif)|Plant_Simple),
                                      ziformula = ~1,
                                      family = gaussian(),
                                      data = fitness.data)


GF_MIX_LIN_intercept_GF_Plot_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                             scale(hete_motif) + 
                                             scale(StrengthIn) + scale(Ratio) + 
                                             (1|G_F) + (1|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                           ziformula = ~1,
                                           family = gaussian(),
                                           data = fitness.data)

GF_MIX_LIN_intercept_Plot_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                          scale(hete_motif) + 
                                          scale(StrengthIn) + scale(Ratio) + 
                                          (1|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                        ziformula = ~1,
                                        family = gaussian(),
                                        data = fitness.data)

#########################################3
# LINEAL--------RD INTERCEPT+SLOPE---------NO-ZI

GF_MIX_LIN_inter_slop_GF <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                      scale(hete_motif) + 
                                      scale(StrengthIn) + scale(Ratio) +
                                      (1+scale(homo_motif)|G_F)+ (1+scale(homo_motif)|Plant_Simple),
                                    family = gaussian(),
                                    data = fitness.data)


GF_MIX_LIN_inter_slop_GF_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                           scale(hete_motif) + 
                                           scale(StrengthIn) + scale(Ratio) + 
                                           (1+scale(homo_motif)|G_F) + (1+scale(homo_motif)|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                         #ziformula = ~1,
                                         family = gaussian(),
                                         data = fitness.data)

GF_MIX_LIN_inter_slop_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                        scale(hete_motif) + 
                                        scale(StrengthIn) + scale(Ratio) + 
                                        (1+scale(homo_motif)|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                      #ziformula = ~1,
                                      family = gaussian(),
                                      data = fitness.data)


# LINEAL--------RD INTER+slope---------ZI

GF_MIX_LIN_inter_slop_GF_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                         scale(hete_motif) + 
                                         scale(StrengthIn) + scale(Ratio) +
                                         (1+scale(homo_motif)|G_F)+ (1+scale(homo_motif)|Plant_Simple),
                                       ziformula = ~1,
                                       family = gaussian(),
                                       data = fitness.data)


GF_MIX_LIN_inter_slop_GF_Plot_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                              scale(hete_motif) + 
                                              scale(StrengthIn) + scale(Ratio) + 
                                              (1+scale(homo_motif)|G_F) + (1+scale(homo_motif)|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                            ziformula = ~1,
                                            family = gaussian(),
                                            data = fitness.data)

GF_MIX_LIN_inter_slop_Plot_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                           scale(hete_motif) + 
                                           scale(StrengthIn) + scale(Ratio) + 
                                           (1+scale(homo_motif)|Plot)+ (1+scale(homo_motif)|Plant_Simple),
                                         ziformula = ~1,
                                         family = gaussian(),
                                         data = fitness.data)


# check residuals ---------------------------------------------------------

# get residuals


res_GF_MIX_NB_slope_GF <- simulateResiduals(fittedModel = GF_MIX_NB_slope_GF, n = 1500)
res_GF_MIX_NB_slope_GF_Plot<- simulateResiduals(fittedModel = GF_MIX_NB_slope_GF_Plot, n = 1500)
res_GF_MIX_NB_slope_Plot<- simulateResiduals(fittedModel = GF_MIX_NB_slope_Plot, n = 1500)
res_GF_MIX_NB_slope_GF_ZI<- simulateResiduals(fittedModel = GF_MIX_NB_slope_GF_ZI, n = 1500)
res_GF_MIX_NB_slope_GF_Plot_ZI<- simulateResiduals(fittedModel = GF_MIX_NB_slope_GF_Plot_ZI, n = 1500)
res_GF_MIX_NB_slope_Plot_ZI<- simulateResiduals(fittedModel = GF_MIX_NB_slope_Plot_ZI, n = 1500)
res_GF_MIX_NB_intercept_GF<- simulateResiduals(fittedModel = GF_MIX_NB_intercept_GF, n = 1500)
res_GF_MIX_NB_intercept_GF_Plot<- simulateResiduals(fittedModel = GF_MIX_NB_intercept_GF_Plot, n = 1500)
res_GF_MIX_NB_intercept_Plot<- simulateResiduals(fittedModel = GF_MIX_NB_intercept_Plot, n = 1500)
res_GF_MIX_NB_intercept_GF_ZI<- simulateResiduals(fittedModel = GF_MIX_NB_intercept_GF_ZI, n = 1500)
res_GF_MIX_NB_intercept_GF_Plot_ZI<- simulateResiduals(fittedModel = GF_MIX_NB_intercept_GF_Plot_ZI, n = 1500)
res_GF_MIX_NB_intercept_Plot_ZI<- simulateResiduals(fittedModel = GF_MIX_NB_intercept_Plot_ZI, n = 1500)
res_GF_MIX_NB_inter_slop_GF<- simulateResiduals(fittedModel = GF_MIX_NB_inter_slop_GF, n = 1500)
res_GF_MIX_NB_inter_slop_GF_Plot<- simulateResiduals(fittedModel = GF_MIX_NB_inter_slop_GF_Plot, n = 1500)
res_GF_MIX_NB_inter_slop_Plot<- simulateResiduals(fittedModel = GF_MIX_NB_inter_slop_Plot, n = 1500)
res_GF_MIX_NB_inter_slop_GF_ZI<- simulateResiduals(fittedModel = GF_MIX_NB_inter_slop_GF_ZI, n = 1500)
res_GF_MIX_NB_inter_slop_GF_Plot_ZI<- simulateResiduals(fittedModel = GF_MIX_NB_inter_slop_GF_Plot_ZI, n = 1500)
res_GF_MIX_NB_inter_slop_Plot_ZI<- simulateResiduals(fittedModel = GF_MIX_NB_inter_slop_Plot_ZI, n = 1500)



res_GF_MIX_LIN_slope_GF <- simulateResiduals(fittedModel = GF_MIX_LIN_slope_GF, n = 1500)
res_GF_MIX_LIN_slope_GF_Plot<- simulateResiduals(fittedModel = GF_MIX_LIN_slope_GF_Plot, n = 1500)
res_GF_MIX_LIN_slope_Plot<- simulateResiduals(fittedModel = GF_MIX_LIN_slope_Plot, n = 1500)
res_GF_MIX_LIN_slope_GF_ZI<- simulateResiduals(fittedModel = GF_MIX_LIN_slope_GF_ZI, n = 1500)
res_GF_MIX_LIN_slope_GF_Plot_ZI<- simulateResiduals(fittedModel = GF_MIX_LIN_slope_GF_Plot_ZI, n = 1500)
res_GF_MIX_LIN_slope_Plot_ZI<- simulateResiduals(fittedModel = GF_MIX_LIN_slope_Plot_ZI, n = 1500)
res_GF_MIX_LIN_intercept_GF<- simulateResiduals(fittedModel = GF_MIX_LIN_intercept_GF, n = 1500)
res_GF_MIX_LIN_intercept_GF_Plot<- simulateResiduals(fittedModel = GF_MIX_LIN_intercept_GF_Plot, n = 1500)
res_GF_MIX_LIN_intercept_Plot<- simulateResiduals(fittedModel = GF_MIX_LIN_intercept_Plot, n = 1500)
res_GF_MIX_LIN_intercept_GF_ZI<- simulateResiduals(fittedModel = GF_MIX_LIN_intercept_GF_ZI, n = 1500)
res_GF_MIX_LIN_intercept_GF_Plot_ZI<- simulateResiduals(fittedModel = GF_MIX_LIN_intercept_GF_Plot_ZI, n = 1500)
res_GF_MIX_LIN_intercept_Plot_ZI<- simulateResiduals(fittedModel = GF_MIX_LIN_intercept_Plot_ZI, n = 1500)
res_GF_MIX_LIN_inter_slop_GF<- simulateResiduals(fittedModel = GF_MIX_LIN_inter_slop_GF, n = 1500)
res_GF_MIX_LIN_inter_slop_GF_Plot<- simulateResiduals(fittedModel = GF_MIX_LIN_inter_slop_GF_Plot, n = 1500)
res_GF_MIX_LIN_inter_slop_Plot<- simulateResiduals(fittedModel = GF_MIX_LIN_inter_slop_Plot, n = 1500)
res_GF_MIX_LIN_inter_slop_GF_ZI<- simulateResiduals(fittedModel = GF_MIX_LIN_inter_slop_GF_ZI, n = 1500)
res_GF_MIX_LIN_inter_slop_GF_Plot_ZI<- simulateResiduals(fittedModel = GF_MIX_LIN_inter_slop_GF_Plot_ZI, n = 1500)
res_GF_MIX_LIN_inter_slop_Plot_ZI<- simulateResiduals(fittedModel = GF_MIX_LIN_inter_slop_Plot_ZI, n = 1500)


# Checking Residuals 

plot(res_GF_MIX_NB_slope_GF)#KS +3 Quant desv
plot(res_GF_MIX_NB_slope_GF_Plot)#KS +3 Quant desv
plot(res_GF_MIX_NB_slope_Plot)#1 Quant desv
plot(res_GF_MIX_NB_slope_GF_ZI)#KS +3 Quant desv
plot(res_GF_MIX_NB_slope_GF_Plot_ZI)#KS +3 Quant desv
plot(res_GF_MIX_NB_slope_Plot_ZI)#1 Quant desv
plot(res_GF_MIX_NB_intercept_GF)#3 Quant desv
plot(res_GF_MIX_NB_intercept_GF_Plot) #KS +2 Quant desv
plot(res_GF_MIX_NB_intercept_Plot)#KS +3 Quant desv
plot(res_GF_MIX_NB_intercept_GF_ZI)#3 Quant desv
plot(res_GF_MIX_NB_intercept_GF_Plot_ZI)#KS +3 Quant desv
plot(res_GF_MIX_NB_intercept_Plot_ZI)#KS +3 Quant desv

plot(res_GF_MIX_NB_inter_slop_GF)#3 Quant desv
plot(res_GF_MIX_NB_inter_slop_GF_Plot) #KS +2 Quant desv
plot(res_GF_MIX_NB_inter_slop_Plot)#KS +3 Quant desv
plot(res_GF_MIX_NB_inter_slop_GF_ZI)#3 Quant desv
plot(res_GF_MIX_NB_inter_slop_GF_Plot_ZI)#KS +2 Quant desv
plot(res_GF_MIX_NB_inter_slop_Plot_ZI)#KS +3 Quant desv

plot(res_GF_MIX_LIN_slope_GF)#KS +3 Quant desv
plot(res_GF_MIX_LIN_slope_GF_Plot)#KS +3 Quant desv
plot(res_GF_MIX_LIN_slope_Plot)
plot(res_GF_MIX_LIN_slope_GF_ZI)#KS +3 Quant desv
plot(res_GF_MIX_LIN_slope_GF_Plot_ZI)#KS +3 Quant desv
plot(res_GF_MIX_LIN_slope_Plot_ZI)
plot(res_GF_MIX_LIN_intercept_GF)#3 Quant desv
plot(res_GF_MIX_LIN_intercept_GF_Plot) #KS +2 Quant desv
plot(res_GF_MIX_LIN_intercept_Plot)#KS +3 Quant desv
plot(res_GF_MIX_LIN_intercept_GF_ZI)#3 Quant desv
plot(res_GF_MIX_LIN_intercept_GF_Plot_ZI)#KS +3 Quant desv
plot(res_GF_MIX_LIN_intercept_Plot_ZI)#KS +3 Quant desv

plot(res_GF_MIX_LIN_inter_slop_GF)#3 Quant desv
plot(res_GF_MIX_LIN_inter_slop_GF_Plot) #KS +2 Quant desv
plot(res_GF_MIX_LIN_inter_slop_Plot)#KS +3 Quant desv
plot(res_GF_MIX_LIN_inter_slop_GF_ZI)#3 Quant desv
plot(res_GF_MIX_LIN_inter_slop_GF_Plot_ZI)#KS +3 Quant desv
plot(res_GF_MIX_LIN_inter_slop_Plot_ZI)#KS +3 Quant desv
# check differences between models with and without ZI, when seed > 0 -----
AIC(
  GF_MIX_NB_slope_GF,
  GF_MIX_NB_slope_GF_Plot,
  GF_MIX_NB_slope_Plot,
  GF_MIX_NB_slope_GF_ZI,
  GF_MIX_NB_slope_GF_Plot_ZI,
  GF_MIX_NB_slope_Plot_ZI,
  GF_MIX_NB_intercept_GF,
  GF_MIX_NB_intercept_GF_Plot,
  GF_MIX_NB_intercept_Plot,
  GF_MIX_NB_intercept_GF_ZI,
  GF_MIX_NB_intercept_GF_Plot_ZI,
  GF_MIX_NB_intercept_Plot_ZI,
  GF_MIX_LIN_slope_GF,
  GF_MIX_LIN_slope_GF_Plot,
  GF_MIX_LIN_slope_Plot,
  GF_MIX_LIN_slope_GF_ZI,
  GF_MIX_LIN_slope_GF_Plot_ZI,
  GF_MIX_LIN_slope_Plot_ZI,
  GF_MIX_LIN_intercept_GF,
  GF_MIX_LIN_intercept_GF_Plot,
  GF_MIX_LIN_intercept_Plot,
  GF_MIX_LIN_intercept_GF_ZI,
  GF_MIX_LIN_intercept_GF_Plot_ZI,
  GF_MIX_LIN_intercept_Plot_ZI)


AIC(
  GF_MIX_NB_inter_slop_GF,
  GF_MIX_NB_inter_slop_GF_Plot,
  GF_MIX_NB_inter_slop_Plot,
  GF_MIX_NB_inter_slop_GF_ZI,
  GF_MIX_NB_inter_slop_GF_Plot_ZI,
  GF_MIX_NB_inter_slop_Plot_ZI)
AIC(
  GF_MIX_LIN_inter_slop_GF,
  GF_MIX_LIN_inter_slop_GF_Plot,
  GF_MIX_LIN_inter_slop_Plot,
  GF_MIX_LIN_inter_slop_GF_ZI,
  GF_MIX_LIN_inter_slop_GF_Plot_ZI,
  GF_MIX_LIN_inter_slop_Plot_ZI)

library(performance)
r2(GF_MIX_NB_slope_GF)
r2(GF_MIX_NB_slope_GF_Plot)
r2(GF_MIX_NB_slope_Plot)
r2(GF_MIX_NB_slope_GF_ZI)
r2(GF_MIX_NB_slope_GF_Plot_ZI)
r2(GF_MIX_NB_slope_Plot_ZI)
r2(GF_MIX_NB_intercept_GF)
r2(GF_MIX_NB_intercept_GF_Plot)
r2(GF_MIX_NB_intercept_Plot)
r2(GF_MIX_NB_intercept_GF_ZI)
r2(GF_MIX_NB_intercept_GF_Plot_ZI)
r2(GF_MIX_NB_intercept_Plot_ZI)

r2(GF_MIX_NB_inter_slop_GF)
r2(GF_MIX_NB_inter_slop_GF_Plot)
r2(GF_MIX_NB_inter_slop_Plot)
r2(GF_MIX_NB_inter_slop_GF_ZI)
r2(GF_MIX_NB_inter_slop_GF_Plot_ZI)
r2(GF_MIX_NB_inter_slop_Plot_ZI)

r2(GF_MIX_LIN_slope_GF)
r2(GF_MIX_LIN_slope_GF_Plot)
r2(GF_MIX_LIN_slope_Plot)
r2(GF_MIX_LIN_slope_GF_ZI)
r2(GF_MIX_LIN_slope_GF_Plot_ZI)
r2(GF_MIX_LIN_slope_Plot_ZI)
r2(GF_MIX_LIN_intercept_GF)
r2(GF_MIX_LIN_intercept_GF_Plot)
r2(GF_MIX_LIN_intercept_Plot)
r2(GF_MIX_LIN_intercept_GF_ZI)
r2(GF_MIX_LIN_intercept_GF_Plot_ZI)
r2(GF_MIX_LIN_intercept_Plot_ZI)

r2(GF_MIX_LIN_inter_slop_GF)
r2(GF_MIX_LIN_inter_slop_GF_Plot)
r2(GF_MIX_LIN_inter_slop_Plot)
r2(GF_MIX_LIN_inter_slop_GF_ZI)
r2(GF_MIX_LIN_inter_slop_GF_Plot_ZI)
r2(GF_MIX_LIN_inter_slop_Plot_ZI)

