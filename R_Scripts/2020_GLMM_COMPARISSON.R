
library(tidyverse)


# read data ---------------------------------------------------------------

fitness_final_aux0 <- read.csv(file = "2020_data_models_phenol_overlap.csv",
                         header = TRUE,
                         stringsAsFactors = FALSE) %>%
  rename(Plant_Simple=Plant) %>% mutate(Seeds_GF = round(Seeds_GF))
head(fitness_final_aux0, 20)

PageRank_results <- read.csv(file = "2020_PageRank_results.csv",
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

G_F_list <- read_csv2("Raw_Data/raw_Pollinators_2020_1.csv") %>%
  dplyr::select(G_F,ID) %>% unique() #rename(ID=ID_Simple)

G_F_list <- bind_rows(G_F_list,tibble(G_F="None",ID="None"))

# Fix "Elateridae","Musca_sp."

G_F_list$G_F[G_F_list$ID=="Elateridae"] <- "Flower_beetles"
G_F_list$G_F[G_F_list$ID=="Musca_sp."] <- "Hoverflies"

G_F_list <- unique(G_F_list)

# Sanity check
G_F_list %>% group_by(ID) %>% count() %>% filter(n>1)

fitness_orig <- fitness_final_aux %>% dplyr::left_join(G_F_list,by = "ID")


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

fitness_LEMA$G_F <- as.factor(fitness_LEMA$G_F)
fitness_PUPA$G_F  <- as.factor(fitness_PUPA$G_F)
fitness_CHFU$G_F  <- as.factor(fitness_CHFU$G_F)


# Data of LEMA plants with R > 1

# x <- fitness_LEMA %>% dplyr::select(Plot,Subplot,G_F,Plant_Simple,Ratio) %>% filter(Ratio>1)
# write_csv(x,"data_LEMA_PLANTS_R.csv")

# scale sometimes gives problems
summary(scale(fitness_LEMA$homo_motif))
summary(scale(fitness_LEMA$hete_motif))
summary(scale(fitness_LEMA$DegreeIn))
summary(scale(fitness_LEMA$Real_PR_Layer))
summary(scale(fitness_LEMA$Ratio))


summary(scale(fitness_PUPA$homo_motif))
summary(scale(fitness_PUPA$hete_motif))
summary(scale(fitness_PUPA$DegreeIn))
summary(scale(fitness_PUPA$Real_PR_Layer))
summary(scale(fitness_PUPA$Ratio))


summary(scale(fitness_CHFU$homo_motif))
summary(scale(fitness_CHFU$hete_motif))
summary(scale(fitness_CHFU$DegreeIn))
summary(scale(fitness_CHFU$Real_PR_Layer))
summary(scale(fitness_CHFU$Ratio))

###################################################
# ALL SPECIES MODELS
###################################################

head(fitness.data)
GF_MIX_LIN_intercept_Plot_Plant <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                             scale(hete_motif) + 
                                             scale(StrengthIn) + scale(Ratio) +
                                             (1|Plot) +(1|Plant_Simple) ,
                                           #ziformula = ~1,
                                           family = gaussian(),
                                           data = fitness.data)

summary(GF_MIX_LIN_intercept_Plot_Plant)


#########################################
# INDIVIDUAL MODEL FOR EACH PLANT SPECIES
#########################################

library(nlme)
library(usdm)

vif(as.data.frame(dplyr::select(fitness_PUPA,homo_motif,hete_motif)))
cor.test(fitness_PUPA$homo_motif, fitness_PUPA$hete_motif, method=c("pearson", "kendall", "spearman"))

vif(as.data.frame(dplyr::select(fitness_LEMA,DegreeIn,homo_motif,hete_motif,Ratio)))
vif(as.data.frame(dplyr::select(fitness_PUPA,DegreeIn,homo_motif,Ratio)))
vif(as.data.frame(dplyr::select(fitness_CHFU,DegreeIn,homo_motif,hete_motif,Ratio)))

vif(as.data.frame(dplyr::select(fitness_LEMA,DegreeIn,Real_PR_Layer,DegreeIn,Ratio)))
vif(as.data.frame(dplyr::select(fitness_LEMA,DegreeIn,homo_motif,hete_motif)))
vif(as.data.frame(dplyr::select(fitness_LEMA,PageRank,homo_motif,hete_motif)))


############################

# NB-------RD INTERCEPT---------NO-ZI


GF_LEMA_NB_intercept_Plot <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                                        scale(hete_motif) + 
                                        scale(StrengthIn) + scale(Ratio) +
                                        (1|Plot),
                                      #ziformula = ~1,
                                      family = nbinom1(),
                                      data = fitness_LEMA)

GF_PUPA_NB_intercept_Plot <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                                        scale(hete_motif) + 
                                        scale(StrengthIn) + scale(Ratio) +
                                        (1|Plot),
                                      #ziformula = ~1,
                                      family = nbinom1(),
                                      data = fitness_PUPA)

GF_CHFU_NB_intercept_Plot <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                        scale(hete_motif) + 
                                        scale(StrengthIn) + scale(Ratio) +
                                        (1|Plot),
                                      #ziformula = ~1,
                                      family = nbinom1(),
                                      data = fitness_CHFU)



#########################################3
# LINEAL--------RD INTERCEPT---------NO-ZI


GF_LEMA_LIN_intercept_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                       scale(hete_motif) + 
                                       scale(StrengthIn) + scale(Ratio) +
                                       (1|Plot),
                                     #ziformula = ~1,
                                     family = gaussian(),
                                     data = fitness_LEMA)

GF_PUPA_LIN_intercept_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                       scale(hete_motif) + 
                                       scale(StrengthIn) + scale(Ratio) +
                                       (1|Plot),
                                     #ziformula = ~1,
                                     family = gaussian(),
                                     data = fitness_PUPA)

GF_CHFU_LIN_intercept_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                       scale(hete_motif) + 
                                       scale(StrengthIn) + scale(Ratio) +
                                       (1|Plot),
                                     #ziformula = ~1,
                                     family = gaussian(),
                                     data = fitness_CHFU)



# NB--------RD SLOPE---------ZI


GF_LEMA_NB_intercept_Plot_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                                           scale(hete_motif) + 
                                           scale(StrengthIn) + scale(Ratio) +
                                           (1|Plot),
                                         ziformula = ~1,
                                         family = nbinom1(),
                                         data = fitness_LEMA)

GF_PUPA_NB_intercept_Plot_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                                           scale(hete_motif) + 
                                           scale(StrengthIn) + scale(Ratio) +
                                           (1|Plot),
                                         ziformula = ~1,
                                         family = nbinom1(),
                                         data = fitness_PUPA)

GF_CHFU_NB_intercept_Plot_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                           scale(hete_motif) + 
                                           scale(StrengthIn) + scale(Ratio) +
                                           (1|Plot),
                                         ziformula = ~1,
                                         family = nbinom1(),
                                         data = fitness_CHFU)


# LINEAL--------RD SLOPE---------ZI


GF_LEMA_LIN_intercept_Plot_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                          scale(hete_motif) + 
                                          scale(StrengthIn) + scale(Ratio) +
                                          (1|Plot),
                                        ziformula = ~1,
                                        family = gaussian(),
                                        data = fitness_LEMA)

GF_PUPA_LIN_intercept_Plot_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                          scale(hete_motif) + 
                                          scale(StrengthIn) + scale(Ratio) +
                                          (1|Plot),
                                        ziformula = ~1,
                                        family = gaussian(),
                                        data = fitness_PUPA)

GF_CHFU_LIN_intercept_Plot_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                          scale(hete_motif) + 
                                          scale(StrengthIn) + scale(Ratio) +
                                          (1|Plot),
                                        ziformula = ~1,
                                        family = gaussian(),
                                        data = fitness_CHFU)

# check residuals ---------------------------------------------------------

# get residuals


res_GF_LEMA_NB_intercept_Plot<- simulateResiduals(fittedModel = GF_LEMA_NB_intercept_Plot, n = 500)
res_GF_PUPA_NB_intercept_Plot<- simulateResiduals(fittedModel = GF_PUPA_NB_intercept_Plot, n = 500)
res_GF_CHFU_NB_intercept_Plot<- simulateResiduals(fittedModel = GF_CHFU_NB_intercept_Plot, n = 500)

res_GF_LEMA_LIN_intercept_Plot<- simulateResiduals(fittedModel = GF_LEMA_LIN_intercept_Plot, n = 500)
res_GF_PUPA_LIN_intercept_Plot<- simulateResiduals(fittedModel = GF_PUPA_LIN_intercept_Plot, n = 500)
res_GF_CHFU_LIN_intercept_Plot<- simulateResiduals(fittedModel = GF_CHFU_LIN_intercept_Plot, n = 500)

res_GF_LEMA_NB_intercept_Plot_ZI<- simulateResiduals(fittedModel = GF_LEMA_NB_intercept_Plot_ZI, n = 500)
res_GF_PUPA_NB_intercept_Plot_ZI<- simulateResiduals(fittedModel = GF_PUPA_NB_intercept_Plot_ZI, n = 500)
res_GF_CHFU_NB_intercept_Plot_ZI<- simulateResiduals(fittedModel = GF_CHFU_NB_intercept_Plot_ZI, n = 500)


res_GF_LEMA_LIN_intercept_Plot_ZI<- simulateResiduals(fittedModel = GF_LEMA_LIN_intercept_Plot_ZI, n = 500)
res_GF_PUPA_LIN_intercept_Plot_ZI<- simulateResiduals(fittedModel = GF_PUPA_LIN_intercept_Plot_ZI, n = 500)
res_GF_CHFU_LIN_intercept_Plot_ZI<- simulateResiduals(fittedModel = GF_CHFU_LIN_intercept_Plot_ZI, n = 500)


# Checking Residuals 
plot(res_GF_LEMA_NB_intercept_Plot)
plot(res_GF_LEMA_LIN_intercept_Plot)#KS +2 Quant desv
plot(res_GF_LEMA_NB_intercept_Plot_ZI)
plot(res_GF_LEMA_LIN_intercept_Plot_ZI)#KS +2 Quant desv


plot(res_GF_PUPA_NB_intercept_Plot)
plot(res_GF_PUPA_LIN_intercept_Plot)#2 Quant desv
plot(res_GF_PUPA_NB_intercept_Plot_ZI)
plot(res_GF_PUPA_LIN_intercept_Plot_ZI)#2 Quant desv

plot(res_GF_CHFU_NB_intercept_Plot)
plot(res_GF_CHFU_LIN_intercept_Plot)#KS + 3 Quant desv
plot(res_GF_CHFU_NB_intercept_Plot_ZI)
plot(res_GF_CHFU_LIN_intercept_Plot_ZI)#KS + 3 Quant desv

summary(GF_LEMA_LIN_intercept_Plot)
summary(GF_PUPA_LIN_intercept_Plot)
summary(GF_CHFU_LIN_intercept_Plot)

# check differences between models with and without ZI, when seed > 0 -----
AIC(
  GF_LEMA_NB_intercept_Plot,
  GF_LEMA_NB_intercept_Plot_ZI,
  GF_LEMA_LIN_intercept_Plot,
  GF_LEMA_LIN_intercept_Plot_ZI)

library(performance)

r2(GF_LEMA_NB_intercept_Plot)
r2(GF_LEMA_NB_intercept_Plot_ZI)
r2(GF_LEMA_LIN_intercept_Plot)
r2(GF_LEMA_LIN_intercept_Plot_ZI)


AIC(
  GF_PUPA_NB_intercept_Plot,
  GF_PUPA_NB_intercept_Plot_ZI,
  GF_PUPA_LIN_intercept_Plot,
  GF_PUPA_LIN_intercept_Plot_ZI)


r2(GF_PUPA_NB_intercept_Plot)
r2(GF_PUPA_NB_intercept_Plot_ZI)
r2(GF_PUPA_LIN_intercept_Plot)
r2(GF_PUPA_LIN_intercept_Plot_ZI)


AIC(
  GF_CHFU_NB_intercept_Plot,
  GF_CHFU_NB_intercept_Plot_ZI,
  GF_CHFU_LIN_intercept_Plot,
  GF_CHFU_LIN_intercept_Plot_ZI)


r2(GF_CHFU_NB_intercept_Plot)
r2(GF_CHFU_NB_intercept_Plot_ZI)
r2(GF_CHFU_LIN_intercept_Plot)
r2(GF_CHFU_LIN_intercept_Plot_ZI)



##########################################
# EFFECTS IN HOMO HETERO
##########################################
#########################################3
# LINEAL--------RD INTERCEPT---------NO-ZI


GF_LEMA_LIN_intercept_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif)*G_F + 
                                        scale(hete_motif)*G_F +
                                        (1|Plot),
                                      #ziformula = ~1,
                                      family = gaussian(),
                                      data = fitness_LEMA)

GF_PUPA_LIN_intercept_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif)*G_F + 
                                        scale(hete_motif)*G_F +
                                        (1|Plot),
                                      #ziformula = ~1,
                                      family = gaussian(),
                                      data = fitness_PUPA)

GF_CHFU_LIN_intercept_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif)*G_F + 
                                        scale(hete_motif)*G_F +
                                        (1|Plot),
                                      #ziformula = ~1,
                                      family = gaussian(),
                                      data = fitness_CHFU)

fitness_CHFU %>% dplyr::select(G_F) %>% group_by(G_F) %>%count()
fitness_LEMA %>% dplyr::select(G_F) %>% group_by(G_F) %>%count()
fitness_PUPA %>% dplyr::select(G_F) %>% group_by(G_F) %>%count()

fitness_LEMA %>% dplyr::select(G_F, visits_GF) %>% group_by(G_F) %>%count(wt=visits_GF)
