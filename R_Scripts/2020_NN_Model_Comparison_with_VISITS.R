
# Compare LM for Caracoles (2020) with alternative models based solely in visits or
# adding visits as a covariable

# How does R^2 vary??


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

################
# LEMA
################

fitness_orig_LEMA <- fitness_orig %>% filter(Plant_Simple=="LEMA")
fitness.data_LEMA <- fitness_orig %>% filter(Plant_Simple=="LEMA")

fitness.data_LEMA$Plot %>% unique()
fitness.data_LEMA$hete_motif %>% unique()

################
# CHFU

fitness_orig_CHFU <- fitness_orig %>% filter(Plant_Simple=="CHFU")
fitness.data_CHFU <- fitness.data %>% filter(Plant_Simple=="CHFU")
fitness.data_CHFU$Plot %>% unique()

# Readjust factor levels
fitness_orig_CHFU$Plot <- as.factor(as.numeric(fitness_orig_CHFU$Plot))
fitness.data_CHFU$Plot <- as.factor(as.numeric(fitness.data_CHFU$Plot))

################
# PUPA

fitness_orig_PUPA <- fitness_orig %>% filter(Plant_Simple=="PUPA")
fitness.data_PUPA <- fitness.data %>% filter(Plant_Simple=="PUPA")

fitness.data_PUPA$Plot %>% unique() 


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

##################
# VISITS MODELS
##################

#################################3
# FINAL MODELS

LEMA_LIN_visits <- lm(log(Seeds_GF) ~ scale(visits_GF)+
                                      Plot,
                                    data = fitness.data_LEMA)
CHFU_LIN_visits <- lm(log(Seeds_GF) ~ scale(visits_GF)+
                                      Plot,
                                    data = fitness.data_CHFU)
PUPA_LIN_visits <- lm(log(Seeds_GF) ~ scale(visits_GF)+
                                      Plot,
                                    data = fitness.data_PUPA)


summary(LEMA_LIN_visits) #R2_adj = 0.4928 (R2_adj_ref = 0.4953)
summary(CHFU_LIN_visits) #R2_adj = 0.167  (R2_adj_ref = 0.1532)
summary(PUPA_LIN_visits) #R2_adj = 0.1266 (R2_adj_ref = 0.1577)

car::vif(LEMA_LIN_visits) # OK
car::vif(CHFU_LIN_visits) # OK
car::vif(PUPA_LIN_visits) # OK


##################################3
# FINAL MODELS

LEMA_LIN_multi_visits <- lm(log(Seeds_GF) ~ scale(visits_GF) + scale(homo_motif) +
                                      scale(hete_motif) + 
                                      scale(StrengthIn) + scale(Ratio)+
                                      Plot,
                                    data = fitness.data_LEMA)
CHFU_LIN_multi_visits <- lm(log(Seeds_GF) ~ scale(visits_GF) + scale(homo_motif) +
                                      scale(hete_motif) + 
                                      scale(StrengthIn) + scale(Ratio)+
                                      Plot,
                                    data = fitness.data_CHFU)
PUPA_LIN_multi_visits <- lm(log(Seeds_GF) ~ scale(visits_GF) + scale(homo_motif) + 
                                      scale(hete_motif) + 
                                      scale(StrengthIn) + scale(Ratio)+
                                      Plot,
                                    data = fitness.data_PUPA)


summary(LEMA_LIN_multi_visits) #R2_adj = 0.4936 (R2_adj_ref = 0.4953)
summary(CHFU_LIN_multi_visits) #R2_adj = 0.1552  (R2_adj_ref = 0.1532)
summary(PUPA_LIN_multi_visits) #R2_adj = 0.1674 (R2_adj_ref = 0.1577)


car::vif(LEMA_LIN_multi_visits) # Homo_motif 2.629673*2.629673 = 6.91518
car::vif(CHFU_LIN_multi_visits) # Homo_motif and visits_GF show collinearity
car::vif(PUPA_LIN_multi_visits) # OK


jtools::summ(LEMA_LIN_multi_visits,confint = TRUE,digits = 3)
jtools::summ(CHFU_LIN_multi_visits,confint = TRUE,digits = 3)
jtools::summ(PUPA_LIN_multi_visits,confint = TRUE,digits = 3)

drop1(LEMA_LIN_multi_visits, test = "F") 
drop1(CHFU_LIN_multi_visits, test = "F")
drop1(PUPA_LIN_multi_visits, test = "F")


############

AIC(LEMA_LIN_intercept_Plot_Plant,LEMA_LIN_visits,LEMA_LIN_multi_visits)
AIC(CHFU_LIN_intercept_Plot_Plant,CHFU_LIN_visits,CHFU_LIN_multi_visits)
AIC(PUPA_LIN_intercept_Plot_Plant,PUPA_LIN_visits,PUPA_LIN_multi_visits)

############################
#Visualization of slopes by using visreg

p1 <- visreg(LEMA_LIN_visits,"visits_GF",xlab="visits",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="L. maroccanus")
p2 <- visreg(CHFU_LIN_visits,"visits_GF",xlab="visits",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="C. fuscatum")
p3 <- visreg(PUPA_LIN_visits,"visits_GF",xlab="visits",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="P. paludosa")


library(gridExtra)
library(gtable)
library(grid)


grid.arrange(
  grobs = list(p2, p1, p3),
  widths = c(1, 1, 1),
  layout_matrix = rbind(c(1, 2, 3)),
  left = textGrob("log(Seeds per individual)", rot = 90, vjust = 1),
  up = textGrob("Covariate: visits", rot = 90, vjust = 1)
)


############
# Visualization covariate model

############################
#Visualization of slopes by using visreg

p1 <- visreg(LEMA_LIN_multi_visits,"homo_motif",xlab="Homo triplet",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="L. maroccanus")
p2 <- visreg(CHFU_LIN_multi_visits,"homo_motif",xlab="Homo triplet",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="C. fuscatum")
p3 <- visreg(PUPA_LIN_multi_visits,"homo_motif",xlab="Homo triplet",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="P. paludosa")


p5 <- visreg(LEMA_LIN_multi_visits,"hete_motif",xlab="Hetero triplet",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="L. maroccanus")
p6 <- visreg(CHFU_LIN_multi_visits,"hete_motif",xlab="Hetero triplet",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="C. fuscatum")
p7 <- visreg(PUPA_LIN_multi_visits,"hete_motif",xlab="Hetero triplet",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="P. paludosa")


p9 <- visreg(LEMA_LIN_multi_visits,"StrengthIn",xlab="Within layer\n centrality",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="L. maroccanus")
p10 <- visreg(CHFU_LIN_multi_visits,"StrengthIn",xlab="Within layer\n centrality",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="C. fuscatum")
p11 <- visreg(PUPA_LIN_multi_visits,"StrengthIn",xlab="Within layer\n centrality",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="P. paludosa")


p13 <- visreg(LEMA_LIN_multi_visits,"Ratio",xlab="Among layer\n centrality ratio",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="L. maroccanus")
p14 <- visreg(CHFU_LIN_multi_visits,"Ratio",xlab="Among layer\n centrality ratio",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="C. fuscatum")
p15 <- visreg(PUPA_LIN_multi_visits,"Ratio",xlab="Among layer\n centrality ratio",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="P. paludosa")


p17 <- visreg(LEMA_LIN_multi_visits,"visits_GF",xlab="visits",ylab=NULL,gg = TRUE,  partial=T, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="L. maroccanus")

p18 <- visreg(CHFU_LIN_multi_visits,"visits_GF",xlab="visits",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="C. fuscatum")
p19 <- visreg(PUPA_LIN_multi_visits,"visits_GF",xlab="visits",ylab=NULL,gg = TRUE, partial=TRUE, rug=FALSE)+
  theme_bw()+#geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="P. paludosa")



grid.arrange(
  grobs = list(p2, p1, p3, p6,p5,p7,p10,p9,
               p11, p14, p13,p15,p18,p17,p19),
  widths = c(1, 1, 1),
  layout_matrix = rbind(c(1, 2, 3),
                        c(4,5, 6),
                        c(7, 8, 9),
                        c(10, 11, 12),
                        c(13, 14, 15)),
  left = textGrob("log(Seeds per individual)", rot = 90, vjust = 1)
)

