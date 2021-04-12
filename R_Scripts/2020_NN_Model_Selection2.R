
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
ggplot(fitness_orig %>% filter(DegreeIn>0),aes(x=Fruit_GF,y=Seeds_GF))+
  geom_point(alpha=0.3)+
  facet_wrap(~Plant_Simple)+
  labs(title = "With visits")

ggplot(fitness_orig %>% filter(DegreeIn==0),aes(x=Fruit_GF,y=Seeds_GF))+
  geom_point(alpha=0.3)+
  facet_wrap(~Plant_Simple)+
  labs(title = "Without visits")

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
###################

LEMA_NB_intercept_Plot_Plant <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                             scale(hete_motif) +
                                             scale(StrengthIn) + scale(Ratio) +
                                             (1|Plot),
                                           ziformula = ~1,
                                           family = nbinom2(),
                                           data = fitness_orig_LEMA)


CHFU_NB_intercept_Plot_Plant <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                             scale(hete_motif) +
                                             scale(StrengthIn) + scale(Ratio) +
                                             (1|Plot),
                                           ziformula = ~1,
                                           family = nbinom2(),
                                           data = fitness_orig_CHFU)

PUPA_NB_intercept_Plot_Plant <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                             scale(hete_motif) +
                                             scale(StrengthIn) + scale(Ratio) +
                                             (1|Plot),
                                           ziformula = ~1,
                                           family = nbinom2(),
                                           data = fitness_orig_PUPA)

summary(LEMA_NB_intercept_Plot_Plant)
summary(CHFU_NB_intercept_Plot_Plant)
summary(PUPA_NB_intercept_Plot_Plant)

res_LEMA_NB_intercept_Plot_Plant <- simulateResiduals(fittedModel = LEMA_NB_intercept_Plot_Plant, n = 500)
res_CHFU_NB_intercept_Plot_Plant <- simulateResiduals(fittedModel = CHFU_NB_intercept_Plot_Plant, n = 500)
res_PUPA_NB_intercept_Plot_Plant <- simulateResiduals(fittedModel = PUPA_NB_intercept_Plot_Plant, n = 500)

testZeroInflation(res_LEMA_NB_intercept_Plot_Plant)
testDispersion(res_LEMA_NB_intercept_Plot_Plant)

testZeroInflation(res_CHFU_NB_intercept_Plot_Plant)
testDispersion(res_CHFU_NB_intercept_Plot_Plant)

testZeroInflation(res_PUPA_NB_intercept_Plot_Plant)
testDispersion(res_PUPA_NB_intercept_Plot_Plant)

plot(res_LEMA_NB_intercept_Plot_Plant)
plot(res_CHFU_NB_intercept_Plot_Plant)
plot(res_PUPA_NB_intercept_Plot_Plant)

performance::r2(LEMA_NB_intercept_Plot_Plant)
performance::r2(CHFU_NB_intercept_Plot_Plant)
performance::r2(PUPA_NB_intercept_Plot_Plant)



coefficients(LEMA_NB_intercept_Plot_Plant)
performance::check_collinearity(LEMA_NB_intercept_Plot_Plant,component = "conditional") # All the ( GVIF^(1/(2*Df)) )^2 < 5 
performance::check_collinearity(CHFU_NB_intercept_Plot_Plant,component = "conditional")
performance::check_collinearity(PUPA_NB_intercept_Plot_Plant,component = "conditional")

############################
jtools::summ(LEMA_LIN_intercept_Plot_Plant,confint = TRUE,digits = 3)
jtools::summ(CHFU_LIN_intercept_Plot_Plant,confint = TRUE,digits = 3)
jtools::summ(PUPA_LIN_intercept_Plot_Plant,confint = TRUE,digits = 3)
jtools::summ(PUPA_LIN_intercept_Plot_Plant_int,confint = TRUE,digits = 3)



###############################
LEMA_NB_visits<- glmmTMB(Seeds_GF ~ scale(visits_GF)+
                                          (1|Plot),
                                        ziformula = ~1,
                                        family = nbinom2(),
                                        data = fitness_orig_LEMA)
performance::r2(LEMA_NB_visits)

########################################
############################
#Visualization of slopes by using visreg

library(ggeffects)
library(scales)


#####################################
dev.off()

par(mfrow = c(4,3),mar=c(4,4,2,1)+0.5)

visreg(CHFU_NB_intercept_Plot_Plant,"homo_motif",xlab="Homo triplet",ylab="Seeds",
       main=expression(italic("C. fuscatum")),scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ homo_motif, data = fitness_orig_CHFU, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

visreg(LEMA_NB_intercept_Plot_Plant,"homo_motif",xlab="Homo triplet",ylab="Seeds",
       main=expression(italic("L. maroccanus")),scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ homo_motif, data = fitness_orig_LEMA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

visreg(PUPA_NB_intercept_Plot_Plant,"homo_motif",xlab="Homo triplet",ylab="Seeds",
       main=expression(italic("P. paludosa")),scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ homo_motif, data = fitness_orig_PUPA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

##############

visreg(CHFU_NB_intercept_Plot_Plant,"hete_motif",xlab="Hetero triplet",ylab="Seeds",
       main=expression(italic("C. fuscatum")),scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ hete_motif, data = fitness_orig_CHFU, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

visreg(LEMA_NB_intercept_Plot_Plant,"hete_motif",xlab="Hetero triplet",ylab="Seeds",
       main=expression(italic("L. maroccanus")),scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ hete_motif, data = fitness_orig_LEMA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 



visreg(PUPA_NB_intercept_Plot_Plant,"hete_motif",xlab="Hetero triplet",ylab="Seeds",
       main=expression(italic("P. paludosa")),scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ hete_motif, data = fitness_orig_PUPA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

##############
##############


visreg(CHFU_NB_intercept_Plot_Plant,"StrengthIn",xlab="Within-layer\n centrality",ylab="Seeds",
       main=expression(italic("C. fuscatum")),scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ StrengthIn, data = fitness_orig_CHFU, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 


visreg(LEMA_NB_intercept_Plot_Plant,"StrengthIn",xlab="Within-layer\n centrality",ylab="Seeds",
       main=expression(italic("L. maroccanus")),scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ StrengthIn, data = fitness_orig_LEMA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 


visreg(PUPA_NB_intercept_Plot_Plant,"StrengthIn",xlab="Within-layer\n centrality",ylab="Seeds",
       main=expression(italic("P. paludosa")),scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ StrengthIn, data = fitness_orig_PUPA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

##############
##############

visreg(CHFU_NB_intercept_Plot_Plant,"Ratio",xlab="Among-layer\n centrality ratio",ylab="Seeds",
       main=expression(italic("C. fuscatum")),scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ Ratio, data = fitness_orig_CHFU, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

visreg(LEMA_NB_intercept_Plot_Plant,"Ratio",xlab="Among-layer\n centrality ratio",ylab="Seeds",
       main=expression(italic("L. maroccanus")),scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ Ratio, data = fitness_orig_LEMA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 



visreg(PUPA_NB_intercept_Plot_Plant,"Ratio",xlab="Among-layer\n centrality ratio",ylab="Seeds",
       main=expression(italic("P. paludosa")),scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ Ratio, data = fitness_orig_PUPA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

# save 600 x 800

dev.off()

# Random intercepts
dev.off()

par(mfrow = c(1,3),mar=c(4,4,2,1)+0.5)
visreg(CHFU_NB_intercept_Plot_Plant,"Plot",xlab="Plot",ylab="Seeds",
       main=expression(italic("C. fuscatum")),scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+

visreg(LEMA_NB_intercept_Plot_Plant,"Plot",xlab="Plot",ylab="Seeds",
       main=expression(italic("L. maroccanus")),scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+


visreg(PUPA_NB_intercept_Plot_Plant,"Plot",xlab="Plot",ylab="Seeds",
       main=expression(italic("P. paludosa")),scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+

###############################################
#save 600 x 250

########################################
# COMPARISON WITH MODELS WITH VISITS
########################################

###################

LEMA_NB_intercept_Plot_Plant_vist <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                          scale(hete_motif) +
                                          scale(StrengthIn) + scale(Ratio) +
                                            scale(visits_GF)+(1|Plot),
                                        ziformula = ~1,
                                        family = nbinom2(),
                                        data = fitness_orig_LEMA)


CHFU_NB_intercept_Plot_Plant_vist <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                          scale(hete_motif) +
                                          scale(StrengthIn) + scale(Ratio) +
                                            scale(visits_GF)+(1|Plot),
                                        ziformula = ~1,
                                        family = nbinom2(),
                                        data = fitness_orig_CHFU)

PUPA_NB_intercept_Plot_Plant_vist <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                          scale(hete_motif) +
                                          scale(StrengthIn) + scale(Ratio) +
                                            scale(visits_GF)+(1|Plot),
                                        ziformula = ~1,
                                        family = nbinom2(),
                                        data = fitness_orig_PUPA)

summary(LEMA_NB_intercept_Plot_Plant_vist)
summary(CHFU_NB_intercept_Plot_Plant_vist)
summary(PUPA_NB_intercept_Plot_Plant_vist)

r2(LEMA_NB_intercept_Plot_Plant_vist)
r2(CHFU_NB_intercept_Plot_Plant_vist)
r2(PUPA_NB_intercept_Plot_Plant_vist)


res_LEMA_NB_intercept_Plot_Plant_vist <- simulateResiduals(fittedModel = LEMA_NB_intercept_Plot_Plant_vist, n = 500)
res_CHFU_NB_intercept_Plot_Plant_vist <- simulateResiduals(fittedModel = CHFU_NB_intercept_Plot_Plant_vist, n = 500)
res_PUPA_NB_intercept_Plot_Plant_vist <- simulateResiduals(fittedModel = PUPA_NB_intercept_Plot_Plant_vist, n = 500)

testZeroInflation(res_LEMA_NB_intercept_Plot_Plant_vist)
testDispersion(res_LEMA_NB_intercept_Plot_Plant_vist)

testZeroInflation(res_CHFU_NB_intercept_Plot_Plant_vist)
testDispersion(res_CHFU_NB_intercept_Plot_Plant_vist)

testZeroInflation(res_PUPA_NB_intercept_Plot_Plant_vist)
testDispersion(res_PUPA_NB_intercept_Plot_Plant_vist)

plot(res_LEMA_NB_intercept_Plot_Plant_vist)
plot(res_CHFU_NB_intercept_Plot_Plant_vist)
plot(res_PUPA_NB_intercept_Plot_Plant_vist)

performance::r2(LEMA_NB_intercept_Plot_Plant_vist)
performance::r2(LEMA_NB_intercept_Plot_Plant)
performance::r2(CHFU_NB_intercept_Plot_Plant_vist)
performance::r2(CHFU_NB_intercept_Plot_Plant)
performance::r2(PUPA_NB_intercept_Plot_Plant_vist)
performance::r2(PUPA_NB_intercept_Plot_Plant)

AIC(LEMA_NB_intercept_Plot_Plant,LEMA_NB_intercept_Plot_Plant_vist)
AIC(CHFU_NB_intercept_Plot_Plant,CHFU_NB_intercept_Plot_Plant_vist)
AIC(PUPA_NB_intercept_Plot_Plant,PUPA_NB_intercept_Plot_Plant_vist)

performance::check_collinearity(LEMA_NB_intercept_Plot_Plant_vist,component = "conditional") # All the ( GVIF^(1/(2*Df)) )^2 < 5 
performance::check_collinearity(CHFU_NB_intercept_Plot_Plant_vist,component = "conditional")
performance::check_collinearity(PUPA_NB_intercept_Plot_Plant_vist,component = "conditional")

########################################
############################
#Visualization of slopes by using visreg

library(ggeffects)
library(scales)


#####################################
dev.off()

par(mfrow = c(5,3),mar=c(3.95,4,1,1))

visreg(CHFU_NB_intercept_Plot_Plant_vist,"homo_motif",xlab="Homo triplet",ylab="Seeds",
       main="C. fuscatum",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ homo_motif, data = fitness_orig_CHFU, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

visreg(LEMA_NB_intercept_Plot_Plant_vist,"homo_motif",xlab="Homo triplet",ylab="Seeds",
       main="L. maroccanus",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ homo_motif, data = fitness_orig_LEMA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 
Seeds
visreg(PUPA_NB_intercept_Plot_Plant_vist,"homo_motif",xlab="Homo triplet",ylab="Seeds",
       main="P. paludosa",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ homo_motif, data = fitness_orig_PUPA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

##############

visreg(CHFU_NB_intercept_Plot_Plant_vist,"hete_motif",xlab="Hetero triplet",ylab="Seeds",
       main="C. fuscatum",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ hete_motif, data = fitness_orig_CHFU, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

visreg(LEMA_NB_intercept_Plot_Plant_vist,"hete_motif",xlab="Hetero triplet",ylab="Seeds",
       main="L. maroccanus",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ hete_motif, data = fitness_orig_LEMA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 



visreg(PUPA_NB_intercept_Plot_Plant_vist,"hete_motif",xlab="Hetero triplet",ylab="Seeds",
       main="P. paludosa",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ hete_motif, data = fitness_orig_PUPA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

##############
##############


visreg(CHFU_NB_intercept_Plot_Plant_vist,"StrengthIn",xlab="Within layer\n centrality",ylab="Seeds",
       main="C. fuscatum",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ StrengthIn, data = fitness_orig_CHFU, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 


visreg(LEMA_NB_intercept_Plot_Plant_vist,"StrengthIn",xlab="Within layer\n centrality",ylab="Seeds",
       main="L. maroccanus",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ StrengthIn, data = fitness_orig_LEMA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 


visreg(PUPA_NB_intercept_Plot_Plant_vist,"StrengthIn",xlab="Within layer\n centrality",ylab="Seeds",
       main="P. paludosa",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ StrengthIn, data = fitness_orig_PUPA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

##############
##############

visreg(CHFU_NB_intercept_Plot_Plant_vist,"Ratio",xlab="Among layer\n centrality ratio",ylab="Seeds",
       main="C. fuscatum",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ Ratio, data = fitness_orig_CHFU, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

visreg(LEMA_NB_intercept_Plot_Plant_vist,"Ratio",xlab="Among layer\n centrality ratio",ylab="Seeds",
       main="L. maroccanus",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ Ratio, data = fitness_orig_LEMA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 



visreg(PUPA_NB_intercept_Plot_Plant_vist,"Ratio",xlab="Among layer\n centrality ratio",ylab="Seeds",
       main="P. paludosa",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ Ratio, data = fitness_orig_PUPA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

##############

visreg(CHFU_NB_intercept_Plot_Plant_vist,"visits_GF",xlab="Visits",ylab="Seeds",
       main="C. fuscatum",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ visits_GF, data = fitness_orig_CHFU, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

visreg(LEMA_NB_intercept_Plot_Plant_vist,"visits_GF",xlab="Visits",ylab="Seeds",
       main="L. maroccanus",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ visits_GF, data = fitness_orig_LEMA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 



visreg(PUPA_NB_intercept_Plot_Plant_vist,"visits_GF",xlab="Visits",ylab="Seeds",
       main="P. paludosa",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ visits_GF, data = fitness_orig_PUPA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

# save
#
#
#
#
#
dev.off()

hist(fitness_orig_LEMA$visits_GF)
hist(fitness_orig_CHFU$visits_GF)
hist(fitness_orig_PUPA$visits_GF)

###############################################
# ONLY VISITS MODELS
###############################################

LEMA_NB_vist <- glmmTMB(Seeds_GF~ scale(visits_GF)+(1|Plot),
                                             ziformula = ~1,
                                             family = nbinom2(),
                                             data = fitness_orig_LEMA)


CHFU_NB_vist <- glmmTMB(Seeds_GF ~ scale(visits_GF)+(1|Plot),
                                             ziformula = ~1,
                                             family = nbinom2(),
                                             data = fitness_orig_CHFU)

PUPA_NB_vist <- glmmTMB(Seeds_GF ~ scale(visits_GF)+(1|Plot),
                                             ziformula = ~1,
                                             family = nbinom2(),
                                             data = fitness_orig_PUPA)


summary(LEMA_NB_vist)
summary(CHFU_NB_vist)
summary(PUPA_NB_vist)

performance::r2(LEMA_NB_vist)
performance::r2(LEMA_NB_intercept_Plot_Plant)
performance::r2(CHFU_NB_vist)
performance::r2(CHFU_NB_intercept_Plot_Plant)
performance::r2(PUPA_NB_vist)
performance::r2(PUPA_NB_intercept_Plot_Plant)


##############
dev.off()

par(mfrow = c(1,3),mar=c(4,4,2,1)+0.5)
visreg(CHFU_NB_vist,"visits_GF",xlab="Visits",ylab="Seeds",
       main="C. fuscatum",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ visits_GF, data = fitness_orig_CHFU, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

visreg(LEMA_NB_vist,"visits_GF",xlab="Visits",ylab="Seeds",
       main="L. maroccanus",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ visits_GF, data = fitness_orig_LEMA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 



visreg(PUPA_NB_vist,"visits_GF",xlab="Visits",ylab="Seeds",
       main="P. paludosa",scale="response", rug=FALSE)#,gg = TRUE, partial=TRUE)#, rug=FALSE)+
points(Seeds_GF ~ visits_GF, data = fitness_orig_PUPA, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5),
       pch = 20) 

dev.off()