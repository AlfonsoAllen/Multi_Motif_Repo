# Fit GLMMs for multi-motif study in Caracoles (2019)

library(tidyverse)
library(MASS)
library(fitdistrplus)
library(lme4)
library(lmerTest)
library(glmmTMB)
#library(bbmle) ## for AICtab
library(DHARMa)
library(performance)
library(visreg)
library(RColorBrewer)
library(nlme)
library(usdm)


#####################################
# Read data ---------------------------------------------------------------

fitness_final_aux <- read.csv(file = "NN_data_models_phenol_overlap.csv",
                         header = TRUE,
                         stringsAsFactors = FALSE) %>%
  mutate(Seeds_GF = round(Seeds_GF))


###############
# Add G_F

G_F_list <- read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv") %>%
  dplyr::select(G_F,ID_Simple) %>% rename(ID=ID_Simple) %>% unique()

G_F_list <- bind_rows(G_F_list,tibble(G_F="None",ID="None"))

# Fix "Odontomyia_sp."

G_F_list$G_F[G_F_list$ID=="Odontomyia_sp."] <- "Hoverflies"
G_F_list <- unique(G_F_list)

# Sanity check
G_F_list %>% group_by(ID) %>% count() %>% filter(n>1)

fitness_orig <- fitness_final_aux %>% dplyr::left_join(G_F_list,by = "ID")


# Turn ID, GF and Plot into factors
fitness_orig$Plot <- as.factor(fitness_orig$Plot)
fitness_orig$ID <- as.factor(fitness_orig$ID)
fitness_orig$G_F <- as.factor(fitness_orig$G_F)


# pre-analysis ------------------------------------------------------------

# remove fitness = 0

fitness.data <- subset(fitness_orig,Seeds_GF > 0) # & !Plant_Simple %in% c("ME"))
#fitness.data <- subset(fitness.data,ID!="None")

fitness.data %>% count(wt=visits_GF)

# class of every column
# fitness.data %>% map_chr(class)

fitness_LEMA <- subset(fitness.data,Plant_Simple == "LEMA") #,G_F != "None")
fitness_PUPA <- subset(fitness.data,Plant_Simple == "PUPA") #,G_F != "None")
fitness_CHFU <- subset(fitness.data,Plant_Simple == "CHFU") #,G_F != "None")

fitness_LEMA$G_F <- as.factor(fitness_LEMA$G_F)
fitness_PUPA$G_F  <- as.factor(fitness_PUPA$G_F)
fitness_CHFU$G_F  <- as.factor(fitness_CHFU$G_F)


# Data of LEMA plants with R > 1

# Number of observations without visits

nrow(fitness_LEMA %>% filter(ID!="None"))
nrow(fitness_PUPA %>% filter(ID!="None"))
nrow(fitness_CHFU %>% filter(ID!="None"))

# scale sometimes gives problems
summary(scale(fitness_LEMA$homo_motif))
summary(scale(fitness_LEMA$hete_motif))
summary(scale(fitness_LEMA$DegreeIn))
summary(scale(fitness_LEMA$Real_PR_Layer))
summary(scale(fitness_LEMA$Ratio))


summary(scale(fitness_PUPA$homo_motif))
summary(scale(fitness_PUPA$hete_motif)) # PUPA has no hete-motifs
summary(scale(fitness_PUPA$DegreeIn))
summary(scale(fitness_PUPA$Real_PR_Layer))
summary(scale(fitness_PUPA$Ratio))


summary(scale(fitness_CHFU$homo_motif))
summary(scale(fitness_CHFU$hete_motif))
summary(scale(fitness_CHFU$DegreeIn))
summary(scale(fitness_CHFU$Real_PR_Layer))
summary(scale(fitness_CHFU$Ratio))


###############################
# ALL PLANT SPECIES MODEL: 2 Centrality Index + Homo_motif + Hetero_motifs
# (WITH and WITHOUT ZERO INFLATION FACTOR)
###############################

#CHECK--------
vif(as.data.frame(dplyr::select(fitness.data,homo_motif,hete_motif,StrengthIn, Ratio)))
cor(as.data.frame(dplyr::select(fitness.data,homo_motif,hete_motif,StrengthIn, Ratio)))


GF_MIX_NB_intercept_Plot_Plant <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                            scale(hete_motif) + 
                                            scale(StrengthIn) + scale(Ratio) +
                                            (1|Plot) +(1|Plant_Simple) ,
                                          #ziformula = ~1,
                                          family = nbinom2(),
                                          data = fitness.data)

GF_MIX_LIN_intercept_Plot_Plant <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                             scale(hete_motif) + 
                                             scale(StrengthIn) + scale(Ratio) +
                                             (1|Plot) +(1|Plant_Simple) ,
                                           #ziformula = ~1,
                                           family = gaussian(),
                                           data = fitness.data)

GF_MIX_LIN_intercept_Plot_Plant0 <- lmer(log(Seeds_GF) ~ scale(homo_motif) +
                                             scale(hete_motif) +
                                             scale(StrengthIn) + scale(Ratio) +
                                             (1|Plot) + (1|Plant_Simple) ,
                                           #ziformula = ~1,
                                           #family = gaussian(),
                                           data = fitness.data)

GF_MIX_NB_intercept_Plot_Plant_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                               scale(hete_motif) + 
                                               scale(StrengthIn) + scale(Ratio) +
                                               (1|Plot) +(1|Plant_Simple) ,
                                             ziformula = ~1,
                                             family = nbinom2(),
                                             data = fitness.data)

GF_MIX_LIN_intercept_Plot_Plant_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                                scale(hete_motif) + 
                                                scale(StrengthIn) + scale(Ratio) +
                                                (1|Plot) +(1|Plant_Simple) ,
                                              ziformula = ~1,
                                              family = gaussian(),
                                              data = fitness.data)

# Summaries

summary(GF_MIX_NB_intercept_Plot_Plant)
summary(GF_MIX_LIN_intercept_Plot_Plant)
summary(GF_MIX_LIN_intercept_Plot_Plant0)
summary(GF_MIX_NB_intercept_Plot_Plant_ZI)
summary(GF_MIX_LIN_intercept_Plot_Plant_ZI)

#############################################
# EXPLORING LMMs hypotheses for GF_MIX_LIN_intercept_Plot_Plant0 model

# linearity : OK
plot(log(fitness.data$Seeds_GF),residuals(GF_MIX_LIN_intercept_Plot_Plant0))

# homoscedasticity: ?? There are 2 stripes on the left

plot(fitted(GF_MIX_LIN_intercept_Plot_Plant0), resid(GF_MIX_LIN_intercept_Plot_Plant0, type = "pearson"))
abline(0,0, col="red")

ggplot(fitness.data)+
  geom_point(aes(x=fitted(GF_MIX_LIN_intercept_Plot_Plant0),
                 y=resid(GF_MIX_LIN_intercept_Plot_Plant0, type = "pearson"),
                 color=Plot))

ggplot(fitness.data)+
  geom_point(aes(x=fitted(GF_MIX_LIN_intercept_Plot_Plant0),
                 y=resid(GF_MIX_LIN_intercept_Plot_Plant0, type = "pearson"),
                 color=Plant_Simple))

fitness.data$RES<- residuals(GF_MIX_LIN_intercept_Plot_Plant0) #extracts the residuals and places them in a new column in our original data table
fitness.data$Abs.Res <-abs(fitness.data$RES) #creates a new column with the absolute value of the residuals
fitness.data$Model.Res2 <- fitness.data$Abs.Res^2 #squares the absolute values of the residuals to provide the more robust estimate

Levene.Model.Plot<- lm(Model.Res2 ~ Plot, data=fitness.data) #ANOVA of the squared residuals
anova(Levene.Model.Plot) #Residuals are different for each plot

Levene.Model.Plant<- lm(Model.Res2 ~ Plant_Simple, data=fitness.data) #ANOVA of the squared residuals
anova(Levene.Model.Plant) #Residuals are different for each plant



# normality of residuals: bufff... left side departures from normality
# When ME is removed, normality appears.
qqnorm(resid(GF_MIX_LIN_intercept_Plot_Plant0)) 
qqline(resid(GF_MIX_LIN_intercept_Plot_Plant0), col = "red") # add a perfect fit line

qplot(color=Plant_Simple,sample = RES, data = fitness.data)
qplot(color=Plot,sample = RES, data = fitness.data)

ggplot(fitness.data, aes(sample=fitness.data$RES))+stat_qq()

# normality of residuals random effects: ~OK
qqnorm(ranef(GF_MIX_LIN_intercept_Plot_Plant0)$Plot[,1] )
qqline(ranef(GF_MIX_LIN_intercept_Plot_Plant0)$Plot[,1], col = "red")

# normality of residuals random effects: CHMI   0.4593920
qqnorm(ranef(GF_MIX_LIN_intercept_Plot_Plant0)$Plant_Simple[,1] )
qqline(ranef(GF_MIX_LIN_intercept_Plot_Plant0)$Plant_Simple[,1], col = "red")

###################################################
# Simulating residuals with Dahrma (glmmTMB models)

res_GF_MIX_NB_intercept_Plot_Plant <- simulateResiduals(fittedModel = GF_MIX_NB_intercept_Plot_Plant, n = 500)
res_GF_MIX_LIN_intercept_Plot_Plant <- simulateResiduals(fittedModel = GF_MIX_LIN_intercept_Plot_Plant, n = 500)
res_GF_MIX_NB_intercept_Plot_Plant_ZI <- simulateResiduals(fittedModel = GF_MIX_NB_intercept_Plot_Plant_ZI, n = 500)
res_GF_MIX_LIN_intercept_Plot_Plant_ZI <- simulateResiduals(fittedModel = GF_MIX_LIN_intercept_Plot_Plant_ZI, n = 500)

# Checking Residuals 
plot(res_GF_MIX_NB_intercept_Plot_Plant) #KS + 3 Quantile deviations
plot(res_GF_MIX_NB_intercept_Plot_Plant_ZI) #KS + 3 Quantile deviations
plot(res_GF_MIX_LIN_intercept_Plot_Plant) #KS + 3 Quantile deviations
plot(res_GF_MIX_LIN_intercept_Plot_Plant_ZI) #KS + 3 Quantile deviations

AIC(
  GF_MIX_NB_intercept_Plot_Plant,
  GF_MIX_NB_intercept_Plot_Plant_ZI,
  GF_MIX_LIN_intercept_Plot_Plant,
  GF_MIX_LIN_intercept_Plot_Plant_ZI)


r2(GF_MIX_NB_intercept_Plot_Plant)
r2(GF_MIX_LIN_intercept_Plot_Plant)
r2(GF_MIX_NB_intercept_Plot_Plant_ZI)
r2(GF_MIX_LIN_intercept_Plot_Plant_ZI)


# Conditional plots

visreg(GF_MIX_LIN_intercept_Plot_Plant, "homo_motif", by="Plant_Simple",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(size=1.5, alpha=0.2, shape=16)+
  facet_wrap(vars(Plant_Simple),nrow = 1,ncol = 5)+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Homo-triplet", y = "Ln(Seeds per individual)",fill=NULL,color=NULL)

visreg(GF_MIX_LIN_intercept_Plot_Plant, "hete_motif", by="Plant_Simple",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(size=1.5, alpha=0.2, shape=16)+
  facet_wrap(vars(Plant_Simple),nrow = 2,ncol = 3)+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Hetero-triplet", y = "Ln(Seeds per individual)",fill=NULL,color=NULL)

visreg(GF_MIX_LIN_intercept_Plot_Plant, "StrengthIn", by="Plant_Simple",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(size=1.5, alpha=0.2, shape=16)+
  facet_wrap(vars(Plant_Simple),nrow = 2,ncol = 3)+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="In-Strength", y = "Ln(Seeds per individual)",fill=NULL,color=NULL)

visreg(GF_MIX_LIN_intercept_Plot_Plant, "Ratio", by="Plant_Simple",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(size=1.5, alpha=0.2, shape=16)+
  facet_wrap(vars(Plant_Simple),nrow = 2,ncol = 3)+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Ratio", y = "Seeds per individual",fill=NULL,color=NULL)

plot_labs <-c(
  "Plot 1",
  "Plot 2",
  "Plot 3",
  "Plot 4",
  "Plot 5",
  "Plot 6",
  "Plot 7",
  "Plot 8",
  "Plot 9"
)
names(plot_labs) <- c(
  '1',
  '2',
  '3',
  '4',
  '5',
  "6",
  '7',
  '8',
  "9"
)

visreg(GF_MIX_LIN_intercept_Plot_Plant, "hete_motif", by="Plot",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(size=1.5, alpha=0.2, shape=16)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Hetero-triplet", y = "Ln(Seeds per individual)",fill=NULL,color=NULL)




####################################################
# EFFECTS OF TYPE OF VISITOR IN HOMO-HETERO MOTIFS
####################################################

fitness.data$type_vist <- NA

fitness.data %>% group_by(G_F) %>% count()

mutualist_list <- c("Bees","Bees/sirphid","Hoverflies")
antagonist_list <- c("Big_beetles","Flower_beetles")
#other_list <- c("Small_flies","Small_beetles","Humbleflies","Beetles","House_flies")

fitness.data$type_vist[fitness.data$G_F %in% mutualist_list] <- "mutualist"
fitness.data$type_vist[fitness.data$G_F %in% antagonist_list] <- "antagonist"
fitness.data$type_vist[!fitness.data$G_F %in% c(mutualist_list,antagonist_list,"None")] <- "other"
#fitness.data$type_vist[fitness.data$G_F=="None"] <- "none"

# CHARACTERIZATION OF VISITOR TYPES PER PLANT SPECIES
fitness.data %>% group_by(type_vist) %>% count()
fitness.data %>% group_by(type_vist,Plant_Simple) %>% count() %>% filter(type_vist=="mutualist") # None and other dominates 
fitness.data %>% group_by(type_vist,Plant_Simple) %>% count() %>% filter(type_vist=="antagonist") # None and other dominates 



ggplot(fitness.data %>% filter(G_F!="None"), aes(fill=type_vist, y=visits_GF, x=Plant_Simple)) + 
  geom_bar(position="stack", stat="identity")+ theme_bw()+
  scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Plant species", y = "Number of visits",fill=NULL)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  theme(legend.position="bottom")

#########################################3
# LINEAL--------RD INTERCEPT---------NO-ZI

# The following models show convergence problems

GF_MIX_NB_type_Plot0 <- glmmTMB((Seeds_GF) ~ scale(homo_motif)*type_vist + 
                                 scale(hete_motif)*type_vist +
                                 (1|Plot)+
                                 (1|Plant_Simple),
                               #ziformula = ~1,
                               family = nbinom2(),
                               data = fitness.data %>%
                                 filter(!is.na(type_vist)))

GF_MIX_NB_type_Plot <-   glmer.nb((Seeds_GF) ~ scale(homo_motif)*type_vist + #Convergence failed
                                  scale(hete_motif)*type_vist +
                                  (1|Plot)+
                                  (1|Plant_Simple),
                                #ziformula = ~1,
                                #family = nbinom2(),
                                data = fitness.data %>%
                                  filter(!is.na(type_vist)))


GF_MIX_LIN_type_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif)*type_vist + 
                                  scale(hete_motif)*type_vist +
                                  (1|Plot)+
                                  (1|Plant_Simple),
                                #ziformula = ~1,
                                family = gaussian(),
                                data = fitness.data %>%
                                  filter(!is.na(type_vist)))

GF_MIX_LIN_type_Plot0 <- lmer(log(Seeds_GF) ~ scale(homo_motif)*type_vist + #boundary (singular)
                                  scale(hete_motif)*type_vist +
                                  (1|Plot)+
                                  (1|Plant_Simple),
                                #ziformula = ~1,
                                #family = gaussian(),
                                data = fitness.data %>%
                                  filter(!is.na(type_vist)))

summary(GF_MIX_LIN_type_Plot) 

GF_MIX_LIN_type_Plot <- lmer(log(Seeds_GF) ~ scale(homo_motif)*type_vist + 
                               scale(hete_motif)*type_vist +
                               (1|Plot)+
                               (1|Plant_Simple),
                             #ziformula = ~1,
                             #family = gaussian(),
                             data = fitness.data %>%
                               filter(!is.na(type_vist)))

summary(GF_MIX_NB_type_Plot0)
summary(GF_MIX_LIN_type_Plot)

res_GF_MIX_NB_type_Plot <- simulateResiduals(fittedModel = GF_MIX_NB_type_Plot0, n = 500)
res_GF_MIX_LIN_type_Plot <- simulateResiduals(fittedModel = GF_MIX_LIN_type_Plot, n = 500)

# Checking Residuals 
plot(res_GF_MIX_NB_type_Plot) #KS + 3 Quantile deviations
plot(res_GF_MIX_LIN_type_Plot)#KS + 3 Quantile deviations

AIC(
  GF_MIX_NB_type_Plot0,
  GF_MIX_LIN_type_Plot)


r2(GF_MIX_NB_type_Plot0)
r2(GF_MIX_LIN_type_Plot)


# Results for lineal model
visreg(GF_MIX_LIN_type_Plot, "homo_motif", by="type_vist",gg = TRUE,
       overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(aes(color=fitness.data%>% filter(!is.na(type_vist))%>% 
                   dplyr::select(Plant_Simple) %>% pull()),
             size=1.5, alpha=0.5, shape=16)+
  facet_wrap(vars(type_vist),nrow = 5,ncol = 3)+
  #scale_color_brewer(palette = 'Paired')+ 
  labs(x ="Homo-triplet", y = "Ln(Seeds per individual)",fill=NULL,color=NULL)+
  theme(legend.position="bottom")

visreg(GF_MIX_LIN_type_Plot, "hete_motif", by="type_vist",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(aes(color=fitness.data%>% filter(!is.na(type_vist))%>% 
                   dplyr::select(Plant_Simple) %>% pull()),
             size=1.5, alpha=0.5, shape=16)+
  #scale_color_brewer(palette = 'Paired')+ 
  facet_wrap(vars(type_vist),nrow = 2,ncol = 3)+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Hetero-triplet", y = "Ln(Seeds per individual)",fill=NULL,color=NULL)+
  theme(legend.position="bottom")


visreg(GF_MIX_LIN_type_Plot, "homo_motif", by="Plant_Simple",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(size=1.5, alpha=0.2, shape=16)+
  facet_wrap(vars(Plant_Simple),nrow = 2,ncol = 3)+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Homo-triplet", y = "Ln(Seeds per individual)",fill=NULL,color=NULL)

visreg(GF_MIX_LIN_type_Plot, "hete_motif", by="Plant_Simple",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(size=1.5, alpha=0.2, shape=16)+
  facet_wrap(vars(Plant_Simple),nrow = 2,ncol = 3)+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Hetero-triplet", y = "Ln(Seeds per individual)",fill=NULL,color=NULL)

plot_labs <-c(
  "Plot 1",
  "Plot 2",
  "Plot 3",
  "Plot 4",
  "Plot 5",
  "Plot 6",
  "Plot 7",
  "Plot 8",
  "Plot 9"
)
names(plot_labs) <- c(
  '1',
  '2',
  '3',
  '4',
  '5',
  "6",
  '7',
  '8',
  "9"
)

visreg(GF_MIX_LIN_type_Plot, "homo_motif", by="Plot",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(size=1.5, alpha=0.2, shape=16)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Homo-triplet", y = "Ln(Seeds per individual)",fill=NULL,color=NULL)


visreg(GF_MIX_LIN_type_Plot, "hete_motif", by="Plot",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(size=1.5, alpha=0.2, shape=16)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Hetero-triplet", y = "Ln(Seeds per individual)",fill=NULL,color=NULL)


#################################
# EXPLORING EFFECT OF A SINGLE TYPE OF VISITOR
#################################

GF_MIX_LIN_mutua_Plot <- lmer(log(Seeds_GF) ~ scale(homo_motif) + 
                               scale(hete_motif) +
                               (1|Plot)+
                               (1|Plant_Simple),
                             #ziformula = ~1,
                             #family = gaussian(),
                             data = fitness.data %>%
                               filter(type_vist=="mutualist"))

GF_MIX_LIN_antag_Plot <- lmer(log(Seeds_GF) ~ scale(homo_motif) + 
                                scale(hete_motif) +
                                (1|Plot),
                              #ziformula = ~1,
                              #family = gaussian(),
                              data = fitness.data %>%
                                filter(type_vist=="antagonist"))

GF_MIX_LIN_other_Plot <- lmer(log(Seeds_GF) ~ scale(homo_motif) + 
                                scale(hete_motif) +
                                (1|Plot),
                              #ziformula = ~1,
                              #family = gaussian(),
                              data = fitness.data %>%
                                filter(type_vist=="other"))

summary(GF_MIX_LIN_mutua_Plot)
summary(GF_MIX_LIN_antag_Plot)
summary(GF_MIX_LIN_other_Plot)


#########################################
# INDIVIDUAL MODEL FOR EACH PLANT SPECIES
#########################################

vif(as.data.frame(dplyr::select(fitness_PUPA,homo_motif,hete_motif)))
cor.test(fitness_PUPA$homo_motif, fitness_PUPA$hete_motif, method=c("pearson", "kendall", "spearman"))

vif(as.data.frame(dplyr::select(fitness_LEMA,DegreeIn,homo_motif,hete_motif,Ratio)))
vif(as.data.frame(dplyr::select(fitness_PUPA,DegreeIn,homo_motif,Ratio)))
vif(as.data.frame(dplyr::select(fitness_CHFU,DegreeIn,homo_motif,hete_motif,Ratio)))
vif(as.data.frame(dplyr::select(fitness_CHFU,DegreeIn,homo_motif,hete_motif,Ratio)))

vif(as.data.frame(dplyr::select(fitness_LEMA,DegreeIn,Real_PR_Layer,DegreeIn,Ratio)))
vif(as.data.frame(dplyr::select(fitness_LEMA,DegreeIn,homo_motif,hete_motif)))
vif(as.data.frame(dplyr::select(fitness_LEMA,Real_PR_Multi,homo_motif,hete_motif,StrengthIn, Ratio)))

############################
#INDIVIDUAL PLANT MODELS
############################

# NB-------RD INTERCEPT---------NO-ZI


GF_LEMA_NB_intercept_Plot <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                       scale(hete_motif) + 
                                       scale(StrengthIn) + scale(Ratio) +
                                       (1|Plot) ,
                                     #ziformula = ~1,
                                     family = nbinom1(),
                                     data = fitness_LEMA)

GF_PUPA_NB_intercept_Plot <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                                       #scale(hete_motif) + 
                                       scale(StrengthIn) + scale(Ratio) +
                                       (1|Plot) ,
                                     #ziformula = ~1,
                                     family = nbinom1(),
                                     data = fitness_PUPA)

GF_CHFU_NB_intercept_Plot <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                       scale(hete_motif) + 
                                       scale(StrengthIn) + scale(Ratio) +
                                       (1|Plot) ,
                                     #ziformula = ~1,
                                     family = nbinom1(),
                                     data = fitness_CHFU)



#########################################3
# LINEAL--------RD INTERCEPT---------NO-ZI


GF_LEMA_LIN_intercept_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                        scale(hete_motif) + 
                                        scale(StrengthIn) + scale(Ratio) +
                                        (1|Plot) ,
                                      #ziformula = ~1,
                                      family = gaussian(),
                                      data = fitness_LEMA)

GF_PUPA_LIN_intercept_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                        #scale(hete_motif) + 
                                        scale(StrengthIn) + scale(Ratio) +
                                        (1|Plot) ,
                                      #ziformula = ~1,
                                      family = gaussian(),
                                      data = fitness_PUPA)

GF_CHFU_LIN_intercept_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                        scale(hete_motif) + 
                                        scale(StrengthIn) + scale(Ratio) +
                                        (1|Plot) ,
                                      #ziformula = ~1,
                                      family = gaussian(),
                                      data = fitness_CHFU)



# NB--------RD SLOPE---------ZI


GF_LEMA_NB_intercept_Plot_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                                          scale(hete_motif) + 
                                          scale(StrengthIn) + scale(Ratio) +
                                          (1|Plot) ,
                                        ziformula = ~1,
                                        family = nbinom1(),
                                        data = fitness_LEMA)

GF_PUPA_NB_intercept_Plot_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) + #Conv problems
                                          #scale(hete_motif) + 
                                          scale(StrengthIn) + scale(Ratio) +
                                          (1|Plot) ,
                                        ziformula = ~1,
                                        family = nbinom1(),
                                        data = fitness_PUPA)

GF_CHFU_NB_intercept_Plot_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                          scale(hete_motif) + 
                                          scale(StrengthIn) + scale(Ratio) +
                                          (1|Plot) ,
                                        ziformula = ~1,
                                        family = nbinom1(),
                                        data = fitness_CHFU)


# LINEAL--------RD SLOPE---------ZI


GF_LEMA_LIN_intercept_Plot_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                           scale(hete_motif) + 
                                           scale(StrengthIn) + scale(Ratio) +
                                           (1|Plot) ,
                                         ziformula = ~1,
                                         family = gaussian(),
                                         data = fitness_LEMA)

GF_PUPA_LIN_intercept_Plot_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + #Conv probl.
                                           #scale(hete_motif) + 
                                           scale(StrengthIn) + scale(Ratio) +
                                           (1|Plot) ,
                                         ziformula = ~1,
                                         family = gaussian(),
                                         data = fitness_PUPA)

GF_CHFU_LIN_intercept_Plot_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                           scale(hete_motif) + 
                                           scale(StrengthIn) + scale(Ratio) +
                                           (1|Plot) ,
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
plot(res_GF_LEMA_NB_intercept_Plot)#KS +3 Quant desv
plot(res_GF_LEMA_LIN_intercept_Plot)#KS +2 Quant desv
plot(res_GF_LEMA_NB_intercept_Plot_ZI)#KS +3 Quant desv
plot(res_GF_LEMA_LIN_intercept_Plot_ZI)#KS +2 Quant desv


plot(res_GF_PUPA_NB_intercept_Plot)
plot(res_GF_PUPA_LIN_intercept_Plot)
plot(res_GF_PUPA_NB_intercept_Plot_ZI)
plot(res_GF_PUPA_LIN_intercept_Plot_ZI)

plot(res_GF_CHFU_NB_intercept_Plot)#3 Quant desv
plot(res_GF_CHFU_LIN_intercept_Plot)#KS + 3 Quant desv
plot(res_GF_CHFU_NB_intercept_Plot_ZI)#3 Quant desv
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

fitness_LEMA$type_vist <- NA
fitness_PUPA$type_vist <- NA
fitness_CHFU$type_vist <- NA

fitness.data %>% group_by(G_F) %>% count()

mutualist_list <- c("Bees","Bees/sirphid","Hoverflies")
antagonist_list <- c("Big_beetles","Flower_beetles")
#other_list <- c("Small_flies","Small_beetles","Humbleflies","Beetles","House_flies")

fitness_LEMA$type_vist[fitness_LEMA$G_F %in% mutualist_list] <- "mutualist"
fitness_LEMA$type_vist[fitness_LEMA$G_F %in% antagonist_list] <- "antagonist"
fitness_LEMA$type_vist[!fitness_LEMA$G_F %in% c(mutualist_list,antagonist_list,c("None"))] <- "other"
#fitness_LEMA$type_vist[is.na(fitness_LEMA$type_vist)] <- "none"

fitness_PUPA$type_vist[fitness_PUPA$G_F %in% mutualist_list] <- "mutualist"
fitness_PUPA$type_vist[fitness_PUPA$G_F %in% antagonist_list] <- "antagonist"
fitness_PUPA$type_vist[!fitness_PUPA$G_F %in% c(mutualist_list,antagonist_list,c("None"))] <- "other"
#fitness_PUPA$type_vist[is.na(fitness_PUPA$type_vist)] <- "none"

fitness_CHFU$type_vist[fitness_CHFU$G_F %in% mutualist_list] <- "mutualist"
fitness_CHFU$type_vist[fitness_CHFU$G_F %in% antagonist_list] <- "antagonist"
fitness_CHFU$type_vist[!fitness_CHFU$G_F %in% c(mutualist_list,antagonist_list,"None")] <- "other"
#fitness_CHFU$type_vist[is.na(fitness_CHFU$type_vist)] <- "none"

fitness_LEMA %>% group_by(type_vist) %>%count(wt=visits_GF)
fitness_PUPA %>% group_by(type_vist) %>%count(wt=visits_GF)
fitness_CHFU %>% group_by(type_vist) %>%count(wt=visits_GF)

#########################################3
# LINEAL--------RD INTERCEPT---------NO-ZI


GF_LEMA_LIN_type_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif)*type_vist + 
                                   scale(hete_motif)*type_vist +
                                   (1|Plot),
                                 #ziformula = ~1,
                                 family = gaussian(),
                                 data = fitness_LEMA %>% filter(!is.na(type_vist)))

# PUPA HAS NOT ENOUGH GROUPS
GF_PUPA_LIN_type_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif)*type_vist + 
                                   #scale(hete_motif)*type_vist +
                                   (1|Plot),
                                 #ziformula = ~1,
                                 family = gaussian(),
                                 data = fitness_PUPA %>% filter(!is.na(type_vist)))

GF_CHFU_LIN_type_Plot <- lmer(log(Seeds_GF) ~ scale(homo_motif)*type_vist + 
                                   scale(hete_motif)*type_vist +
                                   (1|Plot),
                                 #ziformula = ~1,
                                 #family = gaussian(),
                                 data = fitness_CHFU %>% filter(!is.na(type_vist)))

summary(GF_LEMA_LIN_type_Plot)
summary(GF_PUPA_LIN_type_Plot) # NO group variance 
summary(GF_CHFU_LIN_type_Plot) # No significant terms

visreg(GF_LEMA_LIN_type_Plot, "homo_motif", by="type_vist",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(size=1.5, alpha=0.2, shape=16)+
  facet_wrap(vars(type_vist),nrow = 2,ncol = 3)+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Homo-triplet", y = "Ln(Seeds per individual)",fill=NULL,color=NULL)

visreg(GF_LEMA_LIN_type_Plot, "hete_motif", by="type_vist",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(size=1.5, alpha=0.2, shape=16)+
  facet_wrap(vars(type_vist),nrow = 2,ncol = 3)+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Hetero-triplet", y = "Ln(Seeds per individual)",fill=NULL,color=NULL)

visreg(GF_CHFU_LIN_type_Plot, "homo_motif", by="type_vist",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(size=1.5, alpha=0.2, shape=16)+
  facet_wrap(vars(type_vist),nrow = 2,ncol = 3)+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Homo-triplet", y = "Ln(Seeds per individual)",fill=NULL,color=NULL)

visreg(GF_CHFU_LIN_type_Plot, "hete_motif", by="type_vist",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(size=1.5, alpha=0.2, shape=16)+
  facet_wrap(vars(type_vist),nrow = 2,ncol = 3)+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Hetero-triplet", y = "Ln(Seeds per individual)",fill=NULL,color=NULL)


# Get models' Residuals

res_GF_LEMA_LIN_type_Plot <- simulateResiduals(fittedModel = GF_LEMA_LIN_type_Plot, n = 500)
res_GF_CHFU_LIN_type_Plot <- simulateResiduals(fittedModel = GF_CHFU_LIN_type_Plot, n = 500)

# Checking Residuals 
plot(res_GF_LEMA_LIN_type_Plot)# KS + 2 Quantile deviations
plot(res_GF_CHFU_LIN_type_Plot) # KS + 3 Quantile deviations

# AIC
AIC(GF_LEMA_LIN_type_Plot)
AIC(GF_CHFU_LIN_type_Plot)

r2(GF_LEMA_LIN_type_Plot)
r2(GF_CHFU_LIN_type_Plot)


##############################-----------
# ONE MODEL PER PLANT AND TYPE OF VISITOR

GF_LEMA_LIN_mutua_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                    scale(hete_motif) +
                                    (1|Plot),
                                  #ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness_LEMA %>% filter(type_vist=="mutualist"))

GF_PUPA_LIN_mutua_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                    #scale(hete_motif) +
                                    (1|Plot),
                                  #ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness_PUPA %>% filter(type_vist=="mutualist"))

# In the case of CHFU, hete-motifs are equal to zero

GF_CHFU_LIN_mutua_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                    #scale(hete_motif) +
                                    (1|Plot),
                                  #ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness_CHFU %>% filter(type_vist=="mutualist"))


GF_LEMA_LIN_antag_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                    scale(hete_motif) +
                                    (1|Plot),
                                  #ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness_LEMA %>% filter(type_vist=="antagonist"))

# In the case of PUPA, homo and hete motifs are equal for antagonist
# There is only one observation

GF_PUPA_LIN_antag_Plot <- glmmTMB(log(Seeds_GF) ~ 
                                    (1|Plot),
                                  #ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness_PUPA %>% filter(type_vist=="antagonist"))

fitness_PUPA %>% filter(type_vist=="antagonist")


GF_CHFU_LIN_antag_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                    #scale(hete_motif) +
                                    (1|Plot),
                                  #ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness_CHFU %>% filter(type_vist=="antagonist"))

fitness_CHFU %>% filter(type_vist=="antagonist")

GF_LEMA_LIN_other_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                    scale(hete_motif) +
                                    (1|Plot),
                                  #ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness_LEMA %>% filter(type_vist=="other"))

# In the case of PUPA, homo and hete motifs are equal for antagonist

fitness_PUPA %>% filter(type_vist=="other")

GF_PUPA_LIN_other_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                    #scale(hete_motif) +
                                    (1|Plot),
                                  #ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness_PUPA %>% filter(type_vist=="other"))

GF_CHFU_LIN_other_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                    scale(hete_motif) +
                                    (1|Plot),
                                  #ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness_CHFU %>% filter(type_vist=="other"))

summary(GF_LEMA_LIN_mutua_Plot)
summary(GF_LEMA_LIN_antag_Plot)
summary(GF_LEMA_LIN_other_Plot)

summary(GF_PUPA_LIN_mutua_Plot) #NO VARIANCE (GROUPING FACTOR): 2 groups
summary(GF_PUPA_LIN_antag_Plot) #NO VARIANCE (GROUPING FACTOR): 1 groups
summary(GF_PUPA_LIN_other_Plot) #NO VARIANCE (GROUPING FACTOR): 2 groups

summary(GF_CHFU_LIN_mutua_Plot) #NO VARIANCE (GROUPING FACTOR): 3 groups
summary(GF_CHFU_LIN_antag_Plot) #NO VARIANCE (GROUPING FACTOR): 2 groups
summary(GF_CHFU_LIN_other_Plot)

# Get models' Residuals

res_GF_LEMA_LIN_mutua_Plot <- simulateResiduals(fittedModel = GF_LEMA_LIN_mutua_Plot, n = 500)
res_GF_LEMA_LIN_antag_Plot <- simulateResiduals(fittedModel = GF_LEMA_LIN_antag_Plot, n = 500)
res_GF_LEMA_LIN_other_Plot <- simulateResiduals(fittedModel = GF_LEMA_LIN_other_Plot, n = 500)

res_GF_CHFU_LIN_mutua_Plot <- simulateResiduals(fittedModel = GF_CHFU_LIN_mutua_Plot, n = 500)
res_GF_CHFU_LIN_antag_Plot <- simulateResiduals(fittedModel = GF_CHFU_LIN_antag_Plot, n = 500)
res_GF_CHFU_LIN_other_Plot <- simulateResiduals(fittedModel = GF_CHFU_LIN_other_Plot, n = 500)

res_GF_PUPA_LIN_mutua_Plot <- simulateResiduals(fittedModel = GF_PUPA_LIN_mutua_Plot, n = 500)
res_GF_PUPA_LIN_antag_Plot <- simulateResiduals(fittedModel = GF_PUPA_LIN_antag_Plot, n = 500)
res_GF_PUPA_LIN_other_Plot <- simulateResiduals(fittedModel = GF_PUPA_LIN_other_Plot, n = 500)


# Checking Residuals 
plot(res_GF_LEMA_LIN_mutua_Plot)#2 Quant desv
plot(res_GF_LEMA_LIN_antag_Plot)#KS +3 Quant desv
plot(res_GF_LEMA_LIN_other_Plot)#KS +1 Quant desv

plot(res_GF_PUPA_LIN_mutua_Plot)
plot(res_GF_PUPA_LIN_antag_Plot)
plot(res_GF_PUPA_LIN_other_Plot)

plot(res_GF_CHFU_LIN_mutua_Plot)
plot(res_GF_CHFU_LIN_antag_Plot)
plot(res_GF_CHFU_LIN_other_Plot)#KS +3 Quant desv


# check convergence between models -----

AIC(GF_LEMA_LIN_mutua_Plot)
AIC( GF_LEMA_LIN_antag_Plot)
AIC( GF_LEMA_LIN_other_Plot)

AIC(GF_PUPA_LIN_mutua_Plot)
AIC(  GF_PUPA_LIN_antag_Plot)
AIC(  GF_PUPA_LIN_other_Plot)


AIC(  GF_CHFU_LIN_mutua_Plot)
AIC(  GF_CHFU_LIN_antag_Plot)
AIC(  GF_CHFU_LIN_other_Plot)


r2(GF_LEMA_LIN_mutua_Plot)
r2(GF_LEMA_LIN_antag_Plot)
r2(GF_LEMA_LIN_other_Plot)

r2(GF_PUPA_LIN_mutua_Plot)
r2(GF_PUPA_LIN_antag_Plot)
r2(GF_PUPA_LIN_other_Plot)

r2(GF_CHFU_LIN_mutua_Plot)
r2(GF_CHFU_LIN_antag_Plot)
r2(GF_CHFU_LIN_other_Plot)

###########################################
# ONE NEG. BIN. MODEL PER PLANT AND TYPE OF VISITOR 


GF_LEMA_NB_mutua_Plot <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                   scale(hete_motif) +
                                   (1|Plot),
                                 #ziformula = ~1,
                                 family = nbinom2(),
                                 data = fitness_LEMA %>% filter(type_vist=="mutualist"))

GF_PUPA_NB_mutua_Plot <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                   #scale(hete_motif) +
                                   (1|Plot),
                                 #ziformula = ~1,
                                 family = nbinom2(),
                                 data = fitness_PUPA %>% filter(type_vist=="mutualist"))

# In the case of CHFU, hete-motifs are equal to zero

GF_CHFU_NB_mutua_Plot <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                   #scale(hete_motif) +
                                   (1|Plot),
                                 #ziformula = ~1,
                                 family = nbinom2(),
                                 data = fitness_CHFU %>% filter(type_vist=="mutualist"))


GF_LEMA_NB_antag_Plot <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                   scale(hete_motif) +
                                   (1|Plot),
                                 #ziformula = ~1,
                                 family = nbinom2(),
                                 data = fitness_LEMA %>% filter(type_vist=="antagonist"))

# In the case of PUPA, homo and hete motifs are equal for antagonist
# There is only one observation

GF_PUPA_NB_antag_Plot <- glmmTMB((Seeds_GF) ~ 
                                   (1|Plot),
                                 #ziformula = ~1,
                                 family = nbinom2(),
                                 data = fitness_PUPA %>% filter(type_vist=="antagonist"))

fitness_PUPA %>% filter(type_vist=="antagonist")


GF_CHFU_NB_antag_Plot <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                   #scale(hete_motif) +
                                   (1|Plot),
                                 #ziformula = ~1,
                                 family = nbinom2(),
                                 data = fitness_CHFU %>% filter(type_vist=="antagonist"))

fitness_CHFU %>% filter(type_vist=="antagonist")

GF_LEMA_NB_other_Plot <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                   scale(hete_motif) +
                                   (1|Plot),
                                 #ziformula = ~1,
                                 family = nbinom2(),
                                 data = fitness_LEMA %>% filter(type_vist=="other"))

# In the case of PUPA, homo and hete motifs are equal for antagonist

fitness_PUPA %>% filter(type_vist=="other")

GF_PUPA_NB_other_Plot <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                   #scale(hete_motif) +
                                   (1|Plot),
                                 #ziformula = ~1,
                                 family = nbinom2(),
                                 data = fitness_PUPA %>% filter(type_vist=="other"))

GF_CHFU_NB_other_Plot <- glmmTMB((Seeds_GF) ~ scale(homo_motif) + 
                                   scale(hete_motif) +
                                   (1|Plot),
                                 #ziformula = ~1,
                                 family = nbinom2(),
                                 data = fitness_CHFU %>% filter(type_vist=="other"))

summary(GF_LEMA_NB_mutua_Plot)
summary(GF_LEMA_NB_antag_Plot)
summary(GF_LEMA_NB_other_Plot)

summary(GF_PUPA_NB_mutua_Plot) #NO VARIANCE (GROUPING FACTOR): 2 groups
summary(GF_PUPA_NB_antag_Plot) #NO VARIANCE (GROUPING FACTOR): 1 groups
summary(GF_PUPA_NB_other_Plot)

summary(GF_CHFU_NB_mutua_Plot) 
summary(GF_CHFU_NB_antag_Plot) #NO VARIANCE (GROUPING FACTOR): 2 groups
summary(GF_CHFU_NB_other_Plot)

# Get models' Residuals

res_GF_LEMA_NB_mutua_Plot <- simulateResiduals(fittedModel = GF_LEMA_NB_mutua_Plot, n = 500)
res_GF_LEMA_NB_antag_Plot <- simulateResiduals(fittedModel = GF_LEMA_NB_antag_Plot, n = 500)
res_GF_LEMA_NB_other_Plot <- simulateResiduals(fittedModel = GF_LEMA_NB_other_Plot, n = 500)

res_GF_CHFU_NB_mutua_Plot <- simulateResiduals(fittedModel = GF_CHFU_NB_mutua_Plot, n = 500)
res_GF_CHFU_NB_antag_Plot <- simulateResiduals(fittedModel = GF_CHFU_NB_antag_Plot, n = 500)
res_GF_CHFU_NB_other_Plot <- simulateResiduals(fittedModel = GF_CHFU_NB_other_Plot, n = 500)

res_GF_PUPA_NB_mutua_Plot <- simulateResiduals(fittedModel = GF_PUPA_NB_mutua_Plot, n = 500)
res_GF_PUPA_NB_antag_Plot <- simulateResiduals(fittedModel = GF_PUPA_NB_antag_Plot, n = 500)
res_GF_PUPA_NB_other_Plot <- simulateResiduals(fittedModel = GF_PUPA_NB_other_Plot, n = 500)


# Checking Residuals 
plot(res_GF_LEMA_NB_mutua_Plot)#
plot(res_GF_LEMA_NB_antag_Plot)#KS + 3 Quant desv
plot(res_GF_LEMA_NB_other_Plot)#KS +1 Quant desv

plot(res_GF_PUPA_NB_mutua_Plot)
plot(res_GF_PUPA_NB_antag_Plot)
plot(res_GF_PUPA_NB_other_Plot)

plot(res_GF_CHFU_NB_mutua_Plot)
plot(res_GF_CHFU_NB_antag_Plot)
plot(res_GF_CHFU_NB_other_Plot) #KS +3 Quant desv

testDispersion(res_GF_CHFU_NB_other_Plot)

# check convergence between models -----

AIC(GF_LEMA_NB_mutua_Plot)
AIC( GF_LEMA_NB_antag_Plot)
AIC( GF_LEMA_NB_other_Plot)

AIC(GF_PUPA_NB_mutua_Plot)
AIC(  GF_PUPA_NB_antag_Plot)
AIC(  GF_PUPA_NB_other_Plot)


AIC(  GF_CHFU_NB_mutua_Plot)
AIC(  GF_CHFU_NB_antag_Plot)
AIC(  GF_CHFU_NB_other_Plot)


r2(GF_LEMA_NB_mutua_Plot)
r2(GF_LEMA_NB_antag_Plot)
r2(GF_LEMA_NB_other_Plot)

r2(GF_PUPA_NB_mutua_Plot) #Some variance components equal zero.
r2(GF_PUPA_NB_antag_Plot) #Some variance components equal zero.
r2(GF_PUPA_NB_other_Plot)

r2(GF_CHFU_NB_mutua_Plot)
r2(GF_CHFU_NB_antag_Plot) #Some variance components equal zero.
r2(GF_CHFU_NB_other_Plot)
