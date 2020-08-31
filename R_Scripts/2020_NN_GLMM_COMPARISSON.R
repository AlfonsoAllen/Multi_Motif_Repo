
library(tidyverse)

#################################
# read data ---------------------------------------------------------------

fitness_final_aux <- read.csv(file = "2020_NN_data_models_phenol_overlap.csv",
                         header = TRUE,
                         stringsAsFactors = FALSE) %>%
  rename(Plant_Simple=Plant) %>% mutate(Seeds_GF = round(Seeds_GF))

#########################
# Add G_F

G_F_list <- read_csv2("Raw_Data/raw_Pollinators_2020_1.csv") %>%
  dplyr::select(G_F,ID_Simple) %>% unique() %>% rename(ID=ID_Simple)

G_F_list <- bind_rows(G_F_list,tibble(G_F="None",ID="None"))

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

#############################
# pre-analysis ------------------------------------------------------------

# remove fitness = 0

fitness.data <- subset(fitness_orig,Seeds_GF > 0)


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
summary(scale(fitness_PUPA$hete_motif))
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
summary(GF_MIX_NB_intercept_Plot_Plant_ZI)
summary(GF_MIX_LIN_intercept_Plot_Plant_ZI)

# Simulating residuals with Dahrma

res_GF_MIX_NB_intercept_Plot_Plant <- simulateResiduals(fittedModel = GF_MIX_NB_intercept_Plot_Plant, n = 500)
res_GF_MIX_LIN_intercept_Plot_Plant <- simulateResiduals(fittedModel = GF_MIX_LIN_intercept_Plot_Plant, n = 500)
res_GF_MIX_NB_intercept_Plot_Plant_ZI<- simulateResiduals(fittedModel = GF_MIX_NB_intercept_Plot_Plant_ZI, n = 500)
res_GF_MIX_LIN_intercept_Plot_Plant_ZI <- simulateResiduals(fittedModel = GF_MIX_LIN_intercept_Plot_Plant_ZI, n = 500)

# Checking Residuals 
plot(res_GF_MIX_NB_intercept_Plot_Plant)
plot(res_GF_MIX_NB_intercept_Plot_Plant_ZI)
plot(res_GF_MIX_LIN_intercept_Plot_Plant)
plot(res_GF_MIX_LIN_intercept_Plot_Plant_ZI)

AIC(
  GF_MIX_NB_intercept_Plot_Plant,
  GF_MIX_NB_intercept_Plot_Plant_ZI,
  GF_MIX_LIN_intercept_Plot_Plant,
  GF_MIX_LIN_intercept_Plot_Plant_ZI)

library(performance)

r2(GF_MIX_NB_intercept_Plot_Plant)
r2(GF_MIX_LIN_intercept_Plot_Plant)
r2(GF_MIX_NB_intercept_Plot_Plant_ZI)
r2(GF_MIX_LIN_intercept_Plot_Plant_ZI)


library(visreg)
visreg(GF_MIX_LIN_intercept_Plot_Plant)

visreg(GF_MIX_LIN_intercept_Plot_Plant, "homo_motif", by="Plant_Simple",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(size=1.5, alpha=0.2, shape=16)+
  facet_wrap(vars(Plant_Simple),nrow = 2,ncol = 3)+
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

fitness.data %>% group_by(type_vist,Plant_Simple) %>% count() %>% filter(type_vist=="mutualist") # None and other dominates 
fitness.data %>% group_by(type_vist,Plant_Simple) %>% count() %>% filter(type_vist=="antagonist") # None and other dominates 

library(RColorBrewer)


ggplot(fitness.data %>% filter(G_F!="None"), aes(fill=type_vist, y=visits_GF, x=Plant_Simple)) + 
  geom_bar(position="stack", stat="identity")+ theme_bw()+
  scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Plant species", y = "Number of visits",fill=NULL)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  theme(legend.position="bottom")

#########################################3
# LINEAL--------RD INTERCEPT---------NO-ZI

# The following models show convergence problems

GF_MIX_NB_type_Plot <- glmmTMB((Seeds_GF) ~ scale(homo_motif)*type_vist + 
                                   scale(hete_motif)*type_vist +
                                   (0+scale(homo_motif)+scale(hete_motif)|Plot)+
                                   (0+scale(homo_motif)+scale(hete_motif)|Plant_Simple),
                                 #ziformula = ~1,
                                 family = nbinom2(),
                                 data = fitness.data %>%
                                 filter(!is.na(type_vist)))


GF_MIX_LIN_type_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif)*type_vist + 
                                 scale(hete_motif)*type_vist +
                                  (0+scale(homo_motif)+scale(hete_motif)|Plot)+
                                  (0+scale(homo_motif)+scale(hete_motif)|Plant_Simple),
                               #ziformula = ~1,
                               family = gaussian(),
                               data = fitness.data %>%
                                 filter(!is.na(type_vist)))

summary(GF_MIX_NB_type_Plot)
summary(GF_MIX_LIN_type_Plot)

res_GF_MIX_NB_type_Plot <- simulateResiduals(fittedModel = GF_MIX_NB_type_Plot, n = 500)
res_GF_MIX_LIN_type_Plot <- simulateResiduals(fittedModel = GF_MIX_LIN_type_Plot, n = 500)

# Checking Residuals 
plot(res_GF_MIX_NB_type_Plot)
plot(res_GF_MIX_NB_type_Plot)

AIC(
  GF_MIX_NB_type_Plot,
  GF_MIX_LIN_type_Plot)

library(performance)

r2(GF_MIX_NB_type_Plot)
r2(GF_MIX_LIN_type_Plot)


library(visreg)
library(RColorBrewer)


visreg(GF_MIX_LIN_type_Plot, "homo_motif", by="type_vist",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(aes(color=fitness.data%>% filter(!is.na(type_vist))%>% 
                   dplyr::select(Plant_Simple) %>% pull()),
             size=1.5, alpha=0.5, shape=16)+
  facet_wrap(vars(type_vist),nrow = 5,ncol = 3)+
  scale_color_brewer(palette = 'Paired')+ 
  labs(x ="Homo-triplet", y = "Ln(Seeds per individual)",fill=NULL,color=NULL)+
  theme(legend.position="bottom")

visreg(GF_MIX_LIN_type_Plot, "hete_motif", by="type_vist",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(aes(color=fitness.data%>% filter(!is.na(type_vist))%>% 
                   dplyr::select(Plant_Simple) %>% pull()),
             size=1.5, alpha=0.5, shape=16)+
  scale_color_brewer(palette = 'Paired')+ 
  facet_wrap(vars(type_vist),nrow = 2,ncol = 3)+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Hetero-triplet", y = "Ln(Seeds per individual)",fill=NULL,color=NULL)



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

visreg(GF_MIX_LIN_type_Plot, "hete_motif", by="Plot",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(size=1.5, alpha=0.2, shape=16)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Hetero-triplet", y = "Ln(Seeds per individual)",fill=NULL,color=NULL)




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
                                        (0+scale(homo_motif)|Plot) ,
                                      #ziformula = ~1,
                                      family = nbinom1(),
                                      data = fitness_LEMA)

GF_PUPA_NB_intercept_Plot <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                                        scale(hete_motif) + 
                                        scale(StrengthIn) + scale(Ratio) +
                                        (0+scale(homo_motif)|Plot) ,
                                      #ziformula = ~1,
                                      family = nbinom1(),
                                      data = fitness_PUPA)

GF_CHFU_NB_intercept_Plot <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                        scale(hete_motif) + 
                                        scale(StrengthIn) + scale(Ratio) +
                                        (0+scale(homo_motif)|Plot) ,
                                      #ziformula = ~1,
                                      family = nbinom1(),
                                      data = fitness_CHFU)



#########################################3
# LINEAL--------RD INTERCEPT---------NO-ZI


GF_LEMA_LIN_intercept_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                       scale(hete_motif) + 
                                       scale(StrengthIn) + scale(Ratio) +
                                       (0+scale(homo_motif)|Plot) ,
                                     #ziformula = ~1,
                                     family = gaussian(),
                                     data = fitness_LEMA)

GF_PUPA_LIN_intercept_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                       scale(hete_motif) + 
                                       scale(StrengthIn) + scale(Ratio) +
                                       (0+scale(homo_motif)|Plot) ,
                                     #ziformula = ~1,
                                     family = gaussian(),
                                     data = fitness_PUPA)

GF_CHFU_LIN_intercept_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                       scale(hete_motif) + 
                                       scale(StrengthIn) + scale(Ratio) +
                                       (0+scale(homo_motif)|Plot) ,
                                     #ziformula = ~1,
                                     family = gaussian(),
                                     data = fitness_CHFU)



# NB--------RD SLOPE---------ZI


GF_LEMA_NB_intercept_Plot_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                                           scale(hete_motif) + 
                                           scale(StrengthIn) + scale(Ratio) +
                                           (0+scale(homo_motif)|Plot) ,
                                         ziformula = ~1,
                                         family = nbinom1(),
                                         data = fitness_LEMA)

GF_PUPA_NB_intercept_Plot_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) + 
                                           scale(hete_motif) + 
                                           scale(StrengthIn) + scale(Ratio) +
                                           (0+scale(homo_motif)|Plot) ,
                                         ziformula = ~1,
                                         family = nbinom1(),
                                         data = fitness_PUPA)

GF_CHFU_NB_intercept_Plot_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                           scale(hete_motif) + 
                                           scale(StrengthIn) + scale(Ratio) +
                                           (0+scale(homo_motif)|Plot) ,
                                         ziformula = ~1,
                                         family = nbinom1(),
                                         data = fitness_CHFU)


# LINEAL--------RD SLOPE---------ZI


GF_LEMA_LIN_intercept_Plot_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                          scale(hete_motif) + 
                                          scale(StrengthIn) + scale(Ratio) +
                                          (0+scale(homo_motif)|Plot) ,
                                        ziformula = ~1,
                                        family = gaussian(),
                                        data = fitness_LEMA)

GF_PUPA_LIN_intercept_Plot_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                          scale(hete_motif) + 
                                          scale(StrengthIn) + scale(Ratio) +
                                          (0+scale(homo_motif)|Plot) ,
                                        ziformula = ~1,
                                        family = gaussian(),
                                        data = fitness_PUPA)

GF_CHFU_LIN_intercept_Plot_ZI <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                                          scale(hete_motif) + 
                                          scale(StrengthIn) + scale(Ratio) +
                                          (0+scale(homo_motif)|Plot) ,
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

fitness_LEMA %>% group_by(type_vist) %>%count()
fitness_PUPA %>% group_by(type_vist) %>%count()
fitness_CHFU %>% group_by(type_vist) %>%count()

#########################################3
# LINEAL--------RD INTERCEPT---------NO-ZI

# The following models show convergence problems

GF_LEMA_LIN_type_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif)*type_vist + 
                                   scale(hete_motif)*type_vist +
                                   (0+scale(homo_motif)+scale(hete_motif)|Plot),
                                 #ziformula = ~1,
                                 family = gaussian(),
                                 data = fitness_LEMA %>% filter(!is.na(type_vist)))

# PUPA HAS NOT ENOUGH GROUPS
GF_PUPA_LIN_type_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif)*type_vist + 
                                   scale(hete_motif)*type_vist +
                                   (0+scale(homo_motif)+scale(hete_motif)|Plot),
                                 #ziformula = ~1,
                                 family = gaussian(),
                                 data = fitness_PUPA %>% filter(!is.na(type_vist)))

GF_CHFU_LIN_type_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif)*type_vist + 
                                   scale(hete_motif)*type_vist +
                                   (0+scale(homo_motif)+scale(hete_motif)|Plot),
                                 #ziformula = ~1,
                                 family = gaussian(),
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


# Get models' Residuals

res_GF_LEMA_LIN_type_Plot <- simulateResiduals(fittedModel = GF_LEMA_LIN_type_Plot, n = 500)
res_GF_CHFU_LIN_type_Plot <- simulateResiduals(fittedModel = GF_CHFU_LIN_type_Plot, n = 500)

# Checking Residuals 
plot(res_GF_LEMA_LIN_type_Plot)
plot(res_GF_CHFU_LIN_type_Plot) # 3 Quantile deviations

# AIC
AIC(GF_LEMA_LIN_type_Plot)
AIC(GF_CHFU_LIN_type_Plot)

r2(GF_LEMA_LIN_type_Plot)
r2(GF_CHFU_LIN_type_Plot)


##############################-----------
# ONE MODEL PER PLANT AND TYPE OF VISITOR

GF_LEMA_LIN_mutua_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                    scale(hete_motif) +
                                    (0+scale(homo_motif)|Plot),
                                  #ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness_LEMA %>% filter(type_vist=="mutualist"))

GF_PUPA_LIN_mutua_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                    #scale(hete_motif) +
                                    (0+scale(homo_motif)|Plot),
                                  #ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness_PUPA %>% filter(type_vist=="mutualist"))

# In the case of CHFU, hete-motifs are equal to zero

GF_CHFU_LIN_mutua_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                    #scale(hete_motif) +
                                    (0+scale(homo_motif)|Plot),
                                  #ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness_CHFU %>% filter(type_vist=="mutualist"))


GF_LEMA_LIN_antag_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                    scale(hete_motif) +
                                    (0+scale(homo_motif)|Plot),
                                  #ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness_LEMA %>% filter(type_vist=="antagonist"))

# In the case of PUPA, homo and hete motifs are equal for antagonist
# There is only one observation

GF_PUPA_LIN_antag_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                    scale(hete_motif) +
                                    (1|Plot),
                                  #ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness_PUPA %>% filter(type_vist=="antagonist"))

fitness_PUPA %>% filter(type_vist=="antagonist")


GF_CHFU_LIN_antag_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                    scale(hete_motif) +
                                    (0+scale(homo_motif)|Plot),
                                  #ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness_CHFU %>% filter(type_vist=="antagonist"))

fitness_CHFU %>% filter(type_vist=="antagonist")

GF_LEMA_LIN_other_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                    scale(hete_motif) +
                                    (0+scale(homo_motif)|Plot),
                                  #ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness_LEMA %>% filter(type_vist=="other"))

# In the case of PUPA, homo and hete motifs are equal for antagonist

fitness_PUPA %>% filter(type_vist=="other")

GF_PUPA_LIN_other_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                    scale(hete_motif) +
                                    (1|Plot),
                                  #ziformula = ~1,
                                  family = gaussian(),
                                  data = fitness_PUPA %>% filter(type_vist=="other"))

GF_CHFU_LIN_other_Plot <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) + 
                                    scale(hete_motif) +
                                    (0+scale(homo_motif)|Plot),
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
plot(res_GF_LEMA_LIN_mutua_Plot)#KS +2 Quant desv
plot(res_GF_LEMA_LIN_antag_Plot)#2 Quant desv
plot(res_GF_LEMA_LIN_other_Plot)#KS +2 Quant desv

plot(res_GF_PUPA_LIN_mutua_Plot)
plot(res_GF_PUPA_LIN_antag_Plot)
plot(res_GF_PUPA_LIN_other_Plot)

plot(res_GF_CHFU_LIN_mutua_Plot)
plot(res_GF_CHFU_LIN_antag_Plot)
plot(res_GF_CHFU_LIN_other_Plot)


# check convergence between models -----
AIC(
  GF_LEMA_LIN_mutua_Plot,
  GF_LEMA_LIN_antag_Plot,
  GF_LEMA_LIN_other_Plot)

AIC(
  GF_PUPA_LIN_mutua_Plot,
  GF_PUPA_LIN_antag_Plot,
  GF_PUPA_LIN_other_Plot)


AIC(
  GF_CHFU_LIN_mutua_Plot,
  GF_CHFU_LIN_antag_Plot,
  GF_CHFU_LIN_other_Plot)


r2(GF_LEMA_LIN_mutua_Plot)
r2(GF_LEMA_LIN_antag_Plot)
r2(GF_LEMA_LIN_other_Plot)

r2(GF_PUPA_LIN_mutua_Plot)
r2(GF_PUPA_LIN_antag_Plot)
r2(GF_PUPA_LIN_other_Plot)

r2(GF_CHFU_LIN_mutua_Plot)
r2(GF_CHFU_LIN_antag_Plot)
r2(GF_CHFU_LIN_other_Plot)
