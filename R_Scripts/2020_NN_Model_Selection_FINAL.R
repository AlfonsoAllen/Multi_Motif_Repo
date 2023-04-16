
library(tidyverse)
source("R_Scripts/functions.R")

# Load data for models
#fitness.data.GF <- load_data_models_2020_without_agg()
fitness_orig_init <- load_data_models_2020_2() 

# Add consespecific and heterospecific probabilities

use_efficiency <- T

if(use_efficiency != T){
  
  plant_stationary_prob_results_file <- "Processed_data/2020_NN_plant_stationary_prob_results.csv"
  plant_stationary_prob_results_file_UNCOUPLED <- "Processed_data/2020_NN_plant_stationary_prob_results_UNCOUPLED.csv"

  }else{
    
    plant_stationary_prob_results_file <- "Processed_data/2020_NN_plant_stationary_prob_results_efficiency.csv"
    plant_stationary_prob_results_file_UNCOUPLED <- "Processed_data/2020_NN_plant_stationary_prob_results_efficiency_UNCOUPLED.csv"
    
}


Prob_results <- read_csv(plant_stationary_prob_results_file) %>%
  separate(name,sep=" ",c("Subplot","Plant_Simple")) %>% dplyr::select(-type,-layer)

number_plant_nodes <- fitness_orig_init %>% group_by(Plot) %>% count() %>%
  rename(total_number_plant_nodes=n)


Prob_results_UNCOUPLED <- read_csv(plant_stationary_prob_results_file_UNCOUPLED) %>%
  separate(name,sep=" ",c("Subplot","Plant_Simple")) %>% 
  dplyr::select(-type,-layer,- number_plant_nodes_with_visits)

Prob_results$Plot <- as.factor(Prob_results$Plot)
Prob_results_UNCOUPLED$Plot <- as.factor(Prob_results_UNCOUPLED$Plot)

fitness_orig <- fitness_orig_init %>% 
  left_join(Prob_results, by=c("Plot","Subplot","Plant_Simple")) %>%
  left_join(Prob_results_UNCOUPLED, by=c("Plot","Subplot","Plant_Simple")) %>%
  left_join(number_plant_nodes, by="Plot")



# We take into account that each multilayer contain a different amount of plant nodes without visits
# 
# fitness_orig$consp_prob[is.na(fitness_orig$consp_prob)] <- 0
# fitness_orig$heter_prob[is.na(fitness_orig$heter_prob)] <- 0
# 
# fitness_orig$consp_prob_UNCOUPLED[is.na(fitness_orig$consp_prob_UNCOUPLED)] <- 0
# fitness_orig$heter_prob_UNCOUPLED[is.na(fitness_orig$heter_prob_UNCOUPLED)] <- 0
# 
# fitness_orig$consp_prob[fitness_orig$DegreeIn != 0] <-
#   fitness_orig$consp_prob[fitness_orig$DegreeIn != 0]*fitness_orig$number_plant_nodes_with_visits[fitness_orig$DegreeIn != 0]/fitness_orig$total_number_plant_nodes[fitness_orig$DegreeIn != 0]
# fitness_orig$heter_prob[fitness_orig$DegreeIn != 0] <-
#   fitness_orig$heter_prob[fitness_orig$DegreeIn != 0]*fitness_orig$number_plant_nodes_with_visits[fitness_orig$DegreeIn != 0]/fitness_orig$total_number_plant_nodes[fitness_orig$DegreeIn != 0]


fitness_orig <- fitness_orig  %>% mutate(ratio=heter_prob/(consp_prob))
fitness.data <- subset(fitness_orig,Seeds_GF > 0)

fitness_orig %>% filter(type_seed_per_fruit=="Individual",type_fruit=="Individual")
#########

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
fitness.data_CHFU$Plot <- as.factor(fitness.data_CHFU$Plot)
fitness.data_CHFU$Plot %>% unique()

# Readjust factor levels
fitness_orig_CHFU$Plot <- as.factor(as.numeric(fitness_orig_CHFU$Plot))
fitness.data_CHFU$Plot <- as.factor(as.numeric(fitness.data_CHFU$Plot))

################
# PUPA

fitness_orig_PUPA <- fitness_orig %>% filter(Plant_Simple=="PUPA")
fitness.data_PUPA <- fitness.data %>% filter(Plant_Simple=="PUPA")

fitness.data_PUPA$Plot %>% unique() 

# Correlation between degree homospecific motifs

ggplot(fitness_orig %>% filter(DegreeIn>0))+
  geom_point(aes(x=homo_motif,y=DegreeIn,color=Plant_Simple),
             alpha=0.5)#+

# facet_wrap(~Plot)
ggplot(fitness_orig_CHFU %>% filter(DegreeIn>0))+
  geom_point(aes(x=homo_motif,y=DegreeIn,color=Plant_Simple),
             alpha=0.5)#+

ggplot(fitness_orig_LEMA %>% filter(DegreeIn>0))+
  geom_point(aes(x=homo_motif,y=DegreeIn,color=Plant_Simple),
             alpha=0.5)#+

ggplot(fitness_orig_PUPA %>% filter(DegreeIn>0))+
  geom_point(aes(x=homo_motif,y=DegreeIn,color=Plant_Simple),
             alpha=0.5)#+


cor(fitness_orig$homo_motif,fitness_orig$DegreeIn,method = "spearman")

cor(fitness_orig_CHFU$homo_motif,fitness_orig_CHFU$DegreeIn,method = "spearman")

cor(fitness_orig_LEMA$homo_motif,fitness_orig_LEMA$DegreeIn,method = "spearman")

cor(fitness_orig_PUPA$homo_motif,fitness_orig_PUPA$DegreeIn,method = "spearman")

cor(fitness_orig$StrengthIn,fitness_orig$consp_prob,method = "spearman")

cor(fitness_orig_CHFU$StrengthIn,fitness_orig_CHFU$consp_prob,method = "spearman")

cor(fitness_orig_LEMA$StrengthIn,fitness_orig_LEMA$consp_prob,method = "spearman")

cor(fitness_orig_PUPA$StrengthIn,fitness_orig_PUPA$consp_prob,method = "spearman")


cor.test(fitness_orig$consp_prob_UNCOUPLED[fitness_orig$DegreeIn>0],
         fitness_orig$consp_prob[fitness_orig$DegreeIn>0],method = "spearman")

cor.test(fitness_orig_CHFU$consp_prob_UNCOUPLED[fitness_orig_CHFU$DegreeIn>0],
    fitness_orig_CHFU$consp_prob[fitness_orig_CHFU$DegreeIn>0],method = "spearman")

cor.test(fitness_orig_LEMA$consp_prob_UNCOUPLED[fitness_orig_LEMA$DegreeIn>0],
    fitness_orig_LEMA$consp_prob[fitness_orig_LEMA$DegreeIn>0],method = "spearman")

cor.test(fitness_orig_PUPA$consp_prob_UNCOUPLED[fitness_orig_PUPA$DegreeIn>0],
    fitness_orig_PUPA$consp_prob[fitness_orig_PUPA$DegreeIn>0],method = "spearman")


pairwise.wilcox.test(fitness_orig$consp_prob_UNCOUPLED, fitness_orig$Plant_Simple)


cor.test(fitness_orig$consp_prob_UNCOUPLED[fitness_orig$DegreeIn>0],
         fitness_orig$visits_GF[fitness_orig$DegreeIn>0],method = "spearman")
cor.test(fitness_orig$heter_prob[fitness_orig$DegreeIn>0],
         fitness_orig$visits_GF[fitness_orig$DegreeIn>0],method = "spearman")
cor.test(fitness_orig$Seeds_GF[fitness_orig$DegreeIn>0],
         fitness_orig$visits_GF[fitness_orig$DegreeIn>0],method = "spearman")
cor.test(fitness_orig_LEMA$Seeds_GF[fitness_orig_LEMA$DegreeIn>0],
         fitness_orig_LEMA$visits_GF[fitness_orig_LEMA$DegreeIn>0],method = "spearman")
cor.test(fitness_orig_CHFU$Seeds_GF[fitness_orig_CHFU$DegreeIn>0],
         fitness_orig_CHFU$visits_GF[fitness_orig_CHFU$DegreeIn>0],method = "spearman")
cor.test(fitness_orig_PUPA$Seeds_GF[fitness_orig_PUPA$DegreeIn>0],
         fitness_orig_PUPA$visits_GF[fitness_orig_PUPA$DegreeIn>0],method = "spearman")

cor.test(fitness_orig_CHFU$consp_prob_UNCOUPLED[fitness_orig_CHFU$DegreeIn>0],
         fitness_orig_CHFU$consp_prob[fitness_orig_CHFU$DegreeIn>0],method = "spearman")

cor.test(fitness_orig_LEMA$consp_prob_UNCOUPLED[fitness_orig_LEMA$DegreeIn>0],
         fitness_orig_LEMA$consp_prob[fitness_orig_LEMA$DegreeIn>0],method = "spearman")

cor.test(fitness_orig_PUPA$consp_prob_UNCOUPLED[fitness_orig_PUPA$DegreeIn>0],
         fitness_orig_PUPA$consp_prob[fitness_orig_PUPA$DegreeIn>0],method = "spearman")


cor(fitness_orig$StrengthIn[fitness_orig$DegreeIn>0],
    fitness_orig$heter_prob[fitness_orig$DegreeIn>0],method = "spearman")

cor(fitness_orig$StrengthIn[fitness_orig$DegreeIn>0],
    fitness_orig$heter_prob[fitness_orig$DegreeIn>0],method = "spearman")


cor(fitness_orig_CHFU$StrengthIn[fitness_orig_CHFU$DegreeIn>0],
    fitness_orig_CHFU$heter_prob[fitness_orig_CHFU$DegreeIn>0],method = "spearman")

cor(fitness_orig_LEMA$StrengthIn[fitness_orig_LEMA$DegreeIn>0],
    fitness_orig_LEMA$heter_prob[fitness_orig_LEMA$DegreeIn>0],method = "spearman")

cor(fitness_orig_PUPA$StrengthIn[fitness_orig_PUPA$DegreeIn>0],
    fitness_orig_PUPA$heter_prob[fitness_orig_PUPA$DegreeIn>0],method = "spearman")

cor(fitness_orig$homo_motif[fitness_orig$DegreeIn>0],
    fitness_orig$DegreeIn[fitness_orig$DegreeIn>0],method = "spearman")

cor(fitness_orig_CHFU$homo_motif[fitness_orig_CHFU$DegreeIn>0],
    fitness_orig_CHFU$DegreeIn[fitness_orig_CHFU$DegreeIn>0],method = "spearman")

cor(fitness_orig_LEMA$homo_motif[fitness_orig_LEMA$DegreeIn>0],
    fitness_orig_LEMA$DegreeIn[fitness_orig_LEMA$DegreeIn>0],method = "spearman")

cor(fitness_orig_PUPA$homo_motif[fitness_orig_PUPA$DegreeIn>0],
    fitness_orig_PUPA$DegreeIn[fitness_orig_PUPA$DegreeIn>0],method = "spearman")




cor.test((fitness_orig_CHFU$visits_GF),
         fitness_orig_CHFU$Seeds_GF,method = "spearman")
cor.test((fitness_orig_LEMA$visits_GF),
         fitness_orig_LEMA$Seeds_GF,method = "spearman")
cor.test((fitness_orig_PUPA$visits_GF),
         fitness_orig_PUPA$Seeds_GF,method = "spearman")
cor.test((fitness_orig_CHFU$visits_GF[fitness_orig_CHFU$Seeds_GF>0]),
         log(fitness_orig_CHFU$Seeds_GF[fitness_orig_CHFU$Seeds_GF>0]),method = "pearson")
cor.test((fitness_orig_LEMA$visits_GF[fitness_orig_LEMA$Seeds_GF>0]),
         log10(fitness_orig_LEMA$Seeds_GF[fitness_orig_LEMA$Seeds_GF>0]),method = "pearson")
cor.test((fitness_orig_PUPA$visits_GF[fitness_orig_PUPA$Seeds_GF>0]),
         log10(fitness_orig_PUPA$Seeds_GF[fitness_orig_PUPA$Seeds_GF>0]),method = "pearson")
length(fitness_orig_CHFU$visits_GF)
length(fitness_orig_LEMA$visits_GF)
length(fitness_orig_PUPA$visits_GF)
###################

####################################
# NEW MODELS

###################

LEMA_NB_deg_uncoupled <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                          scale(hete_motif) +
                                         scale(consp_prob_UNCOUPLED) +scale(heter_prob),
                                        #ziformula = ~1,
                                        family = nbinom2(),
                                        data = fitness_orig_LEMA%>%ungroup() %>%
                                          filter(DegreeIn>0))

LEMA_NB_deg_uncoupled2 <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED) +scale(heter_prob),
                                 data = fitness_orig_LEMA%>%ungroup() %>%
                                   filter(DegreeIn>0))


CHFU_NB_deg_uncoupled <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                          scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED) + scale(heter_prob),
                                        #ziformula = ~1,
                                        family = nbinom2(),
                                        data = fitness_orig_CHFU%>%ungroup() %>%
                                          filter(DegreeIn>0))

CHFU_NB_deg_uncoupled2 <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                  scale(hete_motif) +
                                  scale(consp_prob_UNCOUPLED) + scale(heter_prob),
                                data = fitness_orig_CHFU%>%ungroup() %>%
                                  filter(DegreeIn>0))

PUPA_NB_deg_uncoupled <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                          scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED) + scale(heter_prob),
                                        #ziformula = ~1,
                                        family = nbinom2(),
                                        data = fitness_orig_PUPA%>%ungroup() %>%
                                          filter(DegreeIn>0))

PUPA_NB_deg_uncoupled2 <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED) + scale(heter_prob),
                                 control = glm.control(maxit = 50),
                                 data = fitness_orig_PUPA %>% ungroup() %>%
                                   filter(DegreeIn>0))


summary(LEMA_NB_deg_uncoupled)
summary(LEMA_NB_deg_uncoupled2)
summary(CHFU_NB_deg_uncoupled)
summary(CHFU_NB_deg_uncoupled2)
summary(PUPA_NB_deg_uncoupled)
summary(PUPA_NB_deg_uncoupled2)


res_LEMA_NB_deg_uncoupled <- simulateResiduals(fittedModel = LEMA_NB_deg_uncoupled, n = 1500)
res_CHFU_NB_deg_uncoupled <- simulateResiduals(fittedModel = CHFU_NB_deg_uncoupled, n = 1500)
res_PUPA_NB_deg_uncoupled <- simulateResiduals(fittedModel = PUPA_NB_deg_uncoupled, n = 1500)
res_LEMA_NB_deg_uncoupled2 <- simulateResiduals(fittedModel = LEMA_NB_deg_uncoupled2, n = 1500)
res_CHFU_NB_deg_uncoupled2 <- simulateResiduals(fittedModel = CHFU_NB_deg_uncoupled2, n = 1500)
res_PUPA_NB_deg_uncoupled2 <- simulateResiduals(fittedModel = PUPA_NB_deg_uncoupled2, n = 1500)

plot(res_LEMA_NB_deg_uncoupled)
plot(res_LEMA_NB_deg_uncoupled2)
plot(res_CHFU_NB_deg_uncoupled)
plot(res_CHFU_NB_deg_uncoupled2)
plot(res_PUPA_NB_deg_uncoupled)
plot(res_PUPA_NB_deg_uncoupled2)

performance::r2(LEMA_NB_deg_uncoupled)
performance::r2(LEMA_NB_deg_uncoupled2)
performance::r2(CHFU_NB_deg_uncoupled)
performance::r2(CHFU_NB_deg_uncoupled2)
performance::r2(PUPA_NB_deg_uncoupled)
performance::r2(PUPA_NB_deg_uncoupled2)

performance::check_collinearity(LEMA_NB_deg_uncoupled,component = "conditional") # All the ( GVIF^(1/(2*Df)) )^2 < 5 
performance::check_collinearity(CHFU_NB_deg_uncoupled,component = "conditional")
performance::check_collinearity(PUPA_NB_deg_uncoupled,component = "conditional")

performance::check_collinearity(LEMA_NB_deg_uncoupled2) # All the ( GVIF^(1/(2*Df)) )^2 < 5 
performance::check_collinearity(CHFU_NB_deg_uncoupled2)
performance::check_collinearity(PUPA_NB_deg_uncoupled2)

#####################################
dev.off()
png("New_Figures/fig6_2.png", width=1961*2, height = 1961*2*800/600, res=300*2)

jtools::plot_summs(CHFU_NB_deg_uncoupled2, LEMA_NB_deg_uncoupled2,PUPA_NB_deg_uncoupled2,inner_ci_level = .9,
                   coefs = c("Homo-triplets" = "scale(homo_motif)", 
                             "Hetero-triplets" = "scale(hete_motif)",
                             "Conspecific pollen\narrival probability\n(uncoupled layers)" = "scale(consp_prob_UNCOUPLED)", 
                             "Heterospecific pollen\narrival probability" = "scale(heter_prob)"),
                   model.names = c("C. Fuscatum", "L. Maroccanus", "P. Paludosa"),legend.title = "Plant sp.")+ 
  theme(legend.text = element_text(face = "italic"))

dev.off()

library(ggtext)
library(scales)

My_Theme = theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  plot.title = element_text(size = 16))

p1 <- jtools::effect_plot(CHFU_NB_deg_uncoupled2, pred = "homo_motif", interval = TRUE, plot.points = TRUE,
                    x.label = "Homo-triplets",
                    y.label = "Seeds",
                    colors = "gray",line.thickness = 1.1,
            jitter = 0.05, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("C. Fuscatum")))+My_Theme
p2 <- jtools::effect_plot(LEMA_NB_deg_uncoupled2, pred = "homo_motif", interval = TRUE, plot.points = TRUE,
                    x.label = "Homo-triplets",
                    y.label = "Seeds",
                    colors = "black",line.thickness = 2.2,
                    jitter = 0.05, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. Maroccanus")))+My_Theme
p3 <- jtools::effect_plot(PUPA_NB_deg_uncoupled2, pred = "homo_motif", interval = TRUE, plot.points = TRUE,
                    x.label = "Homo-triplets",
                    y.label = "Seeds",
                    colors = "black",line.thickness = 2.2,
                    jitter = 0.05, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. Paludosa")))+My_Theme



p4 <- jtools::effect_plot(CHFU_NB_deg_uncoupled2, pred = "hete_motif", interval = TRUE, plot.points = TRUE,
                    x.label = "Hetero-triplets",
                    y.label = "Seeds",
                    colors = "gray",line.thickness = 1,
                    jitter = 0.05, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
labs(title=expression(italic("C. Fuscatum")))+My_Theme
p5 <- jtools::effect_plot(LEMA_NB_deg_uncoupled2, pred = "hete_motif", interval = TRUE, plot.points = TRUE,
                    x.label = "Hetero-triplets",
                    y.label = "Seeds",
                    colors = "gray",line.thickness = 1,
                    jitter = 0.05, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. Maroccanus")))+My_Theme
p6 <- jtools::effect_plot(PUPA_NB_deg_uncoupled2, pred = "hete_motif", interval = TRUE, plot.points = TRUE,
                    x.label = "Hetero-triplets",
                    y.label = "Seeds",
                    colors = "black",line.thickness = 2.2,
                    jitter = 0.05, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. Paludosa")))+My_Theme



p7 <- jtools::effect_plot(CHFU_NB_deg_uncoupled2, pred = consp_prob_UNCOUPLED, interval = TRUE, plot.points = TRUE,
                    x.label = "Conspecific pollen\narrival probability\n(uncoupled layers)",
                    y.label = "Seeds", 
                    colors = "black",line.thickness = 2.2,
                    jitter = 0.0, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("C. Fuscatum")))+My_Theme
p8 <- jtools::effect_plot(LEMA_NB_deg_uncoupled2, pred = consp_prob_UNCOUPLED, interval = TRUE, plot.points = TRUE,
                    x.label = "Conspecific pollen\narrival probability\n(uncoupled layers)",
                    y.label = "Seeds", 
                    colors = "black",line.thickness = 2.2,
                    jitter = 0.0, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. Maroccanus")))+My_Theme
p9 <- jtools::effect_plot(PUPA_NB_deg_uncoupled2, pred = consp_prob_UNCOUPLED, interval = TRUE, plot.points = TRUE,
                    x.label = "Conspecific pollen\narrival probability\n(uncoupled layers)",
                    y.label = "Seeds", 
                    colors = "gray",line.thickness = 1,
                    jitter = 0.0, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. Paludosa")))+My_Theme



p10 <- jtools::effect_plot(CHFU_NB_deg_uncoupled2, pred = heter_prob, interval = TRUE, plot.points = TRUE,
                    x.label = "Heterospecific pollen\narrival probability",
                    y.label = "Seeds", 
                    colors = "gray",line.thickness = 1,
                    jitter = 0, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("C. Fuscatum")))+My_Theme
p11 <- jtools::effect_plot(LEMA_NB_deg_uncoupled2, pred = heter_prob, interval = TRUE, plot.points = TRUE,
                    x.label = "Heterospecific pollen\narrival probability",
                    y.label = "Seeds", 
                    colors = "black",line.thickness = 2.2,
                    jitter = 0, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. Maroccanus")))+My_Theme
p12 <- jtools::effect_plot(PUPA_NB_deg_uncoupled2, pred = heter_prob, interval = TRUE, plot.points = TRUE,
                    x.label = "Heterospecific pollen\narrival probability",
                    y.label = "Seeds", 
                    colors = "black",line.thickness = 2.2,
                    jitter = 0, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. Paludosa")))+My_Theme

library(patchwork)
png("New_Figures/fig6_2.png", width=2000*4, height = 2000*2*2, res=300*2)
(p7+p8+p9)/(p10+p11+p12)/(p1+p2+p3)/(p4+p5+p6)
dev.off()



####################################
# NEW MODELS WITH VISITS
####################################

LEMA_NB_deg_uncoupled_visits <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED) +scale(heter_prob)+
                                     scale(visits_GF),
                                 data = fitness_orig_LEMA%>%ungroup() %>%
                                   filter(DegreeIn>0))

CHFU_NB_deg_uncoupled_visits <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED) + scale(heter_prob)+
                                   scale(visits_GF),
                                 data = fitness_orig_CHFU%>%ungroup() %>%
                                   filter(DegreeIn>0))

PUPA_NB_deg_uncoupled_visits <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED) + scale(heter_prob)+
                                   scale(visits_GF),
                                 control = glm.control(maxit = 50),
                                 data = fitness_orig_PUPA %>% ungroup() %>%
                                   filter(DegreeIn>0))


summary(LEMA_NB_deg_uncoupled_visits)
summary(CHFU_NB_deg_uncoupled_visits)
summary(PUPA_NB_deg_uncoupled_visits)


res_LEMA_NB_deg_uncoupled_visits <- simulateResiduals(fittedModel = LEMA_NB_deg_uncoupled_visits, n = 1500)
res_CHFU_NB_deg_uncoupled_visits <- simulateResiduals(fittedModel = CHFU_NB_deg_uncoupled_visits, n = 1500)
res_PUPA_NB_deg_uncoupled_visits <- simulateResiduals(fittedModel = PUPA_NB_deg_uncoupled_visits, n = 1500)

plot(res_LEMA_NB_deg_uncoupled_visits)
plot(res_CHFU_NB_deg_uncoupled_visits)
plot(res_PUPA_NB_deg_uncoupled_visits)

performance::r2(LEMA_NB_deg_uncoupled_visits)
performance::r2(CHFU_NB_deg_uncoupled_visits)
performance::r2(PUPA_NB_deg_uncoupled_visits)

performance::check_collinearity(LEMA_NB_deg_uncoupled_visits) # All the ( GVIF^(1/(2*Df)) )^2 < 5 
performance::check_collinearity(CHFU_NB_deg_uncoupled_visits)
performance::check_collinearity(PUPA_NB_deg_uncoupled_visits)

#####################################

library(ggtext)
library(scales)

My_Theme = theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  plot.title = element_text(size = 16))

p1_visits <- jtools::effect_plot(CHFU_NB_deg_uncoupled_visits, pred = "homo_motif", interval = TRUE, plot.points = TRUE,
                          x.label = "Homo-triplets",
                          y.label = "Seeds", 
                          colors = "gray",line.thickness = 1,
                          jitter = 0.05, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("C. Fuscatum")))+My_Theme
p2_visits <- jtools::effect_plot(LEMA_NB_deg_uncoupled_visits, pred = "homo_motif", interval = TRUE, plot.points = TRUE,
                          x.label = "Homo-triplets",
                          y.label = "Seeds", 
                          colors = "black",line.thickness = 2.2,
                          jitter = 0.05, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. Maroccanus")))+My_Theme
p3_visits <- jtools::effect_plot(PUPA_NB_deg_uncoupled_visits, pred = "homo_motif", interval = TRUE, plot.points = TRUE,
                          x.label = "Homo-triplets",
                          y.label = "Seeds", 
                          colors = "gray",line.thickness = 1,
                          jitter = 0.05, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. Paludosa")))+My_Theme



p4_visits <- jtools::effect_plot(CHFU_NB_deg_uncoupled_visits, pred = "hete_motif", interval = TRUE, plot.points = TRUE,
                          x.label = "Hetero-triplets",
                          y.label = "Seeds", 
                          colors = "gray",line.thickness = 1,
                          jitter = 0.05, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("C. Fuscatum")))+My_Theme
p5_visits <- jtools::effect_plot(LEMA_NB_deg_uncoupled_visits, pred = "hete_motif", interval = TRUE, plot.points = TRUE,
                          x.label = "Hetero-triplets",
                          y.label = "Seeds", 
                          colors = "gray",line.thickness = 1,
                          jitter = 0.05, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. Maroccanus")))+My_Theme
p6_visits <- jtools::effect_plot(PUPA_NB_deg_uncoupled_visits, pred = "hete_motif", interval = TRUE, plot.points = TRUE,
                          x.label = "Hetero-triplets",
                          y.label = "Seeds", 
                          colors = "gray",line.thickness = 1,
                          jitter = 0.05, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. Paludosa")))+My_Theme



p7_visits <- jtools::effect_plot(CHFU_NB_deg_uncoupled_visits, pred = consp_prob_UNCOUPLED, interval = TRUE, plot.points = TRUE,
                          x.label = "Conspecific pollen\narrival probability\n(uncoupled layers)",
                          y.label = "Seeds",  
                          colors = "black",line.thickness = 2.2,
                          jitter = 0.0, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("C. Fuscatum")))+My_Theme
p8_visits <- jtools::effect_plot(LEMA_NB_deg_uncoupled_visits, pred = consp_prob_UNCOUPLED, interval = TRUE, plot.points = TRUE,
                          x.label = "Conspecific pollen\narrival probability\n(uncoupled layers)",
                          y.label = "Seeds",  
                          colors = "black",line.thickness = 2.2,
                          jitter = 0.0, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. Maroccanus")))+My_Theme
p9_visits <- jtools::effect_plot(PUPA_NB_deg_uncoupled_visits, pred = consp_prob_UNCOUPLED, interval = TRUE, plot.points = TRUE,
                          x.label = "Conspecific pollen\narrival probability\n(uncoupled layers)",
                          y.label = "Seeds",  
                          colors = "gray",line.thickness = 1,
                          jitter = 0.0, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. Paludosa")))+My_Theme



p10_visits <- jtools::effect_plot(CHFU_NB_deg_uncoupled_visits, pred = heter_prob, interval = TRUE, plot.points = TRUE,
                           x.label = "Heterospecific pollen\narrival probability",
                           y.label = "Seeds",  
                           colors = "gray",line.thickness = 1,
                           jitter = 0, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("C. Fuscatum")))+My_Theme
p11_visits <- jtools::effect_plot(LEMA_NB_deg_uncoupled_visits, pred = heter_prob, interval = TRUE, plot.points = TRUE,
                           x.label = "Heterospecific pollen\narrival probability",
                           y.label = "Seeds",  
                           colors = "black",line.thickness = 2.2,
                           jitter = 0, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. Maroccanus")))+My_Theme
p12_visits <- jtools::effect_plot(PUPA_NB_deg_uncoupled_visits, pred = heter_prob, interval = TRUE, plot.points = TRUE,
                           x.label = "Heterospecific pollen\narrival probability",
                           y.label = "Seeds", 
                           colors = "gray", line.thickness = 1,
                           jitter = 0, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. Paludosa")))+My_Theme


p13_visits <- jtools::effect_plot(CHFU_NB_deg_uncoupled_visits, pred = visits_GF, interval = TRUE, plot.points = TRUE,
                                  x.label = "Visits",
                                  y.label = "Seeds",  
                                  colors = "gray",line.thickness = 1,
                                  jitter = 0, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("C. Fuscatum")))+My_Theme
p14_visits <- jtools::effect_plot(LEMA_NB_deg_uncoupled_visits, pred = visits_GF, interval = TRUE, plot.points = TRUE,
                                  x.label = "Visits",
                                  y.label = "Seeds",  
                                  colors = "black",line.thickness = 2.2,
                                  jitter = 0, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. Maroccanus")))+My_Theme
p15_visits <- jtools::effect_plot(PUPA_NB_deg_uncoupled_visits, pred = visits_GF, interval = TRUE, plot.points = TRUE,
                                  x.label = "Visits",
                                  y.label = "Seeds",  
                                  colors = "black",line.thickness = 2.2,
                                  jitter = 0, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. Paludosa")))+My_Theme

library(patchwork)
png("New_Figures/figA16_NEW.png", width=2000*4, height = 2000*5, res=300*2)
(p7_visits+p8_visits+p9_visits)/(p10_visits+p11_visits+p12_visits)/(p1_visits+p2_visits+p3_visits)/(p4_visits+p5_visits+p6_visits)/(p13_visits+p14_visits+p15_visits)
dev.off()


####################################
# ONLY VISITS MODEL
####################################

LEMA_visits <- glm.nb(Seeds_GF ~ scale(visits_GF),
                                       data = fitness_orig_LEMA%>%ungroup() %>%
                                         filter(DegreeIn>0))

CHFU_visits <- glm.nb(Seeds_GF ~ scale(visits_GF),
                                       data = fitness_orig_CHFU%>%ungroup() %>%
                                         filter(DegreeIn>0))

PUPA_visits <- glm.nb(Seeds_GF ~ scale(visits_GF),
                                       control = glm.control(maxit = 50),
                                       data = fitness_orig_PUPA %>% ungroup() %>%
                                         filter(DegreeIn>0))


summary(LEMA_visits)
summary(CHFU_visits)
summary(PUPA_visits)


res_LEMA_visits <- simulateResiduals(fittedModel = LEMA_visits, n = 1500)
res_CHFU_visits <- simulateResiduals(fittedModel = CHFU_visits, n = 1500)
res_PUPA_visits <- simulateResiduals(fittedModel = PUPA_visits, n = 1500)

plot(res_LEMA_visits)
plot(res_CHFU_visits)
plot(res_PUPA_visits)

performance::r2(LEMA_visits)
performance::r2(CHFU_visits)
performance::r2(PUPA_visits)




####################################
# NEW MODELS

###################

LEMA_NB_deg_uncoupled_plot <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED) +scale(heter_prob) + Plot,
                                 #ziformula = ~1,
                                 family = nbinom2(),
                                 data = fitness_orig_LEMA%>%ungroup() %>%
                                   filter(DegreeIn>0))

LEMA_NB_deg_uncoupled2_plot <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED) +scale(heter_prob) + Plot,
                                 data = fitness_orig_LEMA%>%ungroup() %>%
                                   filter(DegreeIn>0))


CHFU_NB_deg_uncoupled_plot <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED) + scale(heter_prob) + Plot,
                                 #ziformula = ~1,
                                 family = nbinom2(),
                                 data = fitness_orig_CHFU%>%ungroup() %>%
                                   filter(DegreeIn>0))

CHFU_NB_deg_uncoupled2_plot <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED) + scale(heter_prob) + Plot,
                                 data = fitness_orig_CHFU%>%ungroup() %>%
                                   filter(DegreeIn>0))

PUPA_NB_deg_uncoupled_plot <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED) + scale(heter_prob) + Plot,
                                 #ziformula = ~1,
                                 family = nbinom2(),
                                 data = fitness_orig_PUPA%>%ungroup() %>%
                                   filter(DegreeIn>0))

PUPA_NB_deg_uncoupled2_plot <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED) + scale(heter_prob) + Plot,
                                 control = glm.control(maxit = 50),
                                 data = fitness_orig_PUPA %>% ungroup() %>%
                                   filter(DegreeIn>0))


summary(LEMA_NB_deg_uncoupled_plot)
summary(LEMA_NB_deg_uncoupled2_plot)
summary(CHFU_NB_deg_uncoupled_plot)
summary(CHFU_NB_deg_uncoupled2_plot)
summary(PUPA_NB_deg_uncoupled_plot)
summary(PUPA_NB_deg_uncoupled2_plot)


res_LEMA_NB_deg_uncoupled_plot <- simulateResiduals(fittedModel = LEMA_NB_deg_uncoupled_plot, n = 1500)
res_CHFU_NB_deg_uncoupled_plot <- simulateResiduals(fittedModel = CHFU_NB_deg_uncoupled_plot, n = 1500)
res_PUPA_NB_deg_uncoupled_plot <- simulateResiduals(fittedModel = PUPA_NB_deg_uncoupled_plot, n = 1500)
res_LEMA_NB_deg_uncoupled2_plot <- simulateResiduals(fittedModel = LEMA_NB_deg_uncoupled2_plot, n = 1500)
res_CHFU_NB_deg_uncoupled2_plot <- simulateResiduals(fittedModel = CHFU_NB_deg_uncoupled2_plot, n = 1500)
res_PUPA_NB_deg_uncoupled2_plot <- simulateResiduals(fittedModel = PUPA_NB_deg_uncoupled2_plot, n = 1500)

plot(res_LEMA_NB_deg_uncoupled_plot)
plot(res_LEMA_NB_deg_uncoupled2_plot)
plot(res_CHFU_NB_deg_uncoupled_plot)
plot(res_CHFU_NB_deg_uncoupled2_plot)
plot(res_PUPA_NB_deg_uncoupled_plot)
plot(res_PUPA_NB_deg_uncoupled2_plot)

performance::r2(LEMA_NB_deg_uncoupled_plot)
performance::r2(LEMA_NB_deg_uncoupled2_plot)
performance::r2(CHFU_NB_deg_uncoupled_plot)
performance::r2(CHFU_NB_deg_uncoupled2_plot)
performance::r2(PUPA_NB_deg_uncoupled_plot)
performance::r2(PUPA_NB_deg_uncoupled2_plot)

performance::check_collinearity(LEMA_NB_deg_uncoupled_plot) # All the ( GVIF^(1/(2*Df)) )^2 < 5 
performance::check_collinearity(CHFU_NB_deg_uncoupled_plot)
performance::check_collinearity(PUPA_NB_deg_uncoupled_plot)


library(scales)


ggplot(fitness_orig_LEMA %>% ungroup() %>%
         filter(DegreeIn>0), aes(x=scale(homo_motif), y= scale(consp_prob_UNCOUPLED)))+
  geom_point()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  geom_smooth(method="lm")+
  labs(title=expression(italic("L. Maroccanus")))

ggplot(fitness_orig_PUPA %>% ungroup() %>%
         filter(DegreeIn>0), aes(x=scale(homo_motif), y= scale(consp_prob_UNCOUPLED)))+
  geom_point()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  geom_smooth(method="lm")+
  labs(title=expression(italic("P. Paludosa")))


ggplot(fitness_orig_LEMA %>% ungroup() %>%
         filter(DegreeIn>0), aes(x=scale(homo_motif), y= scale(consp_prob_UNCOUPLED)))+
  geom_point()+
  geom_smooth(method = "lm", formula = y~x+I(x^2))+
  labs(title=expression(italic("L. Maroccanus")))

ggplot(fitness_orig_PUPA %>% ungroup() %>%
         filter(DegreeIn>0), aes(x=scale(homo_motif), y= scale(consp_prob_UNCOUPLED)))+
  geom_point()+
  geom_smooth(method = "lm", formula = y~x+I(x^2))+
  labs(title=expression(italic("P. Paludosa")))

ggplot(fitness_orig_LEMA %>% ungroup() %>%
         filter(DegreeIn>0), aes(x=scale(homo_motif), y= scale(heter_prob)))+
  geom_point()+
  geom_smooth(method = "lm", formula = y~x+I(x^2))+
  labs(title=expression(italic("L. Maroccanus")))

ggplot(fitness_orig_PUPA %>% ungroup() %>%
         filter(DegreeIn>0), aes(x=scale(homo_motif), y= scale(heter_prob)))+
  geom_point()+
  geom_smooth(method = "lm", formula = y~x+I(x^2))+
  labs(title=expression(italic("P. Paludosa")))


cor.test(fitness.data_CHFU$homo_motif,fitness.data_CHFU$consp_prob_UNCOUPLED, method = "spearman") #NO
cor.test(fitness.data_CHFU$homo_motif,fitness.data_CHFU$heter_prob, method = "spearman") #0.42
cor.test(fitness.data_CHFU$heter_prob,fitness.data_CHFU$consp_prob_UNCOUPLED, method = "spearman") #0.43
cor.test(fitness.data_CHFU$visits_GF,fitness.data_CHFU$consp_prob_UNCOUPLED, method = "spearman") # 0.62
cor.test(fitness.data_CHFU$visits_GF,fitness.data_CHFU$heter_prob, method = "spearman") # 0.52
cor.test(fitness.data_CHFU$visits_GF,fitness.data_CHFU$homo_motif, method = "spearman") # 0.90
cor.test(fitness.data_CHFU$visits_GF,fitness.data_CHFU$hete_motif, method = "spearman") # 0.46

cor.test(fitness.data_LEMA$homo_motif,fitness.data_LEMA$consp_prob_UNCOUPLED, method = "spearman") #0.13
cor.test(fitness.data_LEMA$homo_motif,fitness.data_LEMA$heter_prob, method = "spearman") #0.56
cor.test(fitness.data_LEMA$heter_prob,fitness.data_LEMA$consp_prob_UNCOUPLED, method = "spearman") #0.36
cor.test(fitness.data_LEMA$visits_GF,fitness.data_LEMA$consp_prob_UNCOUPLED, method = "spearman") # 0.60
cor.test(fitness.data_LEMA$visits_GF,fitness.data_LEMA$heter_prob, method = "spearman") # 0.73
cor.test(fitness.data_LEMA$visits_GF,fitness.data_LEMA$homo_motif, method = "spearman") # 0.82
cor.test(fitness.data_LEMA$visits_GF,fitness.data_LEMA$hete_motif, method = "spearman") # 0.30

cor.test(fitness.data_PUPA$homo_motif,fitness.data_PUPA$consp_prob_UNCOUPLED, method = "spearman") # No
cor.test(fitness.data_PUPA$homo_motif,fitness.data_PUPA$heter_prob, method = "spearman") # No
cor.test(fitness.data_PUPA$heter_prob,fitness.data_PUPA$consp_prob_UNCOUPLED, method = "spearman") # NO
cor.test(fitness.data_PUPA$Seeds_GF,fitness.data_PUPA$consp_prob_UNCOUPLED, method = "spearman") #0.26
cor.test(fitness.data_PUPA$Seeds_GF,fitness.data_PUPA$heter_prob, method = "spearman") # NO
cor.test(fitness.data_PUPA$visits_GF,fitness.data_PUPA$consp_prob_UNCOUPLED, method = "spearman") # 0.69
cor.test(fitness.data_PUPA$visits_GF,fitness.data_PUPA$heter_prob, method = "spearman") # NO
cor.test(fitness.data_PUPA$visits_GF,fitness.data_PUPA$homo_motif, method = "spearman") # 0.77
cor.test(fitness.data_PUPA$visits_GF,fitness.data_PUPA$hete_motif, method = "spearman") # 0.69

cor.test(fitness.data_CHFU$hete_motif,fitness.data_CHFU$consp_prob_UNCOUPLED, method = "spearman") # NO
cor.test(fitness.data_CHFU$hete_motif,fitness.data_CHFU$heter_prob, method = "spearman") # NO

cor.test(fitness.data_LEMA$hete_motif,fitness.data_LEMA$consp_prob_UNCOUPLED, method = "spearman") # NO
cor.test(fitness.data_LEMA$hete_motif,fitness.data_LEMA$heter_prob, method = "spearman") # 0.3 Síg

cor.test(fitness.data_PUPA$hete_motif,fitness.data_PUPA$consp_prob_UNCOUPLED, method = "spearman")# NO
cor.test(fitness.data_PUPA$hete_motif,fitness.data_PUPA$heter_prob, method = "spearman") #0.35 Sí

