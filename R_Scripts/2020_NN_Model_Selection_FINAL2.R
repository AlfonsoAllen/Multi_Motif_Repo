
library(tidyverse)
source("R_Scripts/functions.R")

# Load data for models
#fitness.data.GF <- load_data_models_2020_without_agg()
fitness_orig_init <- load_data_models_2020_3() 

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


fitness_orig <- fitness_orig  %>% mutate(ratio=heter_prob/(consp_prob_UNCOUPLED))
fitness.data <- subset(fitness_orig,Seeds_GF > 0)

fitness_orig %>% filter(type_seed_per_fruit=="Individual")
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
ggplot(fitness_orig %>% filter(DegreeIn>0),aes(x=Seeds_GF,y=visits_GF))+
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

###################################
length(fitness_orig_CHFU$visits_GF)
length(fitness_orig_LEMA$visits_GF)
length(fitness_orig_PUPA$visits_GF)
###################

####################################
# NEW MODELS
###################

ALL_NB_deg_uncoupled <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED)+scale(heter_prob) + Plant_Simple,
                                 #ziformula = ~1,
                                 family = nbinom1(),
                                 data = fitness_orig %>% ungroup() %>%
                                   filter(DegreeIn>0, Plant_Simple %in% c("LEMA","CHFU","PUPA")))


LEMA_NB_deg_uncoupled <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                          scale(hete_motif) +
                                         scale(consp_prob_UNCOUPLED)+scale(heter_prob),
                                        #ziformula = ~1,
                                        family = nbinom1(),
                                        data = fitness_orig_LEMA%>%ungroup() %>%
                                          filter(DegreeIn>0))

LEMA_NB_deg_uncoupled2 <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED)+scale(heter_prob) ,
                                 data = fitness_orig_LEMA%>%ungroup() %>%
                                   filter(DegreeIn>0))


CHFU_NB_deg_uncoupled <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                          scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED)+scale(heter_prob),
                                        #ziformula = ~1,
                                        family = nbinom1(),
                                        data = fitness_orig_CHFU%>%ungroup() %>%
                                          filter(DegreeIn>0))

CHFU_NB_deg_uncoupled2 <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                  scale(hete_motif) +
                                  scale(consp_prob_UNCOUPLED)+scale(heter_prob) ,
                                data = fitness_orig_CHFU%>%ungroup() %>%
                                  filter(DegreeIn>0))

PUPA_NB_deg_uncoupled <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                          scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED)+scale(heter_prob),
                                        #ziformula = ~1,
                                        family = nbinom1(),
                                        data = fitness_orig_PUPA%>%ungroup() %>%
                                          filter(DegreeIn>0))

PUPA_NB_deg_uncoupled2 <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED)+scale(heter_prob) ,
                                 control = glm.control(maxit = 50),
                                 data = fitness_orig_PUPA %>% ungroup() %>%
                                   filter(DegreeIn>0))

summary(ALL_NB_deg_uncoupled)
summary(LEMA_NB_deg_uncoupled)
summary(LEMA_NB_deg_uncoupled2)
summary(CHFU_NB_deg_uncoupled)
summary(CHFU_NB_deg_uncoupled2)
summary(PUPA_NB_deg_uncoupled)
summary(PUPA_NB_deg_uncoupled2)

res_ALL_NB_deg_uncoupled <- simulateResiduals(fittedModel = ALL_NB_deg_uncoupled, n = 1500)
res_LEMA_NB_deg_uncoupled <- simulateResiduals(fittedModel = LEMA_NB_deg_uncoupled, n = 1500)
res_CHFU_NB_deg_uncoupled <- simulateResiduals(fittedModel = CHFU_NB_deg_uncoupled, n = 1500)
res_PUPA_NB_deg_uncoupled <- simulateResiduals(fittedModel = PUPA_NB_deg_uncoupled, n = 1500)
res_LEMA_NB_deg_uncoupled2 <- simulateResiduals(fittedModel = LEMA_NB_deg_uncoupled2, n = 1500)
res_CHFU_NB_deg_uncoupled2 <- simulateResiduals(fittedModel = CHFU_NB_deg_uncoupled2, n = 1500)
res_PUPA_NB_deg_uncoupled2 <- simulateResiduals(fittedModel = PUPA_NB_deg_uncoupled2, n = 1500)

plot(res_ALL_NB_deg_uncoupled) #Almost OK
plot(res_LEMA_NB_deg_uncoupled) # OK
plot(res_LEMA_NB_deg_uncoupled2) # OK
plot(res_CHFU_NB_deg_uncoupled) # OK
plot(res_CHFU_NB_deg_uncoupled2) # OK
plot(res_PUPA_NB_deg_uncoupled) # OK
plot(res_PUPA_NB_deg_uncoupled2) # OK

performance::r2(ALL_NB_deg_uncoupled)
performance::r2(LEMA_NB_deg_uncoupled) # 0.52????
performance::r2(LEMA_NB_deg_uncoupled2) # Nagelkerke's R2: 0.272
performance::r2(CHFU_NB_deg_uncoupled)
performance::r2(CHFU_NB_deg_uncoupled2) # Nagelkerke's R2: 0.046
performance::r2(PUPA_NB_deg_uncoupled)
performance::r2(PUPA_NB_deg_uncoupled2) # Nagelkerke's R2: 0.100

performance::check_collinearity(LEMA_NB_deg_uncoupled) # OK 
performance::check_collinearity(CHFU_NB_deg_uncoupled) # OK
performance::check_collinearity(PUPA_NB_deg_uncoupled) # OK

performance::check_collinearity(LEMA_NB_deg_uncoupled2) # OK
performance::check_collinearity(CHFU_NB_deg_uncoupled2)# OK
performance::check_collinearity(PUPA_NB_deg_uncoupled2)# OK

#####################################
dev.off()
png("New_Figures/fig6_3.png", width=1961*2, height = 1961*2*800/600, res=300*2)

jtools::plot_summs(CHFU_NB_deg_uncoupled2, LEMA_NB_deg_uncoupled2,PUPA_NB_deg_uncoupled2,inner_ci_level = .9,
                   coefs = c("Homo-triplets" = "scale(homo_motif)", 
                             "Hetero-triplets" = "scale(hete_motif)",
                             "Conspecific pollen\narrival probability\n(uncoupled layers)" = "scale(consp_prob_UNCOUPLED)", 
                             "Heterospecific pollen\narrival probability" = "scale(heter_prob)"),
                   model.names = c("C. fuscatum", "L. maroccanus", "P. paludosa"),legend.title = "Plant sp.")+ 
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

p1 <- jtools::effect_plot(CHFU_NB_deg_uncoupled, pred = "homo_motif", interval = TRUE, plot.points = TRUE,
                    x.label = "Homo-triplets",
                    y.label = "Seeds/Flower",
                    colors = "gray",line.thickness = 1.1,
            jitter = 0.05, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("C. fuscatum")))+My_Theme
p2 <- jtools::effect_plot(LEMA_NB_deg_uncoupled, pred = "homo_motif", interval = TRUE, plot.points = TRUE,
                    x.label = "Homo-triplets",
                    y.label = "Seeds/Flower",
                    colors = "black",line.thickness = 2.2,
                    jitter = 0.05, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. maroccanus")))+My_Theme
p3 <- jtools::effect_plot(PUPA_NB_deg_uncoupled, pred = "homo_motif", interval = TRUE, plot.points = TRUE,
                    x.label = "Homo-triplets",
                    y.label = "Seeds/Flower",
                    colors = "gray",line.thickness = 1.1,
                    jitter = 0.05, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. paludosa")))+My_Theme



p4 <- jtools::effect_plot(CHFU_NB_deg_uncoupled, pred = "hete_motif", interval = TRUE, plot.points = TRUE,
                    x.label = "Hetero-triplets",
                    y.label = "Seeds/Flower",
                    colors = "gray",line.thickness = 1,
                    jitter = 0.05, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
labs(title=expression(italic("C. fuscatum")))+My_Theme
p5 <- jtools::effect_plot(LEMA_NB_deg_uncoupled, pred = "hete_motif", interval = TRUE, plot.points = TRUE,
                    x.label = "Hetero-triplets",
                    y.label = "Seeds/Flower",
                    colors = "gray",line.thickness = 1,
                    jitter = 0.05, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. maroccanus")))+My_Theme
p6 <- jtools::effect_plot(PUPA_NB_deg_uncoupled, pred = "hete_motif", interval = TRUE, plot.points = TRUE,
                    x.label = "Hetero-triplets",
                    y.label = "Seeds/Flower",
                    colors = "black",line.thickness = 2.2,
                    jitter = 0.05, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. paludosa")))+My_Theme



p7 <- jtools::effect_plot(CHFU_NB_deg_uncoupled, pred = consp_prob_UNCOUPLED, interval = TRUE, plot.points = TRUE,
                    x.label = "Conspecific pollen\narrival probability\n(uncoupled layers)",
                    y.label = "Seeds/Flower", 
                    colors = "gray",line.thickness = 1,
                    jitter = 0.0, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("C. fuscatum")))+My_Theme
p8 <- jtools::effect_plot(LEMA_NB_deg_uncoupled, pred = consp_prob_UNCOUPLED, interval = TRUE, plot.points = TRUE,
                    x.label = "Conspecific pollen\narrival probability\n(uncoupled layers)",
                    y.label = "Seeds/Flower", 
                    colors = "gray",line.thickness = 1,
                    jitter = 0.0, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. maroccanus")))+My_Theme
p9 <- jtools::effect_plot(PUPA_NB_deg_uncoupled, pred = consp_prob_UNCOUPLED, interval = TRUE, plot.points = TRUE,
                    x.label = "Conspecific pollen\narrival probability\n(uncoupled layers)",
                    y.label = "Seeds/Flower", 
                    colors = "gray",line.thickness = 1,
                    jitter = 0.0, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. paludosa")))+My_Theme



p10 <- jtools::effect_plot(CHFU_NB_deg_uncoupled, pred = heter_prob, interval = TRUE, plot.points = TRUE,
                    x.label = "Heterospecific pollen\narrival probability",
                    y.label = "Seeds/Flower", 
                    colors = "gray",line.thickness = 1,
                    jitter = 0, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("C. fuscatum")))+My_Theme
p11 <- jtools::effect_plot(LEMA_NB_deg_uncoupled, pred = heter_prob, interval = TRUE, plot.points = TRUE,
                    x.label = "Heterospecific pollen\narrival probability",
                    y.label = "Seeds/Flower", 
                    colors = "black",line.thickness = 2.2,
                    jitter = 0, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. maroccanus")))+My_Theme
p12 <- jtools::effect_plot(PUPA_NB_deg_uncoupled, pred = heter_prob, interval = TRUE, plot.points = TRUE,
                    x.label = "Heterospecific pollen\narrival probability",
                    y.label = "Seeds/Flower", 
                    colors = "gray",line.thickness = 1,
                    jitter = 0, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. paludosa")))+My_Theme

library(patchwork)
png("New_Figures/fig6_3.png", width=2000*4, height = 2000*2*2, res=300*2)
(p7+p8+p9)/(p10+p11+p12)/(p1+p2+p3)/(p4+p5+p6)
dev.off()



####################################
# NEW MODELS WITH VISITS
####################################

ALL_NB_deg_uncoupled_visits <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                         scale(hete_motif) +
                                         scale(consp_prob_UNCOUPLED)+scale(heter_prob) +
                                         scale(visits_GF),
                                       data = fitness_orig %>% ungroup() %>%
                                         filter(DegreeIn>0, Plant_Simple %in% c("LEMA","CHFU","PUPA")))


LEMA_NB_deg_uncoupled_visits <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED)+scale(heter_prob) +
                                     scale(visits_GF),
                                 data = fitness_orig_LEMA%>%ungroup() %>%
                                   filter(DegreeIn>0))

CHFU_NB_deg_uncoupled_visits <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED)+scale(heter_prob) +
                                   scale(visits_GF),
                                 data = fitness_orig_CHFU%>%ungroup() %>%
                                   filter(DegreeIn>0))

PUPA_NB_deg_uncoupled_visits <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED)+scale(heter_prob) +
                                   scale(visits_GF),
                                 control = glm.control(maxit = 50),
                                 data = fitness_orig_PUPA %>% ungroup() %>%
                                   filter(DegreeIn>0))

summary(ALL_NB_deg_uncoupled_visits)
summary(LEMA_NB_deg_uncoupled_visits)
summary(CHFU_NB_deg_uncoupled_visits)
summary(PUPA_NB_deg_uncoupled_visits)

res_ALL_NB_deg_uncoupled_visits <- simulateResiduals(fittedModel = ALL_NB_deg_uncoupled_visits, n = 1500)
res_LEMA_NB_deg_uncoupled_visits <- simulateResiduals(fittedModel = LEMA_NB_deg_uncoupled_visits, n = 1500)
res_CHFU_NB_deg_uncoupled_visits <- simulateResiduals(fittedModel = CHFU_NB_deg_uncoupled_visits, n = 1500)
res_PUPA_NB_deg_uncoupled_visits <- simulateResiduals(fittedModel = PUPA_NB_deg_uncoupled_visits, n = 1500)

plot(res_ALL_NB_deg_uncoupled_visits) # NOT OK
plot(res_LEMA_NB_deg_uncoupled_visits) # OK
plot(res_CHFU_NB_deg_uncoupled_visits) # OK 
plot(res_PUPA_NB_deg_uncoupled_visits) # Almost OK

performance::r2(ALL_NB_deg_uncoupled_visits) # Nagelkerke's R2: 0.210
performance::r2(LEMA_NB_deg_uncoupled_visits) # Nagelkerke's R2: 0.317
performance::r2(CHFU_NB_deg_uncoupled_visits) # Nagelkerke's R2: 0.046
performance::r2(PUPA_NB_deg_uncoupled_visits) #  Nagelkerke's R2: 0.229

performance::check_collinearity(ALL_NB_deg_uncoupled_visits) #OK
performance::check_collinearity(LEMA_NB_deg_uncoupled_visits) #OK # All the ( GVIF^(1/(2*Df)) )^2 < 5 
performance::check_collinearity(CHFU_NB_deg_uncoupled_visits) #OK
performance::check_collinearity(PUPA_NB_deg_uncoupled_visits) #OK

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
                          y.label = "Seeds/Flower", 
                          colors = "gray",line.thickness = 1,
                          jitter = 0.05, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("C. fuscatum")))+My_Theme
p2_visits <- jtools::effect_plot(LEMA_NB_deg_uncoupled_visits, pred = "homo_motif", interval = TRUE, plot.points = TRUE,
                          x.label = "Homo-triplets",
                          y.label = "Seeds/Flower", 
                          colors = "black",line.thickness = 2.2,
                          jitter = 0.05, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. maroccanus")))+My_Theme
p3_visits <- jtools::effect_plot(PUPA_NB_deg_uncoupled_visits, pred = "homo_motif", interval = TRUE, plot.points = TRUE,
                          x.label = "Homo-triplets",
                          y.label = "Seeds/Flower", 
                          colors = "gray",line.thickness = 1,
                          jitter = 0.05, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. paludosa")))+My_Theme



p4_visits <- jtools::effect_plot(CHFU_NB_deg_uncoupled_visits, pred = "hete_motif", interval = TRUE, plot.points = TRUE,
                          x.label = "Hetero-triplets",
                          y.label = "Seeds/Flower", 
                          colors = "gray",line.thickness = 1,
                          jitter = 0.05, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("C. fuscatum")))+My_Theme
p5_visits <- jtools::effect_plot(LEMA_NB_deg_uncoupled_visits, pred = "hete_motif", interval = TRUE, plot.points = TRUE,
                          x.label = "Hetero-triplets",
                          y.label = "Seeds/Flower", 
                          colors = "gray",line.thickness = 1,
                          jitter = 0.05, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. maroccanus")))+My_Theme
p6_visits <- jtools::effect_plot(PUPA_NB_deg_uncoupled_visits, pred = "hete_motif", interval = TRUE, plot.points = TRUE,
                          x.label = "Hetero-triplets",
                          y.label = "Seeds/Flower", 
                          colors = "gray",line.thickness = 1,
                          jitter = 0.05, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. paludosa")))+My_Theme



p7_visits <- jtools::effect_plot(CHFU_NB_deg_uncoupled_visits, pred = consp_prob_UNCOUPLED, interval = TRUE, plot.points = TRUE,
                          x.label = "Conspecific pollen\narrival probability\n(uncoupled layers)",
                          y.label = "Seeds/Flower",  
                          colors = "gray",line.thickness = 1,
                          jitter = 0.0, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("C. fuscatum")))+My_Theme
p8_visits <- jtools::effect_plot(LEMA_NB_deg_uncoupled_visits, pred = consp_prob_UNCOUPLED, interval = TRUE, plot.points = TRUE,
                          x.label = "Conspecific pollen\narrival probability\n(uncoupled layers)",
                          y.label = "Seeds/Flower",  
                          colors = "gray",line.thickness = 1,
                          jitter = 0.0, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. maroccanus")))+My_Theme
p9_visits <- jtools::effect_plot(PUPA_NB_deg_uncoupled_visits, pred = consp_prob_UNCOUPLED, interval = TRUE, plot.points = TRUE,
                          x.label = "Conspecific pollen\narrival probability\n(uncoupled layers)",
                          y.label = "Seeds/Flower",  
                          colors = "gray",line.thickness = 1,
                          jitter = 0.0, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. paludosa")))+My_Theme



p10_visits <- jtools::effect_plot(CHFU_NB_deg_uncoupled_visits, pred = heter_prob, interval = TRUE, plot.points = TRUE,
                           x.label = "Heterospecific pollen\narrival probability",
                           y.label = "Seeds/Flower",  
                           colors = "gray",line.thickness = 1,
                           jitter = 0, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("C. fuscatum")))+My_Theme
p11_visits <- jtools::effect_plot(LEMA_NB_deg_uncoupled_visits, pred = heter_prob, interval = TRUE, plot.points = TRUE,
                           x.label = "Heterospecific pollen\narrival probability",
                           y.label = "Seeds/Flower",  
                           colors = "black",line.thickness = 2.2,
                           jitter = 0, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. maroccanus")))+My_Theme
p12_visits <- jtools::effect_plot(PUPA_NB_deg_uncoupled_visits, pred = heter_prob, interval = TRUE, plot.points = TRUE,
                           x.label = "Heterospecific pollen\narrival probability",
                           y.label = "Seeds/Flower", 
                           colors = "gray", line.thickness = 1,
                           jitter = 0, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. paludosa")))+My_Theme


p13_visits <- jtools::effect_plot(CHFU_NB_deg_uncoupled_visits, pred = visits_GF, interval = TRUE, plot.points = TRUE,
                                  x.label = "Visits",
                                  y.label = "Seeds/Flower",  
                                  colors = "gray",line.thickness = 1,
                                  jitter = 0, data = fitness_orig_CHFU %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("C. fuscatum")))+My_Theme
p14_visits <- jtools::effect_plot(LEMA_NB_deg_uncoupled_visits, pred = visits_GF, interval = TRUE, plot.points = TRUE,
                                  x.label = "Visits",
                                  y.label = "Seeds/Flower",  
                                  colors = "black",line.thickness = 2.2,
                                  jitter = 0, data = fitness_orig_LEMA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("L. maroccanus")))+My_Theme
p15_visits <- jtools::effect_plot(PUPA_NB_deg_uncoupled_visits, pred = visits_GF, interval = TRUE, plot.points = TRUE,
                                  x.label = "Visits",
                                  y.label = "Seeds/Flower",  
                                  colors = "black",line.thickness = 2.2,
                                  jitter = 0, data = fitness_orig_PUPA %>% ungroup() %>% filter(DegreeIn > 0))+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)))+
  labs(title=expression(italic("P. paludosa")))+My_Theme

library(patchwork)
png("New_Figures/figA16_NEW2.png", width=2000*4, height = 2000*5, res=300*2)
(p7_visits+p8_visits+p9_visits)/(p10_visits+p11_visits+p12_visits)/(p1_visits+p2_visits+p3_visits)/(p4_visits+p5_visits+p6_visits)/(p13_visits+p14_visits+p15_visits)
dev.off()


####################################
# ONLY VISITS MODEL
####################################
ALL_visits <- glm.nb(Seeds_GF ~ scale(visits_GF),
                      data = fitness_orig %>%ungroup() %>%
                        filter(DegreeIn>0, Plant_Simple %in% c("LEMA","CHFU","PUPA")))

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

summary(ALL_visits)
summary(LEMA_visits)
summary(CHFU_visits)
summary(PUPA_visits)

res_ALL_visits <- simulateResiduals(fittedModel = ALL_visits, n = 1500)
res_LEMA_visits <- simulateResiduals(fittedModel = LEMA_visits, n = 1500)
res_CHFU_visits <- simulateResiduals(fittedModel = CHFU_visits, n = 1500)
res_PUPA_visits <- simulateResiduals(fittedModel = PUPA_visits, n = 1500)

plot(res_ALL_visits) # NOT OK
plot(res_LEMA_visits) # OK
plot(res_CHFU_visits) # OK
plot(res_PUPA_visits) # OK

performance::r2(ALL_visits) # Nagelkerke's R2: 0.126
performance::r2(LEMA_visits) # Nagelkerke's R2: 0.060
performance::r2(CHFU_visits) # Nagelkerke's R2: 0.015
performance::r2(PUPA_visits) # Nagelkerke's R2: 0.142




####################################
# MODELS THAT INCLUDE PLOT AS VARIABLE OR RANDOM FACTOR

###################

ALL_NB_deg_uncoupled_plot_RF <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                        scale(hete_motif) +
                                        scale(consp_prob_UNCOUPLED)+scale(heter_prob)  + (1|Plot),
                                      #ziformula = ~1,
                                      family = nbinom1(),
                                      data = fitness_orig%>%ungroup() %>%
                                        filter(DegreeIn>0, Plant_Simple %in% c("LEMA","CHFU","PUPA")))

ALL_NB_deg_uncoupled2_plot <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                        scale(hete_motif) +
                                        scale(consp_prob_UNCOUPLED)+scale(heter_prob)  + Plot,
                                      data = fitness_orig%>%ungroup() %>%
                                        filter(DegreeIn>0, Plant_Simple %in% c("LEMA","CHFU","PUPA")))

LEMA_NB_deg_uncoupled_plot_RF <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED)+scale(heter_prob)  + (1|Plot),
                                 #ziformula = ~1,
                                 family = nbinom1(),
                                 data = fitness_orig_LEMA%>%ungroup() %>%
                                   filter(DegreeIn>0))

LEMA_NB_deg_uncoupled2_plot <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED) + scale(heter_prob)  + Plot,
                                 data = fitness_orig_LEMA%>%ungroup() %>%
                                   filter(DegreeIn>0))


CHFU_NB_deg_uncoupled_plot_RF <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED)+scale(heter_prob)  + (1|Plot),
                                 #ziformula = ~1,
                                 family = nbinom1(),
                                 data = fitness_orig_CHFU%>%ungroup() %>%
                                   filter(DegreeIn>0))

CHFU_NB_deg_uncoupled2_plot <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED)+scale(heter_prob)  + Plot,
                                 data = fitness_orig_CHFU%>%ungroup() %>%
                                   filter(DegreeIn>0))

PUPA_NB_deg_uncoupled_plot_RF <- glmmTMB(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED)+scale(heter_prob)  + (1|Plot),
                                 #ziformula = ~1,
                                 family = nbinom1(),
                                 data = fitness_orig_PUPA%>%ungroup() %>%
                                   filter(DegreeIn>0))

PUPA_NB_deg_uncoupled2_plot <- glm.nb(Seeds_GF ~ scale(homo_motif) +
                                   scale(hete_motif) +
                                   scale(consp_prob_UNCOUPLED)+scale(heter_prob)  + Plot,
                                 control = glm.control(maxit = 50),
                                 data = fitness_orig_PUPA %>% ungroup() %>%
                                   filter(DegreeIn>0))

summary(ALL_NB_deg_uncoupled_plot_RF)
summary(ALL_NB_deg_uncoupled2_plot)
summary(LEMA_NB_deg_uncoupled_plot)
summary(LEMA_NB_deg_uncoupled2_plot)
summary(CHFU_NB_deg_uncoupled_plot_RF)
summary(CHFU_NB_deg_uncoupled2_plot)
summary(PUPA_NB_deg_uncoupled_plot_RF)
summary(PUPA_NB_deg_uncoupled2_plot)


res_ALL_NB_deg_uncoupled_plot_RF <- simulateResiduals(fittedModel = ALL_NB_deg_uncoupled_plot_RF, n = 1500)
res_LEMA_NB_deg_uncoupled_plot_RF <- simulateResiduals(fittedModel = LEMA_NB_deg_uncoupled_plot_RF, n = 1500)
res_CHFU_NB_deg_uncoupled_plot_RF <- simulateResiduals(fittedModel = CHFU_NB_deg_uncoupled_plot_RF, n = 1500)
res_PUPA_NB_deg_uncoupled_plot_RF <- simulateResiduals(fittedModel = PUPA_NB_deg_uncoupled_plot_RF, n = 1500)
res_ALL_NB_deg_uncoupled2_plot <- simulateResiduals(fittedModel = ALL_NB_deg_uncoupled2_plot, n = 1500)
res_LEMA_NB_deg_uncoupled2_plot <- simulateResiduals(fittedModel = LEMA_NB_deg_uncoupled2_plot, n = 1500)
res_CHFU_NB_deg_uncoupled2_plot <- simulateResiduals(fittedModel = CHFU_NB_deg_uncoupled2_plot, n = 1500)
res_PUPA_NB_deg_uncoupled2_plot <- simulateResiduals(fittedModel = PUPA_NB_deg_uncoupled2_plot, n = 1500)

plot(res_ALL_NB_deg_uncoupled_plot_RF) # Almost OK
plot(res_ALL_NB_deg_uncoupled2_plot) # Not OK
plot(res_LEMA_NB_deg_uncoupled_plot_RF) # Almost OK
plot(res_LEMA_NB_deg_uncoupled2_plot) # OK
plot(res_CHFU_NB_deg_uncoupled_plot_RF) # OK
plot(res_CHFU_NB_deg_uncoupled2_plot) # OK
plot(res_PUPA_NB_deg_uncoupled_plot_RF) # OK
plot(res_PUPA_NB_deg_uncoupled2_plot) # Almost OK

performance::r2(ALL_NB_deg_uncoupled_plot_RF)
performance::r2(ALL_NB_deg_uncoupled2_plot)
performance::r2(LEMA_NB_deg_uncoupled_plot_RF)
performance::r2(LEMA_NB_deg_uncoupled2_plot)
performance::r2(CHFU_NB_deg_uncoupled_plot_RF)
performance::r2(CHFU_NB_deg_uncoupled2_plot)
performance::r2(PUPA_NB_deg_uncoupled_plot_RF)
performance::r2(PUPA_NB_deg_uncoupled2_plot)

performance::check_collinearity(ALL_NB_deg_uncoupled_plot_RF) # OK
performance::check_collinearity(LEMA_NB_deg_uncoupled_plot_RF) # OK # OK # All the ( GVIF^(1/(2*Df)) )^2 < 5 
performance::check_collinearity(CHFU_NB_deg_uncoupled_plot_RF) # OK
performance::check_collinearity(PUPA_NB_deg_uncoupled_plot_RF) # OK