
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
fitness_PUPA %>% filter(Plot=="7") %>% dplyr::select(homo_motif) %>% pull() %>%summary()
fitness_PUPA %>% filter(Plot=="8") %>% dplyr::select(homo_motif) %>% pull() %>%summary()

hist(fitness_PUPA %>% filter(Plot=="7") %>% dplyr::select(homo_motif) %>% pull() %>% scale())
hist(fitness_PUPA %>% filter(Plot=="8") %>% dplyr::select(homo_motif) %>% pull()%>% scale(),18)

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

#########################################
# INDIVIDUAL MODEL FOR EACH PLANT SPECIES
#########################################

library(nlme)
library(usdm)

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
                                       #scale(hete_motif) + 
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
                                        #scale(hete_motif) + 
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
                                          #scale(hete_motif) + 
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
                                           #scale(hete_motif) + 
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
#GLMs for PUPA CHFU
##########################################

PUPA_glm1 <- lm(log(Seeds_GF) ~ scale(homo_motif) +
                  scale(StrengthIn) + scale(Ratio),
                #ziformula = ~1,
                #family = quasipoisson,
                data = fitness_PUPA)

PUPA_glm11 <- glm(log(Seeds_GF) ~ scale(homo_motif) +
                    scale(StrengthIn) + scale(Ratio),
                  #ziformula = ~1,
                  family = gaussian(),
                  data = fitness_PUPA)


PUPA_glm2 <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                       scale(StrengthIn) + scale(Ratio),
                     #ziformula = ~1,
                     family = gaussian(),
                     data = fitness_PUPA)

CHFU_glm1 <- lm(log(Seeds_GF) ~ scale(homo_motif) +
                  scale(hete_motif) + 
                  scale(StrengthIn) + scale(Ratio),
                #ziformula = ~1,
                #family = quasipoisson,
                data = fitness_CHFU)


CHFU_glm11 <- glm(log(Seeds_GF) ~ scale(homo_motif) +
                    scale(hete_motif) + 
                    scale(StrengthIn) + scale(Ratio),
                  #ziformula = ~1,
                  family = gaussian(),
                  data = fitness_CHFU)

CHFU_glm2 <- glmmTMB(log(Seeds_GF) ~ scale(homo_motif) +
                       scale(hete_motif) + 
                       scale(StrengthIn) + scale(Ratio),
                     #ziformula = ~1,
                     family = gaussian(),
                     data = fitness_CHFU)

summary(PUPA_glm1)
summary(PUPA_glm2)


summary(CHFU_glm1)
summary(CHFU_glm2)

# glmmTMB residuals show a better behavior

res_PUPA_04 <- simulateResiduals(fittedModel = PUPA_glm1, n = 500)
res_CHFU_04 <- simulateResiduals(fittedModel = CHFU_glm1, n = 500)

res_PUPA_05 <- simulateResiduals(fittedModel = PUPA_glm2, n = 500)
res_CHFU_05 <- simulateResiduals(fittedModel = CHFU_glm2, n = 500)

res_CHFU_06 <- simulateResiduals(fittedModel = m_spamm, n=1500)
# Residuals 

testDispersion(res_PUPA_04)
testDispersion(res_PUPA_05)
testDispersion(res_CHFU_04)
testDispersion(res_CHFU_05)

testZeroInflation(res_PUPA_04)
testZeroInflation(res_PUPA_05)
testZeroInflation(res_CHFU_04)
testZeroInflation(res_CHFU_05)

testSpatialAutocorrelation(res_PUPA_04, x = fitness_PUPA$x, y = fitness_PUPA$y )
testSpatialAutocorrelation(res_CHFU_04, x = fitness_CHFU$x, y = fitness_CHFU$y )
testSpatialAutocorrelation(res_CHFU_06, x = fitness_CHFU$x, y = fitness_CHFU$y )
testSpatialAutocorrelation(res_CHFU_06)


plot(res_PUPA_04)
plot(res_PUPA_05)
plot(res_CHFU_04)
plot(res_CHFU_05)


# more specific plots PUPA
plotResiduals(res_PUPA_04, scale(fitness_PUPA$homo_motif)) #2 strong hete
#plotResiduals(res_PUPA_04, scale(fitness_PUPA$hete_motif)) 
plotResiduals(res_PUPA_04, scale(fitness_PUPA$DegreeIn))#3 strong hete
plotResiduals(res_PUPA_04, scale(fitness_PUPA$Ratio)) #2 strong hete

# more specific plots CHFU
plotResiduals(res_CHFU_04, scale(fitness_CHFU$homo_motif)) #2 strong hete
plotResiduals(res_CHFU_04, scale(fitness_CHFU$hete_motif)) 
plotResiduals(res_CHFU_04, scale(fitness_CHFU$DegreeIn))#3 strong hete
plotResiduals(res_CHFU_04, scale(fitness_CHFU$Ratio)) #2 strong hete




library(car)
Anova(GF_LEMA_LIN_slope_GF_Plot)
Anova(PUPA_glm1)
Anova(CHFU_glm1)

summary(PUPA_glm1)
summary(CHFU_glm1)

library(performance)
performance::r2(GF_LEMA_LIN_slope_GF_Plot) #0.097
performance::r2(PUPA_glm1) #0.085
performance::r2(PUPA_glm11) # 0.128
performance::r2(PUPA_glm2) #0.025
performance::r2(CHFU_glm1) #0.143
performance::r2(CHFU_glm11) #0.19
performance::r2(CHFU_glm2) #0.077

# Conditional plots

library(visreg)
visreg(GF_LEMA_LIN_slope_GF_Plot)#,trans=exp)


visreg(GF_LEMA_LIN_slope_GF_Plot,"homo_motif",xlab="Homo-triplet",ylab="Seeds per individual")

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
visreg(GF_LEMA_LIN_slope_GF_Plot,"homo_motif",xlab="Ho-triplet",ylab="log(Seeds per individual)")
mtext(side=3, "LEMA", outer=F)
visreg(PUPA_glm1,"homo_motif",xlab="Ho-triplet",ylab="log(Seeds per individual)")
mtext(side=3, "PUPA", outer=F)
visreg(CHFU_glm1,"homo_motif",xlab="Ho-triplet",ylab="log(Seeds per individual)")
mtext(side=3, "CHFU", outer=F)
dev.off()

plot.new()
par(mfrow=c(1,3))
visreg(GF_LEMA_LIN_slope_GF_Plot,"hete_motif",xlab="He-triplet",ylab="log(Seeds per individual)")
mtext(side=3, "LEMA", outer=F)
visreg(CHFU_glm1,"hete_motif",xlab="He-triplet",ylab="log(Seeds per individual)")
mtext(side=3, "CHFU", outer=F)
dev.off()

plot.new()
par(mfrow=c(1,3))
visreg(GF_LEMA_LIN_slope_GF_Plot,"DegreeIn",xlab="In-Strength",ylab="log(Seeds per individual)")
mtext(side=3, "LEMA", outer=F)
visreg(PUPA_glm1,"DegreeIn",xlab="In-Strength",ylab="log(Seeds per individual)")
mtext(side=3, "PUPA", outer=F)
visreg(CHFU_glm1,"DegreeIn",xlab="In-Strength",ylab="log(Seeds per individual)")
mtext(side=3, "CHFU", outer=F)
dev.off()

plot.new()
par(mfrow=c(1,3))
visreg(GF_LEMA_LIN_slope_GF_Plot,"Ratio",xlab="Ratio",ylab="log(Seeds per individual)")
mtext(side=3, "LEMA", outer=F)
visreg(PUPA_glm1,"Ratio",xlab="Ratio",ylab="log(Seeds per individual)")
mtext(side=3, "PUPA", outer=F)
visreg(CHFU_glm1,"Ratio",xlab="Ratio",ylab="log(Seeds per individual)")
mtext(side=3, "CHFU", outer=F)
dev.off()



v <- visreg(GF_LEMA_LIN_slope_GF_Plot, "homo_motif", by="G_F")#, trans = exp)
plot(v)
plot(v, overlay=T,trans=exp)

fitness_LEMA %>% group_by(G_F) %>% count()

library(RColorBrewer)
brewer.pal(9, "Paired")


# New facet label names for G_F variable
visitor_labs <-c(
  "Bees",
  "Beetles",
  "Flies",
  "House flies",
  "Hoverflies",
  "No visitors",
  "Small beetles",
  "Small flies",
  "Flower beetles"
)
names(visitor_labs) <- c(
  'Bees',
  'Beetles',
  'Flies',
  'House_flies',
  'Hoverflies',
  "None",
  'Small_beetles',
  'Small_flies',
  "Flower_beetles"
)

visreg(GF_LEMA_LIN_slope_GF_Plot, "homo_motif", by="G_F",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(size=1.5, alpha=0.2, shape=16)+
  facet_wrap(vars(G_F),nrow = 3,ncol = 3,labeller=labeller(G_F= visitor_labs))+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Homo-triplet", y = "log(Seeds per individual)",fill=NULL,color=NULL)




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

visreg(GF_LEMA_LIN_slope_GF_Plot, "homo_motif", by="Plot",gg = TRUE, overlay=F, partial=FALSE, rug=FALSE)+
  theme_bw() +
  #geom_point(aes(color=G_F), size=2, alpha=0.3, shape=16)+
  geom_point(size=1.5, alpha=0.2, shape=16)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Homo-triplet", y = "log(Seeds per individual)",fill=NULL,color=NULL)



p1 <- visreg(GF_LEMA_LIN_slope_GF_Plot,"homo_motif",xlab="Homo triplet",ylab=NULL,gg = TRUE, partial=FALSE, rug=FALSE)+
  theme_bw()+geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="LEMA")
p2 <- visreg(PUPA_glm1,"homo_motif",xlab="Homo triplet",ylab=NULL,gg = TRUE, partial=FALSE, rug=FALSE)+
  theme_bw()+geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="PUPA")
p3 <- visreg(CHFU_glm1,"homo_motif",xlab="Homo triplet",ylab=NULL,gg = TRUE, partial=FALSE, rug=FALSE)+
  theme_bw()+geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="CHFU")

p4 <- visreg(GF_LEMA_LIN_slope_GF_Plot,"hete_motif",xlab="Hetero triplet",ylab=NULL,gg = TRUE, partial=FALSE, rug=FALSE)+
  theme_bw()+geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="LEMA")
p5 <- visreg(CHFU_glm1,"hete_motif",xlab="Hetero triplet",ylab=NULL,gg = TRUE, partial=FALSE, rug=FALSE)+
  theme_bw()+geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="CHFU")


p6 <- visreg(GF_LEMA_LIN_slope_GF_Plot,"DegreeIn",xlab="In-Strength",ylab=NULL,gg = TRUE, partial=FALSE, rug=FALSE)+
  theme_bw()+geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="LEMA")
p7 <- visreg(PUPA_glm1,"DegreeIn",xlab="In-Strength",ylab=NULL,gg = TRUE, partial=FALSE, rug=FALSE)+
  theme_bw()+geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="PUPA")
p8 <- visreg(CHFU_glm1,"DegreeIn",xlab="In-Strength",ylab=NULL,gg = TRUE, partial=FALSE, rug=FALSE)+
  theme_bw()+geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="CHFU")

p9 <- visreg(GF_LEMA_LIN_slope_GF_Plot,"Ratio",xlab="Ratio",ylab=NULL,gg = TRUE, partial=FALSE, rug=FALSE)+
  theme_bw()+geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="LEMA")
p10 <- visreg(PUPA_glm1,"Ratio",xlab="Ratio",ylab=NULL,gg = TRUE, partial=FALSE, rug=FALSE)+
  theme_bw()+geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="PUPA")
p11 <- visreg(CHFU_glm1,"Ratio",xlab="Ratio",ylab=NULL,gg = TRUE, partial=FALSE, rug=FALSE)+
  theme_bw()+geom_point(size=1.5, alpha=0.2, shape=16)+
  labs(title ="CHFU")

library(gridExtra)
library(gtable)
library(grid)


grid.arrange(
  grobs = list(p1, p2, p3,p4, p5,p6,p7,p8,p9,p10,p11),
  widths = c(1, 1, 1),
  layout_matrix = rbind(c(1, 2, 3),
                        c(6, 7, 8),
                        c(9, 10, 11),
                        c(4, NA, 5)),
  left = textGrob("log(Seeds per individual)", rot = 90, vjust = 1)
)
#Save in 800x800


summary(GF_LEMA_LIN_slope_GF_Plot)
