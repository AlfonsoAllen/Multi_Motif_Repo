# CANTIDAD DE SEMILLAS

# he estado viendo data/competition.csv y he comprobado el caso de LEMA, 
# y por ejemplo en 2019 hay 48 semillas por 1 fruto, cuando el año pasado 
# contamos 64 semillas en el caso de LEMA, igual de ahí viene la diferencia

# para 2019 cogí las estimaciones de 2016 porque había muy pocas especies en 2019 para
# las que había estimación de semillas viables (en concreto CHFU,LEMA,SCLA)
# las estimaciones de frutos->semillas viables de 2016 se usaron para 2016,2017 y 2018
# así que las usé también para 2019 en ese fichero


# He comprobado los datos de MEEL y MESU, y MEEL tiene siempre 1 semilla,
# y MESU varía de 1 a 2. Igual se podría hacer 
# una media
# 
# CHMI igual puede tener 0 semillas cuando la abundancia de la planta fue 0. Me explico, 
# como no tenía sentido que un polinizador visitase a una planta  inexistente, decidimos 
# poner abundancias 1 a las que tenían 0, pero el número de frutos y número de semillas 
# se seguían quedando a 0, ya que esa planta cuando se hizo el muestreo de competencia 
# no estaba.


library(tidyverse)


####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################

fitness_data2 <- read_csv("Raw_data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")

fitness2 <- fitness_data2 %>% filter(Year==2019)

fitness2 <- fitness2 %>% select(-Order,-Family,-Superfamily,-ID) %>% rename(ID=ID_Simple) %>%
  mutate(date_raw=as.Date(paste(Day,Month,Year,sep="/"), "%d/%m/%Y"),
         Week=as.numeric(format(date_raw, "%V")))


fitness2$Line <- NA

for (i in 1:nrow(fitness2 )){
  if(fitness2$Plot[i] %in% c(1,2,3)){fitness2$Line[i] <- 1}
  else if(fitness2$Plot[i] %in% c(4,5,6)){fitness2$Line[i] <- 2}
  else{fitness2$Line[i] <- 3}
} 

fitness3 <- fitness2 %>% group_by(Line,Plot,Subplot,Plant_Simple,Week,ID) %>%
  summarize(Seeds_GF = mean(Seed,na.rm = T),
            Fruit_GF = mean(Fruit,na.rm = T),
            visits_GF = sum(Visits,na.rm = T))

#####################################
# Uploading motifs data
#####################################

caracoles_motif <- read_csv("Processed_data/Motifs_WEEK/Caracoles_WEEK_SPECIES.csv")

caracoles_motif <- caracoles_motif %>% separate(Subplot_Plant_Label,c("Subplot",
                                                                      "Plant_Simple")," ")

#####################################
# Merging motifs data and fitness
#####################################

fitness_aux <- caracoles_motif %>% left_join(fitness3, by=c("Plot","Subplot","Plant_Simple","ID","Week"))

#Sanity check: visits from caracoles and fitness3 are equal
fitness_aux %>% filter(visits_GF!=Visits_tot)

#####################################
# COMPETITION
#####################################
competition <- read_csv2("Raw_Data/competition.txt")

competition$ME_iden <- NA
competition$ME_iden[competition$focal=="MEEL"] <- "MEEL" # we add this dummy variable to identify ME
competition$ME_iden[competition$focal=="MESU"] <- "MESU"
competition$focal[competition$focal=="MEEL"] <- "ME"
competition$focal[competition$focal=="MESU"] <- "ME"


#CHMI: 55 semillas planta
#MEEL (ME): 1 en compet. no hay en poll.
#MESU (ME): 2 en compet. 1 en poll.
#CHFU: 76 poll.
#LEMA: 64 poll.
#PUPA: 35 poll.

competition_fil <- competition %>% 
  filter(year==2019,focal %in% c("PUPA","LEMA","CHFU", "CHMI","ME")) %>%
  select(-neighbour,-number) %>% unique() %>% rename(Year=year,Plot=plot,Subplot=subplot,
                                                     Plant_Simple=focal,
                                                     Fruit_GF = fruit,
                                                     Seeds_GF = seed)

fitness <- competition_fil %>% full_join(fitness_aux,by=c("Plot","Subplot",
                                                          "Plant_Simple"))

#According to fitness dataframe, poll. dataset does not contain MEEL observations

for (i in 1:nrow(fitness)){
  
  if (is.na(fitness$Fruit_GF.y[i]) | is.nan(fitness$Fruit_GF.y[i])){
    fitness$Fruit_GF.y[i] <- fitness$Fruit_GF.x[i]
    if (fitness$Plant_Simple[i]=="CHMI"){fitness$Seeds_GF.y[i] <- 55*fitness$Fruit_GF.y[i]}
    else if(fitness$Plant_Simple[i]=="PUPA"){fitness$Seeds_GF.y[i] <- 35*fitness$Fruit_GF.y[i]}
    else if(fitness$Plant_Simple[i]=="LEMA"){fitness$Seeds_GF.y[i] <- 64*fitness$Fruit_GF.y[i]}
    else if(fitness$Plant_Simple[i]=="CHFU"){fitness$Seeds_GF.y[i] <- 76*fitness$Fruit_GF.y[i]}
    else{fitness$Seeds_GF.y[i] <- fitness$Seeds_GF.x[i]}
    #CHMI: 55 semillas planta
    #MEEL (ME): 1 en compet. no hay en poll.
    #MESU (ME): 2 en compet. 1 en poll.
    #CHFU: 76 poll.
    #LEMA: 64 poll.
    #PUPA: 35 poll.
    }
}

# Test number of fruits: There are 12 differences between compet. dataset and poll. dataset
fitness %>% select(Plot,Subplot,Plant_Simple, Fruit_GF.x,Fruit_GF.y)%>%
  filter(Fruit_GF.x!=Fruit_GF.y) %>% group_by(Plot,Subplot,Plant_Simple,Fruit_GF.x,Fruit_GF.y) %>%
  count()

# We use fruit and seed data from pollinator dataset

fitness_final <- fitness %>% select(Plot,Subplot,Plant_Simple,Seeds_GF.y,
                                    Fruit_GF.y,visits_GF,ID,homo_motif,
                                    hete_motif) %>%
  rename(Seeds_GF = Seeds_GF.y, Fruit_GF = Fruit_GF.y)

fitness_final$visits_GF[is.na(fitness_final$visits_GF)] <- 0
fitness_final$homo_motif[is.na(fitness_final$homo_motif)] <- 0
fitness_final$hete_motif[is.na(fitness_final$hete_motif)] <- 0
fitness_final$ID[is.na(fitness_final$ID)] <- "None"


str(fitness_final)

#9     E3      LEMA       Seeds_GF=NA     Fruit_GF=NA    visits_GF=3 Melyridae  0   0
fitness_final <- fitness_final %>% filter(!is.na(Seeds_GF))


#############################################
# ADDING CENTRALITY MEASSURES
#############################################

for (i in 1:9){
  
  file_i = paste0("C:/Users/USUARIO/Desktop/Multi_Motif_Repo/Processed_data/Muxviz_Pheno_Overlap/centrality_table_Plot",i)
  Centrality_i <- read_delim(file_i,";")
  
  # Remove data for eigenvectors
  Centrality_i <- Centrality_i %>% filter(Layer=="1-Multi") %>% mutate(Plot=i) %>%
    select(-Layer,-Node,-Eigenvector) %>% separate(Label,c("Subplot","Plant_Simple")," ")
  
  
  if(i==1){centrality <- Centrality_i}else{centrality <- bind_rows(centrality,Centrality_i)}
  
}

fitness_final <- fitness_final %>% left_join(centrality, by=c("Plot","Subplot","Plant_Simple"))

# Removing centrality NAs
fitness_final[is.na(fitness_final)] <- 0

fitness_final$Plot <- as.factor(fitness_final$Plot)
fitness_final$Subplot <- as.factor(fitness_final$Subplot)
fitness_final$ID <- as.factor(fitness_final$ID)
fitness_final$Plant_Simple <- as.factor(fitness_final$Plant_Simple)


write_csv(fitness_final,"data_models_phenol_overlap.csv")

#############################################
# EXPLORING THE DATA
##############################################

mean(fitness_final$Seeds_GF) # calculate mean: 277.9159
var(fitness_final$Seeds_GF) # calculate variance 57469.86

mean(fitness_final$Fruit_GF) # calculate mean: 19.79
var(fitness_final$Fruit_GF) # calculate variance 4251.57

# Distribution plot for seeds

ggplot(fitness_final, aes(x = Seeds_GF)) +geom_histogram(binwidth = 50)+
  geom_density(aes(y= 50 * ..count..),alpha = .2, fill = "#FF6666")

# Distribution plot for Fruit

ggplot(fitness_final, aes(x = Fruit_GF)) +geom_histogram(binwidth = 1)+
  geom_density(aes(y= 1 * ..count..),alpha = .2, fill = "#FF6666")

# Data LEMA

fitness_LEMA <- fitness_final %>% filter(Plant_Simple=="LEMA")

mean(fitness_LEMA$Seeds_GF) # calculate mean: 212.26
var(fitness_LEMA$Seeds_GF) # calculate variance 23183.75

ggplot(fitness_LEMA, aes(x = Seeds_GF)) +geom_histogram(binwidth = 10)+
  geom_density(aes(y= 10 * ..count..),alpha = .2, fill = "#FF6666")

 
fitness_CHFU <- fitness_final %>% filter(Plant_Simple=="CHFU")

mean(fitness_CHFU$Seeds_GF) # calculate mean: 340.636
var(fitness_CHFU$Seeds_GF) # calculate variance 69094.83

ggplot(fitness_CHFU, aes(x = Seeds_GF)) +geom_histogram(binwidth = 10)+
  geom_density(aes(y= 10 * ..count..),alpha = .2, fill = "#FF6666")

fitness_PUPA <- fitness_final %>% filter(Plant_Simple=="PUPA")

mean(fitness_PUPA$Seeds_GF) # calculate mean: 380.047
var(fitness_PUPA$Seeds_GF) # calculate variance 68470.24

ggplot(fitness_PUPA, aes(x = Seeds_GF)) +geom_histogram(binwidth = 10)+
  geom_density(aes(y= 10 * ..count..),alpha = .2, fill = "#FF6666")


fitness_ME <- fitness_final %>% filter(Plant_Simple=="ME")

mean(fitness_ME$Seeds_GF) # calculate mean: 302.008
var(fitness_ME$Seeds_GF) # calculate variance 88171.09

ggplot(fitness_ME, aes(x = Seeds_GF)) +geom_histogram(binwidth = 10)+
  geom_density(aes(y= 10 * ..count..),alpha = .2, fill = "#FF6666")

fitness_CHMI <- fitness_final %>% filter(Plant_Simple=="CHMI")

mean(fitness_CHMI$Seeds_GF) # calculate mean: 495
var(fitness_CHMI$Seeds_GF) # calculate variance 667315

ggplot(fitness_CHMI, aes(x = Seeds_GF)) +geom_histogram(binwidth = 10)+
  geom_density(aes(y= 10 * ..count..),alpha = .2, fill = "#FF6666")

############################
library(vcd)

# Poissoness and ng_binomialness
# If the distribution fits the data, the plot should show a straight line
# the open points show the observed count metameters; the filled points show the
# confidence interval centers, and the dashed lines show the conf_level confidence intervals for
# each point

distplot(fitness_final$Seeds_GF, type = "nbinomial")
distplot(fitness_final$Seeds_GF, type = "poisson")

distplot(fitness_final$Fruit_GF, type = "nbinomial")
distplot(fitness_final$Fruit_GF, type = "poisson")

distplot(fitness_LEMA$Seeds_GF, type = "nbinomial")
distplot(fitness_CHFU$Seeds_GF, type = "nbinomial")
distplot(fitness_PUPA$Seeds_GF, type = "nbinomial")
distplot(fitness_ME$Seeds_GF, type = "nbinomial")
distplot(fitness_CHMI$Seeds_GF, type = "nbinomial")

distplot(fitness_LEMA$Seeds_GF, type = "poisson")
distplot(fitness_CHFU$Seeds_GF, type = "poisson")
distplot(fitness_PUPA$Seeds_GF, type = "poisson")
distplot(fitness_ME$Seeds_GF, type = "poisson")
distplot(fitness_CHMI$Seeds_GF, type = "poisson")

distplot(fitness_LEMA$Fruit_GF, type = "nbinomial") #OK
distplot(fitness_CHFU$Fruit_GF, type = "nbinomial") #OK
distplot(fitness_PUPA$Fruit_GF, type = "nbinomial") #NO
distplot(fitness_ME$Fruit_GF, type = "nbinomial") #NO
distplot(fitness_CHMI$Fruit_GF, type = "nbinomial") #NO

distplot(fitness_LEMA$Fruit_GF, type = "poisson") #?
distplot(fitness_CHFU$Fruit_GF, type = "poisson") #?
distplot(fitness_PUPA$Fruit_GF, type = "poisson")
distplot(fitness_ME$Fruit_GF, type = "poisson")
distplot(fitness_CHMI$Fruit_GF, type = "poisson")

############################
# Since mean and variance are quite different, 
# we fit negative binomials to avoid overdispersion
############################

########################################
# Intercept model, Random efect: Plot
########################################

library(MASS)
library(fitdistrplus)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(bbmle) ## for AICtab
## cosmetic

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# LEMA --------
#Poisson
m1.poisson_LEMA <- glmer(Seeds_GF ~ (1|Plot), family = poisson,
            data =fitness_LEMA)
summary(m1.poisson_LEMA)

# Negative binomial
m1.nbinom_LEMA <- glmer.nb(Seeds_GF ~ (1|Plot),
                    data =fitness_LEMA)
summary(m1.nbinom_LEMA)

#Check overdispersion
overdisp_fun(m1.poisson_LEMA) #Yes, there is (as expected) 
overdisp_fun(m1.nbinom_LEMA)

# Zero inflated
m1.nbinom_LEMA2 <- glmmTMB(Seeds_GF ~ +(1|Plot),family=nbinom2(),
                           ziformula= ~ 1,
                           data =fitness_LEMA)

summary(m1.nbinom_LEMA2) #ZI term is significant

# PUPA ------
fitness_PUPA %>% group_by(Plot) %>% count()

# Poisson
m1.poisson_PUPA <- glmer(Seeds_GF ~ (1|Plot), family = poisson,
                           data =fitness_PUPA)
# Negative binomial
m1.nbinom_PUPA <- glmer.nb(Seeds_GF ~ +(1|Plot),
                           data =fitness_PUPA)# Warnings : Very few plots (random factors)
summary(m1.nbinom_PUPA)

# Zero inflation
m1.nbinom_PUPA2 <- glmmTMB(Seeds_GF ~ +(1|Plot),family=nbinom2(),
                           ziformula= ~ 1,
                           data =fitness_PUPA)
summary(m1.nbinom_PUPA2)

# CHFU--------
fitness_CHFU %>% group_by(Plot) %>% count()

# Negative binomial
m1.nbinom_CHFU <- glmer.nb(Seeds_GF ~ (1|Plot),
                           data =fitness_CHFU)
summary(m1.nbinom_CHFU)

# Negative binomial with zero inflation term
m1.nbinom_CHFU2 <- glmmTMB(Seeds_GF ~ (1|Plot),
                           ziformula = ~ 1,
                           family = nbinom2(),
                           data = fitness_CHFU)
summary(m1.nbinom_CHFU2)


# ME-----
fitness_ME %>% group_by(Plot) %>% count()

# Negative binomial
m1.nbinom_ME <- glmer.nb(Seeds_GF ~ (1|Plot),
                           data =fitness_ME)
summary(m1.nbinom_ME)

# Negative binomial with zero inflation term
m1.nbinom_ME2 <- glmmTMB(Seeds_GF ~ (1|Plot),
                         ziformula = ~ 1,
                         family = nbinom2(),
                         data =fitness_ME)
summary(m1.nbinom_ME2)

# CHMI----
fitness_CHMI %>% group_by(ID) %>% count()
fitness_CHMI %>% group_by(Plot) %>% count()

fitness_CHMI %>% group_by(Plot,ID) %>% count()

# Negative Binomial
m1.nbinom_CHMI <- glmer.nb(Seeds_GF ~ (1|Plot),
                         data =fitness_CHMI) #Convergence Warnings: Very few plots (random effects)
summary(m1.nbinom_CHMI)

# Negative Binomial 
m1.nbinom_CHMI_b <- glmer.nb(Seeds_GF ~ (1|Plot:Subplot),
                           data =fitness_CHMI) #Convergence Warnings: Very few plots (random effects)
summary(m1.nbinom_CHMI_b)

# Negative binomial with zero inflation term
m1.nbinom_CHMI2 <- glmmTMB(Seeds_GF ~ (1|Plot),
                           ziformula = ~ 1,
                           family = nbinom2(),
                           data = fitness_CHMI)
summary(m1.nbinom_CHMI2)

########################################
# NEGATIVE BINOMIAL
########################################

# In glmmTMB
# nbinom1 (also called quasi-poisson) variance = µ * phi
# where µ is the mean and phi is the over-dispersion parameter
# Negative binomial distribution: linear parameterization
# 
# nbinom2 (the default negative binomial in most packages) variance = µ(1+µ/k)
# also written µ + (µ^2)/k Negative binomial distribution: quadratic parameterization


####################
# LEMA
######################

# CORRELATIONS
# A rule of thumb for interpreting the variance inflation factor:
# 1 = not correlated.
# Between 1 and 5 = moderately correlated.
# Greater than 5 = highly correlated.

library(usdm)

vif(as.data.frame(dplyr::select(fitness_LEMA,StrengthIn,homo_motif,hete_motif,Hub,Katz)))
vif(as.data.frame(dplyr::select(fitness_LEMA,StrengthIn,PageRank,DegreeIn)))
vif(as.data.frame(dplyr::select(fitness_LEMA,DegreeIn,homo_motif,hete_motif)))
vif(as.data.frame(dplyr::select(fitness_LEMA,PageRank,homo_motif,hete_motif)))

fitness_LEMA %>% group_by(Plot) %>% count()

library(DHARMa)
#Seeds_GF ~ DegreeIn + homo_motif + ID + (1|Plot)
m2.nbinom_LEMA_ZI <- glmmTMB(Seeds_GF ~ ID + scale(DegreeIn) + scale(homo_motif) + scale(hete_motif) +(1|Plot),
                           ziformula= ~ 1,
                           family = nbinom2(),
                           data = fitness_LEMA)
summary(m2.nbinom_LEMA_ZI)
# nbinom2 returns an overdispersion parameter (usually denoted θ or k); in contrast to most other
# families, larger θ corresponds to a lower variance which is µ(1 + µ/θ).
# Model Overdispersion = 5.27
var_real <- var(fitness_LEMA$Seeds_GF)
mu_real <- mean(fitness_LEMA$Seeds_GF)
# Real overdispersion
mu_real/(-1+var_real/mu_real) # 1.96

m2.nbinom_LEMA <- glmmTMB(Seeds_GF ~ ID + scale(DegreeIn) + scale(homo_motif) + scale(hete_motif) +(1|Plot),
                             ziformula= ~ 0,
                             family = nbinom2(),
                             data = fitness_LEMA)

AIC(m2.nbinom_LEMA_ZI,m2.nbinom_LEMA) #AIC for ZI_model is smaller than the alternative model

library(DHARMa)

# get residuals
simulationOutput <- simulateResiduals(fittedModel = m2.nbinom_LEMA_ZI, n = 250)
# Checking Residuals 
testDispersion(m2.nbinom)
# Check zero inflation
testZeroInflation(simulationOutput)
plot(simulationOutput)

# Simulation outliers (data points that are outside the range of simulated values)
# are highlighted as red stars. These points should be carefully interpreted, because
# we actually don’t know “how much” these values deviate from the model expectation.

plotResiduals(simulationOutput, fitness_LEMA$homo_motif)
plotResiduals(simulationOutput, fitness_LEMA$hete_motif)
plotResiduals(simulationOutput, fitness_LEMA$ID)
plotResiduals(simulationOutput, fitness_LEMA$DegreeIn)

##########################################
# CHFU
##########################################

vif(as.data.frame(dplyr::select(fitness_CHFU,StrengthIn,homo_motif,hete_motif,Hub,Katz)))
vif(as.data.frame(dplyr::select(fitness_CHFU,DegreeIn, homo_motif)))
vif(as.data.frame(dplyr::select(fitness_CHFU,PageRank, homo_motif,hete_motif)))

fitness_CHFU %>% group_by(Plot) %>% count()

x <- fitness_CHFU %>% group_by(homo_motif,ID) %>% count()

# Seeds_GF ~ ID + homo_motif + hete_motif + (1|Plot)
m2.nbinom_CHFU_ZI <- glmmTMB(Seeds_GF ~ scale(homo_motif) + scale(hete_motif) +(1|Plot),
                     ziformula= ~ 1,
                     family= nbinom2(),
                     data = fitness_CHFU)

summary(m2.nbinom_CHFU_ZI)

# Model Overdispersion = 3.62
var_real <- var(fitness_CHFU$Seeds_GF)
mu_real <- mean(fitness_CHFU$Seeds_GF)
# Real overdispersion
mu_real/(-1+var_real/mu_real) # 1.96


m2.nbinom_CHFU <- glmmTMB(Seeds_GF ~ scale(homo_motif) + scale(hete_motif) +(1|Plot),
                             ziformula= ~ 0,
                             family= nbinom2(),
                             data = fitness_CHFU)

AIC(m2.nbinom_CHFU_ZI,m2.nbinom_CHFU) #AIC for ZI_model is smaller than the alternative model

# get residuals
simulationOutput <- simulateResiduals(fittedModel = m2.nbinom_CHFU_ZI, n = 250)
testDispersion(simulationOutput)
plot(simulationOutput)
testZeroInflation(simulationOutput)

plotResiduals(simulationOutput, fitness_CHFU$homo_motif)
plotResiduals(simulationOutput, fitness_CHFU$hete_motif)
plotResiduals(simulationOutput, fitness_CHFU$ID)
plotResiduals(simulationOutput, fitness_CHFU$DegreeIn)


##########################################
# PUPA
##########################################

# Every pupa has hete_motif = 0

vif(as.data.frame(dplyr::select(fitness_PUPA,StrengthIn,homo_motif,hete_motif,Hub,Katz)))
vif(as.data.frame(dplyr::select(fitness_PUPA,DegreeIn, homo_motif)))
vif(as.data.frame(dplyr::select(fitness_PUPA,PageRank, homo_motif)))

x <- fitness_PUPA %>% group_by(ID,Plot) %>% count()


m2.nbinom_PUPA_ZI <- glmmTMB(Seeds_GF ~ scale(DegreeIn) + scale(homo_motif) + (1|Plot),
                     ziformula= ~ 1,
                     family= nbinom2(),
                     data = fitness_PUPA)
summary(m2.nbinom_PUPA_ZI)

# Model Overdispersion = 5.33
var_real <- var(fitness_PUPA$Seeds_GF)
mu_real <- mean(fitness_PUPA$Seeds_GF)
# Real overdispersion
mu_real/(-1+var_real/mu_real) # 2.12


m2.nbinom_PUPA <- glmmTMB(Seeds_GF ~ scale(DegreeIn) + scale(homo_motif) + (1|Plot),
                             ziformula= ~ 0,
                             family= nbinom2(),
                             data = fitness_PUPA)

AIC(m2.nbinom_PUPA_ZI,m2.nbinom_PUPA) #AIC for ZI_model is smaller than the alternative model


# get residuals
simulationOutput <- simulateResiduals(fittedModel = m2.nbinom_PUPA_ZI, n = 250)
testDispersion(simulationOutput)
plot(simulationOutput)
testZeroInflation(simulationOutput)

plotResiduals(simulationOutput, fitness_PUPA$homo_motif)
plotResiduals(simulationOutput, fitness_PUPA$hete_motif) #Only Zero hete_motifs
plotResiduals(simulationOutput, fitness_PUPA$ID)
plotResiduals(simulationOutput, fitness_PUPA$DegreeIn)


##########################################
# ME
##########################################

# Every pupa has hete_motif = 0

vif(as.data.frame(dplyr::select(fitness_ME,StrengthIn,homo_motif,hete_motif,Hub,Katz)))
vif(as.data.frame(dplyr::select(fitness_ME,DegreeIn, homo_motif,hete_motif)))
vif(as.data.frame(dplyr::select(fitness_ME,PageRank, homo_motif,hete_motif)))
vif(as.data.frame(dplyr::select(fitness_ME,StrengthIn, homo_motif,hete_motif)))

fitness_ME %>% group_by(Plot) %>% count()


m2.nbinom_ME_ZI <- glmmTMB(Seeds_GF ~ scale(DegreeIn) + scale(homo_motif) + scale(hete_motif) + (1|Plot),
                     ziformula= ~ 1,
                     family= nbinom2(),
                     data = fitness_ME)
summary(m2.nbinom_ME_ZI)

# Model Overdispersion = 3.02
var_real <- var(fitness_ME$Seeds_GF)
mu_real <- mean(fitness_ME$Seeds_GF)
# Real overdispersion
mu_real/(-1+var_real/mu_real) # 1.03


m2.nbinom_ME <- glmmTMB(Seeds_GF ~ scale(DegreeIn) + scale(homo_motif) + scale(hete_motif) + (1|Plot),
                           ziformula= ~ 0,
                           family= nbinom2(),
                           data = fitness_ME)

AIC(m2.nbinom_ME_ZI,m2.nbinom_ME) #AIC for ZI_model is smaller than the alternative model

# get residuals
simulationOutput <- simulateResiduals(fittedModel = m2.nbinom_ME_ZI, n = 2500)
testDispersion(simulationOutput)
plot(simulationOutput)
testZeroInflation(simulationOutput)

plotResiduals(simulationOutput, fitness_ME$homo_motif)
plotResiduals(simulationOutput, fitness_ME$hete_motif)
plotResiduals(simulationOutput, fitness_ME$ID)
plotResiduals(simulationOutput, fitness_ME$DegreeIn)

##########################################
# CHMI
##########################################

# Every pupa has hete_motif = 0

vif(as.data.frame(dplyr::select(fitness_CHMI,StrengthIn,homo_motif,hete_motif,Hub,Katz)))
vif(as.data.frame(dplyr::select(fitness_CHMI,DegreeIn, homo_motif,hete_motif)))
vif(as.data.frame(dplyr::select(fitness_CHMI,StrengthIn, homo_motif,hete_motif)))

fitness_CHMI %>% group_by(Plot) %>% count()


m2.nbinom_CHFI_ZI <- glmmTMB(Seeds_GF ~ (DegreeIn) + (1|Plot/ID),
                     ziformula= ~ 1,
                     family= nbinom2(),
                     data = fitness_CHMI)
summary(m2.nbinom_CHFI_ZI)


m2.nbinom_CHFI <- glmmTMB(Seeds_GF ~ (DegreeIn) + (1|Plot/ID),
                          ziformula= ~ 0,
                          family= nbinom2(),
                          data = fitness_CHMI)

AIC(m2.nbinom_CHFI_ZI,m2.nbinom_CHFI)


# get residuals
simulationOutput <- simulateResiduals(fittedModel = m2.nbinom_CHFI, n = 2500)
testDispersion(simulationOutput)
plot(simulationOutput)
testZeroInflation(simulationOutput)
