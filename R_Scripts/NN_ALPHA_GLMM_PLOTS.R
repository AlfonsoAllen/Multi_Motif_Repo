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

fitness_data2 <- read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID_RAPE.csv")

fitness2 <- fitness_data2 %>% filter(Year==2019)

fitness2 <- fitness2 %>% dplyr::select(-Order,-Family,-Superfamily,-ID) %>% rename(ID=ID_Simple) %>%
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
competition <- read_csv("Raw_Data/competition_RAPE.csv")

competition$ME_iden <- NA
competition$ME_iden[competition$focal=="MEEL"] <- "MEEL" # we add this dummy variable to identify ME
competition$ME_iden[competition$focal=="MESU"] <- "MESU"


# Rename MEEL and MESU -> ME
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
  unique() %>% rename(Year=year,Plot=plot,Subplot=subplot,
                      Plant_Simple=focal,
                      Fruit_GF = fruit,
                      Seeds_GF = seed)

fitness <- competition_fil %>% full_join(fitness_aux,by=c("Plot","Subplot",
                                                          "Plant_Simple"))


# There are pollinator observations in ME without ME_iden. No competition data is registered
# at those sites. MESU and MEEL have zero abundances at those places.
# Assuming that ME <> MESu, we set ME_iden

fitness$ME_iden[fitness$Plant_Simple=="ME"& is.na(fitness$ME_iden)] <- "MESU"

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
fitness %>% dplyr::select(Plot,Subplot,Plant_Simple, Fruit_GF.x,Fruit_GF.y,ME_iden)%>%
  filter(Fruit_GF.x!=Fruit_GF.y) %>% group_by(Plot,Subplot,Plant_Simple,Fruit_GF.x,Fruit_GF.y,ME_iden) %>%
  count()

# We use fruit and seed data from pollinator dataset

fitness_final <- fitness %>% dplyr::select(Plot,Subplot,Plant_Simple,Seeds_GF.y,
                                    Fruit_GF.y,visits_GF,ID,homo_motif,
                                    hete_motif,ME_iden) %>%
  rename(Seeds_GF = Seeds_GF.y, Fruit_GF = Fruit_GF.y)

# Removing NAs from motifs, animals and visits
fitness_final$visits_GF[is.na(fitness_final$visits_GF)] <- 0
fitness_final$homo_motif[is.na(fitness_final$homo_motif)] <- 0
fitness_final$hete_motif[is.na(fitness_final$hete_motif)] <- 0
fitness_final$ID[is.na(fitness_final$ID)] <- "None"

#9     E3      LEMA       Seeds_GF=NA     Fruit_GF=NA    visits_GF=3 Melyridae  0   0
fitness_final <- fitness_final %>% filter(!is.na(Seeds_GF))

#############################################
# ADDING CENTRALITY MEASSURES
#############################################


centrality <- read_csv("NN_ALPHA_PageRank_results.csv") %>%
  separate(species,c("Subplot","Plant_Simple")," ")



fitness_final <- fitness_final %>% left_join(centrality, by=c("Plot","Subplot","Plant_Simple"))

# Sanity checks

fitness_final %>% filter(ID=="None",!is.na(Real_PR_Multi)) 
fitness_final %>% filter(ID!="None",is.na(Real_PR_Multi)) 


fitness_final <- fitness_final %>% filter(!is.na(Seeds_GF))

# Removing centrality NAs
fitness_final$Real_PR_Multi[is.na(fitness_final$Real_PR_Multi)] <- 1
fitness_final$Real_PR_Layer[is.na(fitness_final$Real_PR_Layer)] <- 1
fitness_final$Ratio[is.na(fitness_final$Ratio)] <- 1
fitness_final$ME_iden[is.na(fitness_final$ME_iden)] <- "0"
fitness_final[is.na(fitness_final)] <- 0

write_csv(fitness_final,"NN_ALPHA_data_models_phenol_overlap.csv")
