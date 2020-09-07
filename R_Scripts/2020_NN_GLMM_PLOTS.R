
library(tidyverse)


####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################

fitness_data2 <- read_csv2("Raw_Data/raw_Pollinators_2020_1.csv")

fitness2 <- fitness_data2 %>% filter(!is.na(Plant),Plant!="0",Subplot!="OUT",Plant!="Ground")

fitness2 <- fitness2 %>% dplyr::select(-Order,-Family,-Genus,-ID) %>% rename(ID=ID_Simple) %>%
  mutate(date_raw=as.Date(paste(Day,Month,Year,sep="/"), "%d/%m/%Y"),
         Week=as.numeric(format(date_raw, "%V")))


fitness2$Line <- NA

for (i in 1:nrow(fitness2 )){
  if(fitness2$Plot[i] %in% c(1,2,3)){fitness2$Line[i] <- 1}
  else if(fitness2$Plot[i] %in% c(4,5,6)){fitness2$Line[i] <- 2}
  else{fitness2$Line[i] <- 3}
} 


#######################
# ADD SEED DATA
#######################
# Aclaración de María:
# LEMA tiene dos fechas diferentes de recogida de frutos, porque se alargó la temporada 
# Pero la media está hecha teniendo en cuenta las dos fechas juntas 

seed_data_raw <- read_csv2("Raw_Data/Fitness_2020.csv")

seed_data <- seed_data_raw %>% 
  dplyr::select(Plot,Subplot,Plant,`Seeds/Fruit`,`Fruits/Panicle`) %>% unique() %>%
  rename(Seed=`Seeds/Fruit`,Fruit=`Fruits/Panicle`) %>% group_by(Plot,Subplot,Plant) %>%
  summarize(Seed = mean(Seed,na.rm = T),
            Fruit = mean(Fruit,na.rm = T)
            )
                                                                   
fitness2_seeds <- fitness2 %>%
  left_join(seed_data,by=c("Plot","Subplot","Plant"))

fitness2_seeds %>% filter(ID=="Small bee 1_8_6") %>% dplyr::select(Plot,Subplot,Plant,Seed,Fruit)

fitness3_aux <- fitness2_seeds %>% group_by(Line,Plot,Subplot,Plant,Week,ID) %>%
  summarize(Seeds_GF = mean(Seed,na.rm = T),
            Fruit_GF = mean(Fruit,na.rm = T),
            visits_GF = sum(Visits,na.rm = T))

# Intentamos suplir los datos de semillas y visitas que faltan con los datos de competición

#####################

competition <- read_csv2("Raw_Data/competition_caracoles2020.csv") %>%
  filter(!is.na(Seeds))

competition_fil <- competition %>% 
  filter(Sp.Focal %in% c("CETE","CHFU","CHMI","LEMA","PUPA","SPRU")) %>%
  dplyr::select(Year,Plot,Subplot,Sp.Focal,Fruits,Seeds) %>% 
  mutate(Fruit_GF=1,Seeds_GF=Seeds/Fruits) %>%
  unique() %>% rename(Plant=Sp.Focal)

# Corregimos las repeticiones

competition_fil <- competition_fil %>% group_by(Year,Plot,Subplot,Plant) %>%
  summarize(Seeds_GF = mean(Seeds_GF,na.rm = T),
            Fruit_GF = mean(Fruit_GF,na.rm = T)
  )

# Sanity check
competition_fil %>% group_by(Year,Plot,Subplot,Plant) %>% count() %>% filter(n>1)

for (i in 1:nrow(fitness3_aux)){
  
  if (is.nan(fitness3_aux$Seeds_GF[i])){
    data_comp_aux <- competition_fil %>%
      filter(Plot==fitness3_aux$Plot[i],
             Subplot==fitness3_aux$Subplot[i],
             Plant==fitness3_aux$Plant[i])
    
    if(nrow(data_comp_aux)){
      fitness3_aux$Seeds_GF[i] <- data_comp_aux$Seeds_GF
      fitness3_aux$Fruit_GF[i] <- 1
    }
    
    
  }
  
}



######################


fitness3 <- fitness3_aux# %>% filter(!is.nan(Seeds_GF))

# Sanity check

fitness3 %>% dplyr::select(Plot,Subplot,Plant,ID,Week) %>% group_by(Plot,Subplot,Plant,ID,Week)%>%
  count() %>% filter(n>1)

#####################################
# Uploading motifs data
#####################################

caracoles_motif <- read_csv("Processed_data/Motifs_WEEK/2020_Caracoles_WEEK_SPECIES.csv")

caracoles_motif <- caracoles_motif %>% separate(Subplot_Plant_Label,c("Subplot",
                                                                      "Plant")," ")
#unique(caracoles_motif$ID[grep(" ",caracoles_motif$ID,ignore.case = T)])
#####################################
# Merging motifs data and fitness
#####################################

fitness_aux <- fitness3 %>% left_join(caracoles_motif, by=c("Plot","Subplot","Plant","ID","Week"))

fitness_aux %>% filter(ID=="Small bee 1_8_6")
fitness_aux %>% filter(ID=="Beetle black with red end")
fitness_aux %>% filter(ID=="Black small melyridae")

#Sanity check: visits from caracoles and fitness3 are equal
fitness_aux %>% filter(visits_GF!=Visits_tot)


##################################
# Agregating by week
##################################

fitness_aux <- fitness_aux %>% group_by(Line,Plot,Subplot,Plant,ID) %>%
  summarise(
    Seeds_GF = mean(Seeds_GF,na.rm = T),
    Fruit_GF = mean(Fruit_GF,na.rm = T),
    visits_GF = sum(visits_GF,na.rm = T),
    Visits_tot = sum(Visits_tot,na.rm = T),
    homo_motif = sum(homo_motif,na.rm = T),
    hete_motif = sum(hete_motif,na.rm = T)
  )

# sanity chek

fitness_aux %>% dplyr::select(Plot,Subplot,Plant,ID) %>% group_by(Plot,Subplot,Plant,ID)%>%
  count() %>% filter(n>1)

#####################################
# COMPETITION
#####################################

fitness_aux %>% group_by(Plant) %>% count()

# En 2020 no hay MEEL

competition <- read_csv2("Raw_Data/competition_caracoles2020.csv") %>%
  filter(!is.na(Seeds))



competition_fil <- competition %>% 
  filter(Sp.Focal %in% c("CETE","CHFU","CHMI","LEMA","PUPA","SPRU")) %>%
  dplyr::select(Year,Plot,Subplot,Sp.Focal,Fruits,Seeds) %>% 
  mutate(Fruits_unit=1,Seeds_unit=Seeds/Fruits) %>%
  unique() %>% rename(Plant=Sp.Focal,
                      Fruit_GF = Fruits_unit,
                      Seeds_GF = Seeds_unit)

# There are plants with several measures

competition_fil %>% dplyr::select(Plot,Subplot,Plant) %>% group_by(Plot,Subplot,Plant)%>%
  count() %>% filter(n>1)


competition_fil <- competition_fil %>% group_by(Year,Plot,Subplot,Plant) %>%
  summarize(Seeds_GF = mean(Seeds_GF,na.rm = T),
            Fruit_GF = mean(Fruit_GF,na.rm = T),
            Seeds = mean(Seeds,na.rm = T),
            Fruits = mean(Fruits,na.rm = T)
            )

# Sanity check

competition_fil %>% dplyr::select(Plot,Subplot,Plant) %>% group_by(Plot,Subplot,Plant)%>%
  count() %>% filter(n>1)

fitness <- fitness_aux %>% full_join(competition_fil,by=c("Plot","Subplot",
                                                          "Plant"))

# Sanity check

fitness %>% filter(Plot==1,Subplot=="F2") #Están las observaciones sin semillas

fitness %>% dplyr::select(Plot,Subplot,Plant,ID) %>% group_by(Plot,Subplot,Plant,ID)%>%
  count() %>% filter(n>1)

#No todas las plantas de competencia tienen datos de fitness que no son promedio

# x <- competition_fil %>% full_join(seed_data,by=c("Plot","Subplot",
#                                                     "Plant"))

for (i in 1:nrow(fitness)){
  
  if (is.na(fitness$Fruit_GF.y[i]) | is.nan(fitness$Fruit_GF.y[i])){
    fitness$Fruit_GF.y[i] <- fitness$Fruits[i]
    fitness$Seeds_GF.y[i] <- fitness$Seeds[i]
  }else{
    fitness$Fruit_GF.y[i] <- fitness$Fruits[i]
    fitness$Seeds_GF.y[i] <- fitness$Seeds_GF.y[i]*fitness$Fruits[i]
    }
}


# We use fruit and seed data from pollinator dataset

fitness_final <- fitness %>% dplyr::select(Plot,Subplot,Plant,Seeds_GF.y,
                                    Fruit_GF.y,visits_GF,ID,homo_motif,
                                    hete_motif) %>%
  rename(Seeds_GF = Seeds_GF.y, Fruit_GF = Fruit_GF.y)


# Removing NAs from motifs, animals and visits
fitness_final$visits_GF[is.na(fitness_final$visits_GF)] <- 0
fitness_final$homo_motif[is.na(fitness_final$homo_motif)] <- 0
fitness_final$hete_motif[is.na(fitness_final$hete_motif)] <- 0
fitness_final$ID[is.na(fitness_final$ID)] <- "None"

# Sanity check

fitness_final %>% filter(Plot==1,Subplot=="F2") 

fitness_final %>% dplyr::select(Plot,Subplot,Plant,ID) %>% group_by(Plot,Subplot,Plant,ID)%>%
  count() %>% filter(n>1)

#9     E3      LEMA       Seeds_GF=NA     Fruit_GF=NA    visits_GF=3 Melyridae  0   0

#############################################
# ADDING NN CENTRALITY MEASSURES
#############################################


centrality <- read_csv("2020_NN_PageRank_results.csv") %>%
  separate(species,c("Subplot","Plant")," ")



fitness_final <- fitness_final %>% left_join(centrality, by=c("Plot","Subplot","Plant"))

# Sanity checks

fitness_final %>% filter(ID=="None",!is.na(Real_PR_Multi)) 
fitness_final %>% filter(ID!="None",is.na(Real_PR_Multi)) 
fitness_final %>% filter(Plot==1,Subplot=="F2") 

fitness_final <- fitness_final %>% filter(!is.na(Seeds_GF))

# Removing centrality NAs
fitness_final$Real_PR_Multi[is.na(fitness_final$Real_PR_Multi)] <- 1
fitness_final$Real_PR_Layer[is.na(fitness_final$Real_PR_Layer)] <- 1
fitness_final$Ratio[is.na(fitness_final$Ratio)] <- 1
fitness_final[is.na(fitness_final)] <- 0

fitness_final <- fitness_final %>% ungroup() %>% dplyr::select(-Line)

write_csv(fitness_final,"2020_NN_data_models_phenol_overlap.csv")
