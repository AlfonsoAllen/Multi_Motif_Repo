
library(tidyverse)


####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################

fitness_data2 <- read_csv2("Raw_Data/final_Pollinators_2020.csv")
# Remove points from ID names
fitness_data2$ID <- sub("\\.", "", fitness_data2$ID)
fitness_data2$ID_Simple <- sub("\\.", "", fitness_data2$ID_Simple)

# filter tabanidae
fitness_data2 <- fitness_data2 %>% filter(ID != "Tabanidae")


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

# Plant species with visits
fitness2$Plant %>% unique()

fitness2 %>% filter(Plot==6) %>% select(Subplot,ID) %>% unique()
#######################
# ADD SEED DATA
#######################
# Aclaración de María:
# LEMA tiene dos fechas diferentes de recogida de frutos, porque se alargó la temporada 
# Pero la media está hecha teniendo en cuenta las dos fechas juntas 

seed_data_raw <- read_csv2("Raw_Data/Fitness_2020.csv")

seed_data_raw %>% group_by(Plot, Subplot, Plant) %>% count() %>% filter(n==1)
seed_data_raw %>% group_by(Plot, Subplot, Plant) %>% count() %>% filter(n==2)
seed_data_raw %>% group_by(Plot, Subplot, Plant) %>% count() %>% filter(n==3)
seed_data_raw %>% group_by(Plot, Subplot, Plant) %>% count() %>% filter(n>3)
seed_data_raw %>% group_by(Plot, Subplot, Plant) %>% count() %>% ungroup() %>% 
  select(n) %>% pull() %>% mean()
seed_data_raw %>% group_by(Plot, Subplot, Plant) %>% count() %>% ungroup() %>% 
  select(n) %>% pull() %>% sd()

seed_data_raw %>% group_by(Plot, Subplot) %>% count() %>% filter(n==1)
seed_data_raw %>% group_by(Plot, Subplot) %>% count() %>% filter(n==2)
seed_data_raw %>% group_by(Plot, Subplot) %>% count() %>% filter(n==3)
seed_data_raw %>% group_by(Plot, Subplot) %>% count() %>% filter(n==4)
seed_data_raw %>% group_by(Plot, Subplot) %>% count() %>% filter(n==5)
seed_data_raw %>% group_by(Plot, Subplot) %>% count() %>% filter(n>5)
seed_data_raw %>% group_by(Plot, Subplot) %>% count() %>% ungroup() %>% 
  select(n) %>% pull() %>% mean()
seed_data_raw %>% group_by(Plot, Subplot) %>% count() %>% ungroup() %>% 
  select(n) %>% pull() %>% sd()

seed_data <- seed_data_raw %>% 
  dplyr::select(Plot,Subplot,Plant,`Seeds/Fruit`,`Fruits/Panicle`) %>% unique() %>%
  rename(Seed=`Seeds/Fruit`,Fruit=`Fruits/Panicle`) %>% group_by(Plot,Subplot,Plant) %>%
  summarize(Seed = mean(Seed,na.rm = T),
            Fruit = mean(Fruit,na.rm = T)
            ) %>% arrange(Plot,Subplot)
                         

fitness2_seeds_aux <- fitness2 %>%
  left_join(seed_data,by=c("Plot","Subplot","Plant"))

# There are plants with Plot and subplot set to NA. We add to those their mean values

fitness2_seeds_with_ind_data <- fitness2_seeds_aux %>% filter(!is.na(Seed))
fitness2_seeds_without_ind_data <- fitness2_seeds_aux %>% filter(is.na(Seed)) %>%
  dplyr::select(-Seed,-Fruit)

mean_seed <- seed_data_raw %>% 
  dplyr::select(Plant,`Mean Seeds/Fruit`,`Mean Fruits/panicle`) %>% unique() %>%
  rename(Seed=`Mean Seeds/Fruit`,Fruit=`Mean Fruits/panicle`)

mean_seed_plot <- seed_data_raw %>% 
  dplyr::select(Plot,Subplot,Plant,`Seeds/Fruit`,`Fruits/Panicle`) %>% unique() %>%
  rename(Seed=`Seeds/Fruit`,Fruit=`Fruits/Panicle`) %>% group_by(Plot,Plant) %>%
  summarize(Seed = mean(Seed,na.rm = T),
            Fruit = mean(Fruit,na.rm = T)
  ) %>% arrange(Plot) # %>% filter(Plant=="LEMA")


# fitness2_seeds_without_ind_data <- fitness2_seeds_without_ind_data %>%
#   left_join(mean_seed,by="Plant")


fitness2_seeds_without_ind_data <- fitness2_seeds_without_ind_data %>%
  left_join(mean_seed_plot,by=c("Plot","Plant"))

sum(is.na(fitness2_seeds_without_ind_data$Seed))
sum(fitness2_seeds_without_ind_data$Seed==0)

# Merge data with individual values and mean values, respectively

fitness2_seeds <- bind_rows(fitness2_seeds_with_ind_data,fitness2_seeds_without_ind_data)

# Sanity check
fitness2_seeds %>% filter(is.na(Seed))


fitness3_aux <- fitness2_seeds %>% group_by(Line,Plot,Subplot,Plant,Week,ID) %>%
  summarize(Seeds_GF = mean(Seed,na.rm = T),
            Fruit_GF = mean(Fruit,na.rm = T),
            visits_GF = sum(Visits,na.rm = T))

fitness3 <- fitness3_aux %>% filter(!is.nan(Seeds_GF))

fitness3$Plant %>% unique()

# Sanity check

fitness3 %>% dplyr::select(Plot,Subplot,Plant,ID,Week) %>% group_by(Plot,Subplot,Plant,ID,Week)%>%
  count() %>% filter(n>1)

#####################################
# Uploading motifs data
#####################################

caracoles_motif <- read_csv("Processed_data/Motifs_WEEK/2020_Caracoles_WEEK_SPECIES.csv")

caracoles_motif <- caracoles_motif %>% separate(Subplot_Plant_Label,c("Subplot",
                                                                      "Plant")," ")

#####################################
# Merging motifs data and fitness
#####################################

fitness_aux <- fitness3 %>% left_join(caracoles_motif, by=c("Plot","Subplot","Plant","ID","Week"))

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

# sanity check

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
  filter(Sp.Focal %in% c("BEMA","CETE","CHFU","CHMI","LEMA","MESU","PUPA",
                         "SCLA","SOAS","SPRU")) %>%
  dplyr::select(Year,Plot,Subplot,Sp.Focal,Fruits,Seeds) %>% 
  mutate(Fruits_unit=1,Seeds_unit=Seeds/Fruits) %>%
  unique() %>% rename(Plant=Sp.Focal,
                      Fruit_GF = Fruits_unit,
                      Seeds_GF = Seeds_unit)

# In this database they are using mean numbers of seeds per fruit
competition_fil %>% select(Plant,Seeds_GF) %>% unique() 

#########################################
#########################################
# # If possible we used the individual values of seeds
competition_fil_ind <- competition_fil %>%
  left_join(seed_data_raw, by = c("Plot", "Subplot", "Plant"))

competition_fil_ind$Seeds_GF[!is.na(competition_fil_ind$`Seeds/Fruit`)] <-
  competition_fil_ind$`Seeds/Fruit`[!is.na(competition_fil_ind$`Seeds/Fruit`)]

competition_fil_ind$Seeds[!is.na(competition_fil_ind$`Seeds/Fruit`)] <-
  competition_fil_ind$Seeds_GF[!is.na(competition_fil_ind$`Seeds/Fruit`)] *
  competition_fil_ind$Fruits[!is.na(competition_fil_ind$`Seeds/Fruit`)]

competition_fil <- competition_fil_ind[,1:8] %>% rename(Year = Year.x)

########################################
########################################
# There are plants with several measures

competition_fil %>% dplyr::select(Plot,Subplot,Plant) %>% group_by(Plot,Subplot,Plant)%>%
  count() %>% filter(n>1)


competition_fil <- competition_fil %>% group_by(Year,Plot,Subplot,Plant) %>%
  summarize(Seeds_GF = mean(Seeds_GF,na.rm = T),
            Fruit_GF = mean(Fruit_GF,na.rm = T),
            Seeds = mean(Seeds,na.rm = T),
            Fruits = mean(Fruits,na.rm = T)
            )

# Sanity checks
competition_fil$Plant %>% unique()

competition_fil %>% dplyr::select(Plot,Subplot,Plant) %>% group_by(Plot,Subplot,Plant)%>%
  count() %>% filter(n>1)


# fitness <- competition_fil %>% full_join(fitness_aux,by=c("Plot","Subplot",
#                                                           "Plant"))

fitness <- fitness_aux %>% full_join(competition_fil,by=c("Plot","Subplot",
                                                          "Plant"))

# Sanity check

fitness %>% dplyr::select(Plot,Subplot,Plant,ID) %>% group_by(Plot,Subplot,Plant,ID)%>%
  count() %>% filter(n>1)

#No todas las plantas de competencia tienen datos de fitness que no son promedio

# x <- competition_fil %>% full_join(seed_data,by=c("Plot","Subplot",
#                                                     "Plant"))

x <- fitness %>% filter(is.na(Fruits))

# Some focals (81) that received visits do not have information
# about the total amount of fruits
# we will use the average number of fruit of the species in a given plot
# If there are no data for a given plot, we will use data 
# the average number of fruit of the species in Caracoles

mean_fruits <- competition_fil %>% ungroup() %>% 
  dplyr::select(Plot,Plant,Fruits) %>% group_by(Plot,Plant)%>%
  dplyr::summarise_all(mean)

mean_fruits  %>% ungroup() %>% group_by(Plot)%>% count() #Notice that not all the
#plants are present in each subplot

# Estimating total seeds and total fruits for focals
for (i in 1:nrow(fitness)){
  
  # If there are no data for a given subplot, we will use data 
  # the average number of fruit of the species in the plot
  if (is.na(fitness$Fruits[i]) | is.nan(fitness$Fruits[i])){
  
    
    fitness$Fruits[i] <- competition_fil %>% ungroup() %>%
    filter(Plot==fitness$Plot[i],
           Plant== fitness$Plant[i]) %>%
      dplyr::select(Fruits) %>% pull() %>% mean()
    
  }
  
  # If there are no data for a plot, we will use data 
  # the average number of fruit of the species in Caracoles
  if (is.na(fitness$Fruits[i]) | is.nan(fitness$Fruits[i])){
    
    
    fitness$Fruits[i] <- competition_fil %>% ungroup() %>%
      filter(Plant== fitness$Plant[i]) %>%
      dplyr::select(Fruits) %>% pull() %>% mean()
    
  }
  
  if (is.na(fitness$Fruit_GF.x[i]) | is.nan(fitness$Fruit_GF.x[i])){ #Data from competition
    fitness$Fruit_GF.x[i] <- fitness$Fruits[i]
    fitness$Seeds_GF.x[i] <- fitness$Seeds[i]
  }else{ #Data from pollination surveys
    fitness$Fruit_GF.x[i] <- fitness$Fruits[i]
    fitness$Seeds_GF.x[i] <- fitness$Seeds_GF.x[i]*fitness$Fruits[i]
    }
}


# We use fruit and seed data from pollinator dataset

fitness_final <- fitness %>% ungroup() %>% 
  dplyr::select(Plot,Subplot,Plant,Seeds_GF.x,
                Fruit_GF.x,visits_GF,ID,homo_motif,hete_motif) %>%
  rename(Seeds_GF = Seeds_GF.x, Fruit_GF = Fruit_GF.x)


# Removing NAs from motifs, animals and visits
fitness_final$visits_GF[is.na(fitness_final$visits_GF)] <- 0
fitness_final$homo_motif[is.na(fitness_final$homo_motif)] <- 0
fitness_final$hete_motif[is.na(fitness_final$hete_motif)] <- 0
fitness_final$ID[is.na(fitness_final$ID)] <- "None"

# Sanity check

fitness_final %>% dplyr::select(Plot,Subplot,Plant,ID) %>% group_by(Plot,Subplot,Plant,ID)%>%
  count() %>% filter(n>1)

#9     E3      LEMA       Seeds_GF=NA     Fruit_GF=NA    visits_GF=3 Melyridae  0   0
fitness_final <- fitness_final %>% filter(!is.na(Seeds_GF))

#############################################
# ADDING CENTRALITY MEASSURES
#############################################

centrality <- read_csv("Processed_data/2020_NN_NEW_PageRank_results.csv") 

countSpaces <- function(s) { sapply(gregexpr(" ", s), function(p) { sum(p>=0) } ) }

centrality <- transform(centrality, baz = countSpaces(species))

centrality %>% filter(baz>1) %>% select(species) %>% unique()


centrality <- centrality %>%
  separate(species,c("Subplot","Plant")," ") %>% select(-baz)

fitness_final <- fitness_final %>% left_join(centrality, by=c("Plot","Subplot","Plant"))

# Sanity checks

fitness_final %>% filter(ID=="None",!is.na(Real_PR_Multi)) 
fitness_final %>% filter(ID!="None",is.na(Real_PR_Multi)) 


fitness_final <- fitness_final %>% filter(!is.na(Seeds_GF))

# Removing centrality NAs
fitness_final$Real_PR_Multi[is.na(fitness_final$Real_PR_Multi)] <- 0
fitness_final$Real_PR_Layer[is.na(fitness_final$Real_PR_Layer)] <- 0
fitness_final$Ratio[is.na(fitness_final$Ratio)] <- 1
fitness_final[is.na(fitness_final)] <- 0


write_csv(fitness_final,"Processed_data/2020_NN_NEW_data_models_phenol_overlap.csv")

