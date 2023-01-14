
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
  dplyr::select(Plot,Subplot,Plant,`Seeds/Fruit`) %>% unique() %>%
  rename(Seed=`Seeds/Fruit`) %>% group_by(Plot,Subplot,Plant) %>%
  summarize(Seed = mean(Seed,na.rm = T)) %>% arrange(Plot,Subplot)


fitness2_seeds_aux <- fitness2 %>%
  left_join(seed_data,by=c("Plot","Subplot","Plant"))

# There are plants with Plot and subplot set to NA. We add to those their mean values

fitness2_seeds_with_ind_data <- fitness2_seeds_aux %>% filter(!is.na(Seed))
fitness2_seeds_without_ind_data <- fitness2_seeds_aux %>% filter(is.na(Seed)) %>%
  dplyr::select(-Seed)

mean_seed_plot <- seed_data_raw %>% 
  dplyr::select(Plot,Subplot,Plant,`Seeds/Fruit`) %>% #unique() %>%
  rename(Seed=`Seeds/Fruit`) %>% group_by(Plot,Plant) %>%
  summarize(Seed = mean(Seed,na.rm = T)) %>% arrange(Plot)


# fitness2_seeds_without_ind_data <- fitness2_seeds_without_ind_data %>%
#   left_join(mean_seed,by="Plant")


fitness2_seeds_without_ind_data <- fitness2_seeds_without_ind_data %>%
  left_join(mean_seed_plot,by=c("Plot","Plant"))

# There are 24 entries without seed per fruit information and 3 entries with 0 seeds

sum(is.na(fitness2_seeds_without_ind_data$Seed))
sum(fitness2_seeds_without_ind_data$Seed[!is.na(fitness2_seeds_without_ind_data$Seed)]==0)

fitness2_seeds_without_ind_data[is.na(fitness2_seeds_without_ind_data$Seed),]


mean_seeds_caracoles <- seed_data_raw %>% rename(Seed=`Seeds/Fruit`) %>%
  dplyr::select(Plot,Subplot,Plant,Seed) %>% #unique() %>%
  group_by(Plant) %>%
  summarize(Seed = mean(Seed,na.rm = T)) %>% arrange(Plant)

fitness2_seeds_with_mean_Caracoles_data <- fitness2_seeds_without_ind_data %>% 
  filter(is.na(Seed)) %>%
  dplyr::select(-Seed) %>% left_join(mean_seeds_caracoles, by = "Plant")

fitness2_seeds_with_mean_plot_data <- fitness2_seeds_without_ind_data %>% 
  filter(!is.na(Seed))




# Merge data with individual values and mean values, respectively

fitness2_seeds <- bind_rows(fitness2_seeds_with_ind_data %>% mutate(type_seed_per_fruit = "Individual"),
                            fitness2_seeds_with_mean_plot_data %>% mutate(type_seed_per_fruit = "Plot"),
                            fitness2_seeds_with_mean_Caracoles_data %>% mutate(type_seed_per_fruit = "Caracoles"))

# Sanity check
fitness2_seeds %>% filter(is.na(Seed))


fitness3_aux <- fitness2_seeds %>% group_by(Line,Plot,Subplot,Plant,Week,ID,type_seed_per_fruit) %>%
  summarize(Seeds_GF = mean(Seed,na.rm = T),
            visits_GF = sum(Visits,na.rm = T))

fitness3 <- fitness3_aux %>% filter(!is.nan(Seeds_GF))

fitness3$Plant %>% unique()

# Sanity check

fitness3 %>% dplyr::select(Plot,Subplot,Plant,ID,Week,type_seed_per_fruit) %>% group_by(Plot,Subplot,Plant,ID,Week,type_seed_per_fruit)%>%
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

fitness_aux_motifs <- fitness3 %>% left_join(caracoles_motif, by=c("Plot","Subplot","Plant","ID","Week"))

#Sanity check: visits from caracoles and fitness3 are equal
fitness_aux_motifs %>% filter(visits_GF!=Visits_tot)

##################################
# Agregating by week
##################################

fitness_aux <- fitness_aux_motifs %>% group_by(Line,Plot,Subplot,Plant,ID,type_seed_per_fruit) %>%
  summarise(
    Seeds_GF = mean(Seeds_GF,na.rm = T),
    visits_GF = sum(visits_GF,na.rm = T),
    Visits_tot = sum(Visits_tot,na.rm = T),
    homo_motif = sum(homo_motif,na.rm = T),
    hete_motif = sum(hete_motif,na.rm = T)
  )

# sanity check

fitness_aux %>% dplyr::select(Plot,Subplot,Plant,ID,type_seed_per_fruit) %>% group_by(Plot,Subplot,Plant,ID,type_seed_per_fruit)%>%
  count() %>% filter(n>1)

#####################################
# COMPETITION
#####################################
# We use these data to obtain the total number of fruits per plant and then
# estimate the total number of seeds

fitness_aux %>% group_by(Plant) %>% count()

# En 2020 no hay MEEL

competition <- read_csv2("Raw_Data/competition_caracoles2020.csv") %>%
  filter(!is.na(Seeds))


competition_fil <- competition %>% 
  filter(Sp.Focal %in% c("BEMA","CETE","CHFU","CHMI","LEMA","MESU","PUPA",
                         "SCLA","SOAS","SPRU")) %>%
  dplyr::select(Year,Plot,Subplot,Sp.Focal,Fruits,Seeds) %>%
  unique() %>% dplyr::select(-Seeds) %>%
  group_by(Year,Plot,Subplot,Sp.Focal) %>%
  summarise(
    Fruit_GF = mean(Fruits,na.rm = T)
  ) %>% ungroup() %>% rename(Plant = Sp.Focal)



#########################################
#########################################



fitness_4 <- fitness_aux %>% full_join(competition_fil,by=c("Plot","Subplot",
                                                          "Plant"))

fitness_4$Year <- 2020

fitness_4_fruit_with_ind_data <- fitness_4 %>% filter(!is.na(Fruit_GF))
fitness_4_fruit_without_ind_data <- fitness_4 %>% filter(is.na(Fruit_GF)) %>%
  dplyr::select(-Fruit_GF)

mean_fruit_plot <- competition_fil %>% 
  dplyr::select(Plot,Subplot,Plant,Fruit_GF) %>% #unique() %>%
  group_by(Plot,Plant) %>%
  summarize(Fruit_GF = mean(Fruit_GF,na.rm = T)) %>% arrange(Plot)


fitness_4_fruit_with_mean_plot_data <- fitness_4_fruit_without_ind_data %>%
  left_join(mean_fruit_plot,by=c("Plot","Plant"))

# There are 19 entries without fruit information. All are LEMA in PLOT 6

sum(is.na(fitness_4_fruit_with_mean_plot_data$Fruit_GF))
sum(fitness_4_fruit_with_mean_plot_data$Fruit_GF[!is.na(fitness_4_fruit_with_mean_plot_data$Fruit_GF)]==0)

fitness_4_fruit_with_mean_plot_data[is.na(fitness_4_fruit_with_mean_plot_data$Fruit_GF),]

mean_fruit_caracoles <- competition_fil %>% 
  dplyr::select(Plot,Subplot,Plant,Fruit_GF) %>% #unique() %>%
  group_by(Plant) %>%
  summarize(Fruit_GF = mean(Fruit_GF,na.rm = T)) %>% arrange(Plant)

fitness_4_fruit_with_mean_caracoles_data <- fitness_4_fruit_with_mean_plot_data %>% 
  filter(is.na(Fruit_GF)) %>%
  dplyr::select(-Fruit_GF) %>% left_join(mean_fruit_caracoles, by = "Plant")

fitness_4_fruit_with_mean_plot_data <- fitness_4_fruit_with_mean_plot_data %>% 
  filter(!is.na(Fruit_GF))

fitness_5 <- bind_rows(fitness_4_fruit_with_ind_data %>% mutate(type_fruit = "Individual"),
                       fitness_4_fruit_with_mean_caracoles_data %>% mutate(type_fruit = "Caracoles"),
                       fitness_4_fruit_with_mean_plot_data  %>% mutate(type_fruit = "Plot")) %>%
  dplyr::select(-Visits_tot) %>%
  mutate(seed_per_fruit = Seeds_GF,
    Seeds_GF = seed_per_fruit * Fruit_GF)

fitness_5 %>% filter(is.na(Seeds_GF))
fitness_5 %>% filter(is.na(Fruit_GF))

for (i in 1:nrow(fitness_5 )){
  if(fitness_5$Plot[i] %in% c(1,2,3)){fitness_5$Line[i] <- 1}
  else if(fitness_5$Plot[i] %in% c(4,5,6)){fitness_5$Line[i] <- 2}
  else{fitness_5$Line[i] <- 3}
} 


fitness_5_with_seeds <- fitness_5 %>% filter(!is.na(Seeds_GF))

fitness_5_no_seed_info <- fitness_5 %>% filter(is.na(Seeds_GF)) %>%
  select(-Seeds_GF,-seed_per_fruit) %>%
  left_join(seed_data, by = c("Plot", "Subplot","Plant")) %>%
  rename(seed_per_fruit = Seed) %>%
  mutate(Seeds_GF = Fruit_GF * seed_per_fruit)

fitness_5_with_seeds_2 <- fitness_5_no_seed_info %>% filter(!is.na(Seeds_GF))
fitness_5_without_seeds_2 <- fitness_5_no_seed_info %>% filter(is.na(Seeds_GF))


fitness_5_mean_seeds <- fitness_5_without_seeds_2 %>%
  select(-Seeds_GF,-seed_per_fruit) %>%
  left_join(mean_seed_plot,by=c("Plot","Plant"))%>%
  rename(seed_per_fruit = Seed) %>%
  mutate(Seeds_GF = Fruit_GF * seed_per_fruit)

# There are 415 entries without seed per fruit information and 23 entries with 0 seeds

sum(is.na(fitness_5_mean_seeds$Seeds_GF))
sum(fitness_5_mean_seeds$Seeds_GF[!is.na(fitness_5_mean_seeds$Seeds_GF)]==0)

fitness_5_mean_seeds[is.na(fitness_5_mean_seeds$Seeds_GF),]

fitness_5_seeds_with_mean_Caracoles_data <- fitness_5_mean_seeds %>% 
  filter(is.na(Seeds_GF)) %>%
  dplyr::select(-Seeds_GF,-seed_per_fruit) %>% 
  left_join(mean_seeds_caracoles, by = "Plant")%>%
  rename(seed_per_fruit = Seed) %>%
  mutate(Seeds_GF = Fruit_GF * seed_per_fruit)

fitness_5_seeds_with_mean_plot_data <- fitness_5_mean_seeds %>% 
  filter(!is.na(Seeds_GF))


# Merge data with individual values and mean values, respectively

fitness_5_final <- bind_rows(fitness_5_with_seeds %>% mutate(type_seed_per_fruit = "Individual"),
                       fitness_5_with_seeds_2  %>% mutate(type_seed_per_fruit = "Individual"),
                       fitness_5_seeds_with_mean_plot_data  %>% mutate(type_seed_per_fruit = "Plot"),
                       fitness_5_seeds_with_mean_Caracoles_data %>% mutate(type_seed_per_fruit = "Caracoles"))


fitness_5_final %>% filter(is.na(Seeds_GF))
fitness_5_final %>% filter(is.na(Fruit_GF))
  
########################################
########################################


# Removing NAs from motifs, animals and visits
fitness_5_final$visits_GF[is.na(fitness_5_final$visits_GF)] <- 0
fitness_5_final$homo_motif[is.na(fitness_5_final$homo_motif)] <- 0
fitness_5_final$hete_motif[is.na(fitness_5_final$hete_motif)] <- 0
fitness_5_final$ID[is.na(fitness_5_final$ID)] <- "None"

# Sanity check

fitness_5_final %>% dplyr::select(Plot,Subplot,Plant,ID) %>% group_by(Plot,Subplot,Plant,ID)%>%
  count() %>% filter(n>1)

fitness_6 <- fitness_5_final %>% filter(!is.na(Seeds_GF))

#############################################
# ADDING CENTRALITY MEASSURES
#############################################

centrality <- read_csv("Processed_data/2020_NN_NEW_PageRank_results.csv") 

countSpaces <- function(s) { sapply(gregexpr(" ", s), function(p) { sum(p>=0) } ) }

centrality <- transform(centrality, baz = countSpaces(species))

centrality %>% filter(baz>1) %>% select(species) %>% unique()


centrality <- centrality %>%
  separate(species,c("Subplot","Plant")," ") %>% select(-baz)

fitness_final <- fitness_6 %>% left_join(centrality, by=c("Plot","Subplot","Plant"))

# Sanity checks

fitness_final %>% filter(ID=="None",!is.na(Real_PR_Multi)) 
fitness_final %>% filter(ID!="None",is.na(Real_PR_Multi)) 


fitness_final <- fitness_final %>% filter(!is.na(Seeds_GF))

# Removing centrality NAs
fitness_final$Real_PR_Multi[is.na(fitness_final$Real_PR_Multi)] <- 0
fitness_final$Real_PR_Layer[is.na(fitness_final$Real_PR_Layer)] <- 0
fitness_final$Ratio[is.na(fitness_final$Ratio)] <- 1
fitness_final[is.na(fitness_final)] <- 0


write_csv(fitness_final,"Processed_data/2020_NN_NEW_data_models_phenol_overlap_2.csv")

