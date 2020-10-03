# We will add to the previous data of floral visitors and competition 
# data for Rape

library(tidyverse)

# Old data

old_data <- read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")

# Vistis RAPE

Raw_RAPE <- read_csv2("Raw_Data/Metadata_Pollinators_2019_2016_bueno.csv") %>%
  as_tibble() %>%
  filter(Year == 2019,
         Plant_Simple=="RAPE",
         Plot !="OUT", Subplot != "OUT",
         Order %in% c("Coleoptera", "Diptera", "Hymenoptera", "Lepidoptera" ),
         !G_F %in% c("Ants","Mosquitoes"),
         !ID_Simple %in% c("Coccinella_septempunctata","Larva",
                           "Chrysididae","Diplazon_sp.") ) %>%
  select(Day,Month,Year,Plot,Subplot,Group,Order,Superfamily,Family,Species,
  ID,ID_Simple,G_F,Visits,Plant_Simple)

#names(read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv"))
# NAME OF FIELDS

# [1] "Day"          "Month"        "Year"         "Plot"         "Subplot"     
# [6] "Group"        "Order"        "Superfamily"  "Family"       "Species"     
# [11] "ID"           "ID_Simple"    "G_F"          "Visits"       "Plant_Simple"
# [16] "Fruit"        "Seed"  


# Fruit and seed headers are missing

##################
# ADDING SEED AND FRUIT

# DAVID: es un problema no resuelto. 
# en mis análisis (archivo competencia.csv en /data) 
# creo recordar que cogí el valor que tuviera mayor número de semillas. 
# Pregunté en su momento y no supieron decirme a qué se debían esas 
# duplicidades
# puedes usar ese archivo (espera, se llama competition.csv, 
#                          jaja, estoy hablando de memoria)
# ahí deberían estar los datos que yo limpié para 2015:2019

seeds_fruits_2019 <- read_csv2("Raw_Data/competition.csv") %>%
  rename(Plot=plot,Subplot=subplot,
         Plant_Simple=focal,Fruit=fruit,Seed=seed,Year=year) %>%
  filter(Year==2019) %>%
  select(-neighbour,-number) %>% unique()

#Sanity check
x <- seeds_fruits_2019  %>%
  filter(Plant_Simple %in% c("LEMA","CHFU","ME","PUPA", "CHMI", "RAPE")) %>%
  group_by(Plot, Subplot, Plant_Simple) %>% 
  count() %>% filter(n>1)

seeds_fruits_2019  %>%
  group_by(Plant_Simple) %>% 
  count()

seeds_fruits_2019  %>%
  filter(Plant_Simple=="RAPE") #NO FRUITS FOR RAPE in competition.csv

# There are also duplicities. Example:

seeds_fruits_2019 %>% filter(Plot==2, Subplot=="F5", Plant_Simple=="CHMI")

# Data Objective
Raw_Maria <- read_csv2("Raw_Data/Data_Caracoles_pol_plants_2019.txt")
Raw_Maria %>% filter(Plot==2, Subplot=="F5", Plant_Simple=="CHMI")

# Apparently there is no need to modify competitions records. 
# Focal with potential duplicities does not receive floral visitors

#############################
# We look for seeds and fruits in Seeds_and_fruits_caracoles_simplified_2019.csv 

read_csv2("Raw_Data/Seeds_and_fruits_caracoles_simplified_2019.csv")%>%
  group_by(Plant_Simple) %>% 
  count() 

# No matches for RAPE

##############################

# We look for seeds and fruits in competition_2019.csv

read_csv2("Raw_Data/competition_2019.csv") %>%
  select(-competitor,-number) %>% unique() %>%
  group_by(focal) %>%
  count()

# THERE ARE RAPES

RAPES_com <- read_csv2("Raw_Data/competition_2019.csv")%>%
  select(-competitor,-number) %>% unique() %>% 
  filter(focal=="RAPE") # There are rapes with visitors

read_csv2("Raw_Data/competition_2019.csv")%>%
  select(-competitor,-number,-day,-month) %>% unique() %>% 
  filter(focal=="RAPE") %>% 
  group_by(plot,subplot) %>% count() %>% filter(n>1)
# There are no duplicated RAPES

#######################
# We save and prepare the results for RAPE
seeds_fruits_2019_RAPE <- read_csv2("Raw_Data/competition_2019.csv")%>%
  select(-competitor,-number,-day,-month) %>% unique() %>% 
  filter(focal=="RAPE") %>%
  rename(Plot=plot,Subplot=subplot,Plant_Simple=focal,
         Fruit=fruit,Year=year)%>%
  mutate(Seed=8*Fruit)

# Info Oscar: Alfonso, no sé los datos que habrá, 
# RAPE produce 8 semillas por fruto, 
# lo digo por si solo salen frutos

final_seeds_fruits_2019 <- bind_rows(seeds_fruits_2019,
                                     seeds_fruits_2019_RAPE)

###########################

RAPE_Seed <- Raw_RAPE %>% 
  left_join(final_seeds_fruits_2019,
            by=c("Year","Plot","Subplot","Plant_Simple"))

# There are 2 rape focals without seed data


###########################

new_data <- bind_rows(old_data,RAPE_Seed)

###########################
# Sanity check

new_data %>% filter(Seed == 0 | is.na(Seed)) # 127 entries without seeds

new_data %>% filter(Visits == 0) # Every row refers to observed visits

#############################

write_csv(new_data,"Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID_RAPE.csv")


final_seeds_fruits_2019 <- final_seeds_fruits_2019 %>%
  rename(plot=Plot,subplot=Subplot,
         focal=Plant_Simple,fruit=Fruit,seed=Seed,year=Year)


write_csv(final_seeds_fruits_2019,"Raw_Data/competition_RAPE.csv")
