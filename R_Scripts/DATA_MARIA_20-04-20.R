library(tidyverse)

# Load raw database
#Raw_with_date <- read.table("Raw_Data/Metadata_Pollinators_2019_2016_bueno.txt",
#                            header=T, sep=";")
Raw_with_date <- read_csv2("Raw_Data/Metadata_Pollinators_2019_2016_bueno.txt")
str(Raw_with_date)

Raw_with_date <- as_tibble(Raw_with_date)

# Filter data
Raw_with_date_f <- filter(Raw_with_date,Year == 2019,
                        Plant_Simple %in% c("LEMA","CHFU","ME","PUPA", "CHMI"),
                        Plot !="OUT", Subplot != "OUT",
                        Order %in% c("Coleoptera", "Diptera", "Hymenoptera", "Lepidoptera" ),
                        !G_F %in% c("Ants","Mosquitoes"),
                        !ID_Simple %in% c("Coccinella_septempunctata","Larva",
                                         "Chrysididae","Diplazon_sp.") )




# Remove columns
Raw_with_date_f <- Raw_with_date_f %>% select(-Sex,-ID2,-Jar,-YESNO_jar,-Time,-Obs.comments,
                                          -Weather,-Input_by,-X24)

# Group visits by location, pollinator taxonomy and plant label
Raw_group <- Raw_with_date_f %>% group_by(Plot, Subplot, Plant_Simple,Order,Family, Species, ID, G_F,
                             ID_Simple, Plant_Simple) %>% 
  summarise (num.visits = sum(Visits))

# Data Objective
Raw_Maria <- read_csv2("Raw_Data/Data_Caracoles_pol_plants_2019.txt")


###############################################################################
# Sanity check
nrow(Raw_Maria)==nrow(Raw_group)

# Comparison between Raw_group y Raw_Maria

a1 <- Raw_group  %>% select(Plot, Subplot, Plant_Simple,Order,Family, Species, ID, G_F,
                                                ID_Simple, Plant_Simple,num.visits)
a2 <- Raw_Maria %>% select(Plot, Subplot, Plant_Simple,Order,Family, Species, ID, G_F,
                                                ID_Simple, Plant_Simple,num.visits)
anti_join(a2,a1)
difference <- setdiff(a1,a2)

# NOTE: We can leave the additional entries as morphospecies

###############################################################

Raw_Maria_aux <- Raw_Maria %>% select(Plot, Subplot, Plant_Simple,Order,Family, Species, ID, G_F,
                           ID_Simple, Fruit,Seed)

Raw_with_date_final <- Raw_with_date_f %>%
  left_join(Raw_Maria_aux, by = c("Plot", "Subplot", "Plant_Simple","Order","Family",
                                  "Species", "ID", "G_F",
                                  "ID_Simple"))

###########################
# Sanity check

Raw_Maria_aux %>% filter(Seed == 0) # 83 entries without seeds

Raw_with_date_final %>% filter(Visits == 0) # Every row refers to observed visits

#############################

write_csv(Raw_with_date_final,"Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")



read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")%>% 
  filter(Year==2019)%>% group_by(G_F) %>%
  count(wt=Visits) 

read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")%>% 
  filter(Year==2019)%>% group_by(G_F,Order,Family, Species) %>%
  count(wt=Visits) 

read_csv2("Raw_Data/Metadata_Pollinators_2019_2016_bueno_MARIA_Sept_2020.csv") %>%
  filter(Year == 2019,
         Plant_Simple %in% c("LEMA","CHFU","ME","PUPA", "CHMI"),
         Plot !="OUT", Subplot != "OUT",
         Order %in% c("Coleoptera", "Diptera", "Hymenoptera", "Lepidoptera" ),
         !G_F %in% c("Ants","Mosquitoes"),
         !ID_Simple %in% c("Coccinella_septempunctata","Larva",
                           "Chrysididae","Diplazon_sp.") ) %>%
  group_by(G_F) %>%
  count(wt=Visits) 
  
