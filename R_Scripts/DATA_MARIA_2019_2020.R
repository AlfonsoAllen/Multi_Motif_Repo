library(tidyverse)

# Load raw database

data_19 <- read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")
data_20 <- read_csv2("Raw_Data/Pollinators_2020_unpocomaslimpia.csv") %>% 
  filter(!is.na(Plant),Plant!="0",Subplot!="OUT",Plant!="Ground")


data_20 %>% group_by(Plant) %>% count(wt=`#visits`)

ID_simple_G_F_19 <- data_19 %>% select(ID_Simple,G_F) %>% unique()

ID_20 <- data_20 %>% select(ID) %>% unique()

ID_20$ID[ID_20$ID %in% ID_simple_G_F_19$ID_Simple]

new_ID <- ID_20$ID[!ID_20$ID %in% ID_simple_G_F_19$ID_Simple]

new_IDs_visits <- data_20 %>% filter(ID %in% new_ID) %>% group_by(ID) %>%count(wt=`#visits`)
