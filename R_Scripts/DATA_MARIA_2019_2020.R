library(tidyverse)

# Load raw database

data_19 <- read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")
data_20 <- read_csv2("Raw_Data/raw_Pollinators_2020_1.csv") %>% 
  filter(!is.na(Plant),Plant!="0",Subplot!="OUT",Plant!="Ground")


data_20 %>% group_by(Plant) %>% count(wt=Visits)

ID_simple_G_F_19 <- data_19 %>% select(ID_Simple,G_F) %>% unique()

ID_simple_G_F_20_visits <- data_20 %>% select(ID_Simple,G_F,Visits)

ID_simple_G_F_20 <- ID_simple_G_F_20_visits %>% 
  filter(!ID_Simple %in% ID_simple_G_F_19$ID_Simple) %>% group_by(ID_Simple,G_F) %>%
  count(wt=Visits) %>% rename(cumulated_visits=n)

ID_simple_G_F_20 <- ID_simple_G_F_20_visits %>% 
  filter(!G_F %in% ID_simple_G_F_19$G_F) %>% group_by(ID_Simple,G_F) %>%
  count(wt=Visits) %>% rename(cumulated_visits=n)

data_20 %>% filter(ID_Simple %in% c("1/29/5","2_28_5"))

write_csv(ID_simple_G_F_20,"new_IDs_2020.csv")
