# load libraries
library(tidyverse)


####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################

metadata <- read_csv2("Metadata_Pollinators_Plants_fruits_abundancesmodified.csv")

metadata_19 <- metadata %>% filter(Year==2019,Subplot!="OUT",!is.na(G_F)) %>% #QUITAR OUTS
  select(Plot,Subplot,G_F,Plant_Simple,Visits)

metadata_19 %>% group_by(Plant_Simple)%>%count(wt=Visits)

metadata_19 %>% filter(Plant_Simple=="HOMA")

seeds_caracoles_2019 <- read_csv2("Seeds_caracoles_2019.csv")
seeds_caracoles_2019 %>% group_by(PLANT)%>%count()

abundances_2019 <- read_csv2("abundances_2019.csv")
seeds_caracoles_2019 %>% group_by(PLANT)%>%count()