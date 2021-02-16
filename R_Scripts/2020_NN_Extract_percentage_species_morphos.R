# Percentage de species and morfos

# load libraries
library(tidyverse)
library(bipartite)
library(matlib)

####################################################################
# Loadind Plant-pollinator dataset (Caracoles): visits, abundances, seeds
####################################################################
fitness_data2 <- read_csv2("Raw_Data/final_Pollinators_2020.csv")

# Inspect ID, ID_Simple

fitness_data2$ID %>% unique() %>% sort()
fitness_data2$ID_Simple %>% unique() %>% sort()

fitness_data2$ID[grep(" ",fitness_data2$ID)] # No labels contain spaces
fitness_data2$ID_Simple[grep(" ",fitness_data2$ID_Simple)]# No labels contain spaces

fitness_data2$ID[grep("\\.",fitness_data2$ID)] # Labels contain dots
fitness_data2$ID_Simple[grep("\\.",fitness_data2$ID_Simple)]# Labels contain dots

# Remove points from ID names
fitness_data2$ID <- sub("\\.", "", fitness_data2$ID)
fitness_data2$ID_Simple <- sub("\\.", "", fitness_data2$ID_Simple)

# filter tabanidae
fitness_data2 <- fitness_data2 %>% filter(ID != "Tabanidae")

# Filtering & relabeling
fitness2 <- fitness_data2 %>% filter(!is.na(Plant),
                                     Plant!="0",
                                     Subplot!="OUT",
                                     Plant!="Ground")


fitness <- fitness2 %>% group_by(Plot,Subplot,Plant,ID_Simple) %>%
  count(wt=Visits) %>% rename(Visits_tot = n,ID=ID_Simple)

fitness <- fitness %>% mutate(Subplot_Plant_Label=paste(Subplot,Plant,sep=" "))


fitness_ID <- fitness %>% group_by(ID) %>% count()

# Amount of morphos
n_morphos <- sum(fitness_ID$n[grepl("_sp", fitness_ID$ID, fixed = TRUE)])
# Amount of species
n_species <- sum(fitness_ID$n[grepl("_", fitness_ID$ID, fixed = TRUE)]) - n_morphos
#Amount family
n_family <- sum(fitness_ID$n) - n_morphos - n_species

#percentage of morphos
round(100 * n_morphos / (n_morphos + n_species + n_family), 2)

#percentage of species
round(100 * n_species / (n_morphos + n_species + n_family), 2)

#percentage of family
round(100 * n_family / (n_morphos + n_species + n_family), 2)
