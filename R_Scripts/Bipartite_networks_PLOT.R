# load libraries
library(tidyverse)
library(bipartite)
library(matlib)

####################################################################
# Loadind Plant-pollinator dataset (Caracoles): visits, abundances, seeds
####################################################################
fitness_data2 <- read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")

# Filtering & relabeling
fitness2 <- fitness_data2 %>% filter(Year==2019)


fitness <- fitness2 %>% group_by(Plot,Subplot,Plant_Simple,ID_Simple) %>%
  count(wt=Visits) %>% rename(Visits_tot = n, ID=ID_Simple)

fitness <- fitness %>% mutate(Subplot_Plant_Label=paste(Subplot,Plant_Simple,sep=" "))%>%
  filter(!Plant_Simple %in% c("Lysimachia_arvensis","HOMA"))


fitness %>% group_by(ID) %>% count()

###########################################
#Plants-interactions: Generating a bipartite network for each plot
###########################################

for (i in 1:9){ #i: Plot

fitness_data <- fitness %>% filter(Plot==i)

testdata_19 <-   data.frame(higher = fitness_data$ID,
                            lower = fitness_data$Subplot_Plant_Label,
                            webID = fitness_data$Plot,
                            freq = fitness_data$Visits_tot)

list_incid_matrix_19 <- frame2webs(testdata_19,type.out="list")

incid_matrix_i <- list_incid_matrix_19[[1]] 

write.table(incid_matrix_i, file = paste0("Plot",i,".txt"),
            row.names=FALSE, col.names = F)
}
