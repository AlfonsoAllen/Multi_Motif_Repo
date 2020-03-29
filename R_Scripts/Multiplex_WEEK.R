# load libraries
library(tidyverse)
library(bipartite)
library(matlib)

####################################################################
# Loadind Plant-pollinator dataset (Caracoles): visits, abundances, seeds
####################################################################

fitness <- read_csv("Processed_data/Motifs_WEEK/Caracoles_WEEK.csv")

fitness <- fitness %>% separate(Subplot_Plant_Label,c("Subplot","Plant_Simple"),
                                          sep = " ",remove = FALSE) %>%
  filter(!Plant_Simple %in% c("Lysimachia_arvensis","HOMA"))

fitness %>% group_by(G_F) %>% count()

###########################################
#Plants-interactions: Generating a bipartite network for each plot
###########################################

for (i in 1:9){ #i: Plot

fitness_data <- fitness %>% filter(Plot==i)

testdata_19 <-   data.frame(higher = fitness_data$G_F,
                            lower = fitness_data$Subplot_Plant_Label,
                            webID = fitness_data$Week,
                            freq = fitness_data$Visits_tot)

list_incid_matrix_19 <- frame2webs(testdata_19,type.out="list")

  for (j in 1:length(list_incid_matrix_19)){

    plant_i <- j
    print(names(list_incid_matrix_19)[plant_i])
    
    incid_matrix_i <- list_incid_matrix_19[[plant_i]] 
    
    layer_csv <- paste("Processed_data/Multilayer_WEEK/Plot_",i,"_layer_WEEK_",names(list_incid_matrix_19)[plant_i],".csv",sep="")
    write.csv(incid_matrix_i,layer_csv)
  }
}

fitness %>% group_by(Week) %>%count()


###########################################
#Plants-interactions: Generating a bipartite network for CARACOLES
###########################################

fitness_C <- fitness %>% unite("Plot_Sub_Plant", c("Plot","Subplot_Plant_Label"), sep=" ", remove = FALSE)

  
testdata_19 <-   data.frame(higher = fitness_C$G_F,
                            lower = fitness_C$Plot_Sub_Plant,
                            webID = fitness_C$Week,
                            freq = fitness_C$Visits_tot)
  
list_incid_matrix_19 <- frame2webs(testdata_19,type.out="list")
  
for (j in 1:length(list_incid_matrix_19)){
    
    plant_i <- j
    print(names(list_incid_matrix_19)[plant_i])
    
    incid_matrix_i <- list_incid_matrix_19[[plant_i]] 
    
    layer_csv <- paste("Processed_data/Multilayer/Caracoles_layer_WEEK_",names(list_incid_matrix_19)[plant_i],".csv",sep="")
    write.csv(incid_matrix_i,layer_csv)
}

