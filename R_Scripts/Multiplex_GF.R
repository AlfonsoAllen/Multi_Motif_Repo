# load libraries
library(tidyverse)
library(bipartite)
library(matlib)

####################################################################
# Loadind Plant-pollinator dataset (Caracoles): visits, abundances, seeds
####################################################################

fitness_data2 <- read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")

# Filtering & relabeling
fitness2 <- fitness_data2 %>% filter(Year==2019,!Subplot == "OUT" & !is.na(G_F))


fitness <- fitness2 %>% group_by(Plot,Subplot,Plant_Simple,G_F) %>%
  count(wt=Visits) %>% rename(Visits_tot = n)

fitness <- fitness %>% mutate(Subplot_Plant_Label=paste(Subplot,Plant_Simple,sep=" "))%>%
  filter(!Plant_Simple %in% c("Lysimachia_arvensis","HOMA"))

fitness %>% group_by(G_F) %>% count()

###########################################
#Plants-interactions: Generating a bipartite network for each plot
###########################################

for (i in 1:9){ #i: Plot

fitness_data <- fitness %>% filter(Plot==i)

testdata_19 <-   data.frame(higher = fitness_data$G_F,
                            lower = fitness_data$Subplot_Plant_Label,
                            webID = fitness_data$Plant_Simple,
                            freq = fitness_data$Visits_tot)

list_incid_matrix_19 <- frame2webs(testdata_19,type.out="list")

  for (j in 1:length(list_incid_matrix_19)){

    plant_i <- j
    print(names(list_incid_matrix_19)[plant_i])
    
    incid_matrix_i <- list_incid_matrix_19[[plant_i]] 
    
    layer_csv <- paste("Processed_data/Multilayer/Plot_",i,"_layer_",names(list_incid_matrix_19)[plant_i],".csv",sep="")
    write.csv(incid_matrix_i,layer_csv)
  }
}

fitness %>% group_by(Plant_Simple) %>%count()


###########################################
#Plants-interactions: Generating a bipartite network for CARACOLES
###########################################

fitness_C <- fitness %>% unite("Plot_Sub_Plant", c("Plot","Subplot_Plant_Label"), sep=" ", remove = FALSE)

  
testdata_19 <-   data.frame(higher = fitness_C$G_F,
                            lower = fitness_C$Plot_Sub_Plant,
                            webID = fitness_C$Plant_Simple,
                            freq = fitness_C$Visits_tot)
  
list_incid_matrix_19 <- frame2webs(testdata_19,type.out="list")
  
for (j in 1:length(list_incid_matrix_19)){
    
    plant_i <- j
    print(names(list_incid_matrix_19)[plant_i])
    
    incid_matrix_i <- list_incid_matrix_19[[plant_i]] 
    
    layer_csv <- paste("Processed_data/Multilayer/Caracoles_layer_",names(list_incid_matrix_19)[plant_i],".csv",sep="")
    write.csv(incid_matrix_i,layer_csv)
}
