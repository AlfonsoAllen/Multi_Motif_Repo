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


fitness <- fitness2 %>% group_by(Plot,Subplot,Plant_Simple,G_F) %>%
  count(wt=Visits) %>% rename(Visits_tot = n, ID=G_F)

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
                            webID = fitness_data$Plant_Simple,
                            freq = fitness_data$Visits_tot)

list_incid_matrix_19 <- frame2webs(testdata_19,type.out="list")

  for (j in 1:length(list_incid_matrix_19)){

    plant_i <- j
    print(names(list_incid_matrix_19)[plant_i])
    
    incid_matrix_i <- list_incid_matrix_19[[plant_i]] 
    
    layer_csv <- paste("Processed_data/Multilayer_GF/Plot_",i,"_layer_",names(list_incid_matrix_19)[plant_i],".csv",sep="")
    write.csv(incid_matrix_i,layer_csv)
  }
}

fitness %>% group_by(Plant_Simple) %>%count()

###########################################
#Plants-interactions: Generating a bipartite network for each line
###########################################

fitness_L <- fitness %>% unite("Plot_Sub_Plant", c("Plot","Subplot_Plant_Label"), sep=" ", remove = FALSE)

fitness_L$Line <- NA

for (i in 1:nrow(fitness_L)){
  if(fitness_L$Plot[i] %in% c(1,2,3)){fitness_L$Line[i] <- 1}
  else if(fitness_L$Plot[i] %in% c(4,5,6)){fitness_L$Line[i] <- 2}
  else{fitness_L$Line[i] <- 3}
}

for (i in 1:3){ #i: Line
  
  fitness_data <- fitness_L %>% filter(Line==i)
  
  testdata_19 <-   data.frame(higher = fitness_data$ID,
                              lower = fitness_data$Plot_Sub_Plant,
                              webID = fitness_data$Plant_Simple,
                              freq = fitness_data$Visits_tot)
  
  list_incid_matrix_19 <- frame2webs(testdata_19,type.out="list")
  
  for (j in 1:length(list_incid_matrix_19)){
    
    plant_i <- j
    print(names(list_incid_matrix_19)[plant_i])
    
    incid_matrix_i <- list_incid_matrix_19[[plant_i]] 
    
    layer_csv <- paste("Processed_data/Multilayer_GF/Line_",i,"_layer_",names(list_incid_matrix_19)[plant_i],".csv",sep="")
    write.csv(incid_matrix_i,layer_csv)
  }
}

fitness_L %>% group_by(Plant_Simple) %>%count()


#######

###########################################
#Plants-interactions: Generating a bipartite network for CARACOLES
###########################################

fitness_C <- fitness %>% unite("Plot_Sub_Plant", c("Plot","Subplot_Plant_Label"), sep=" ", remove = FALSE)

  
testdata_19 <-   data.frame(higher = fitness_C$ID,
                            lower = fitness_C$Plot_Sub_Plant,
                            webID = fitness_C$Plant_Simple,
                            freq = fitness_C$Visits_tot)
  
list_incid_matrix_19 <- frame2webs(testdata_19,type.out="list")
  
for (j in 1:length(list_incid_matrix_19)){
    
    plant_i <- j
    print(names(list_incid_matrix_19)[plant_i])
    
    incid_matrix_i <- list_incid_matrix_19[[plant_i]] 
    
    layer_csv <- paste("Processed_data/Multilayer_GF/Caracoles_layer_",names(list_incid_matrix_19)[plant_i],".csv",sep="")
    write.csv(incid_matrix_i,layer_csv)
}
