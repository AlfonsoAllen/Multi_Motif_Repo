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


fitness %>% group_by(ID) %>% count()

###########################################
#Plants-interactions: Generating a bipartite network for each plot
###########################################

for (i in 1:9){ #i: Plot

  fitness_data <- fitness %>% filter(Plot==i)
  
  testdata_19 <-   data.frame(higher = fitness_data$ID,
                              lower = fitness_data$Subplot_Plant_Label,
                              webID = fitness_data$Plant,
                              freq = fitness_data$Visits_tot)
  
  list_incid_matrix_19 <- frame2webs(testdata_19,type.out="list")
  
    for (j in 1:length(list_incid_matrix_19)){
  
      plant_i <- j
      print(names(list_incid_matrix_19)[plant_i])
      
      incid_matrix_i <- list_incid_matrix_19[[plant_i]] 
      
      layer_csv <- paste("Processed_data/Multilayer_Species/2020_Plot_",i,"_layer_",names(list_incid_matrix_19)[plant_i],".csv",sep="")
      write.csv(incid_matrix_i,layer_csv)
  }
}

