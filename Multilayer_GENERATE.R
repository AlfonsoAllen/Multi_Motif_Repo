# load libraries
library(tidyverse)
library(bipartite)
library(matlib)

####################################################################
# Loadind Plant-pollinator dataset (Caracoles): visits, abundances, seeds
####################################################################

fitness_data <- read_csv2("FV_2019_FINAL_visits_abundances_seeds.csv")

fitness_data %>% group_by(G_F) %>% count()

fitness <- fitness_data %>% group_by(Plot,Subplot,Plant_Simple,G_F,num.plantas,Fruit, Seed) %>%
  count(wt=Visits) %>% rename(Visits_tot = n)

#####################################
# Filtering & relabeling
#####################################

fitness <- fitness %>% filter(!Subplot == "OUT" & !is.na(G_F))

fitness <- fitness %>% mutate(Seeds_tot = num.plantas*Seed,
                              weight = Visits_tot/num.plantas)


#####################################
# Data for individuals
#####################################

for (i in 1:nrow(fitness)){
  for (j in 1:fitness$num.plantas[i]){
    
    if (i==1 & j == 1){
      
      fitness_ind <- fitness[i,] %>% mutate(plant_ind=j)
      
    } else {
      aux <- fitness[i,] %>% mutate(plant_ind=j)
      fitness_ind <- fitness_ind %>% bind_rows(aux)
      
    }
  }
}

fitness_ind$Plant_Label <- paste(fitness_ind$Subplot,
                             fitness_ind$Plant_Simple,
                             fitness_ind$plant_ind,
                             sep = " ")


###########################################
#Plants-interactions: Generating a bipartite network for each plot
###########################################

for (i in 1:9){ #i: Plot

fitness_ind_plot <- fitness_ind %>% filter(Plot==i)

testdata_19 <-   data.frame(higher = fitness_ind_plot$G_F,
                            lower = fitness_ind_plot$Plant_Label,
                            webID = fitness_ind_plot$Plant_Simple,
                            freq = fitness_ind_plot$weight)

list_incid_matrix_19 <- frame2webs(testdata_19,type.out="list")

  for (j in 1:length(list_incid_matrix_19)){
  
    plant_i <- j
    print(names(list_incid_matrix_19)[plant_i])
    
    incid_matrix_i <- list_incid_matrix_19[[plant_i]] 
    
    
    layer_csv <- paste("Multilayer/Plot_",i,"_layer_",names(list_incid_matrix_19)[plant_i],".csv",sep="")
    write.csv(incid_matrix_i,layer_csv)
  }
}

fitness %>% group_by(Plant_Simple) %>%count()


###########################################
#Plants-interactions: Generating a bipartite network for CARACOLES
###########################################

  
testdata_19 <-   data.frame(higher = fitness_ind$G_F,
                            lower = fitness_ind$Plant_Label,
                            webID = fitness_ind$Plant_Simple,
                            freq = fitness_ind$weight)
  
list_incid_matrix_19 <- frame2webs(testdata_19,type.out="list")
  
for (j in 1:length(list_incid_matrix_19)){
    
    plant_i <- j
    print(names(list_incid_matrix_19)[plant_i])
    
    incid_matrix_i <- list_incid_matrix_19[[plant_i]] 
    
    layer_csv <- paste("Multilayer/Caracoles_layer_",names(list_incid_matrix_19)[plant_i],".csv",sep="")
    write.csv(incid_matrix_i,layer_csv)
}
