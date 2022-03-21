# load libraries
library(tidyverse)
library(bipartite)
library(igraph)

source("R_Scripts/functions.R")

####################################################################
# Loadind example dataset

data_simulation_raw <- read_csv("Processed_data/Data_simulation/data_changing_individuals_V2.csv")

list_plots <- data_simulation_raw$Plot %>% unique()

data_changing_individuals_MOTIFS <- NULL

for(i.plot in 1:length(list_plots)){
  
  data_simulation <- data_simulation_raw %>% filter(Plot %in% list_plots[i.plot])
  
  aggregate_total <- NULL
  
  for (week_i in sort(unique(data_simulation$Week))){
    
    print("WEEK")
    print(week_i)
    
    example_week_i <- data_simulation %>% filter(Week==week_i)
    
    # Aggregate visits by week 
    
    example_week_i <- example_week_i %>% group_by(Plot,Subplot,Plant,ID) %>%
      count(wt=Visits) %>% rename(Visits_tot = n)
    
    list_visitors_week_i <- example_week_i$ID %>% unique()
    
    if((nrow(example_week_i)>length(list_visitors_week_i))&
       (length(list_visitors_week_i)>1)){
      
      example_week_i$Subplot_Plant_Label <- paste(example_week_i$Subplot,example_week_i$Plant,sep = " ")
      
      aggregate_week_i <- example_week_i %>% ungroup() %>% 
        select(Plot,ID,Subplot_Plant_Label,Visits_tot)
      
      
      aggregate_week_i <- homo_hete_motifs(aggregate_week_i)
      aggregate_week_i <- aggregate_week_i %>% mutate(Week=week_i)
      
    }else if(length(list_visitors_week_i)==1){
      
      example_week_i$Subplot_Plant_Label <- paste(example_week_i$Subplot,
                                                  example_week_i$Plant,sep = " ")
      
      example_week_i$Week <- week_i
      example_week_i$homo_motif <- 0
      example_week_i$hete_motif <- 0
      
      for(i.plant_ind in 1:nrow(example_week_i)){
        
        number_homo_ind <- example_week_i %>% filter(Plant==example_week_i$Plant[i.plant_ind]) %>%
          nrow()
        
        number_hete_ind <- example_week_i %>% filter(Plant!=example_week_i$Plant[i.plant_ind]) %>%
          nrow()
        
        example_week_i$homo_motif[i.plant_ind] <- number_homo_ind-1
        example_week_i$hete_motif[i.plant_ind] <- number_hete_ind
        
      }
      
      aggregate_week_i <- example_week_i[,c("Plot","ID","Subplot_Plant_Label","Visits_tot",
                                            "homo_motif","hete_motif","Week")]
      
    }else{
      example_week_i$Subplot_Plant_Label <- paste(example_week_i$Subplot,
                                                  example_week_i$Plant,sep = " ")
      example_week_i$Week <- week_i
      example_week_i$homo_motif <- 0
      example_week_i$hete_motif <- 0
      aggregate_week_i <- example_week_i[,c("Plot","ID","Subplot_Plant_Label","Visits_tot",
                                            "homo_motif","hete_motif","Week")]
    }
    
    
    
    aggregate_total <-  bind_rows(aggregate_total, aggregate_week_i) 
    
  }
  
  aggregate_total_plot <- aggregate_total %>% 
    separate(Subplot_Plant_Label,c("Subplot","Plant"), " ")
  
  data_changing_individuals_MOTIFS <- bind_rows(data_changing_individuals_MOTIFS,
                                               aggregate_total_plot)
  
  # Commented for security reasons
  # write_csv(data_changing_individuals_MOTIFS,
  #           "Processed_data/Data_simulation/data_changing_individuals_MOTIFS_V2.csv")
  
}



