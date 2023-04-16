
# We supposed that efficiency is equal to 1

# Load relevant libraries
library(attempt)
library(igraph)
library(tidyverse)
library(expm)
source("R_scripts/functions.R")


dir_ini <- getwd() # Register working directory

data_simulation_raw <- read_csv("Processed_data/Data_simulation/data_changing_interlinks_V2.csv")

# Aggregate visits per plant individual, week and floral visitor
example <- data_simulation_raw %>% group_by(Plot,Subplot,Plant,ID,Week) %>% count(wt = Visits) %>%
  rename(Visits=n)
  
list_sites <- example$Plot %>% unique() # List of networks IDs

stationary_probabilities_consp_heter_final <- NULL # Variable to save probability metrics

for (Plot_i in list_sites){
  
  print(Plot_i)
  
  ########################################################
  # CREATE (NODE-COLORED) MULTILAYER (DIRECTED) EDGE LIST
  ########################################################
  
  plot_edge_list <- example %>% ungroup() %>% filter(Plot==Plot_i) %>%
    group_by(Plot,Subplot,Plant,ID) %>% count(wt=Visits) %>% ungroup() %>%
    mutate(from = paste0(Subplot," ",Plant)) %>% 
    rename(to = ID, weight = n, species = Plant) %>%
    dplyr::select(from, to, weight, species)
  
  
  # Extract multilayer info
  
  pollinators <- sort(unique(plot_edge_list$to)) 
  plants <- sort(unique(plot_edge_list$from))
  layer_plant <- sort(unique(plot_edge_list$species))
  
  # Sanity check
  intersect(pollinators, plants)
  
  A <- length(pollinators) # Number of pollinators
  P <- length(plants) # Number of plants
  S <- A + P
  
  # Create a table with node metadata
  physical_nodes <- tibble(node_id=1:S,
                           type=c(rep('plant',P),rep('pollinator',A)),
                           species=c(plants,pollinators))
  layer_metadata <- tibble(layer_id=1:length(layer_plant), layer_name=layer_plant)
  
  
  # Create a list that contains directed and weighted intra- and inter-links
  
  S_edge_list <- create_weighted_link_list(plot_edge_list)
  
  # Extract conspecific and heterospecif probabilities for Plot_i
  
  #######################################
  # CREATE RANDOM-WALK TRANSITION MATRIX
  #######################################
  
  
  NN_edge_list <- S_edge_list %>% filter(weight>0)
  
  for (i in 1:nrow(NN_edge_list)){
    
    if(NN_edge_list$node_from[i] %in% pollinators){
      NN_edge_list$node_from[i] <- paste0(NN_edge_list$node_from[i]," ",NN_edge_list$layer_from[i])
    }
    
    if(NN_edge_list$node_to[i] %in% pollinators){
      NN_edge_list$node_to[i] <- paste0(NN_edge_list$node_to[i]," ",NN_edge_list$layer_to[i])
    }
    
  }
  
  NN_edge_list_final <- NN_edge_list %>% select(node_from,node_to,weight) %>%
    rename(from = node_from, to = node_to)
  
  NN_nodes_data <- tibble(name = NN_edge_list_final$from %>% unique(),
                          type = NA,
                          layer = NA)
  
  for (i in 1:nrow(NN_nodes_data)) {
    
    split_name_i <- strsplit(NN_nodes_data$name[i], " ")
    sp_i <- split_name_i[[1]][1]
    if(grepl("ind",sp_i,ignore.case = T)==T){
      NN_nodes_data$type[i] <- "plant"
    }else{
      NN_nodes_data$type[i] <- "pollinator"
    }
    NN_nodes_data$layer[i] <- split_name_i[[1]][2]
  }
  
  NN_nodes <- nrow(NN_nodes_data)
  
  supra_adj_matrix <- matrix(rep(0,NN_nodes*NN_nodes),nrow = NN_nodes, ncol=NN_nodes) 
  
  colnames(supra_adj_matrix) <- NN_nodes_data$name
  rownames(supra_adj_matrix) <- NN_nodes_data$name
  
  for (i in 1:nrow(NN_nodes_data)) {
    
    for(j in 1:nrow(NN_nodes_data)){
      
      sp_from <- NN_nodes_data$name[i]
      sp_to <- NN_nodes_data$name[j]
      position_from <- which(NN_nodes_data$name == NN_nodes_data$name[i])
      position_to <- which(NN_nodes_data$name == NN_nodes_data$name[j])
      
      aux_NN_link_info <- NN_edge_list_final %>%
        filter(from == sp_from, to == sp_to)
      
      if(nrow(aux_NN_link_info)>0){
        
        supra_adj_matrix[position_from,position_to] <- as.numeric(aux_NN_link_info$weight)
        
      }
      
      
    }
    
  }
  
  transition_matrix <- supra_adj_matrix
  
  for (i in 1:nrow(supra_adj_matrix)) {
    
    sum_supra_adj_row <- sum(supra_adj_matrix[i,])
    
    transition_matrix[i,] <- supra_adj_matrix[i,]/sum_supra_adj_row
    
  }
  
  transpose_transition_matrix <- t(transition_matrix)
  
  stationary_transition_matrix <- 0.5*(transpose_transition_matrix %^% 1e9 + transpose_transition_matrix %^% (1e9+1))
  
  # Delete insect columns
  
  position_pollinator_nodes <- which(NN_nodes_data$type == "pollinator")
  stationary_transition_matrix_plant_aux <- 
    stationary_transition_matrix[,-position_pollinator_nodes]
  
  stationary_transition_matrix_plant <- 
    stationary_transition_matrix_plant_aux[-position_pollinator_nodes,]
  
  stationary_probabilities <- tibble(name=colnames(stationary_transition_matrix_plant),
                                     Plot = Plot_i,
                                     consp_prob_UNCOUPLED = NA,
                                     heter_prob_UNCOUPLED = NA)
  
  stationary_probabilities_consp_heter <- stationary_probabilities %>%
    left_join(NN_nodes_data, by = "name")
  
  for (i in 1:nrow(stationary_probabilities_consp_heter)) {
    
    layer_i <- stationary_probabilities_consp_heter$layer[i]
    consp_rows <- which(stationary_probabilities_consp_heter$layer == layer_i)
    
    prob_find_pollen_in_i_from_a_plant <- sum(stationary_transition_matrix_plant[i,])
    
    stationary_probabilities_consp_heter$consp_prob_UNCOUPLED[i] <- 
      sum(stationary_transition_matrix_plant[i,consp_rows])
    
    stationary_probabilities_consp_heter$heter_prob_UNCOUPLED[i] <- prob_find_pollen_in_i_from_a_plant -
      stationary_probabilities_consp_heter$consp_prob_UNCOUPLED[i]
    
  }
  
  # Since pollinators do not produce pollen grains and assuming that each patch
  # produces pollen with equal probability, then:
  
  stationary_probabilities_consp_heter$number_plant_nodes_with_visits <- nrow(stationary_probabilities_consp_heter)
  
  stationary_probabilities_consp_heter$consp_prob_UNCOUPLED <- 
    stationary_probabilities_consp_heter$consp_prob_UNCOUPLED/nrow(stationary_probabilities_consp_heter)
  stationary_probabilities_consp_heter$heter_prob_UNCOUPLED <- 
    stationary_probabilities_consp_heter$heter_prob_UNCOUPLED/nrow(stationary_probabilities_consp_heter)
    
  # Update final variables with Plot_i results
  
  stationary_probabilities_consp_heter_final <- bind_rows(stationary_probabilities_consp_heter_final,
                                                          stationary_probabilities_consp_heter)

  
  setwd(dir_ini)
  
  # Commented for security reasons
  write_csv(stationary_probabilities_consp_heter_final,"Processed_data/Data_simulation/data_changing_interlinks_PROBABILITIES_UNCOUPLED.csv")
  

}

sessionInfo()
