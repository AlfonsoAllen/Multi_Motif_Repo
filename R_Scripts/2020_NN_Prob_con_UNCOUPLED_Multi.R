
# Load relevant libraries

library(tidyverse)
library(igraph)
library(expm)
source("R_scripts/functions.R")



#Access layers files
dir_ini <- getwd()

#Load data on pollinator visits
pollination <- read_csv2("Raw_Data/final_Pollinators_2020.csv")
# Remove points from ID names
pollination$ID <- sub("\\.", "", pollination$ID)
pollination$ID_Simple <- sub("\\.", "", pollination$ID_Simple)

# filter tabanidae
pollination<- pollination %>% filter(ID != "Tabanidae")


pollination$Line <- NA

for (i in 1:nrow(pollination)){
  if(pollination$Plot[i] %in% c(1,2,3)){pollination$Line[i] <- 1}
  else if(pollination$Plot[i] %in% c(4,5,6)){pollination$Line[i] <- 2}
  else{pollination$Line[i] <- 3}
}


unique(pollination$ID_Simple[grep(" ",pollination$ID_Simple,ignore.case = T)])
#No labels with spaces -> Good!


for (Plot_i in 1:9){
  
  ##########################
  #ESTIMATE PHENOLOGY
  ##########################
  
  #Filter pollination data
  pollination_20_i <- pollination %>% filter(Year==2020,!is.na(Plant),Plant!="0",Subplot!="OUT",Plant!="Ground")
  
  pollination_20_i <- pollination_20_i %>% select(Day,Month,Year,Line,Plot,Subplot,Plant,ID_Simple,Visits) %>%
    mutate(date_raw=as.Date(paste(Day,Month,Year,sep="/"), "%d/%m/%Y"),
           Week=as.numeric(format(date_raw, "%V"))) %>%
    rename(ID=ID_Simple)
  
  
  ###########################
  # CREATE MULTILAYER FOR Plot_i
  ###########################
  
  folder_base <- paste(dir_ini,"/Processed_data/Multilayer_Species/",sep="")
  
  files_base <- list.files(folder_base)
  
  setwd(folder_base)
  
  # Extract layer files for Plot_i
  
  list_files_field_level <- files_base[grepl(paste("Plot_",Plot_i,sep = ""), files_base) &
                                         grepl("2020", files_base) ]
  
  # Extract edge_list for each layer
  for (i in 1:length(list_files_field_level)){
    
    # Extract the incidence matrix
    inc_matrix <- read.csv(list_files_field_level[i], header=T, row.names=1)
    
    # Create a graph for each layer
    g_i <- graph_from_incidence_matrix(inc_matrix, directed = FALSE, weighted = T)
    
    # Get the edge_list from the graph and add plant (layer) information
    plant <- strsplit(list_files_field_level[i],".csv")
    plant <- strsplit(plant[[1]][1],".layer_")
    plant <- plant[[1]][2]
    
    g_i_edge_list <- as_tibble(igraph::as_data_frame(g_i, 'edges')) %>% mutate(species=plant)
    
    if (i==1){
      plot_edge_list <- g_i_edge_list
    }
    else{
      plot_edge_list <- plot_edge_list %>% bind_rows(g_i_edge_list)
    }
  }
  
  # Extract multilayer info
  
  pollinators <- sort(unique(plot_edge_list$to)) 
  plants <- sort(unique(plot_edge_list$from))
  layer_plant <- sort(unique(plot_edge_list$species))
  intersect(pollinators, plants)
  A <- length(pollinators) # Number of pollinators
  P <- length(plants) # Number of plants
  S <- A + P
  
  # Create a table with node metadata
  physical_nodes <- tibble(node_id=1:S,
                           type=c(rep('plant',P),rep('pollinator',A)),
                           species=c(plants,pollinators))
  layer_metadata <- tibble(layer_id=1:length(layer_plant), layer_name=layer_plant)
  
  # Replace the node names with node_ids
  
  Plot_edgelist_complete <- tibble(layer_from=plot_edge_list$species,
                                   node_from=plot_edge_list$from,
                                   layer_to=plot_edge_list$species,
                                   node_to=plot_edge_list$to,
                                   weight=plot_edge_list$weight)
  
  ##########
  plant_strength <- Plot_edgelist_complete %>% group_by(layer_from,node_from) %>% 
    count(wt = weight) %>% rename(strength = n)
  
  pollinator_strength <- Plot_edgelist_complete %>% group_by(layer_from,node_to) %>% 
    count(wt = weight) %>% rename(strength = n)
  ##########
  
  #Create the scaled directed list (previous list was meant to be undirected)
  
  #From plant to pollinator
  
  S_Links_Plant_Poll <- Plot_edgelist_complete %>% left_join(plant_strength,
                                                             by=c("layer_from","node_from")) %>%
    mutate(weight=weight/strength) %>% select(-strength)
  
  S_Links_Poll_Plant <- Plot_edgelist_complete %>% left_join(pollinator_strength,
                                                             by=c("layer_from","node_to")) %>%
    
    mutate(weight=weight/strength) %>% select(-strength) %>%
    rename(node_from=node_to,node_to=node_from)
  
  
  S_edge_list <- bind_rows(S_Links_Plant_Poll,S_Links_Poll_Plant)
  
  S_edge_list_i <- S_edge_list %>% mutate(Plot=Plot_i)
  
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
    if(nchar(sp_i)==2){
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
  
  if (Plot_i==1){stationary_prob_final <- stationary_probabilities_consp_heter}else{stationary_prob_final <- bind_rows(stationary_prob_final,stationary_probabilities_consp_heter)}
  
  #############################
  
  setwd(dir_ini)

  
}

write_csv(stationary_prob_final,"Processed_data/2020_NN_plant_stationary_prob_results_UNCOUPLED.csv")
