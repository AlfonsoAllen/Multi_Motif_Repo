# Load relevant libraries
library(infomapecology)
# Infomap installation guide: https://github.com/Ecological-Complexity-Lab/infomap_ecology_package
library(attempt)
library(igraph)
library(bipartite)
library(tidyverse)
library(magrittr)
library(ggalluvial)
source("R_scripts/functions.R")
source("R_scripts/run_infomap_monolayer2.R")


#Access layers files
dir_ini <- getwd()

#Load data on pollinator visits
data_simulation_raw <- read_csv("Processed_data/Data_simulation/data_increasing_mix_V2.csv")

# Aggregate visits per plant individual, week and floral visitor
example <- data_simulation_raw %>% group_by(Plot,Subplot,Plant,ID,Week) %>% count(wt = Visits) %>%
  rename(Visits=n)

list_sites <- example$Plot %>% unique() # List of networks IDs

modularity_info_final <- NULL # Variable to save probability metrics
  
  
for (Plot_i in list_sites){
  
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
  
  #######################################
  # CREATE NETWORK OF NETWORKS LIST
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
  
  nn_nodes <- NN_edge_list_final$from %>% unique()
  
  NN_node_ID <- tibble(node_id=as.numeric(1:length(nn_nodes)),species=nn_nodes)
  
  NN_edge_list_final_ID <- NN_edge_list_final %>% rename(species=from) %>%
    left_join(NN_node_ID,by="species") %>% select(-species) %>% rename(from=node_id,species=to) %>%
    left_join(NN_node_ID,by="species") %>% 
    select(-species) %>% rename(to=node_id) %>% select(from,to,weight)
  
  #######################################
  #Running Infomap
  #######################################
  
  #Running Infomap
  #Setting folder with infomap.exe
  folder_info <- paste(dir_ini,"/R_Scripts",sep="")
  
  # Check Infomap is running
  setwd(folder_info)
  check_infomap() # Make sure file can be run correctly. Should return TRUE
  
  
  #############################
  # NN Modules
  
  NN_node_ID2 <- NN_node_ID %>% rename(node_name=species)
  
  NN_edge_list_final_ID$from <- as.character(NN_edge_list_final_ID$from)
  NN_edge_list_final_ID$to <- as.character(NN_edge_list_final_ID$to)
  
  Plot_NN <- create_monolayer_object(x = NN_edge_list_final_ID, directed = T, bipartite = F, node_metadata = NN_node_ID2)
  
  Plot_NN2 <- Plot_NN
  # Plot_NN2$nodes$node_name <- as.numeric(Plot_NN2$nodes$node_name)
  
  # Run Infomap
  modules_relax_rate_NN <- run_infomap_monolayer(Plot_NN2, flow_model = 'directed', silent=T,trials=1000, two_level=T, seed=200952)
  
  # Extract information
  module_info_i <- tibble(Plot=Plot_i,
                        L=modules_relax_rate_NN$L,
                        m=modules_relax_rate_NN$m)
  
  
  # Update final variables with Plot_i results
  
  modularity_info_final <- bind_rows(modularity_info_final,
                                     module_info_i)
  
  
  setwd(dir_ini)
  
  # Commented for security reasons
  write_csv(modularity_info_final,"Processed_data/Data_simulation/data_changing_increasing_mix_MODULARITY.csv")

}

