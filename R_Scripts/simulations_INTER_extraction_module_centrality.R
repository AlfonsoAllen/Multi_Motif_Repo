# Load relevant libraries
library(infomapecology)
# Infomap installation guide:
# https://github.com/Ecological-Complexity-Lab/infomap_ecology_package
library(attempt)
library(igraph)
library(tidyverse)
source("R_scripts/functions.R")
source("R_scripts/run_infomap_monolayer2.R") # function to parse infomap's tree files


dir_ini <- getwd() # Register working directory

data_simulation_raw <- read_csv("Processed_data/Data_simulation/data_changing_interlinks_V2.csv")

# Aggregate visits per plant individual, week and floral visitor
example <- data_simulation_raw %>% group_by(Plot,Subplot,Plant,ID,Week) %>% count(wt = Visits) %>%
  rename(Visits=n)
  
list_sites <- example$Plot %>% unique() # List of networks IDs

centrality_final <- NULL # Variable to save centrality metrics
plot_modules_NN_final <- NULL # Variable to save module partitions

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
  
  # Extract modules for Plot_i from the network of networks representation with infomap
  
  plot_modules_NN_i <- extract_modules_NN_infomapec(S_edge_list)
  
  # Extract centrality metrics for Plot_i
  
  centrality_i <- centrality_metrics_NN(S_edge_list)
    
  # Update final variables with Plot_i results
  
  plot_modules_NN_final <- bind_rows(plot_modules_NN_final,plot_modules_NN_i)
  centrality_final <- bind_rows(centrality_final,centrality_i)
  
  setwd(dir_ini)
  
  # Commented for security reasons
  write_csv(plot_modules_NN_final,"Processed_data/Data_simulation/data_changing_interlinks_MODULES_V2.csv")
  write_csv(centrality_final,"Processed_data/Data_simulation/data_changing_interlinks_CENTRALITY_V2.csv")

}

sessionInfo()
