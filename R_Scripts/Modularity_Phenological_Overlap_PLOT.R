# Load relevant libraries
library(infomapecology)
# Infomap installation guide: https://github.com/Ecological-Complexity-Lab/infomap_ecology_package
library(attempt)
library(igraph)
library(bipartite)
library(tidyverse)
library(magrittr)
library(ggalluvial)

#Access layers files
dir_ini <- getwd()


#Load data on pollinator visits
pollination <- read_csv("Raw_data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")

for (Plot_i in 1:9){

##########################
#ESTIMATE PHENOLOGY
##########################

#Filter pollination data
pollination_19_i <- pollination %>% filter(Year==2019,Subplot!="OUT",Plot==Plot_i,!is.na(ID))

pollination_19_i <- pollination_19_i %>% select(Day,Month,Year,Plot,Subplot,Plant_Simple,ID,Visits) %>%
  mutate(date_raw=as.Date(paste(Day,Month,Year,sep="/"), "%d/%m/%Y"),
         Week=as.numeric(format(date_raw, "%V")))

#------
# FUNCTION: plant_pheno_overlap 
# Returns the number of weeks of phenological overlap between plant1 and plant2 over
# divided by the total duration of plant 1 phenology

plant_pheno_overlap <- function(plant1,plant2,pollination_19_i) {
  
  plants <- sort(unique(pollination_19_i$Plant_Simple))
  if(sum(c(plant1,plant2) %in% plants)<2){
    print(paste("Error: At least one plant does not belong to plot ",
                Plot_i," experimental phenology",sep=""))
    }else{
    pollination_plant1 <- pollination_19_i %>% filter(Plant_Simple==plant1)
    #Weeks in which plant 1 recieves visits
    plant_1_weeks <- min(pollination_plant1$Week):max(pollination_plant1$Week)
    
    pollination_plant2 <- pollination_19_i %>% filter(Plant_Simple==plant2)
    #Weeks in which plant 2 recieves visits
    plant_2_weeks <- min(pollination_plant2$Week):max(pollination_plant2$Week)
    
    #number of weeks with phenological overlap between plant1 and plant2
    overlap <- sum(plant_1_weeks %in% plant_2_weeks)
    
    #total_number of weeks in the experimental phenology of plot i
    total_weeks_1 <-  length(plant_1_weeks)
    
    return(overlap/total_weeks_1)
    #return(overlap)
    
  }  
}



###########################
# CREATE MULTILAYER FOR Plot_i
###########################

folder_base <- paste(dir_ini,"/Processed_data/Multilayer_Species/",sep="")

files_base <- list.files(folder_base)

setwd(folder_base)

# Extract layer files for Plot_i

list_files_field_level <- files_base[grepl(paste("Plot_",Plot_i,sep = ""), files_base)]

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



pollinators <- sort(unique(plot_edge_list$to)) 
plants <- sort(unique(plot_edge_list$from))
layer_plant <- sort(unique(plot_edge_list$species))
intersect(pollinators, plants)
A <- length(pollinators) # Number of pollinators
P <- length(plants) # Number of plants
S <- A+P

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

Plot_edgelist_complete %>% filter(node_from == "B5 LEMA") # total strength should be 11 (check!)

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

###############
# To create the inter-links we rely on the previous Plot_edgelist_complete
# Here we can extract information on interlayer connections

for (i in 1:length(pollinators)){
  
  polinator_edges <- Plot_edgelist_complete %>% filter(node_to==pollinators[i])
  polinator_layers <- unique(polinator_edges$layer_to)
  #print (polinator_layers)
  if (length(polinator_layers)>1){
    combination_layers <- t(combn(polinator_layers, 2))
    for (j in 1:nrow(combination_layers)){
      #For undirected networks
      # interlink_i<- tibble(layer_from=combination_layers[j,1],
      #                      node_from=pollinators[i],
      #                      layer_to=combination_layers[j,2],
      #                      node_to=pollinators[i],
      #                      weight=1)
      
      #For directed networks
      interlink_i<- tibble(layer_from=c(combination_layers[j,1],combination_layers[j,2]),
                           node_from=c(pollinators[i],pollinators[i]),
                           layer_to=c(combination_layers[j,2],combination_layers[j,1]),
                           node_to=c(pollinators[i],pollinators[i]),
                           weight=c(plant_pheno_overlap(combination_layers[j,1],
                                                        combination_layers[j,2],
                                                        pollination_19_i),
                                    plant_pheno_overlap(combination_layers[j,2],
                                                        combination_layers[j,1],
                                                        pollination_19_i)))
      
      #For undirected
      #Plot_edgelist_complete <- bind_rows(Plot_edgelist_complete,interlink_i)
      
      #For directed
      S_edge_list <- bind_rows(S_edge_list,interlink_i)
    }
  }
}




# Replace the node names with node_ids
S_edge_list_ID <- 
  S_edge_list %>% 
  left_join(physical_nodes, by=c('node_from' = 'species')) %>%  # Join for pollinators
  left_join(physical_nodes, by=c('node_to' = 'species')) %>%  # Join for plants
  select(-node_from, -node_to) %>% 
  select(layer_from, node_from=node_id.x, layer_to, node_to=node_id.y, weight) %>% 
  left_join(layer_metadata, by=c('layer_from' = 'layer_name')) %>%  # Join for plants
  left_join(layer_metadata, by=c('layer_to' = 'layer_name')) %>%  # Join for plants
  select(-layer_from, -layer_to) %>% 
  select(layer_from=layer_id.x, node_from, layer_to=layer_id.y, node_to, weight)


#######################################
#Running Infomap
#######################################

#Running Infomap
#Setting folder with infomap.exe
folder_info <- paste(dir_ini,"/R_Scripts",sep="")

# Check Infomap is running
setwd(folder_info)
check_infomap() # Make sure file can be run correctly. Should return TRUE

# Prepare data
#Plot_multilayer <- create_multilayer_object(extended = Plot_edgelist_complete_ids, nodes = physical_nodes, intra_output_extended = T, inter_output_extended = T)
Plot_multilayer <- create_multilayer_object(extended = S_edge_list_ID, nodes = physical_nodes, intra_output_extended = T, inter_output_extended = T)


# Run Infomap
#modules_relax_rate <- run_infomap_multilayer(Plot_multilayer, relax = F, silent = T, flow_model = 'undirected', trials = 250, seed = 497294, temporal_network = F)
modules_relax_rate <- run_infomap_multilayer(Plot_multilayer, relax = F, silent = T, flow_model = 'directed', trials = 1000, seed = 497294, temporal_network = F)

# Extract information
plot_modules_i <- modules_relax_rate$modules %>% left_join(layer_metadata,by="layer_id")

plot_modules_i$Plot <- Plot_i

setwd(dir_ini)

write_csv(plot_modules_i,paste0("Processed_data/Modularity_Pheno_Overlap/Modularity_Plot",
                                Plot_i,".csv")  )

}

