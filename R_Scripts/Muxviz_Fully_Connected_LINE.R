
library(attempt)
library(igraph)
library(bipartite)
library(tidyverse)

#Access layers files
dir_ini <- getwd()

folder_base <- paste(dir_ini,"/Processed_data/Multilayer_Species/",sep="")

files_base <- list.files(folder_base)

setwd(folder_base)

for (Line_i in 1:3){

# Extract layer files for Line_i

list_files_field_level <- files_base[grepl(paste("Line_",Line_i,sep = ""), files_base)]

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
                           weight=c(1,1))
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

# # Replace the node names with node_ids
# Plot_edgelist_complete_ids <- 
#   Plot_edgelist_complete %>% 
#   left_join(physical_nodes, by=c('node_from' = 'species')) %>%  # Join for pollinators
#   left_join(physical_nodes, by=c('node_to' = 'species')) %>%  # Join for plants
#   select(-node_from, -node_to) %>% 
#   select(layer_from, node_from=node_id.x, layer_to, node_to=node_id.y, weight) %>% 
#   left_join(layer_metadata, by=c('layer_from' = 'layer_name')) %>%  # Join for plants
#   left_join(layer_metadata, by=c('layer_to' = 'layer_name')) %>%  # Join for plants
#   select(-layer_from, -layer_to) %>% 
#   select(layer_from=layer_id.x, node_from, layer_to=layer_id.y, node_to, weight)

#########################################################################################
###########################################
#Saving MuxViz Files for posterior analysis
###########################################

folder_muxviz_root <- paste(dir_ini,"/Processed_data/Muxviz_Fully_Connected/",sep="")
newfolder <- paste("Line_",Line_i,"/", sep="")
folder_muxviz <- paste0(folder_muxviz_root, newfolder)
dir.create(folder_muxviz)

general_multilayer_layout <- physical_nodes %>% rename(nodeID=node_id,nodeLabel=species)%>%
  select(nodeID,nodeLabel) 

mutate(general_multilayer_layout,
       nodeLabel=str_replace(general_multilayer_layout$nodeLabel," ", "_"))

write_delim(general_multilayer_layout,
            paste(folder_muxviz,"general_multilayer_layout_Line",Line_i,".txt",sep=""),
            delim = " ")

general_multilayer_layers <- layer_metadata %>% rename(layerID=layer_id,layerLabel=layer_name)

write_delim(general_multilayer_layers,
            paste(folder_muxviz,"general_multilayer_layers_Line",Line_i,".txt",sep=""),
            delim = " ")

general_multilayer <- S_edge_list_ID %>%
  select(node_from,layer_from,node_to,layer_to,weight) %>%
  arrange(node_from,layer_from,node_to,layer_to)

write_delim(general_multilayer,
            paste(folder_muxviz,"general_multilayer_Line",Line_i,".edges",sep=""),
            col_names=FALSE,
            delim = " ")


write_delim(as.data.frame(
  paste(paste(folder_muxviz,"general_multilayer_Line",Line_i,".edges",sep=""),
        paste(folder_muxviz,"general_multilayer_layers_Line",Line_i,".txt",sep=""),
        paste(folder_muxviz,"general_multilayer_layout_Line",Line_i,".txt",sep=""),sep=";")),
  paste(folder_muxviz,"general_multilayer_config_Line",Line_i,".txt",sep=""),
  col_names=FALSE,
  quote_escape = "none",
  delim = ";")

write.table(as.data.frame(
  paste(paste(folder_muxviz,"general_multilayer_Line",Line_i,".edges",sep=""),
        paste(folder_muxviz,"general_multilayer_layers_Line",Line_i,".txt",sep=""),
        paste(folder_muxviz,"general_multilayer_layout_Line",Line_i,".txt",sep=""),sep=";")),
  paste(folder_muxviz,"general_multilayer_config_Line",Line_i,".txt",sep=""),
  row.names=FALSE,col.names = FALSE,sep="", quote = FALSE)

#################################################################################
}

setwd(dir_ini)
