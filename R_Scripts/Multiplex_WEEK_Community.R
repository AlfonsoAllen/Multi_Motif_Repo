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

folder_base <- paste(dir_ini,"/Processed_data/Multilayer_WEEK",sep="")

files_base <- list.files(folder_base)

setwd(folder_base)

# Select the plot number to build up the multilayer
Plot_i <- 1

# Extract layer files for Plot_i

list_files_field_level <- files_base[grepl(paste("Plot_",Plot_i,sep = ""), files_base)]

# Extract edge_list for each layer
for (i in 1:length(list_files_field_level)){

  # Extract the incidence matrix
  inc_matrix <- read.csv(list_files_field_level[i], header=T, row.names=1)
  
  # Create a graph for each layer
  g_i <- graph_from_incidence_matrix(inc_matrix, directed = FALSE, weighted = T)
  
  # Get the edge_list from the graph and add week (layer) information
  week <- strsplit(list_files_field_level[i],".csv")
  week <- strsplit(week[[1]][1],".WEEK_")
  week <- as.numeric(week[[1]][2])
  
  g_i_edge_list <- as_tibble(igraph::as_data_frame(g_i, 'edges')) %>% mutate(Week=week)

  if (i==1){
    plot_edge_list <- g_i_edge_list
  }
  else{
    plot_edge_list <- plot_edge_list %>% bind_rows(g_i_edge_list)
  }
}



pollinators <- sort(unique(plot_edge_list$to)) 
plants <- sort(unique(plot_edge_list$from))
layer_week <- sort(unique(plot_edge_list$Week))
intersect(pollinators, plants)
A <- length(pollinators) # Number of pollinators
P <- length(plants) # Number of plants
S <- A+P

# Create a table with node metadata
physical_nodes <- tibble(node_id=1:S,
                         type=c(rep('plant',P),rep('pollinator',A)),
                         species=c(plants,pollinators))
layer_metadata <- tibble(layer_id=1:length(layer_week), layer_name=layer_week)

# Replace the node names with node_ids

Plot_edgelist_complete <- tibble(layer_from=plot_edge_list$Week,
                                        node_from=plot_edge_list$from,
                                        layer_to=plot_edge_list$Week,
                                        node_to=plot_edge_list$to,
                                        weight=plot_edge_list$weight)


# Replace the node names with node_ids
Plot_edgelist_complete_ids <- 
  Plot_edgelist_complete %>% 
  left_join(physical_nodes, by=c('node_from' = 'species')) %>%  # Join for pollinators
  left_join(physical_nodes, by=c('node_to' = 'species')) %>%  # Join for plants
  select(-node_from, -node_to) %>% 
  select(layer_from, node_from=node_id.x, layer_to, node_to=node_id.y, weight) %>% 
  left_join(layer_metadata, by=c('layer_from' = 'layer_name')) %>%  # Join for plants
  left_join(layer_metadata, by=c('layer_to' = 'layer_name')) %>%  # Join for plants
  select(-layer_from, -layer_to) %>% 
  select(layer_from=layer_id.x, node_from, layer_to=layer_id.y, node_to, weight)


# Write to infomap
# Write *Intra format
write_lines('*Intra', 'Plot_multilayer.txt')
Plot_edgelist_complete_ids %>% 
  filter(layer_from==layer_to) %>% 
  select(layer_from, node_from, node_to, weight) %>% 
  write_delim('Plot_multilayer.txt', append=T, col_names = F) #Append=T do not overwrite the
#previous information-->'*Intra'



#Running Infomap
#Setting folder with infomap.exe
folder_info <- paste(dir_ini,"/R_Scripts",sep="")

# Check Infomap is running
setwd(folder_info)
check_infomap() # Make sure file can be run correctly. Should return TRUE

# Prepare data
Plot_multilayer <- create_multilayer_object(extended = Plot_edgelist_complete_ids, nodes = physical_nodes, intra_output_extended = F, inter_output_extended = F)
# Ignore interlayer edges
Plot_multilayer$inter <- NULL

# Run Infomap
modules_relax_rate <- run_infomap_multilayer(Plot_multilayer, relax = T, silent = T, trials = 50, seed = 497294, multilayer_relax_rate = 0.15, multilayer_relax_limit_up = 1, multilayer_relax_limit_down = 0, temporal_network = T)

# Extract information
plot_modules <- modules_relax_rate$modules %>% left_join(layer_metadata,by="layer_id")


# Number of species in a combination of layer-module
plot_modules %>% group_by(module, layer_id) %>% count() %>% 
  ggplot()+
  geom_point(aes(layer_id, module, size=n, color=as.factor(module)))+
  theme_bw()+
  theme(legend.position = 'none', axis.text.x = element_text(angle = -90))


# Alluvial plots
is_lodes_form(plot_modules,
              key = layer_id, value = module, id = node_id, silent = TRUE) # Check if conforms to alluvial type data
ggplot(plot_modules,
       aes(x=layer_name, stratum=as.factor(module), 
           alluvium=node_id, label=as.factor(module), fill=as.factor(module))) + 
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  geom_text(stat = "stratum", size = 3) +
  labs(x='Week', y='Number of species')+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color='black'), axis.text.x = element_text(angle = -90))

# Alluvial plots for plants

plot_modules_plants <- plot_modules %>% filter(type=="plant")

ggplot(plot_modules_plants,
       aes(x=layer_name, stratum=as.factor(module), 
           alluvium=node_id, label=as.factor(module), fill=as.factor(module))) + 
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  geom_text(stat = "stratum", size = 3) +
  labs(x='Week', y='Number of (plant species,subplots)')+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color='black'), axis.text.x = element_text(angle = -90))

# Alluvial plots for pollinators

plot_modules_pollinators <- plot_modules %>% filter(type!="plant")

ggplot(plot_modules_pollinators,
       aes(x=layer_name, stratum=as.factor(module), 
           alluvium=node_id, label=as.factor(module), fill=as.factor(module))) + 
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  geom_text(stat = "stratum", size = 3) +
  labs(x='Week', y='Number of pollinator species')+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color='black'), axis.text.x = element_text(angle = -90))



# Flexibility of species 
plot_modules %>%
  left_join(physical_nodes) %>% 
  group_by(type, node_id) %>% summarise(n_modules=n_distinct(module)) %>% 
  ggplot()+geom_histogram(aes(x=n_modules, fill=type), position = 'dodge')+
  theme_bw()+scale_fill_manual(values = c('dark green','orange'))+labs(x='Number of modules')
