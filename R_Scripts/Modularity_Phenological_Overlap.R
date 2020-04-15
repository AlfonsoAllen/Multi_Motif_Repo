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


# Select the plot number to build up the multilayer
Plot_i <- 8

##########################
#ESTIMATE PHENOLOGY
##########################

#Load data on pollinator visits
pollination <- read_csv("Raw_data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")

#Filter data
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

folder_muxviz_root <- paste(dir_ini,"/Processed_data/Muxviz_Pheno_Overlap/",sep="")
newfolder <- paste("Plot_",Plot_i,"/", sep="")
folder_muxviz <- paste0(folder_muxviz_root, newfolder)
dir.create(folder_muxviz)

general_multilayer_layout <- physical_nodes %>% rename(nodeID=node_id,nodeLabel=species)%>%
  select(nodeID,nodeLabel) 

mutate(general_multilayer_layout,
       nodeLabel=str_replace(general_multilayer_layout$nodeLabel," ", "_"))

write_delim(general_multilayer_layout,
            paste(folder_muxviz,"general_multilayer_layout_Plot",Plot_i,".txt",sep=""),
            delim = " ")

general_multilayer_layers <- layer_metadata %>% rename(layerID=layer_id,layerLabel=layer_name)

write_delim(general_multilayer_layers,
            paste(folder_muxviz,"general_multilayer_layers_Plot",Plot_i,".txt",sep=""),
            delim = " ")

general_multilayer <- S_edge_list_ID %>%
  select(node_from,layer_from,node_to,layer_to,weight) %>%
  arrange(node_from,layer_from,node_to,layer_to)

write_delim(general_multilayer,
            paste(folder_muxviz,"general_multilayer_Plot",Plot_i,".edges",sep=""),
            col_names=FALSE,
            delim = " ")


write_delim(as.data.frame(
  paste(paste(folder_muxviz,"general_multilayer_Plot",Plot_i,".edges",sep=""),
        paste(folder_muxviz,"general_multilayer_layers_Plot",Plot_i,".txt",sep=""),
        paste(folder_muxviz,"general_multilayer_layout_Plot",Plot_i,".txt",sep=""),sep=";")),
  paste(folder_muxviz,"general_multilayer_config_Plot",Plot_i,".txt",sep=""),
  col_names=FALSE,
  quote_escape = "none",
  delim = ";")

write.table(as.data.frame(
  paste(paste(folder_muxviz,"general_multilayer_Plot",Plot_i,".edges",sep=""),
        paste(folder_muxviz,"general_multilayer_layers_Plot",Plot_i,".txt",sep=""),
        paste(folder_muxviz,"general_multilayer_layout_Plot",Plot_i,".txt",sep=""),sep=";")),
  paste(folder_muxviz,"general_multilayer_config_Plot",Plot_i,".txt",sep=""),
  row.names=FALSE,col.names = FALSE,sep="", quote = FALSE)

#################################################################################
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
plot_modules <- modules_relax_rate$modules %>% left_join(layer_metadata,by="layer_id")


# Number of species in a combination of layer-module
plot_modules %>% group_by(module, layer_name) %>% count() %>% 
  ggplot()+
  geom_point(aes(layer_name, module, size=n, color=as.factor(module)))+
  theme_bw()+
  theme(legend.position = 'none', axis.text.x = element_text(angle = 0))

plot_modules %>% filter(module==6)


# Alluvial plots
is_lodes_form(plot_modules,
              key = layer_id, value = module, id = node_id, silent = TRUE) # Check if conforms to alluvial type data
ggplot(plot_modules,
       aes(x=layer_id, stratum=as.factor(module), 
           alluvium=node_id, label=as.factor(module), fill=as.factor(module))) + 
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  geom_text(stat = "stratum", size = 3) +
  labs(x='Plant Layer', y='Number of species',
       title = paste("Plot ",Plot_i,': Modules across plant layers (Fully connected ensemble)',sep=""))+
  scale_x_discrete(limits=as.character(layer_metadata$layer_name))+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color='black'), axis.text.x = element_text(angle = 0))

# Alluvial plots for plants #No fluexes between modules should appear: Check 

plot_modules_plants <- plot_modules %>% filter(type=="plant")

ggplot(plot_modules_plants,
       aes(x=layer_id, stratum=as.factor(module), 
           alluvium=node_id, label=as.factor(module), fill=as.factor(module))) + 
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  geom_text(stat = "stratum", size = 3) +
  labs(x='Plant Layer', y='Number of (plant species,subplots)')+
  scale_x_discrete(limits=as.character(layer_metadata$layer_name))+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color='black'), axis.text.x = element_text(angle = 0))

# Alluvial plots for pollinators

plot_modules_pollinators <- plot_modules %>% filter(type!="plant")

ggplot(plot_modules_pollinators,
       aes(x=layer_id, stratum=as.factor(module), 
           alluvium=node_id, label=as.factor(module), fill=as.factor(module))) + 
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  geom_text(stat = "stratum", size = 3) +
  labs(x='Plant Layer', y='Number of pollinator species')+
  scale_x_discrete(limits=as.character(layer_metadata$layer_name))+
  theme_bw()+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color='black'), axis.text.x = element_text(angle = 0))

plot_modules_pollinators <- plot_modules_pollinators %>%arrange(module)

# Flexibility of species 
plot_modules %>%
  left_join(physical_nodes) %>% 
  group_by(type, node_id) %>% summarise(n_modules=n_distinct(module)) %>% 
  ggplot()+geom_histogram(aes(x=n_modules, fill=type), position = 'dodge')+
  theme_bw()+scale_fill_manual(values = c('dark green','orange'))+labs(x='Number of modules')


setwd(dir_ini)

  
######
library(circlize)
# https://jokergoo.github.io/circlize_book/book/the-chorddiagram-function.html

#colormaps
library(viridis)           # Load
library(RColorBrewer)


test_mod0 <- plot_modules %>% filter(type!="plant")%>%
  group_by(layer_name,module) %>% count() %>%
  select(layer_name,module,n) %>% rename(from =layer_name,to=module,value=n)


c <- length(unique(test_mod0$from))



grid.col0 <- setNames(c(magma(c),
                        rep("grey",length(unique(test_mod0$to)))),
                      union(test_mod0$from, test_mod0$to))

par(cex = .7, mar = c(0, 0, 0, 0))
chordDiagram(test_mod0, preAllocateTracks = 1,
             order=union(test_mod0$from, rev(test_mod0$to)),
             big.gap = 10,grid.col=grid.col0,
             symmetric=FALSE,
             link.arr.type = "big.arrow",
             #transparency of colours
             transparency = 0.2,
             scale = TRUE
)


circos.clear()

###########



test_mod1 <- plot_modules %>% filter(type!="plant")%>%
  group_by(species,module) %>% count() %>%
  select(species,module,n) %>% rename(from =species,to=module,value=n)

a <- 9
b <- 9
c <- length(unique(test_mod1$from))- a-b



grid.col1 <- setNames(c(brewer.pal(n = a, name = "Set1"),
                        brewer.pal(n = b, name = "Pastel1"),
                        magma(c),
                        rep("grey",length(unique(test_mod1$to)))),
                      union(test_mod1$from, test_mod1$to))

par(cex = .7, mar = c(0, 0, 0, 0))
chordDiagram(test_mod1,annotationTrack = "grid", preAllocateTracks = 1,
             order=union(test_mod1$from, rev(test_mod1$to)),
             big.gap = 10,grid.col=grid.col1,
             symmetric=FALSE,
             link.arr.type = "big.arrow",
             #transparency of colours
             transparency = 0.2,
             scale = TRUE
)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

circos.clear()

###########


test_mod2 <- plot_modules %>% filter(type!="plant")%>%
  group_by(species,layer_name) %>% count() %>%
  select(species,layer_name,n) %>% rename(from =species,to=layer_name,value=n)

a <- 9
b <- 9
c <- length(unique(test_mod2$from))- a-b



grid.col2 <- setNames(c(brewer.pal(n = a, name = "Set1"),
                        brewer.pal(n = b, name = "Pastel1"),
                        magma(c),
                        rep("grey",length(unique(test_mod2$to)))),
                      union(test_mod2$from, test_mod2$to))

par(cex = .7, mar = c(0, 0, 0, 0))
chordDiagram(test_mod2,annotationTrack = "grid", preAllocateTracks = 1,
             order=union(test_mod2$from, rev(test_mod2$to)),
             big.gap = 10,grid.col=grid.col2,
             symmetric=FALSE,
             link.arr.type = "big.arrow",
             #transparency of colours
             transparency = 0.2,
             scale = TRUE
)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

circos.clear()


###################
#https://www.r-graph-gallery.com/309-intro-to-hierarchical-edge-bundling.html

layout(matrix(1:2, 1, 2)) 
par(cex = .7, mar = c(0, 0, 0, 0))
chordDiagram(test_mod1,annotationTrack = "grid", preAllocateTracks = 1,
             order=union(test_mod1$from, rev(test_mod1$to)),
             big.gap = 10,grid.col=grid.col1,
             symmetric=FALSE,
             link.arr.type = "big.arrow",
             #transparency of colours
             transparency = 0.2,
             scale = TRUE
)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
#title("Modules (Plot 8)", cex = 1.0)
circos.clear()

#########

fitness_data2 <- read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")

# Filtering & relabeling
fitness2 <- fitness_data2 %>% filter(Year==2019,!Subplot == "OUT" & !is.na(ID))
test3 <- fitness2 %>% select(G_F,ID) %>% rename(from=ID)
test_mod3 <- test_mod1 %>% left_join(unique(test3),by="from")  
test_mod3 <- test_mod3 %>% ungroup() %>% select(-to) %>% rename(to=G_F) %>%select(from,to,value)%>%
  mutate(to=paste(to,"_GF",sep=""))
test_mod3$value <- 1

a <- 9
b <- 9
c <- length(unique(test_mod3$from))- a-b



grid.col3 <- setNames(c(brewer.pal(n = a, name = "Set1"),
                        brewer.pal(n = b, name = "Pastel1"),
                        magma(c),
                        rep("grey",length(unique(test_mod3$to)))),
                      union(test_mod3$from, test_mod3$to))

par(cex = .7, mar = c(0, 0, 0, 0))
chordDiagram(test_mod3,annotationTrack = "grid", preAllocateTracks = 1,
             order=union(test_mod3$from, rev(test_mod3$to)),
             big.gap = 10,grid.col=grid.col3,
             symmetric=FALSE,
             link.arr.type = "big.arrow",
             #transparency of colours
             transparency = 0.2,
             scale = TRUE
)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
#title("Functional Groups (Plot 8)", cex = 1.0)
circos.clear()

