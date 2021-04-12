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
pollination <- read_csv2("Raw_Data/final_Pollinators_2020.csv")
# Remove points from ID names
pollination$ID <- sub("\\.", "", pollination$ID)
pollination$ID_Simple <- sub("\\.", "", pollination$ID_Simple)

# filter tabanidae
pollination <- pollination %>% filter(ID != "Tabanidae")


pollination$Line <- NA

for (i in 1:nrow(pollination)){
  if(pollination$Plot[i] %in% c(1,2,3)){pollination$Line[i] <- 1}
  else if(pollination$Plot[i] %in% c(4,5,6)){pollination$Line[i] <- 2}
  else{pollination$Line[i] <- 3}
}


unique(pollination$ID_Simple[grep(" ",pollination$ID_Simple,ignore.case = T)])
#No labels with spaces -> Good!



cz_edge_module_info_IN <- NULL
cz_edge_module_info_OUT <- NULL

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
  # Fix those pollinator names that were modified during the extraction process
  
  plot_edge_list$to[plot_edge_list$to=="Lassioglosum_.immunitum"] <- "Lassioglosum_ immunitum"
  
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
        
        #For directed networks
        interlink_i<- tibble(layer_from=c(combination_layers[j,1],combination_layers[j,2]),
                             node_from=c(pollinators[i],pollinators[i]),
                             layer_to=c(combination_layers[j,2],combination_layers[j,1]),
                             node_to=c(pollinators[i],pollinators[i]),
                             weight=c(plant_pheno_overlap(combination_layers[j,1],
                                                          combination_layers[j,2],
                                                          pollination_20_i),
                                      plant_pheno_overlap(combination_layers[j,2],
                                                          combination_layers[j,1],
                                                          pollination_20_i)))
 
        #For directed
        S_edge_list <- bind_rows(S_edge_list,interlink_i)
      }
    }
  }
  
  S_edge_list_i <- S_edge_list %>% mutate(Plot=Plot_i)
  
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
  # CREATE NETWORK OF NETWORKS LIST
  #######################################
  
  NN_edge_list <- S_edge_list 
  
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
  
  Plot_NN <- create_monolayer_object(x = NN_edge_list_final_ID, directed = T, bipartite = F, node_metadata = NN_node_ID2)
  
  Plot_NN2 <- Plot_NN
  Plot_NN2$nodes$node_name <- as.numeric(Plot_NN2$nodes$node_name)
  # Run Infomap
  #modules_relax_rate <- run_infomap_multilayer(Plot_multilayer, relax = F, silent = T, flow_model = 'undirected', trials = 250, seed = 497294, temporal_network = F)
  modules_relax_rate_NN <- run_infomap_monolayer2(Plot_NN2, flow_model = 'directed', silent=T,trials=1000, two_level=T, seed=200952)
  
  # Extract information
  plot_modules_i <- modules_relax_rate_NN$modules %>% left_join(NN_node_ID,by="node_id") %>%
    mutate(layer_name=species,module=module_level1) %>% separate(layer_name,c("Org","layer_name")," ") %>%
    left_join(physical_nodes,by="species")
    
  plot_modules_i$type[is.na(plot_modules_i$type)] <- "pollinator"
  
  plot_modules_i$Plot <- Plot_i
  
  setwd(dir_ini)
  
  modules_from <- plot_modules_i %>% dplyr::select(Plot,module,species,layer_name,type) %>%
    rename(module_from=module,node_from=species,layer_from=layer_name,type_from=type)
  modules_to <- modules_from %>%
    rename(module_to=module_from,node_to=node_from,layer_to=layer_from,type_to=type_from)
  
  edge_module_info_i <- NN_edge_list %>% 
    left_join(modules_from, by = c("node_from","layer_from")) %>%
    left_join(modules_to, by = c("node_to","layer_to","Plot"))
  
  cz_edge_module_info_i_IN <- cz_values_function_IN(edge_module_info_i)
  cz_edge_module_info_i_OUT <- cz_values_function_OUT(edge_module_info_i)
  
  cz_edge_module_info_IN <- bind_rows(cz_edge_module_info_IN,cz_edge_module_info_i_IN)
  cz_edge_module_info_OUT <- bind_rows(cz_edge_module_info_OUT,cz_edge_module_info_i_OUT)

}

# Save data: Commented for security reasons
# write_csv(cz_edge_module_info_IN,"Processed_data/cz_edge_module_info_IN.csv")
# write_csv(cz_edge_module_info_OUT,"Processed_data/cz_edge_module_info_OUT.csv")

cz_edge_module_info_IN <- read_csv("Processed_data/cz_edge_module_info_IN.csv")
cz_edge_module_info_OUT <- read_csv("Processed_data/cz_edge_module_info_OUT.csv")

#######################
#UPLOADING THRESHOLDS
cz_edge_module_info_IN_NULL <- read_csv("Processed_data/Data_2020_cz_NULL/2020_NN_results_test_IN_cz_NULL_9.csv")
cz_edge_module_info_OUT_NULL <- read_csv("Processed_data/Data_2020_cz_NULL/2020_NN_results_test_OUT_cz_NULL_9.csv")

z_in_threshold <- quantile(cz_edge_module_info_IN_NULL$z, probs = 0.95, na.rm = T)
c_in_threshold <- quantile(cz_edge_module_info_IN_NULL$c, probs = 0.95, na.rm = T)

z_out_threshold <- quantile(cz_edge_module_info_OUT_NULL$z, probs = 0.95, na.rm = T)
c_out_threshold <- quantile(cz_edge_module_info_OUT_NULL$c, probs = 0.95, na.rm = T)


#################
# PREPARING GRAPH


plot_labs <-c(
  "Plot 1",
  "Plot 2",
  "Plot 3",
  "Plot 4",
  "Plot 5",
  "Plot 6",
  "Plot 7",
  "Plot 8",
  "Plot 9"
)
names(plot_labs) <- c(
  '1',
  '2',
  '3',
  '4',
  '5',
  "6",
  '7',
  '8',
  "9"
)


# read Functional Group ---------------------------------------------------------------

# Add G_F

G_F_list_to <- read_csv2("Raw_Data/final_Pollinators_2020.csv") %>%
  filter(ID != "Tabanidae") %>%
  dplyr::select(G_F,ID_Simple) %>% rename(node_to=ID_Simple) %>% unique() 

# Remove points from ID names
G_F_list_to$node_to <- sub("\\.", "", G_F_list_to$node_to)

G_F_list_from <- G_F_list_to %>% rename(node_from=node_to)


G_F_list_from <- unique(G_F_list_from)


###########

GF_cz_edge_module_info_IN <- cz_edge_module_info_IN %>% 
  mutate(node_to_aux=node_to) %>% separate(node_to,c("node_to","aux")," ") %>%
  dplyr::left_join(G_F_list_to,by = "node_to") %>% mutate(node_to=node_to_aux) %>%
  dplyr::select(-node_to_aux,-aux)

GF_cz_edge_module_info_OUT <- cz_edge_module_info_OUT %>%
  mutate(node_from_aux=node_from) %>% separate(node_from,c("node_from","aux")," ") %>%
  dplyr::left_join(G_F_list_from,by = "node_from") %>% mutate(node_from=node_from_aux) %>%
  dplyr::select(-node_from_aux,-aux)


GF_cz_edge_module_info_IN$G_F[is.na(GF_cz_edge_module_info_IN$G_F)] <- GF_cz_edge_module_info_IN$layer_to[is.na(GF_cz_edge_module_info_IN$G_F)]
GF_cz_edge_module_info_OUT$G_F[is.na(GF_cz_edge_module_info_OUT$G_F)] <- GF_cz_edge_module_info_OUT$layer_from[is.na(GF_cz_edge_module_info_OUT$G_F)]

GF_cz_edge_module_info_IN$G_F[GF_cz_edge_module_info_IN$G_F=="Flower_beetles"] <- "Flower beetles"
GF_cz_edge_module_info_IN$G_F[GF_cz_edge_module_info_IN$G_F=="House_flies"] <- "House flies"
GF_cz_edge_module_info_IN$G_F[GF_cz_edge_module_info_IN$G_F=="Small_beetles"] <- "Small beetles"
GF_cz_edge_module_info_IN$G_F[GF_cz_edge_module_info_IN$G_F=="Big_beetles"] <- "Big beetles"
GF_cz_edge_module_info_IN$G_F[GF_cz_edge_module_info_IN$G_F=="Small_flies"] <- "Small flies"
GF_cz_edge_module_info_IN$G_F[GF_cz_edge_module_info_IN$G_F=="Solitary_bees"] <- "Solitary bees"
GF_cz_edge_module_info_IN$G_F[GF_cz_edge_module_info_IN$G_F=="Wasp"] <- "Wasps"

GF_cz_edge_module_info_OUT$G_F[GF_cz_edge_module_info_OUT$G_F=="Flower_beetles"] <- "Flower beetles"
GF_cz_edge_module_info_OUT$G_F[GF_cz_edge_module_info_OUT$G_F=="House_flies"] <- "House flies"
GF_cz_edge_module_info_OUT$G_F[GF_cz_edge_module_info_OUT$G_F=="Small_beetles"] <- "Small beetles"
GF_cz_edge_module_info_OUT$G_F[GF_cz_edge_module_info_OUT$G_F=="Big_beetles"] <- "Big beetles"
GF_cz_edge_module_info_OUT$G_F[GF_cz_edge_module_info_OUT$G_F=="Small_flies"] <- "Small flies"
GF_cz_edge_module_info_OUT$G_F[GF_cz_edge_module_info_OUT$G_F=="Solitary_bees"] <- "Solitary bees"
GF_cz_edge_module_info_OUT$G_F[GF_cz_edge_module_info_OUT$G_F=="Wasp"] <- "Wasps"

#OUT-STRENGTH-----------------


library(RColorBrewer)

GF_cz_edge_module_info_OUT$G_F %>% unique()

GF_cz_edge_module_info_OUT_exp <- GF_cz_edge_module_info_OUT
GF_cz_edge_module_info_OUT_exp$G_F[GF_cz_edge_module_info_OUT_exp$G_F == "BEMA"] <- "B. macrocarpa"
GF_cz_edge_module_info_OUT_exp$G_F[GF_cz_edge_module_info_OUT_exp$G_F == "CETE"] <- "C. tenuiflorum"
GF_cz_edge_module_info_OUT_exp$G_F[GF_cz_edge_module_info_OUT_exp$G_F == "CHFU"] <- "C. fuscatum"
GF_cz_edge_module_info_OUT_exp$G_F[GF_cz_edge_module_info_OUT_exp$G_F == "CHMI"] <- "C. mixtum"
GF_cz_edge_module_info_OUT_exp$G_F[GF_cz_edge_module_info_OUT_exp$G_F == "LEMA"] <- "L. maroccanus"
GF_cz_edge_module_info_OUT_exp$G_F[GF_cz_edge_module_info_OUT_exp$G_F == "MESU"] <- "M. sulcatus"
GF_cz_edge_module_info_OUT_exp$G_F[GF_cz_edge_module_info_OUT_exp$G_F == "PUPA"] <- "P. paludosa"
GF_cz_edge_module_info_OUT_exp$G_F[GF_cz_edge_module_info_OUT_exp$G_F == "SCLA"] <- "S. laciniata"
GF_cz_edge_module_info_OUT_exp$G_F[GF_cz_edge_module_info_OUT_exp$G_F == "SOAS"] <- "S. asper"
GF_cz_edge_module_info_OUT_exp$G_F[GF_cz_edge_module_info_OUT_exp$G_F == "SPRU"] <- "S. rubra"


ggplot(GF_cz_edge_module_info_OUT_exp%>% filter(G_F %in% c("B. macrocarpa","C. tenuiflorum",
                                                       "C. fuscatum","C. mixtum",
                                                       "L. maroccanus","M. sulcatus",
                                                       "P. paludosa","S. laciniata",
                                                       "S. asper","S. rubra")))+
  geom_point(aes(x=c,y=z,color=G_F),position = "jitter",alpha=.5)+
  scale_color_brewer(palette = 'Paired')+
  geom_vline(xintercept = 1.0)+
  geom_hline(yintercept = 1.62)+
  theme_bw()+
  #scale_colour_brewer(palette = 'Paired')+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  labs(x ="Among module out-strength, out-c", y = "Within module in-strength, out-z",color=NULL)+ 
  theme(legend.position="bottom")+
  theme(legend.text = element_text(face = "italic"))



ggplot(GF_cz_edge_module_info_OUT_exp%>% filter(!G_F %in% c("B. macrocarpa","C. tenuiflorum",
                                                            "C. fuscatum","C. mixtum",
                                                            "L. maroccanus","M. sulcatus",
                                                            "P. paludosa","S. laciniata",
                                                            "S. asper","S. rubra")))+
  geom_point(aes(x=c,y=z,color=G_F),position = "jitter",alpha=.5)+
  scale_color_brewer(palette = 'Paired')+
  geom_vline(xintercept = 1.0)+
  geom_hline(yintercept = 1.62)+
  theme_bw()+
  scale_colour_brewer(palette = 'Paired')+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  labs(x ="Among module out-strength, out-c", y = "Within module out-strength, out-z",color=NULL)+ theme(legend.position="bottom")

total_num_plants <- GF_cz_edge_module_info_OUT %>% filter(G_F %in% c("LEMA","CETE",
                                                 "CHFU","PUPA",
                                                 "BEMA","CHMI",
                                                 "SPRU","MESU",
                                                 "SCLA","SOAS")) %>% nrow() #441 plants
total_num_poll <- GF_cz_edge_module_info_OUT %>% filter(!G_F %in% c("LEMA","CETE",
                                                  "CHFU","PUPA",
                                                  "BEMA","CHMI",
                                                  "SPRU","MESU",
                                                  "SCLA","SOAS")) %>% nrow() #206 visitors 

#HUBS
GF_cz_edge_module_info_OUT %>% filter(G_F %in% c("LEMA","CETE",
                                                 "CHFU","PUPA",
                                                 "BEMA","CHMI",
                                                 "SPRU","MESU",
                                                 "SCLA","SOAS"),z>z_out_threshold,c==1)%>%
  nrow()/total_num_plants #0 plants
GF_cz_edge_module_info_OUT %>% filter(!G_F %in% c("LEMA","CETE",
                                                  "CHFU","PUPA",
                                                  "BEMA","CHMI",
                                                  "SPRU","MESU",
                                                  "SCLA","SOAS"),z>z_out_threshold,c==1)%>%
  nrow()/total_num_poll #5 animal

GF_cz_edge_module_info_OUT %>% filter(!G_F %in% c("LEMA","CETE",
                                                  "CHFU","PUPA",
                                                  "BEMA","CHMI",
                                                  "SPRU","MESU",
                                                  "SCLA","SOAS"),z>z_out_threshold,c==1)%>%
  group_by(node_from, layer_from,Plot,G_F) %>% count()

#Perip
GF_cz_edge_module_info_OUT %>% filter(G_F %in% c("LEMA","CETE",
                                                 "CHFU","PUPA",
                                                 "BEMA","CHMI",
                                                 "SPRU","MESU",
                                                 "SCLA","SOAS"),z<z_out_threshold,c<1)%>%
  nrow()/total_num_plants #10 plants
GF_cz_edge_module_info_OUT %>% filter(G_F %in% c("LEMA","CETE",
                                                 "CHFU","PUPA",
                                                 "BEMA","CHMI",
                                                 "SPRU","MESU",
                                                 "SCLA","SOAS"),z<z_out_threshold,c<1)%>%
  group_by(layer_from,Plot) %>% count() %>% arrange(Plot) #10 plants


GF_cz_edge_module_info_OUT %>% filter(!G_F %in% c("LEMA","CETE",
                                                  "CHFU","PUPA",
                                                  "BEMA","CHMI",
                                                  "SPRU","MESU",
                                                  "SCLA","SOAS"),z<z_out_threshold,c<1)%>%
  nrow()/total_num_poll #3 animal


#connectors
GF_cz_edge_module_info_OUT %>% filter(G_F %in% c("LEMA","CETE",
                                                 "CHFU","PUPA",
                                                 "BEMA","CHMI",
                                                 "SPRU","MESU",
                                                 "SCLA","SOAS"),
                                      z<z_out_threshold,c==1)%>%
  nrow()/total_num_plants #383 plants
GF_cz_edge_module_info_OUT %>% filter(!G_F %in% c("LEMA","CETE",
                                                  "CHFU","PUPA",
                                                  "BEMA","CHMI",
                                                  "SPRU","MESU",
                                                  "SCLA","SOAS"),
                                     z<z_out_threshold,c==1)%>%
  nrow()/total_num_poll #169 animal


# Module hubs
GF_cz_edge_module_info_OUT %>% filter(G_F %in% c("LEMA","CETE",
                                                 "CHFU","PUPA",
                                                 "BEMA","CHMI",
                                                 "SPRU","MESU",
                                                 "SCLA","SOAS"),
                                      z>z_out_threshold,c<1)%>%
  nrow()/total_num_plants #383 plants
GF_cz_edge_module_info_OUT %>% filter(!G_F %in% c("LEMA","CETE",
                                                  "CHFU","PUPA",
                                                  "BEMA","CHMI",
                                                  "SPRU","MESU",
                                                  "SCLA","SOAS"),
                                      z>z_out_threshold,c<1)%>%
  nrow()/total_num_poll #169 animal

GF_cz_edge_module_info_OUT%>% filter(!G_F %in% c("LEMA","CETE",
                                                 "CHFU","PUPA",
                                                 "BEMA","CHMI",
                                                 "SPRU","MESU",
                                                 "SCLA","SOAS"),z>z_out_threshold,c<1) %>% 
  select(node_from,layer_from,Plot) %>% arrange(node_from) %>%
  unique()

# Network HUBS
hubs_OUT <- GF_cz_edge_module_info_OUT%>% filter(!G_F %in% c("LEMA","CETE",
                                                 "CHFU","PUPA",
                                                 "BEMA","CHMI",
                                                 "SPRU","MESU",
                                                 "SCLA","SOAS"),z>z_out_threshold,c==1) %>% 
  select(node_from,layer_from,Plot) %>% arrange(node_from) %>%
  unique()



#IN-STRENGTH-----------------


library(RColorBrewer)
GF_cz_edge_module_info_IN_exp <- GF_cz_edge_module_info_IN
GF_cz_edge_module_info_IN_exp$G_F[GF_cz_edge_module_info_IN_exp$G_F == "BEMA"] <- "B. macrocarpa"
GF_cz_edge_module_info_IN_exp$G_F[GF_cz_edge_module_info_IN_exp$G_F == "CETE"] <- "C. tenuiflorum"
GF_cz_edge_module_info_IN_exp$G_F[GF_cz_edge_module_info_IN_exp$G_F == "CHFU"] <- "C. fuscatum"
GF_cz_edge_module_info_IN_exp$G_F[GF_cz_edge_module_info_IN_exp$G_F == "CHMI"] <- "C. mixtum"
GF_cz_edge_module_info_IN_exp$G_F[GF_cz_edge_module_info_IN_exp$G_F == "LEMA"] <- "L. maroccanus"
GF_cz_edge_module_info_IN_exp$G_F[GF_cz_edge_module_info_IN_exp$G_F == "MESU"] <- "M. sulcatus"
GF_cz_edge_module_info_IN_exp$G_F[GF_cz_edge_module_info_IN_exp$G_F == "PUPA"] <- "P. paludosa"
GF_cz_edge_module_info_IN_exp$G_F[GF_cz_edge_module_info_IN_exp$G_F == "SCLA"] <- "S. laciniata"
GF_cz_edge_module_info_IN_exp$G_F[GF_cz_edge_module_info_IN_exp$G_F == "SOAS"] <- "S. asper"
GF_cz_edge_module_info_IN_exp$G_F[GF_cz_edge_module_info_IN_exp$G_F == "SPRU"] <- "S. rubra"


ggplot(GF_cz_edge_module_info_IN_exp%>% filter(G_F %in% c("B. macrocarpa","C. tenuiflorum",
                                                          "C. fuscatum","C. mixtum",
                                                          "L. maroccanus","M. sulcatus",
                                                          "P. paludosa","S. laciniata",
                                                          "S. asper","S. rubra")))+
  geom_point(aes(x=c,y=z,color=G_F),position = "jitter",alpha=.5)+
  geom_vline(xintercept = c_in_threshold)+
  geom_hline(yintercept = z_in_threshold)+
  theme_bw()+
  scale_colour_brewer(palette = 'Paired')+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  labs(x ="Among module in-strength, in-c", y = "Within module in-strength, in-z",color=NULL)+ 
  theme(legend.position="bottom")+
  theme(legend.text = element_text(face = "italic"))



ggplot(GF_cz_edge_module_info_IN%>% filter(!G_F %in% c("LEMA","CETE",
                                                       "CHFU","PUPA",
                                                       "BEMA","CHMI",
                                                       "SPRU","MESU",
                                                       "SCLA","SOAS")))+
  geom_point(aes(x=c,y=z,color=G_F),position = "jitter",alpha=.5)+
  geom_vline(xintercept = c_in_threshold)+
  geom_hline(yintercept = z_in_threshold)+
  theme_bw()+
  scale_colour_brewer(palette = 'Paired')+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  labs(x ="Among module in-strength, in-c", y = "Within module in-strength, in-z",color=NULL)+ theme(legend.position="bottom")

GF_cz_edge_module_info_IN %>% filter(G_F %in% c("LEMA","CETE",
                                                "CHFU","PUPA",
                                                "BEMA","CHMI",
                                                "SPRU","MESU",
                                                "SCLA","SOAS")) #441 plants

total_num_plants <- GF_cz_edge_module_info_IN %>% filter(G_F %in% c("LEMA","CETE",
                                                "CHFU","PUPA",
                                                "BEMA","CHMI",
                                                "SPRU","MESU",
                                                "SCLA","SOAS")) %>% nrow()

GF_cz_edge_module_info_IN %>% filter(!G_F %in% c("LEMA","CETE",
                                                 "CHFU","PUPA",
                                                 "BEMA","CHMI",
                                                 "SPRU","MESU",
                                                 "SCLA","SOAS")) #206 visitors 

total_num_poll <- GF_cz_edge_module_info_IN %>% filter(!G_F %in% c("LEMA","CETE",
                                                                   "CHFU","PUPA",
                                                                   "BEMA","CHMI",
                                                                   "SPRU","MESU",
                                                                   "SCLA","SOAS")) %>%
  nrow()

#HUBS
GF_cz_edge_module_info_IN %>% filter(G_F %in% c("LEMA","CETE",
                                                "CHFU","PUPA",
                                                "BEMA","CHMI",
                                                "SPRU","MESU",
                                                "SCLA","SOAS"),z>z_in_threshold) #0 plants
GF_cz_edge_module_info_IN %>% filter(!G_F %in% c("LEMA","CETE",
                                                 "CHFU","PUPA",
                                                 "BEMA","CHMI",
                                                 "SPRU","MESU",
                                                 "SCLA","SOAS"),z>z_in_threshold & c==1) #17 animal

GF_cz_edge_module_info_IN %>% filter(G_F %in% c("LEMA","CETE",
                                                "CHFU","PUPA",
                                                "BEMA","CHMI",
                                                "SPRU","MESU",
                                                "SCLA","SOAS"),z>z_in_threshold & c==1) %>%
  nrow()/total_num_plants
GF_cz_edge_module_info_IN %>% filter(!G_F %in% c("LEMA","CETE",
                                                 "CHFU","PUPA",
                                                 "BEMA","CHMI",
                                                 "SPRU","MESU",
                                                 "SCLA","SOAS"),z>z_in_threshold & c==1) %>%
  nrow()/total_num_poll

GF_cz_edge_module_info_IN %>% filter(!G_F %in% c("LEMA","CETE",
                                                 "CHFU","PUPA",
                                                 "BEMA","CHMI",
                                                 "SPRU","MESU",
                                                 "SCLA","SOAS"),z>z_in_threshold & c==1)%>%
  group_by(node_to, layer_to,Plot,G_F) %>% count()

#Perip
GF_cz_edge_module_info_IN %>% filter(G_F %in% c("LEMA","CETE",
                                                "CHFU","PUPA",
                                                "BEMA","CHMI",
                                                "SPRU","MESU",
                                                "SCLA","SOAS"),z<z_in_threshold & c<1) %>%
  nrow()/total_num_plants
GF_cz_edge_module_info_IN %>% filter(G_F %in% c("LEMA","CETE",
                                                "CHFU","PUPA",
                                                "BEMA","CHMI",
                                                "SPRU","MESU",
                                                "SCLA","SOAS"),z<z_in_threshold & c<1)%>%
  group_by(layer_to,Plot) %>% count() %>% arrange(Plot) #39 plants

GF_cz_edge_module_info_IN %>% filter(!G_F %in% c("LEMA","CETE",
                                                 "CHFU","PUPA",
                                                 "BEMA","CHMI",
                                                 "SPRU","MESU",
                                                 "SCLA","SOAS"),z<z_in_threshold & c<1) %>%
  nrow()/total_num_poll


# Connectors
GF_cz_edge_module_info_IN %>% filter(G_F %in% c("LEMA","CETE",
                                                "CHFU","PUPA",
                                                "BEMA","CHMI",
                                                "SPRU","MESU",
                                                "SCLA","SOAS"),
                                     z<z_in_threshold & c==1) %>%
  nrow()/total_num_plants
GF_cz_edge_module_info_IN %>% filter(!G_F %in% c("LEMA","CETE",
                                                 "CHFU","PUPA",
                                                 "BEMA","CHMI",
                                                 "SPRU","MESU",
                                                 "SCLA","SOAS"),
                                     z<z_in_threshold & c==1)%>%
  nrow()/total_num_poll #169 animal


GF_cz_edge_module_info_IN %>% filter(G_F %in% c("LEMA","CETE",
                                                "CHFU","PUPA",
                                                "BEMA","CHMI",
                                                "SPRU","MESU",
                                                "SCLA","SOAS"),z>2.6) %>%
  select(node_to,Plot) %>%
  unique()

GF_cz_edge_module_info_IN%>% filter(G_F %in% c("LEMA","CETE",
                                               "CHFU","PUPA",
                                               "BEMA","CHMI",
                                               "SPRU","MESU",
                                               "SCLA","SOAS"),c<0.6) %>% 
  select(layer_to,Plot) %>%
  group_by(layer_to,Plot) %>% count() 

# Module hubs
GF_cz_edge_module_info_IN %>% filter(G_F %in% c("LEMA","CETE",
                                                "CHFU","PUPA",
                                                "BEMA","CHMI",
                                                "SPRU","MESU",
                                                "SCLA","SOAS"),
                                     z>z_in_threshold & c<1) %>%
  nrow()/total_num_plants
GF_cz_edge_module_info_IN %>% filter(!G_F %in% c("LEMA","CETE",
                                                 "CHFU","PUPA",
                                                 "BEMA","CHMI",
                                                 "SPRU","MESU",
                                                 "SCLA","SOAS"),
                                     z>z_in_threshold & c<1)%>%
  nrow()/total_num_poll #169 animal

GF_cz_edge_module_info_IN %>% filter(!G_F %in% c("LEMA","CETE",
                                                 "CHFU","PUPA",
                                                 "BEMA","CHMI",
                                                 "SPRU","MESU",
                                                 "SCLA","SOAS"),
                                     z>z_in_threshold & c<1) %>%
  select(node_to,layer_to,Plot) %>%
  unique() %>% arrange(node_to)

#Network HUBS
GF_cz_edge_module_info_IN%>% filter(!G_F %in% c("LEMA","CETE",
                                                "CHFU","PUPA",
                                                "BEMA","CHMI",
                                                "SPRU","MESU",
                                                "SCLA","SOAS"),z>2.6) %>%
  select(node_to,layer_to,Plot) %>%
  unique() %>% arrange(node_to)

GF_cz_edge_module_info_IN%>% filter(!G_F %in% c("LEMA","CETE",
                                                "CHFU","PUPA",
                                                "BEMA","CHMI",
                                                "SPRU","MESU",
                                                "SCLA","SOAS"),z<2.6,c>0.6) %>%
  select(node_to,layer_to,Plot) %>%
  unique()

HUBS_pol_IN <- GF_cz_edge_module_info_IN %>% filter(G_F %in% c("LEMA","CETE",
                                                              "CHFU","PUPA",
                                                              "BEMA","CHMI",
                                                              "SPRU","MESU",
                                                              "SCLA","SOAS"),z>2.6) %>% 
  select(node_to,layer_to,Plot) %>%
  unique()


GF_cz_edge_module_info_IN%>% filter(!G_F %in% c("LEMA","CETE",
                                                "CHFU","PUPA",
                                                "BEMA","CHMI",
                                                "SPRU","MESU",
                                                "SCLA","SOAS"),c<0.6) %>% 
  select(node_to,layer_to,Plot) %>%
  unique()

HUBS_pol_IN_2<- HUBS_pol_IN %>% rename(node=node_to,layer=layer_to)
HUBS_pol_OUT_2 <- HUBS_pol_OUT %>%rename(node=node_from,layer=layer_from)

#filter rows in df1 that are not in df2
anti_join(HUBS_pol_IN_2,HUBS_pol_OUT_2)

#filter rows in df1 that are also in df2
semi_join(HUBS_pol_IN_2,HUBS_pol_OUT_2)

