
#-----
##################################################################
# FUNCTION (NAME): motifs_extraction
# INPUT (1) -> visit_list: it contains (Plot,ID,Subplot_Plant_Label(e.g. A1 LEMA),Visits_tot)
# OUTPUT (1) -> list with motifs: it contains
# (plot_id,number_nodes,number_plants,plant_1,subplot_plant_1,
# plant_2,subplot_plant_2,poll_1,poll_2
##################################################################

motifs_extraction <- function(visit_list) {
  
  ##############################################################
  # GENERATE A LIST OF BIPARTITE NETWORKS (1 NETWORK PER PLOT)
  ##############################################################
  
  testdata_19 <-   data.frame(higher = visit_list$ID,
                              lower = visit_list$Subplot_Plant_Label,
                              webID = visit_list$Plot,
                              freq = visit_list$Visits_tot)
  
  list_incid_matrix_19 <- frame2webs(testdata_19,type.out="list")
  
  names_plot <- names(list_incid_matrix_19)
  
  first_motif_encuounter <- F
  
  for (plot_i in 1:length(list_incid_matrix_19)){
    
    
    #print(names_plot[plot_i])
    
    # Incidence matrix for each plot network
    
    incid_matrix_i <- list_incid_matrix_19[[plot_i ]] 
    
    graph_i <- graph_from_incidence_matrix(incid_matrix_i, weighted = T, directed = F)
    
    # motiv in igraph (triplet is equal to a two-path graph, i.e., a star graph with 2 nodes)
    
    pattern <- make_star(3, mode = "undirected")
    
    iso <- subgraph_isomorphisms(pattern, graph_i)      # takes a while
    
    # WARNING subgraph_isomorphisms works with directed graphs and
    # at the end of the day it duplicates the motifs in our list
    # check manual: https://igraph.org/r/doc/subgraph_isomorphisms.html
    
    motifs <- lapply(iso, function (x) { induced_subgraph(graph_i, x) })
    
    
    ##############################################################
    # TRIPLETS ANALYSIS
    ##############################################################
    
    #print(length(motifs))
    
    if (length(motifs)>0){
      
      if (first_motif_encuounter==F){
        first_motif_encuounter <- T
        fist_plot_i <- plot_i}
      
      # We collect triplets information in "motif_3" tibble
      
      tbl_colnames <- c("plot_id","number_nodes","number_plants","plant_1","subplot_plant_1",
                        "plant_2","subplot_plant_2","poll_1","poll_2")
      
      motif_3 <- as_tibble(data.frame(matrix(nrow=length(motifs),ncol=length(tbl_colnames))))
      colnames(motif_3) <- tbl_colnames
      
      motif_3$plot_id <- as.character(names_plot[plot_i])
      
      for (i in 1:length(motifs)){
        
        lista_names_i <- V(motifs[[i]])$name
        motif_3$number_nodes[i] <- length(lista_names_i)
        number_plants <- 0
        
        list_plants <- c()
        list_plants_loc <- c() 
        list_poll <- c()
        
        for (j in 1:length(lista_names_i)) {
          
          if (lista_names_i[j] %in%  visit_list$Subplot_Plant_Label){
            
            location_plant <- strsplit(lista_names_i[j]," ")
            list_plants_loc <- c(list_plants_loc,c(location_plant[[1]][1]))
            list_plants <- c(list_plants,c(location_plant[[1]][2]))
            number_plants <- number_plants + 1
            
          }else{
            
            list_poll <- c(list_poll,c(lista_names_i[j]))
          }
          
        }
        
        if (number_plants==1) {
          
          motif_3$number_plants[i] <- number_plants
          motif_3$plant_1[i] <- list_plants[1]
          motif_3$subplot_plant_1[i] <- list_plants_loc[1]
          motif_3$poll_1[i] <- list_poll[1]
          motif_3$poll_2[i] <- list_poll[2]
          
          iden_plant <- paste(list_plants_loc[1],list_plants[1],sep=" ")
          
          
        } else {
          
          motif_3$number_plants[i] <- number_plants
          motif_3$plant_1[i] <- list_plants[1]
          motif_3$subplot_plant_1[i] <- list_plants_loc[1]
          motif_3$plant_2[i] <- list_plants[2]
          motif_3$subplot_plant_2[i] <- list_plants_loc[2]
          motif_3$poll_1[i] <- list_poll[1]
          
        }
      }
      print(paste("Triplets_plot",names_plot[plot_i],sep=""))
      
      if (plot_i==fist_plot_i){
        motif_3 <- unique(motif_3) #to avoid duplicities due to subgraph isomorphisms	
        motif_3_list <- motif_3
      }else{
        motif_3 <- unique(motif_3)
        motif_3_list <- motif_3_list %>% bind_rows(motif_3) 
      }
    }
    
  }
  return(motif_3_list)
}

#--------
##################################################################
# FUNCTION (NAME): homo_hete_motifs
# INPUT (1) -> visit_list: it contains (Plot,ID,Subplot_Plant_Label(e.g. A1 LEMA),Visits_tot)
# OUTPUT (1) -> visit_list with 2 new columns: homo_motif, hete_motif
##################################################################

homo_hete_motifs <- function(visit_list) {
  
  motif_3_Carac <- motifs_extraction(visit_list)
  
  motif_3_Carac <- motif_3_Carac %>%
    mutate(descript_plant_1=NA,descript_plant_2=NA,Same_plant = NA)
  
  for (i in 1:nrow(motif_3_Carac)){
    
    motif_3_Carac$descript_plant_1[i] <- paste(motif_3_Carac$plot_id[i],
                                               motif_3_Carac$subplot_plant_1[i],
                                               motif_3_Carac$plant_1[i],
                                               sep = " ")
    
    if (!is.na(motif_3_Carac$plant_2[i])){
      motif_3_Carac$descript_plant_2[i] <- paste(motif_3_Carac$plot_id[i],
                                                 motif_3_Carac$subplot_plant_2[i],
                                                 motif_3_Carac$plant_2[i],
                                                 sep = " ")
    }
    if (motif_3_Carac$number_plants[i]==2 && 
        motif_3_Carac$plant_2[i] == motif_3_Carac$plant_1[i] 
    ){
      motif_3_Carac$Same_plant[i]<-TRUE}
    else if(motif_3_Carac$number_plants[i]==2 && 
            motif_3_Carac$plant_2[i] != motif_3_Carac$plant_1[i]){
      
      motif_3_Carac$Same_plant[i]<-FALSE}
    
    else{motif_3_Carac$Same_plant[i]<-NA}
  }
  
  output_funct <- visit_list %>% mutate(homo_motif=NA,hete_motif=NA)
  
  # Homo_Motifs: same plant species in other subplots
  
  for (i in 1:nrow(output_funct)){
    
    descrip <-  paste(output_funct$Plot[i],output_funct$Subplot_Plant_Label[i],sep = " ")
    
    homo_motif <- motif_3_Carac %>% filter(Same_plant & poll_1 == output_funct$ID[i] &
                                             (motif_3_Carac$descript_plant_1 == descrip|motif_3_Carac$descript_plant_2==descrip))
    num_homo_motif <- sum(homo_motif$Same_plant,na.rm = TRUE)
    
    hete_motif <- motif_3_Carac %>% filter(!Same_plant & number_plants==2 & poll_1 == output_funct$ID[i] &
                                             (motif_3_Carac$descript_plant_1== descrip|motif_3_Carac$descript_plant_2==descrip))
    num_hete_motif <- sum(!hete_motif$Same_plant,na.rm = TRUE)
    
    output_funct$homo_motif[i] <- num_homo_motif
    output_funct$hete_motif[i] <- num_hete_motif
  }
  return(output_funct)
}


#------
# FUNCTION: plant_pheno_overlap 
# Returns the number of weeks of phenological overlap between plant1 and plant2 over
# divided by the total duration of plant 1 phenology

plant_pheno_overlap <- function(plant1,plant2,pollination_20_i) {
  
  plantsx <- sort(unique(pollination_20_i$Plant))
  if(sum(c(plant1,plant2) %in% plantsx)<2){
    print(paste("Error: At least one plant does not belong to plot ",
                Plot_i," experimental phenology",sep=""))
  }else{
    pollination_plant1 <- pollination_20_i %>% filter(Plant==plant1)
    #Weeks in which plant 1 recieves visits
    plant_1_weeks <- min(pollination_plant1$Week):max(pollination_plant1$Week)
    
    pollination_plant2 <- pollination_20_i %>% filter(Plant==plant2)
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

#------

create_weighted_link_list <- function(plot_edge_list){
  
  # Replace the node names with node_ids
  
  Plot_edgelist_complete <- tibble(layer_from=plot_edge_list$species,
                                   node_from=plot_edge_list$from,
                                   layer_to=plot_edge_list$species,
                                   node_to=plot_edge_list$to,
                                   weight=plot_edge_list$weight)
  
  # Obtain plant ind. total strength
  plant_strength <- Plot_edgelist_complete %>% group_by(layer_from,node_from) %>% 
    count(wt = weight) %>% rename(strength = n)
  
  # Obtain poll. sp. total strength
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
    print(pollinators[i])
    polinator_edges <- Plot_edgelist_complete %>% filter(node_to==pollinators[i])
    polinator_layers <- unique(polinator_edges$layer_to)
    
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
                                                          example),
                                      plant_pheno_overlap(combination_layers[j,2],
                                                          combination_layers[j,1],
                                                          example)))
        
        
        #For directed
        S_edge_list <- bind_rows(S_edge_list,interlink_i)
      }
    }
  }
  
  S_edge_list_i <- S_edge_list %>% mutate(Plot=Plot_i)
  
  interlinks_Plot_i <- S_edge_list %>% filter(node_from==node_to)
  
  if(nrow(interlinks_Plot_i)==0){ #If there are no interlinks, we create a dummy one
    
    list_possible_layers <- layer_metadata$layer_name[layer_metadata$layer_name!=S_edge_list_i$layer_from[1]]
    
    new_row <- tibble(
      layer_from=list_possible_layers,
      node_from=S_edge_list_i$node_to[1],
      layer_to=S_edge_list_i$layer_to[1],
      node_to=S_edge_list_i$node_to[1],
      weight=0.0,
      Plot=S_edge_list_i$Plot[1],
    )
    
    S_edge_list <- bind_rows(S_edge_list,new_row)
    
  }
  
  return(S_edge_list)
  
}

#------

extract_modules_NN_infomapec <- function(S_edge_list){
  # Transforms S_edge_list to a network of network list
  
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
  
  Plot_NN <- create_monolayer_object(x = NN_edge_list_final_ID, directed = T, bipartite = F, node_metadata = NN_node_ID2)
  
  Plot_NN2 <- Plot_NN
  Plot_NN2$nodes$node_name <- as.numeric(Plot_NN2$nodes$node_name)
  # Run Infomap
  modules_relax_rate_NN <- run_infomap_monolayer2(Plot_NN2, flow_model = 'directed', silent=T,trials=1000, two_level=T, seed=200952)
  
  # Extract information
  plot_modules_NN_i_aux <- modules_relax_rate_NN$modules %>% 
    dplyr::select(node_id,module_level1) %>% rename(module=module_level1) %>%
    left_join(NN_node_ID,by="node_id")
  
  
  plot_modules_NN_i_aux$Plot <- Plot_i
  
  plot_modules_NN_i <- plot_modules_NN_i_aux %>% separate(species,c("species","layer_name")," ") %>% 
    left_join(physical_nodes,by=c("species")) %>% dplyr::select(-node_id.y) %>% 
    rename(node_id=node_id.x)
  
  # Add type to plants
  
  for (i in 1:nrow(plot_modules_NN_i)){
    if(is.na(plot_modules_NN_i$type[i])){
      plot_modules_NN_i$type[i] <- "plant"
      plot_modules_NN_i$species[i] <- paste0(plot_modules_NN_i$species[i]," ",plot_modules_NN_i$layer_name[i])
    }
  }
  
  return(plot_modules_NN_i)
}

#------

centrality_metrics_NN <- function(S_edge_list){
  
  # Create Network of Networks edge list
  
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
  
  
  # Prepare igraph representation to get the centrality metrics
  
  graph_Plot_i <- igraph::graph_from_edgelist(as.matrix(NN_edge_list_final[,1:2]), directed = TRUE)
  
  E(graph_Plot_i)$weight <- pull(NN_edge_list_final[,3])
  
  # Sanity check: pull nodes and edge weights
  igraph::get.data.frame(graph_Plot_i)
  ##############
  # PAgeRAnk Multilayer
  
  page_rank_i <- igraph::page_rank(graph_Plot_i,
                                   directed = TRUE, damping = 0.85,
                                   personalized = NULL, weights = NULL, options = NULL)
  
  
  page_rank_i_nodes <-
    data.frame(species=names(page_rank_i[[1]]),
               Real_PR_Multi=page_rank_i[[1]], row.names=NULL)
  
  # Sanity check: 
  sum(page_rank_i_nodes$Real_PR_Multi) # 1, Correct
  
  #################
  # PageRank Layer:We remove the interlinks
  
  NN_edge_list_final_filter <- NN_edge_list %>% filter(layer_from==layer_to) %>%
    select(node_from,node_to,weight) %>%
    rename(from = node_from, to = node_to)
  
  graph_Plot_i_fil <- igraph::graph_from_edgelist(as.matrix(NN_edge_list_final_filter[,1:2]), directed = TRUE)
  
  E(graph_Plot_i_fil)$weight <- pull(NN_edge_list_final_filter[,3])
  
  # Sanity check: pull nodes and edge weights
  igraph::get.data.frame(graph_Plot_i_fil)
  
  # PageRAnk Layer (we remove the interlinks)
  
  page_rank_i <- igraph::page_rank(graph_Plot_i_fil,
                                   directed = TRUE, damping = 0.85,
                                   personalized = NULL, weights = NULL, options = NULL)
  
  
  page_rank_layer_i_nodes <-
    data.frame(species=names(page_rank_i[[1]]),
               Real_PR_Layer=page_rank_i[[1]], row.names=NULL)
  
  # Sanity check: 
  sum(page_rank_layer_i_nodes$Real_PR_Layer) # 1, Correct
  
  ######################
  # Strength + Degree
  StrengthIn_i <- igraph::strength(graph_Plot_i, mode = c("in"),
                                   loops = TRUE,weights = NULL)
  
  StrengthIn_i_nodes <- data.frame(species=names(StrengthIn_i), StrengthIn=StrengthIn_i, row.names=NULL) 
  
  StrengthOut_i <- igraph::strength(graph_Plot_i, mode = c("out"),
                                    loops = TRUE,weights = NULL)
  
  StrengthOut_i_nodes <- data.frame(species=names(StrengthOut_i), StrengthOut=StrengthOut_i, row.names=NULL) 
  
  
  degreeIn_i <- igraph::degree(graph_Plot_i, mode = c("in"),
                               loops = TRUE)
  
  degreeIn_i_nodes <- data.frame(species=names(degreeIn_i), DegreeIn=degreeIn_i, row.names=NULL) 
  
  degreeOut_i <- igraph::degree(graph_Plot_i, mode = c("out"),
                                loops = TRUE)
  
  degreeOut_i_nodes <- data.frame(species=names(degreeOut_i), DegreeOut=degreeOut_i, row.names=NULL) 
  
  
  centrality_i <- page_rank_i_nodes %>%
    left_join(page_rank_layer_i_nodes,by="species") %>%
    left_join(StrengthIn_i_nodes,by="species") %>%
    left_join(StrengthOut_i_nodes,by="species") %>%
    left_join(degreeIn_i_nodes,by="species") %>%
    left_join(degreeOut_i_nodes,by="species") %>%
    mutate(Delta=Real_PR_Multi-Real_PR_Layer,Ratio=Real_PR_Multi/Real_PR_Layer,Plot=Plot_i)
  
  
  return(centrality_i)

}



#------

# This function generates a ramdom list of visitis between animals and plants
# The IDs of animals are those of the original edge/visit list
# The IDs of focal plants are selected from those of the original edge/visit list
# The total amount of links per animal and plant species is equal to that of the original list

random_plot_edge_list_plants <- function(plot_edge_list,complete_random_selection){
  
  total_visits_anim_plant <- plot_edge_list %>% group_by(to,species) %>% count(wt=weight)
  
  random_plot_edge_list <- NULL
  
  for(i in 1:nrow(total_visits_anim_plant)){
    
    filter_visits <- plot_edge_list %>% filter(species==total_visits_anim_plant$species[i])
    
    candidate_plants <- filter_visits %>%
      group_by(from) %>% count(wt=weight) %>%
      mutate(percentage = n/sum(filter_visits$weight))
    
    if (complete_random_selection == T){
      visited_plants <- sample(candidate_plants$from, total_visits_anim_plant$n[i], replace=TRUE)
    }else{
      visited_plants <- sample(candidate_plants$from, total_visits_anim_plant$n[i], replace=TRUE, prob = candidate_plants$percentage)
    }
    
    
    
    visited_plants2 <- as_tibble(as.data.frame(table(visited_plants))) %>% rename(from=visited_plants,weight=Freq)
    
    visited_plants2$from <- as.character(visited_plants2$from)
    visited_plants2$to <- total_visits_anim_plant$to[i]
    visited_plants2$species <- total_visits_anim_plant$species[i]
    
    visited_plants2 <- visited_plants2 %>% dplyr::select(from,to,weight,species)
    
    random_plot_edge_list <- bind_rows(random_plot_edge_list,visited_plants2)
    
  }
  
  return(random_plot_edge_list)
  
}
#----
# This function generates a ramdom list of visitis between animals and plants
# The IDs of animals are those of the original edge/visit list
# The IDs of focal plants are selected from those of the original edge/visit list
# The total amount of links per animal and plant species is equal to that of the original list

random_plot_edge_list_plants2 <- function(plot_edge_list,complete_random_selection){
  
  total_visits_anim_plant <- plot_edge_list %>% group_by(to) %>% count(wt=weight)
  
  random_plot_edge_list <- NULL
  
  for(i in 1:nrow(total_visits_anim_plant)){
    
    candidate_plants <- plot_edge_list %>%
      group_by(from) %>% count(wt=weight) %>%
      mutate(percentage = n/sum(plot_edge_list$weight))
    
    if (complete_random_selection == T){
      visited_plants <- sample(candidate_plants$from, total_visits_anim_plant$n[i], replace=TRUE)
    }else{
      visited_plants <- sample(candidate_plants$from, total_visits_anim_plant$n[i], replace=TRUE, prob = candidate_plants$percentage)
    }
    
    visited_plants2 <- as_tibble(as.data.frame(table(visited_plants))) %>% rename(from=visited_plants,weight=Freq)
    
    visited_plants2$from <- as.character(visited_plants2$from)
    visited_plants2$to <- total_visits_anim_plant$to[i]
    
    species <- visited_plants2 %>% separate(from,c("Subplot","species")," ") %>% 
      dplyr::select(species)
    
    
    visited_plants2$species <-  pull(species)
    
    visited_plants2 <- visited_plants2 %>% dplyr::select(from,to,weight,species)
    
    random_plot_edge_list <- bind_rows(random_plot_edge_list,visited_plants2)
    
  }
  
  return(random_plot_edge_list)
  
}

#-------
# This function generates a ramdom list of visitis between animals and plants
# The IDs of animals are those of the original edge/visit list
# The IDs of focal plants are selected from those of the original edge/visit list
# The total amount of links per animal and plant species is equal to that of the original list

random_visits_WEEKS <- function(visit_list_week){
  
  visit_list_week_aux <- visit_list_week %>% mutate(aux=Subplot_Plant_Label) %>%
    separate(aux,c("Subplot","Plant")," ") %>% select(-Subplot)
  
  total_visits_anim_plant <- visit_list_week_aux %>% group_by(Plot,ID,Plant) %>% count(wt=Visits_tot)
  
  random_visit_list_week <- NULL
  
  for(i in 1:nrow(total_visits_anim_plant)){
    
    filter_visits <- visit_list_week_aux %>% filter(Plant==total_visits_anim_plant$Plant[i],
                                                    Plot==total_visits_anim_plant$Plot[i])
    
    candidate_plants <- unique(filter_visits$Subplot_Plant_Label)
    
    visited_plants <- sample(candidate_plants, total_visits_anim_plant$n[i], replace=TRUE)
    
    
    visited_plants2 <- as_tibble(as.data.frame(table(visited_plants))) %>% 
      rename(Subplot_Plant_Label=visited_plants,Visits_tot=Freq)
    
    visited_plants2$Subplot_Plant_Label <- as.character(visited_plants2$Subplot_Plant_Label)
    visited_plants2$Plot <- total_visits_anim_plant$Plot[i]
    visited_plants2$ID <- total_visits_anim_plant$ID[i]
    
    visited_plants2 <- visited_plants2 %>% dplyr::select(Plot,ID,Subplot_Plant_Label,Visits_tot)
    
    random_visit_list_week <- bind_rows(random_visit_list_week,visited_plants2)
    
  }
  
  return(random_visit_list_week)
  
}

#-------
# c = 1 - sum( (k.it/k.i)^2) # among-module connectivity = participation coefficient P in Guimerà & Amaral
# 
# z = (k.is - ks.bar) / SD.ks # within-module degree
# 
# k.is = number of links of i to other species in its own module s
# ks.bar = average k.is of all species in module s
# SD.ks = standard deviation of k.is of all species in module s
# k.it = number of links of species i to module t
# k.i = degree of species i
# 
# Note that for any species alone (in its level) in a module the z-value will be NaN,
# since then SD.ks is 0. This is a limitation of the way the z-value is defined
# (in multiples of degree/strength standard deviations).
# 
# Olesen et al. (2006) give critical c and z values of 0.62 and 2.6, respectively.
# Species exceeding these values are deemed connectors or hubs of a network.
# The justification of these thresholds remains unclear to me. They may also not apply
# for the quantitative version.

cz_values_function_IN <- function(edge_module_info_i){
  
  k.is <- edge_module_info_i %>% filter(module_to==module_from) %>% 
    group_by(node_to,layer_to,module_from,Plot) %>% count(wt=weight) %>% rename(k.is = n)
  ks.bar<- k.is %>% group_by(module_from,Plot) %>% summarise(ks.bar=mean(k.is,na.rm = T))
  SD.ks <- k.is %>% group_by(module_from,Plot) %>% summarise(SD.ks=sd(k.is,na.rm = T))
  
  k.it <- edge_module_info_i %>% filter(module_to!=module_from) %>% 
    group_by(node_to,layer_to,module_from,Plot) %>% count(wt=weight) %>% rename(k.it = n)
  
  k.i <- edge_module_info_i %>% group_by(node_to,layer_to,Plot) %>% count(wt=weight) %>% rename(k.i = n)
  
  z_edge_module_info_i <- k.is %>% 
    left_join(ks.bar,by=c("module_from","Plot")) %>% 
    left_join(SD.ks,by=c("module_from","Plot"))
  
  c_edge_module_info_i <- k.it %>% 
    left_join(k.i,by=c("node_to","layer_to","Plot")) %>%
    mutate(aux_c =  (k.it/k.i)^2)
  
  c <- c_edge_module_info_i %>% group_by(node_to,layer_to,Plot) %>% 
    count(wt=aux_c) %>% rename(sum_aux_c = n) %>% mutate(c=1-sum_aux_c)
  
  
  z_edge_module_info_i <- z_edge_module_info_i %>%
    mutate(z = (k.is - ks.bar) / SD.ks)
  
  cz_edge_module_info_i <- z_edge_module_info_i %>% 
    left_join(k.i,by=c("node_to","layer_to","Plot")) %>%
    left_join(c,by=c("node_to","layer_to","Plot"))
  
  cz_edge_module_info_i$sum_aux_c[is.na(cz_edge_module_info_i$sum_aux_c)] <- 0
  cz_edge_module_info_i$c[is.na(cz_edge_module_info_i$c)] <- 1
  
  return(cz_edge_module_info_i)
}

#-----
# c = 1 - sum( (k.it/k.i)^2) # among-module connectivity = participation coefficient P in Guimerà & Amaral
# 
# z = (k.is - ks.bar) / SD.ks # within-module degree
# 
# k.is = number of links of i to other species in its own module s
# ks.bar = average k.is of all species in module s
# SD.ks = standard deviation of k.is of all species in module s
# k.it = number of links of species i to module t
# k.i = degree of species i
# 
# Note that for any species alone (in its level) in a module the z-value will be NaN,
# since then SD.ks is 0. This is a limitation of the way the z-value is defined
# (in multiples of degree/strength standard deviations).
# 
# Olesen et al. (2006) give critical c and z values of 0.62 and 2.6, respectively.
# Species exceeding these values are deemed connectors or hubs of a network.
# The justification of these thresholds remains unclear to me. They may also not apply
# for the quantitative version.

cz_values_function_OUT <- function(edge_module_info_i){
  
  k.is <- edge_module_info_i %>% filter(module_to==module_from) %>% 
    group_by(node_from,layer_from,module_to,Plot) %>% count(wt=weight) %>% rename(k.is = n)
  ks.bar<- k.is %>% group_by(module_to,Plot) %>% summarise(ks.bar=mean(k.is,na.rm = T))
  SD.ks <- k.is %>% group_by(module_to,Plot) %>% summarise(SD.ks=sd(k.is,na.rm = T))
  
  k.it <- edge_module_info_i %>% filter(module_to!=module_from) %>% 
    group_by(node_from,layer_from,module_to,Plot) %>% count(wt=weight) %>% rename(k.it = n)
  
  k.i <- edge_module_info_i %>% group_by(node_from,layer_from,Plot) %>% count(wt=weight) %>% rename(k.i = n)
  
  z_edge_module_info_i <- k.is %>% 
    left_join(ks.bar,by=c("module_to","Plot")) %>% 
    left_join(SD.ks,by=c("module_to","Plot"))
  
  c_edge_module_info_i <- k.it %>% 
    left_join(k.i,by=c("node_from","layer_from","Plot")) %>%
    mutate(aux_c =  (k.it/k.i)^2)
  
  c <- c_edge_module_info_i %>% group_by(node_from,layer_from,Plot) %>% 
    count(wt=aux_c) %>% rename(sum_aux_c = n) %>% mutate(c=1-sum_aux_c)
  
  
  z_edge_module_info_i <- z_edge_module_info_i %>%
    mutate(z = (k.is - ks.bar) / SD.ks)
  
  cz_edge_module_info_i <- z_edge_module_info_i %>% 
    left_join(k.i,by=c("node_from","layer_from","Plot")) %>%
    left_join(c,by=c("node_from","layer_from","Plot"))
  
  cz_edge_module_info_i$sum_aux_c[is.na(cz_edge_module_info_i$sum_aux_c)] <- 0
  cz_edge_module_info_i$c[is.na(cz_edge_module_info_i$c)] <- 1
  return(cz_edge_module_info_i)
}


#----
load_data_models_2020_3 <- function(path_data_models_overlap="Processed_data/2020_NN_NEW_data_models_phenol_overlap_3.csv",
                                    path_raw_data_poll="Raw_Data/final_Pollinators_2020.csv"){
  
  fitness_final_aux <- read.csv(file = path_data_models_overlap,
                                header = TRUE,
                                stringsAsFactors = FALSE) %>%
    rename(Plant_Simple=Plant) %>% mutate(Seeds_GF = round(Seeds_GF))
  
  fitness_final_aux %>% group_by(Plant_Simple) %>% count()
  fitness_final_aux %>% filter(ID != "None") %>% group_by(Plant_Simple) %>% count()
  
  #########################
  # Add G_F
  
  G_F_list <- read_csv2(path_raw_data_poll) %>%
    filter(ID != "Tabanidae") %>%
    dplyr::select(G_F,ID_Simple) %>% unique() %>% rename(ID=ID_Simple)
  
  # Remove points from ID names
  G_F_list$ID <- sub("\\.", "", G_F_list$ID)
  
  G_F_list <- bind_rows(G_F_list,tibble(G_F="None",ID="None"))
  
  G_F_list <- unique(G_F_list)
  
  G_F_list$G_F %>% unique() %>% sort()
  
  # Sanity check
  G_F_list %>% group_by(ID) %>% count() %>% filter(n>1)
  
  
  fitness_orig1 <- fitness_final_aux %>% dplyr::left_join(G_F_list,by = "ID")
  
  
  fitness_orig <- fitness_orig1 %>% #filter(Plant_Simple %in% c("CHFU","LEMA","PUPA",
    #                          "CETE","CHMI","SPRU")) %>% 
    group_by(Plot,Subplot,Plant_Simple,type_seed_per_fruit) %>%
    summarise(Seeds_GF=mean(Seeds_GF),
              visits_GF=sum(visits_GF),
              homo_motif=sum(homo_motif),
              hete_motif=sum(hete_motif),
              Real_PR_Multi=mean(Real_PR_Multi),
              Real_PR_Layer=mean(Real_PR_Layer),
              StrengthIn=mean(StrengthIn),
              StrengthOut=mean(StrengthOut),
              DegreeIn=mean(DegreeIn),
              DegreeOut=mean(DegreeOut),
              Delta=mean(Delta),
              Ratio=mean(Ratio))
  
  
  # Turn ID, GF and Plot into factors
  fitness_orig$Plot <- as.factor(fitness_orig$Plot)
  fitness_orig$Plant_Simple <- as.factor(fitness_orig$Plant_Simple)
  #fitness_orig$ID <- as.factor(fitness_orig$ID)
  #fitness_orig$G_F <- as.factor(fitness_orig$G_F)
  
  return(fitness_orig)
  
  
}

load_data_models_2020_2 <- function(path_data_models_overlap="Processed_data/2020_NN_NEW_data_models_phenol_overlap_2.csv",
                           path_raw_data_poll="Raw_Data/final_Pollinators_2020.csv"){
  
  fitness_final_aux <- read.csv(file = path_data_models_overlap,
                                header = TRUE,
                                stringsAsFactors = FALSE) %>%
    rename(Plant_Simple=Plant) %>% mutate(Seeds_GF = round(Seeds_GF))
  
  fitness_final_aux %>% group_by(Plant_Simple) %>% count()
  fitness_final_aux %>% filter(ID != "None") %>% group_by(Plant_Simple) %>% count()
  
  #########################
  # Add G_F
  
  G_F_list <- read_csv2(path_raw_data_poll) %>%
    filter(ID != "Tabanidae") %>%
    dplyr::select(G_F,ID_Simple) %>% unique() %>% rename(ID=ID_Simple)
  
  # Remove points from ID names
  G_F_list$ID <- sub("\\.", "", G_F_list$ID)
  
  G_F_list <- bind_rows(G_F_list,tibble(G_F="None",ID="None"))
  
  G_F_list <- unique(G_F_list)
  
  G_F_list$G_F %>% unique() %>% sort()
  
  # Sanity check
  G_F_list %>% group_by(ID) %>% count() %>% filter(n>1)
  
  
  fitness_orig1 <- fitness_final_aux %>% dplyr::left_join(G_F_list,by = "ID")
  
  
  fitness_orig <- fitness_orig1 %>% #filter(Plant_Simple %in% c("CHFU","LEMA","PUPA",
    #                          "CETE","CHMI","SPRU")) %>% 
    group_by(Plot,Subplot,Plant_Simple,type_fruit,type_seed_per_fruit) %>%
    summarise(Seeds_GF=mean(Seeds_GF),
              Fruit_GF=mean(Fruit_GF),
              visits_GF=sum(visits_GF),
              homo_motif=sum(homo_motif),
              hete_motif=sum(hete_motif),
              Real_PR_Multi=mean(Real_PR_Multi),
              Real_PR_Layer=mean(Real_PR_Layer),
              StrengthIn=mean(StrengthIn),
              StrengthOut=mean(StrengthOut),
              DegreeIn=mean(DegreeIn),
              DegreeOut=mean(DegreeOut),
              Delta=mean(Delta),
              Ratio=mean(Ratio))
  
  
  # Turn ID, GF and Plot into factors
  fitness_orig$Plot <- as.factor(fitness_orig$Plot)
  fitness_orig$Plant_Simple <- as.factor(fitness_orig$Plant_Simple)
  #fitness_orig$ID <- as.factor(fitness_orig$ID)
  #fitness_orig$G_F <- as.factor(fitness_orig$G_F)
  
  return(fitness_orig)


}

load_data_models_2020 <- function(path_data_models_overlap="Processed_data/2020_NN_NEW_data_models_phenol_overlap.csv",
                                  path_raw_data_poll="Raw_Data/final_Pollinators_2020.csv"){
  
  fitness_final_aux <- read.csv(file = "Processed_data/2020_NN_NEW_data_models_phenol_overlap.csv",
                                header = TRUE,
                                stringsAsFactors = FALSE) %>%
    rename(Plant_Simple=Plant) %>% mutate(Seeds_GF = round(Seeds_GF))
  
  fitness_final_aux %>% group_by(Plant_Simple) %>% count()
  fitness_final_aux %>% filter(ID != "None") %>% group_by(Plant_Simple) %>% count()
  
  #########################
  # Add G_F
  
  G_F_list <- read_csv2("Raw_Data/final_Pollinators_2020.csv") %>%
    filter(ID != "Tabanidae") %>%
    dplyr::select(G_F,ID_Simple) %>% unique() %>% rename(ID=ID_Simple)
  
  # Remove points from ID names
  G_F_list$ID <- sub("\\.", "", G_F_list$ID)
  
  G_F_list <- bind_rows(G_F_list,tibble(G_F="None",ID="None"))
  
  G_F_list <- unique(G_F_list)
  
  G_F_list$G_F %>% unique() %>% sort()
  
  # Sanity check
  G_F_list %>% group_by(ID) %>% count() %>% filter(n>1)
  
  
  fitness_orig1 <- fitness_final_aux %>% dplyr::left_join(G_F_list,by = "ID")
  
  ########################
  # We should aggregate
  
  list_none <- fitness_orig1 %>% filter(ID == "None") %>%
    dplyr::select(Plot,Subplot,Plant_Simple)
  
  list_none %>% group_by(Plot,Subplot,Plant_Simple) %>% count() %>% filter(n>1)
  list_none$non_non=NA
  
  for(i in 1:nrow(list_none)){
    list_non_none_i <- fitness_orig1 %>% filter(Plot==list_none$Plot[i],
                                                Subplot==list_none$Subplot[i],
                                                Plant_Simple==list_none$Plant_Simple[i],
                                                ID != "None")
    list_none$non_non[i] <- nrow(list_non_none_i)
  }
  list_none %>% filter(non_non>1)
  
  fitness_orig <- fitness_orig1 %>% #filter(Plant_Simple %in% c("CHFU","LEMA","PUPA",
    #                          "CETE","CHMI","SPRU")) %>% 
    group_by(Plot,Subplot,Plant_Simple) %>%
    summarise(Seeds_GF=mean(Seeds_GF),
              Fruit_GF=mean(Fruit_GF),
              visits_GF=sum(visits_GF),
              homo_motif=sum(homo_motif),
              hete_motif=sum(hete_motif),
              Real_PR_Multi=mean(Real_PR_Multi),
              Real_PR_Layer=mean(Real_PR_Layer),
              StrengthIn=mean(StrengthIn),
              StrengthOut=mean(StrengthOut),
              DegreeIn=mean(DegreeIn),
              DegreeOut=mean(DegreeOut),
              Delta=mean(Delta),
              Ratio=mean(Ratio))
  
  
  # Turn ID, GF and Plot into factors
  fitness_orig$Plot <- as.factor(fitness_orig$Plot)
  fitness_orig$Plant_Simple <- as.factor(fitness_orig$Plant_Simple)
  #fitness_orig$ID <- as.factor(fitness_orig$ID)
  #fitness_orig$G_F <- as.factor(fitness_orig$G_F)
  
  return(fitness_orig)
  
  
}


# Data without aggregation

load_data_models_2020_without_agg <- function(path_data_models_overlap="Processed_data/2020_NN_NEW_data_models_phenol_overlap.csv",
                                  path_raw_data_poll="Raw_Data/final_Pollinators_2020.csv"){
  
  fitness_final_aux <- read.csv(file = "Processed_data/2020_NN_NEW_data_models_phenol_overlap.csv",
                                header = TRUE,
                                stringsAsFactors = FALSE) %>%
    rename(Plant_Simple=Plant) %>% mutate(Seeds_GF = round(Seeds_GF))
  
  fitness_final_aux %>% group_by(Plant_Simple) %>% count()
  fitness_final_aux %>% filter(ID != "None") %>% group_by(Plant_Simple) %>% count()
  
  #########################
  # Add G_F
  
  G_F_list <- read_csv2("Raw_Data/final_Pollinators_2020.csv") %>%
    filter(ID != "Tabanidae") %>%
    dplyr::select(G_F,ID_Simple) %>% unique() %>% rename(ID=ID_Simple)
  
  # Remove points from ID names
  G_F_list$ID <- sub("\\.", "", G_F_list$ID)
  
  G_F_list <- bind_rows(G_F_list,tibble(G_F="None",ID="None"))
  
  G_F_list <- unique(G_F_list)
  
  G_F_list$G_F %>% unique() %>% sort()
  
  # Sanity check
  G_F_list %>% group_by(ID) %>% count() %>% filter(n>1)
  
  
  fitness_orig1 <- fitness_final_aux %>% dplyr::left_join(G_F_list,by = "ID")

  return(fitness_orig1)
}
#############
# Corecctions for PageRank

corrections_pagerank <- function(){
  
  fitness_orig <- load_data_models_2020()
  
  fitness.data <- subset(fitness_orig,Seeds_GF > 0)
  
  ##################
  
  correction <- tibble(Plot=1:9)
  
  correction$PR_isolated <- NA
  correction$PR_extra_VS_PR <- NA
  
  for(Plot_i in 1:9){
    
    # Multilayer----
    graph_Plot_i_all_links <- readRDS(file = paste0("Processed_data/NN_networks/Plot_",Plot_i,"_NN_intra_inter.rds"))
    page_rank_i <- igraph::page_rank(graph_Plot_i_all_links,
                                     directed = TRUE, damping = 0.85,
                                     personalized = NULL, weights = NULL, options = NULL)
    
    page_rank_i_nodes_all_links <-
      data.frame(species=names(page_rank_i[[1]]),
                 Real_PR_Multi=page_rank_i[[1]], row.names=NULL)
    
    # Add extra nodes----
    
    # Number of extra nodes
    
    isolated_nodes <- fitness.data %>% filter(Plot==as.character(Plot_i),DegreeIn==0) %>%
      mutate(isolated_nodes = paste0(Plant_Simple," ",Subplot)) %>% ungroup() %>%
      dplyr::select(isolated_nodes) %>% nrow()
    
    graph_Plot_i_all_links_extra_nodes <- graph_Plot_i_all_links %>%
      igraph::add_vertices(isolated_nodes, color = "red")
    
    page_rank_i_extra_nodes <- igraph::page_rank(graph_Plot_i_all_links_extra_nodes,
                                                 directed = TRUE, damping = 0.85,
                                                 personalized = NULL, weights = NULL, options = NULL)
    
    page_rank_i_nodes_all_links_extra_nodes <-
      data.frame(species=names(page_rank_i_extra_nodes[[1]]),
                 Real_PR_Multi_extra=page_rank_i_extra_nodes[[1]], row.names=NULL)
    
    PR_isolated <- page_rank_i_nodes_all_links_extra_nodes$Real_PR_Multi_extra[
      is.na(page_rank_i_nodes_all_links_extra_nodes$species)] %>% unique()
    
    # Compare
    page_rank_i_nodes_all_links <- page_rank_i_nodes_all_links %>%
      left_join(page_rank_i_nodes_all_links_extra_nodes,by="species") %>%
      mutate(PR_extra_VS_PR = Real_PR_Multi_extra/Real_PR_Multi)
    
    correction$PR_isolated[Plot_i] <- PR_isolated
    correction$PR_extra_VS_PR[Plot_i] <- 
      unique(round(page_rank_i_nodes_all_links$PR_extra_VS_PR,8))
    
  }

  return(correction)
  
}

corrections_pagerank_2 <- function(){
  
  fitness_orig <- load_data_models_2020_2()
  
  fitness.data <- subset(fitness_orig,Seeds_GF > 0)
  
  ##################
  
  correction <- tibble(Plot=1:9)
  
  correction$PR_isolated <- NA
  correction$PR_extra_VS_PR <- NA
  
  for(Plot_i in 1:9){
    
    # Multilayer----
    graph_Plot_i_all_links <- readRDS(file = paste0("Processed_data/NN_networks/Plot_",Plot_i,"_NN_intra_inter.rds"))
    page_rank_i <- igraph::page_rank(graph_Plot_i_all_links,
                                     directed = TRUE, damping = 0.85,
                                     personalized = NULL, weights = NULL, options = NULL)
    
    page_rank_i_nodes_all_links <-
      data.frame(species=names(page_rank_i[[1]]),
                 Real_PR_Multi=page_rank_i[[1]], row.names=NULL)
    
    # Add extra nodes----
    
    # Number of extra nodes
    
    isolated_nodes <- fitness.data %>% filter(Plot==as.character(Plot_i),DegreeIn==0) %>%
      mutate(isolated_nodes = paste0(Plant_Simple," ",Subplot)) %>% ungroup() %>%
      dplyr::select(isolated_nodes) %>% nrow()
    
    graph_Plot_i_all_links_extra_nodes <- graph_Plot_i_all_links %>%
      igraph::add_vertices(isolated_nodes, color = "red")
    
    page_rank_i_extra_nodes <- igraph::page_rank(graph_Plot_i_all_links_extra_nodes,
                                                 directed = TRUE, damping = 0.85,
                                                 personalized = NULL, weights = NULL, options = NULL)
    
    page_rank_i_nodes_all_links_extra_nodes <-
      data.frame(species=names(page_rank_i_extra_nodes[[1]]),
                 Real_PR_Multi_extra=page_rank_i_extra_nodes[[1]], row.names=NULL)
    
    PR_isolated <- page_rank_i_nodes_all_links_extra_nodes$Real_PR_Multi_extra[
      is.na(page_rank_i_nodes_all_links_extra_nodes$species)] %>% unique()
    
    # Compare
    page_rank_i_nodes_all_links <- page_rank_i_nodes_all_links %>%
      left_join(page_rank_i_nodes_all_links_extra_nodes,by="species") %>%
      mutate(PR_extra_VS_PR = Real_PR_Multi_extra/Real_PR_Multi)
    
    correction$PR_isolated[Plot_i] <- PR_isolated
    correction$PR_extra_VS_PR[Plot_i] <- 
      unique(round(page_rank_i_nodes_all_links$PR_extra_VS_PR,8))
    
  }
  
  return(correction)
  
}

