

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

