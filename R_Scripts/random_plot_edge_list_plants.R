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