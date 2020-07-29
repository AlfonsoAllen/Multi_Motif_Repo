# This function generates a ramdom list of visitis between animals and plants
# The IDs of animals are randomly selected and the probability of being selected for a given org. is
# proportional to the total number of visits of such organism.
# The total amount of links per plant individual is equal to that of the original list

random_plot_edge_list_plant_organism2 <- function(plot_edge_list){
  
  
  prob_visits_anim <- plot_edge_list %>% group_by(to) %>% count(wt=weight) %>%
    mutate(percentage = n/sum(plot_edge_list$weight))
  
  prob_visits_plants <- plot_edge_list %>%
    group_by(from) %>% count(wt=weight) %>%
    mutate(percentage = n/sum(plot_edge_list$weight))
  
  
  total_visits <- sum(plot_edge_list$weight)

  visitors <- sample(prob_visits_anim$to,size=total_visits,replace = T)#,prob = prob_visits_anim$percentage)
  plants <- sample(prob_visits_plants$from,size=total_visits,replace = T)#,prob = prob_visits_plants$percentage)  
  
  plants_aux <- tibble(from=plants) %>% mutate(new_c=from) %>%
    separate(new_c,c("sub","species")," ") %>% 
    select(from,species)
  
  while (length(unique(plants_aux$species))<2){
    
    plants <- sample(prob_visits_plants$from,size=total_visits,replace = T)#,prob = prob_visits_plants$percentage)  
    
    plants_aux <- tibble(from=plants) %>% mutate(new_c=from) %>%
      separate(new_c,c("sub","species")," ") %>% 
      select(from,species)
  }
  
    
  random_plot_edge_list_aux <- tibble(from=plants,to=visitors)
  random_plot_edge_list_aux$weight <- 1
    
  random_plot_edge_list_aux <- random_plot_edge_list_aux %>% group_by(from,to) %>% count(wt=weight) %>%
      rename(weight=n)
    
  random_plot_edge_list <- random_plot_edge_list_aux %>% mutate(new_c=from) %>% 
    separate(new_c,c("sub","species")," ") %>% 
    select(from,to,weight,species)

  return(random_plot_edge_list)
  
}