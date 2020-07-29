# This function generates a ramdom list of visitis between animals and plants
# The IDs of animals are randomly selected and the probability of being selected for a given org. is
# proportional to the total number of visits of such organism.
# The total amount of links per plant individual is equal to that of the original list

random_plot_edge_list_plant_organism <- function(plot_edge_list){
  
  visits_plant <- plot_edge_list %>% group_by(species) %>% count(wt=weight)
  
  random_plot_edge_list <- NULL
  
  for(i in 1:nrow(visits_plant)){
  
    plot_edge_list_fil <- plot_edge_list %>% filter(species==visits_plant$species[i])
    
    prob_visits_anim <- plot_edge_list_fil %>% group_by(to) %>% count(wt=weight) %>%
      mutate(percentage = n/sum(plot_edge_list$weight))
    
    prob_visits_plants <- plot_edge_list_fil %>%
      group_by(from) %>% count(wt=weight) %>%
      mutate(percentage = n/sum(plot_edge_list_fil$weight))
    
      
    visitors <- sample(prob_visits_anim$to,size=visits_plant$n[i],replace = T,prob = prob_visits_anim$percentage)
    plants <- sample(prob_visits_plants$from,size=visits_plant$n[i],replace = T,prob = prob_visits_plants$percentage)  
    
    random_plot_edge_list_aux <- tibble(from=plants,to=visitors)
    random_plot_edge_list_aux$weight <- 1
    
    random_plot_edge_list_aux <- random_plot_edge_list_aux %>% group_by(from,to) %>% count(wt=weight) %>%
      rename(weight=n)
    
    random_plot_edge_list_aux <- random_plot_edge_list_aux %>% mutate(new_c=from) %>% 
      separate(new_c,c("sub","species")," ") %>% 
      select(from,to,weight,species)
    
    random_plot_edge_list <- bind_rows(random_plot_edge_list,random_plot_edge_list_aux)
    
    }

  return(random_plot_edge_list)
  
}