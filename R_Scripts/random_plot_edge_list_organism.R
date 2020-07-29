# This function generates a ramdom list of visitis between animals and plants
# The IDs of animals are randomly selected and the probability of being selected for a given org. is
# proportional to the total number of visits of such organism.
# The total amount of links per plant individual is equal to that of the original list

random_plot_edge_list_organism <- function(plot_edge_list){

  prob_visits_anim <- plot_edge_list %>% group_by(to) %>% count(wt=weight) %>%
    mutate(percentage = n/sum(plot_edge_list$weight))
  
  focal_links <- plot_edge_list %>% group_by(from,species) %>% count(wt=weight)
  
  random_plot_edge_list <- NULL
  
  for(i in 1:nrow(focal_links)){
    
    visitors <- sample(prob_visits_anim$to,size=focal_links$n[i],replace = T,prob = prob_visits_anim$percentage)
    
    visitors2 <- as_tibble(as.data.frame(table(visitors))) %>% rename(to=visitors,weight=Freq)
    
    visitors2$to <- as.character(visitors2$to)
    visitors2$from <- focal_links$from[i]
    visitors2$species <- focal_links$species[i]
    
    visitors2 <- visitors2 %>% dplyr::select(from,to,weight,species)
    
    random_plot_edge_list <- bind_rows(random_plot_edge_list,visitors2)
    
  }

  return(random_plot_edge_list)
  
}