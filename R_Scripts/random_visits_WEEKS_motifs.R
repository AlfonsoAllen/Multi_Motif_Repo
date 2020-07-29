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