
rand_vect <- function(N, M, sd = 2, pos.only = TRUE) {
  vec <- rnorm(N, M/N, sd)
  if (abs(sum(vec)) < 0.01) vec <- vec + 1
  vec <- round(vec / sum(vec) * M)
  deviation <- M - sum(vec)
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
  }
  if (pos.only) while (any(vec < 0)) {
    negs <- vec < 0
    pos  <- vec > 0
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
  }
  vec
}

intralinks_index <- function(number_ind_plants,number_polinator_sp,number_intralinks){
  
  row_index_aux <- c(sample(1:number_ind_plants),
                     sample(1:number_ind_plants, number_intralinks-number_ind_plants, replace=T))
  
  row_index <- sample(row_index_aux)
  
  col_index_aux <- c(sample(1:number_polinator_sp),
                     sample(1:number_polinator_sp, number_intralinks-number_polinator_sp, replace=T))
  
  col_index <- sample(col_index_aux)
  
  
  df_index <- tibble(row_index = row_index,col_index=col_index)
  
  return(df_index)
  
}

intralink_generator <- function(number_ind_plants,
                                number_polinator_sp,
                                number_intralinks){
  
  possible_intralinks_indeces <- intralinks_index(number_ind_plants,
                                                  number_polinator_sp,
                                                  number_intralinks)
  
  df_index_repeated <- possible_intralinks_indeces[duplicated(possible_intralinks_indeces),]
  df_index_repeated
  
  while(nrow(df_index_repeated)>0){
    possible_intralinks_indeces <- intralinks_index(number_ind_plants,
                                                    number_polinator_sp,
                                                    number_intralinks)
    
    df_index_repeated <- possible_intralinks_indeces[duplicated(possible_intralinks_indeces),]
  }
  
  interaction_matrix <- matrix(rep(0,number_ind_plants*number_polinator_sp),
                               nrow = number_ind_plants,ncol = number_polinator_sp)
  
  for(i in 1:nrow(possible_intralinks_indeces)){
    
    interaction_matrix[possible_intralinks_indeces$row_index[i],
                       possible_intralinks_indeces$col_index[i]] <- 1
  }
  
  return(interaction_matrix)
  
}

interlink_corrector <- function(number_plant_individuals,
                                number_plant_individuals_per_plant_sp,
                                number_pollinator_sp,
                                number_intralinks){
  
  interaction_matrix <- intralink_generator(number_plant_individuals,
                                            number_pollinator_sp,number_intralinks)
  
  # Sanity check
  sum(interaction_matrix>0)==number_intralinks
  
  # assign names to species and individuals in our interaction matrix
  
  plant_individual_names <- c(paste0("plant_sp_1 ",
                                     paste0("ind_",1:number_plant_individuals_per_plant_sp[1])),
                              paste0("plant_sp_2 ",
                                     paste0("ind_",1:number_plant_individuals_per_plant_sp[2])),
                              paste0("plant_sp_3 ",
                                     paste0("ind_",1:number_plant_individuals_per_plant_sp[3])),
                              paste0("plant_sp_4 ",
                                     paste0("ind_",1:number_plant_individuals_per_plant_sp[4])),
                              paste0("plant_sp_5 ",
                                     paste0("ind_",1:number_plant_individuals_per_plant_sp[5])))
  
  rownames(interaction_matrix) <- plant_individual_names
  poll_names <- paste0("pol_sp_",1:number_pollinator_sp)
  colnames(interaction_matrix) <- poll_names
  
  # adjust_interlinks
  
  plant_sp_1 <- 1:number_plant_individuals_per_plant_sp[1]
  plant_sp_2 <- max(plant_sp_1)+1:number_plant_individuals_per_plant_sp[2]
  plant_sp_3 <- max(plant_sp_2)+1:number_plant_individuals_per_plant_sp[3]
  plant_sp_4 <- max(plant_sp_3)+1:number_plant_individuals_per_plant_sp[4]
  plant_sp_5 <- max(plant_sp_4)+1:number_plant_individuals_per_plant_sp[5]
  
  poll_plant_sp_1_aux <- interaction_matrix[plant_sp_1,]
  poll_plant_sp_2_aux <- interaction_matrix[plant_sp_2,]
  poll_plant_sp_3_aux <- interaction_matrix[plant_sp_3,]
  poll_plant_sp_4_aux <- interaction_matrix[plant_sp_4,]
  poll_plant_sp_5_aux <- interaction_matrix[plant_sp_5,]
  
  poll_plant_sp_1 <- poll_names[colSums(poll_plant_sp_1_aux)>0]
  poll_plant_sp_2 <- poll_names[colSums(poll_plant_sp_2_aux)>0]
  poll_plant_sp_3 <- poll_names[colSums(poll_plant_sp_3_aux)>0]
  poll_plant_sp_4 <- poll_names[colSums(poll_plant_sp_4_aux)>0]
  poll_plant_sp_5 <- poll_names[colSums(poll_plant_sp_5_aux)>0]
  
  poll_plant_sp_12 <- poll_plant_sp_1[poll_plant_sp_1 %in% poll_plant_sp_2]
  poll_plant_sp_13 <- poll_plant_sp_1[poll_plant_sp_1 %in% poll_plant_sp_3]
  poll_plant_sp_14 <- poll_plant_sp_1[poll_plant_sp_1 %in% poll_plant_sp_4]
  poll_plant_sp_15 <- poll_plant_sp_1[poll_plant_sp_1 %in% poll_plant_sp_5]
  poll_plant_sp_23 <- poll_plant_sp_2[poll_plant_sp_2 %in% poll_plant_sp_3]
  poll_plant_sp_24 <- poll_plant_sp_2[poll_plant_sp_2 %in% poll_plant_sp_4]
  poll_plant_sp_25 <- poll_plant_sp_2[poll_plant_sp_2 %in% poll_plant_sp_5]
  poll_plant_sp_34 <- poll_plant_sp_3[poll_plant_sp_3 %in% poll_plant_sp_4]
  poll_plant_sp_35 <- poll_plant_sp_3[poll_plant_sp_3 %in% poll_plant_sp_5]
  poll_plant_sp_45 <- poll_plant_sp_4[poll_plant_sp_4 %in% poll_plant_sp_5]
  
  interlinks <- length(poll_plant_sp_12)+length(poll_plant_sp_13)+length(poll_plant_sp_14)+
    length(poll_plant_sp_15)+length(poll_plant_sp_23)+length(poll_plant_sp_24)+
    length(poll_plant_sp_25)+length(poll_plant_sp_34)+length(poll_plant_sp_35)+
    length(poll_plant_sp_45)
  
  return(list(interaction_matrix,interlinks))
  
}

interaction_matrix_link_generator <- function(number_plant_individuals,
                                              number_plant_individuals_per_plant_sp,
                                              number_pollinator_sp,
                                              number_intralinks){
  
  interaction_matrix_data <- interlink_corrector(number_plant_individuals,
                                                 number_plant_individuals_per_plant_sp,
                                                 number_pollinator_sp,
                                                 number_intralinks)
  
  interaction_matrix <- interaction_matrix_data[[1]]
  interlinks <- interaction_matrix_data[[2]]
  
  while(interlinks != number_interlinks){
    
    interaction_matrix_data <- interlink_corrector(number_plant_individuals,
                                                   number_plant_individuals_per_plant_sp,
                                                   number_pollinator_sp,
                                                   number_intralinks)
    
    interaction_matrix <- interaction_matrix_data[[1]]
    interlinks <- interaction_matrix_data[[2]]
    
  }
  
  return(interaction_matrix)
  
}

weighted_interaction_matrix_generator <- function(number_plant_individuals,
                                                  number_plant_individuals_per_plant_sp,
                                                  number_pollinator_sp,
                                                  number_intralinks,
                                                  number_interlinks,
                                                  average_intralink_strenght){
  
  interaction_matrix_links <- interaction_matrix_link_generator(number_plant_individuals,
                                                                number_plant_individuals_per_plant_sp,
                                                                number_pollinator_sp,
                                                                number_intralinks)
  
  #Sanity_check
  # cat(paste0("the number of intralinks in link matrix is correct: ",
  #            sum(interaction_matrix_links>0)==number_intralinks,"\n"))
  # 
  # binary_interlinks <- interlink_tester(interaction_matrix_links,
  #                                         number_plant_individuals_per_plant_sp)
  # 
  # cat(paste0("the number of intralinks in link matrix is correct: ",
  #            binary_interlinks == number_interlinks,"\n"))
  
  
  # We add strength
  total_instrength <- number_intralinks * average_intralink_strenght
  
  instrength_plant_individuals <- rand_vect(number_plant_individuals,
                                            total_instrength,sd = 3)
  
  # The instrength of all plant individuals must be greather than zero
  while(sum(instrength_plant_individuals==0)>0){
    
    instrength_plant_individuals <- rand_vect(number_plant_individuals,
                                              total_instrength,sd = 3)
  }
  
  # The degree of all plant individuals must be greather than their instrength
  degree_plant_individuals <- rowSums(interaction_matrix_links) %>% as.numeric()
  
  while(sum(degree_plant_individuals > instrength_plant_individuals)>0){
    instrength_plant_individuals <- sample(instrength_plant_individuals)
  }
  
  # Sanity checks
  # cat(paste0("instrength_plant_ind >= degree: ",
  #            sum(degree_plant_individuals > instrength_plant_individuals)==0,
  #            "\n"))
  
  weighted_interaction_matrix <- interaction_matrix_links
  
  for(i in 1:nrow(weighted_interaction_matrix)){
    
    plant_ind_i_pollinators <- which(weighted_interaction_matrix[i,] > 0) %>% as.numeric()
    
    if(degree_plant_individuals[i]==1){
      
      weighted_interaction_matrix[i,plant_ind_i_pollinators] <- instrength_plant_individuals[i]
      
      #cat(paste0(i,": ",instrength_plant_individuals[i],
      #           "; pollinators: ",plant_ind_i_pollinators,"\n"))
      
    }else if(degree_plant_individuals[i]==instrength_plant_individuals[i]){
      
      for(i.pol in plant_ind_i_pollinators){
        
        weighted_interaction_matrix[i,i.pol] <- 1
        
      }
      
    }else{
      
      weights <- rand_vect(degree_plant_individuals[i],instrength_plant_individuals[i])
      while(sum(weights==0)>0){
        weights <- rand_vect(degree_plant_individuals[i],instrength_plant_individuals[i])
      }
      
      #cat(paste0(i,": ",weights,"; pollinators: ",plant_ind_i_pollinators,"\n"))
      
      for(i.pol in 1:length(plant_ind_i_pollinators)){
        weighted_interaction_matrix[i,plant_ind_i_pollinators[i.pol]] <- weights[i.pol]
      }
      
      
    }
    
  }
  
  #Sanity_checks
  
  # weighted_interlinks <- interlink_tester(weighted_interaction_matrix,
  #                                     number_plant_individuals_per_plant_sp)
  # 
  # cat(paste0("the number of interlinks in weighted interaction matrix is correct: ",
  #            weighted_interlinks == number_interlinks,"\n"))
  # 
  # cat(paste0("the number of intralinks in weighted interaction matrix is correct: ",
  #            sum(weighted_interaction_matrix>0) == number_intralinks,"\n"))
  
  return(weighted_interaction_matrix)
  
}

interlink_tester <- function(weighted_interaction_matrix,
                             number_plant_individuals_per_plant_sp){
  
  # assign names to species and individuals in our interaction matrix
  
  poll_names <- colnames(weighted_interaction_matrix)
  
  # adjust_interlinks
  
  plant_sp_1 <- 1:number_plant_individuals_per_plant_sp[1]
  plant_sp_2 <- max(plant_sp_1)+1:number_plant_individuals_per_plant_sp[2]
  plant_sp_3 <- max(plant_sp_2)+1:number_plant_individuals_per_plant_sp[3]
  plant_sp_4 <- max(plant_sp_3)+1:number_plant_individuals_per_plant_sp[4]
  plant_sp_5 <- max(plant_sp_4)+1:number_plant_individuals_per_plant_sp[5]
  
  poll_plant_sp_1_aux <- weighted_interaction_matrix[plant_sp_1,]
  poll_plant_sp_2_aux <- weighted_interaction_matrix[plant_sp_2,]
  poll_plant_sp_3_aux <- weighted_interaction_matrix[plant_sp_3,]
  poll_plant_sp_4_aux <- weighted_interaction_matrix[plant_sp_4,]
  poll_plant_sp_5_aux <- weighted_interaction_matrix[plant_sp_5,]
  
  poll_plant_sp_1 <- poll_names[colSums(poll_plant_sp_1_aux)>0]
  poll_plant_sp_2 <- poll_names[colSums(poll_plant_sp_2_aux)>0]
  poll_plant_sp_3 <- poll_names[colSums(poll_plant_sp_3_aux)>0]
  poll_plant_sp_4 <- poll_names[colSums(poll_plant_sp_4_aux)>0]
  poll_plant_sp_5 <- poll_names[colSums(poll_plant_sp_5_aux)>0]
  
  poll_plant_sp_12 <- poll_plant_sp_1[poll_plant_sp_1 %in% poll_plant_sp_2]
  poll_plant_sp_13 <- poll_plant_sp_1[poll_plant_sp_1 %in% poll_plant_sp_3]
  poll_plant_sp_14 <- poll_plant_sp_1[poll_plant_sp_1 %in% poll_plant_sp_4]
  poll_plant_sp_15 <- poll_plant_sp_1[poll_plant_sp_1 %in% poll_plant_sp_5]
  poll_plant_sp_23 <- poll_plant_sp_2[poll_plant_sp_2 %in% poll_plant_sp_3]
  poll_plant_sp_24 <- poll_plant_sp_2[poll_plant_sp_2 %in% poll_plant_sp_4]
  poll_plant_sp_25 <- poll_plant_sp_2[poll_plant_sp_2 %in% poll_plant_sp_5]
  poll_plant_sp_34 <- poll_plant_sp_3[poll_plant_sp_3 %in% poll_plant_sp_4]
  poll_plant_sp_35 <- poll_plant_sp_3[poll_plant_sp_3 %in% poll_plant_sp_5]
  poll_plant_sp_45 <- poll_plant_sp_4[poll_plant_sp_4 %in% poll_plant_sp_5]
  
  interlinks <- length(poll_plant_sp_12)+length(poll_plant_sp_13)+length(poll_plant_sp_14)+
    length(poll_plant_sp_15)+length(poll_plant_sp_23)+length(poll_plant_sp_24)+
    length(poll_plant_sp_25)+length(poll_plant_sp_34)+length(poll_plant_sp_35)+
    length(poll_plant_sp_45)
  
  return(interlinks)
  
}

generate_interaction_list_without_phen <- function(number_plant_individuals,
                                                   number_plant_individuals_per_plant_sp,
                                                   number_pollinator_sp,
                                                   number_intralinks,
                                                   number_interlinks,
                                                   average_intralink_strenght){
  
  number_pollinator_sp <- number_intralinks - number_interlinks
  
  number_plant_individuals_per_plant_sp <- rand_vect(number_plant_sp,
                                                     number_plant_individuals,sd = 3)
  
  
  weighted_interaction_matrix <- weighted_interaction_matrix_generator(number_plant_individuals,
                                                                       number_plant_individuals_per_plant_sp,
                                                                       number_pollinator_sp,
                                                                       number_intralinks,
                                                                       number_interlinks,
                                                                       average_intralink_strenght)
  
  real_interlinks <- interlink_tester(weighted_interaction_matrix,
                                      number_plant_individuals_per_plant_sp)
  # Sanity check
  test1 <- all(sum(rowSums(weighted_interaction_matrix))==number_intralinks*average_intralink_strenght)
  test2 <- sum(weighted_interaction_matrix > 0) == number_intralinks
  test3 <- number_interlinks == real_interlinks
  
  cat(paste0("All sanity checks OK: ",all(test1,test2,test3),"\n"))
  
  # Turn matrix into a list of interactions
  
  list_interactions <- NULL
  
  plant_ind_codes <- rownames(weighted_interaction_matrix)
  pollinator_codes <- colnames(weighted_interaction_matrix)
  
  for(i in 1:nrow(weighted_interaction_matrix)){
    
    ind_plant_pollinators <- pollinator_codes[which(weighted_interaction_matrix[i,]>0)]
    ind_plant_int_weights <- weighted_interaction_matrix[i,which(weighted_interaction_matrix[i,]>0)]
    
    list_interactions_ind_i <- tibble(Plant=plant_ind_codes[i],
                                      ID=ind_plant_pollinators,
                                      Visits=ind_plant_int_weights) %>%
      separate(Plant,c("Plant","Subplot")," ")
    
    list_interactions <- bind_rows(list_interactions,list_interactions_ind_i)
    
  }
  
  return(list_interactions)
  
}

aggregate_phen2interaction_list <- function(interaction_list_plant_sp_i_aux, 
                                            mean_phen,
                                            sd_phen){
  
  interactions_with_phen <- NULL
  
  for(i in 1:nrow(interaction_list_plant_sp_i_aux)){
    
    visits_i <- interaction_list_plant_sp_i_aux$Visits[i]
    int_ind_i <- interaction_list_plant_sp_i_aux[i,]
    int_ind_i$Visits <- 1
    
    int_ind_i_phen <- NULL
    
    for(j in 1:visits_i){
      
      int_ind_i_phen <- bind_rows(int_ind_i_phen,int_ind_i)
    }
    
    int_ind_i_phen$Week <- rnorm(visits_i,mean_phen,sd_phen) %>% round(digits = 0)
    
    interactions_with_phen <- bind_rows(interactions_with_phen,int_ind_i_phen)
    
  }
  
  interactions_with_phen_result <- interactions_with_phen %>% group_by(Plant,Subplot,ID,Week) %>%
    count() %>% rename(Visits = n)
  
  return(interactions_with_phen_result)
  
}

generate_interaction_list_with_phen <- function(number_plant_individuals,
                                                number_plant_individuals_per_plant_sp,
                                                number_pollinator_sp, number_intralinks,
                                                number_interlinks,average_intralink_strenght,
                                                mean_phenology_plant_sp_1,
                                                sd_phenology_plant_sp_1,
                                                mean_phenology_plant_sp_2,
                                                sd_phenology_plant_sp_2,
                                                mean_phenology_plant_sp_3,
                                                sd_phenology_plant_sp_3,
                                                mean_phenology_plant_sp_4,
                                                sd_phenology_plant_sp_4,
                                                mean_phenology_plant_sp_5,
                                                sd_phenology_plant_sp_5){
  
  
  if(!number_plant_individuals %in% (number_intralinks - number_interlinks):number_intralinks){

    cat(paste0("The total number of plant individuals should be between ",
               number_intralinks - number_interlinks," and ",number_intralinks,"\n"))

  }else{
    
    interaction_list_without_phen <- generate_interaction_list_without_phen(number_plant_individuals,
                                                                            number_plant_individuals_per_plant_sp,
                                                                            number_pollinator_sp,
                                                                            number_intralinks,
                                                                            number_interlinks,
                                                                            average_intralink_strenght)
    
    interaction_list_plant_sp_1_aux <- interaction_list_without_phen %>% filter(Plant == "plant_sp_1")
    interaction_list_plant_sp_2_aux <- interaction_list_without_phen %>% filter(Plant == "plant_sp_2")
    interaction_list_plant_sp_3_aux <- interaction_list_without_phen %>% filter(Plant == "plant_sp_3")
    interaction_list_plant_sp_4_aux <- interaction_list_without_phen %>% filter(Plant == "plant_sp_4")
    interaction_list_plant_sp_5_aux <- interaction_list_without_phen %>% filter(Plant == "plant_sp_5")
    
    
    interaction_plant_sp_1 <- aggregate_phen2interaction_list(interaction_list_plant_sp_1_aux, 
                                                              mean_phenology_plant_sp_1,
                                                              sd_phenology_plant_sp_1)
    
    interaction_plant_sp_2 <- aggregate_phen2interaction_list(interaction_list_plant_sp_2_aux, 
                                                              mean_phenology_plant_sp_2,
                                                              sd_phenology_plant_sp_2)
    
    interaction_plant_sp_3 <- aggregate_phen2interaction_list(interaction_list_plant_sp_3_aux, 
                                                              mean_phenology_plant_sp_3,
                                                              sd_phenology_plant_sp_3)
    
    interaction_plant_sp_4 <- aggregate_phen2interaction_list(interaction_list_plant_sp_4_aux, 
                                                              mean_phenology_plant_sp_4,
                                                              sd_phenology_plant_sp_4)
    
    interaction_plant_sp_5 <- aggregate_phen2interaction_list(interaction_list_plant_sp_5_aux, 
                                                              mean_phenology_plant_sp_5,
                                                              sd_phenology_plant_sp_5)
    
    interaction_list_with_phen <- bind_rows(interaction_plant_sp_1,interaction_plant_sp_2,
                                            interaction_plant_sp_3,interaction_plant_sp_4,
                                            interaction_plant_sp_5)
    
    interaction_list_with_phen <- interaction_list_with_phen[,c("Subplot","Plant","ID",
                                                                "Visits","Week")]
    
    return(interaction_list_with_phen)
    
  }
  
  
}


