
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

generate_plant_pollinator_interlink_matrix <- function(number_interlinks,number_pollinator_sp){
  
  # Possible pairs Plant 12, 13, 14, 15, 23, 24, 25, 34, 35, 45
  
  number_interlinks_per_pollinator_sp <- rand_vect(number_pollinator_sp, number_interlinks)
  
  possible_number_interlinks_per_pollinator_sp <- c(0,1)
  
  for(x in 2:50){
    possible_number_interlinks_per_pollinator_sp <- c(possible_number_interlinks_per_pollinator_sp,
                                                      ncol(combn(1:x,2)))
  }
  
  while(!all(number_interlinks_per_pollinator_sp %in% 
             possible_number_interlinks_per_pollinator_sp)){
    
    number_interlinks_per_pollinator_sp <- rand_vect(number_pollinator_sp, number_interlinks)
  
  }
  
  name_plant_sp <- paste0("plant_sp_",1:5)
  name_poll_sp <- paste0("pol_sp_",1:number_pollinator_sp)
  
  plant_pollinator_interlink_matrix <- matrix(rep(0,5*number_pollinator_sp),nrow = 5,
                                              ncol = number_pollinator_sp)
  rownames(plant_pollinator_interlink_matrix) <- name_plant_sp
  colnames(plant_pollinator_interlink_matrix) <- name_poll_sp
  
  for(pol_i in 1:number_pollinator_sp){
    
    if(number_interlinks_per_pollinator_sp[pol_i]==1){
      
      plant_sp_interlinks_pol_i <- sample(1:5,size = 2)
      plant_pollinator_interlink_matrix[plant_sp_interlinks_pol_i,pol_i] <- 1
      
    }else if(number_interlinks_per_pollinator_sp[pol_i]==3){
      
      plant_sp_interlinks_pol_i <- sample(1:5,size = 3)
      plant_pollinator_interlink_matrix[plant_sp_interlinks_pol_i,pol_i] <- 1
      
    }else if(number_interlinks_per_pollinator_sp[pol_i]==6){
      
      plant_sp_interlinks_pol_i <- sample(1:5,size = 4)
      plant_pollinator_interlink_matrix[plant_sp_interlinks_pol_i,pol_i] <- 1
      
    }else if(number_interlinks_per_pollinator_sp[pol_i]==10){
      
      plant_pollinator_interlink_matrix[,pol_i] <- 1
    }
  }
  
  return(plant_pollinator_interlink_matrix)
  
}

plant_pollinator_interlink_tester <- function(plant_pollinator_interlink_matrix,
                                              number_pollinator_sp){
  
  # assign names to species and individuals in our interaction matrix
  
  poll_names <- colnames(plant_pollinator_interlink_matrix)
  
  poll_plant_sp_1_aux <- plant_pollinator_interlink_matrix[1,]
  poll_plant_sp_2_aux <- plant_pollinator_interlink_matrix[2,]
  poll_plant_sp_3_aux <- plant_pollinator_interlink_matrix[3,]
  poll_plant_sp_4_aux <- plant_pollinator_interlink_matrix[4,]
  poll_plant_sp_5_aux <- plant_pollinator_interlink_matrix[5,]
  
  poll_plant_sp_1 <- poll_names[poll_plant_sp_1_aux>0]
  poll_plant_sp_2 <- poll_names[poll_plant_sp_2_aux>0]
  poll_plant_sp_3 <- poll_names[poll_plant_sp_3_aux>0]
  poll_plant_sp_4 <- poll_names[poll_plant_sp_4_aux>0]
  poll_plant_sp_5 <- poll_names[poll_plant_sp_5_aux>0]
  
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

generate_plant_individual_names <- function(number_plant_individuals_per_plant_sp){
  
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
  
  return(plant_individual_names)
  
}

generate_empty_interaction_matrix <- function(number_plant_individuals,number_pollinator_sp,
                                              number_plant_individuals_per_plant_sp){
  
  # Create an interaction matrix
  
  interaction_matrix <- matrix(rep(0,number_plant_individuals*number_pollinator_sp),
                               nrow = number_plant_individuals, ncol = number_pollinator_sp)
  
  # assign names to species and individuals in our interaction matrix
  
  plant_individual_names <- generate_plant_individual_names(number_plant_individuals_per_plant_sp)
  rownames(interaction_matrix) <- plant_individual_names
  poll_names <- paste0("pol_sp_",1:number_pollinator_sp)
  colnames(interaction_matrix) <- poll_names
  
  return(interaction_matrix)
  
}

generate_number_plant_individuals_per_plant_sp <- function(number_plant_individuals){
  
  number_plant_individuals_per_plant_sp <- rand_vect(5,number_plant_individuals)
  
  while(!all(number_plant_individuals_per_plant_sp>0)){
    number_plant_individuals_per_plant_sp <- rand_vect(5,number_plant_individuals)
  }
  
  return(number_plant_individuals_per_plant_sp)
  
}

generate_number_intralinks_per_pollinator <- function(number_pollinator_sp, number_intralinks){
  
  number_intralinks_per_pollinator <- rand_vect(number_pollinator_sp,number_intralinks)
  
  while(!all(number_intralinks_per_pollinator>0)){
    number_intralinks_per_pollinator <- rand_vect(number_pollinator_sp,number_intralinks)
  }
  
  return(number_intralinks_per_pollinator)
  
}

generate_list_plant_individuals <- function(number_plant_individuals_per_plant_sp){
  
  plant_sp_1 <- 1:number_plant_individuals_per_plant_sp[1]
  plant_sp_2 <- max(plant_sp_1)+1:number_plant_individuals_per_plant_sp[2]
  plant_sp_3 <- max(plant_sp_2)+1:number_plant_individuals_per_plant_sp[3]
  plant_sp_4 <- max(plant_sp_3)+1:number_plant_individuals_per_plant_sp[4]
  plant_sp_5 <- max(plant_sp_4)+1:number_plant_individuals_per_plant_sp[5]
  
  list_plant_individuals <- list(plant_sp_1,plant_sp_2,plant_sp_3,plant_sp_4,plant_sp_5)
  
  return(list_plant_individuals)
  
}

init_interaction_matrix <- function(number_pollinator_sp,
                                     number_plant_individuals,
                                     number_intralinks,
                                     number_interlinks){
  
  plant_pollinator_interlink_matrix <- generate_plant_pollinator_interlink_matrix(number_interlinks,
                                                                                  number_pollinator_sp)
  
  # Sanity check
  # number_interlinks == plant_pollinator_interlink_tester(plant_pollinator_interlink_matrix,
  #                                                        number_pollinator_sp)
  # 
  
  number_plant_individuals_per_plant_sp <- 
    generate_number_plant_individuals_per_plant_sp(number_plant_individuals)
  
  list_plant_individuals <- generate_list_plant_individuals(number_plant_individuals_per_plant_sp)
  
  number_intralinks_per_pollinator <- generate_number_intralinks_per_pollinator(number_pollinator_sp,
                                                                                number_intralinks)
  # Generate an empty interaction matrix
  
  interaction_matrix <- generate_empty_interaction_matrix(number_plant_individuals,
                                                          number_pollinator_sp,
                                                          number_plant_individuals_per_plant_sp)
  
  return(list(interaction_matrix, number_plant_individuals_per_plant_sp,
              list_plant_individuals, number_intralinks_per_pollinator,
              plant_pollinator_interlink_matrix))
  
}


generate_interaction_matrix <- function(number_pollinator_sp,
                                        number_plant_individuals,
                                        number_intralinks,
                                        number_interlinks){
  
  
  data_interaction_matrix <- init_interaction_matrix(number_pollinator_sp,
                                                     number_plant_individuals,
                                                     number_intralinks,
                                                     number_interlinks)
  
  interaction_matrix <- data_interaction_matrix[[1]]
  number_plant_individuals_per_plant_sp <- data_interaction_matrix[[2]]
  list_plant_individuals <- data_interaction_matrix[[3]]
  number_intralinks_per_pollinator <- data_interaction_matrix[[4]]
  plant_pollinator_interlink_matrix <- data_interaction_matrix[[5]]
  
  cont <- 1
  
  while(all(rowSums(interaction_matrix)>0)==FALSE){
    
    # cat(cont)
    #############
    #We try to generate a proper interaction matrix 1000 times with the info in data_interaction_matrix
    # If we fail we run again data_interaction_matrix (new data from init_interaction_matrix)
    
    if(cont<1000){ 
      
      cont=cont+1
      interaction_matrix <- 0*interaction_matrix
      
    }else{
      
      cont = 1
      data_interaction_matrix <- init_interaction_matrix(number_pollinator_sp,
                                                         number_plant_individuals,
                                                         number_intralinks,
                                                         number_interlinks)
      
      interaction_matrix <- data_interaction_matrix[[1]]
      number_plant_individuals_per_plant_sp <- data_interaction_matrix[[2]]
      list_plant_individuals <- data_interaction_matrix[[3]]
      number_intralinks_per_pollinator <- data_interaction_matrix[[4]]
      plant_pollinator_interlink_matrix <- data_interaction_matrix[[5]]
      
    }
    
    # Adjust number of intralinks
    
    for(pol_i in 1:number_pollinator_sp){
      
      data_for_intralinks_for_each_plant_sp <- 
        generate_intralinks_for_each_plant_sp(plant_pollinator_interlink_matrix,
                                              number_intralinks_per_pollinator, pol_i)
      
      plant_sps <- data_for_intralinks_for_each_plant_sp[[1]]
      intralinks_for_each_plant_sp <- data_for_intralinks_for_each_plant_sp[[2]]
      
      # number of individuals should be larger than the number of intralinks
      
      cont2 <- 1
      
      while((cont2<100) & (test_intralinks_for_each_plant_sp(plant_sps,
                                                             list_plant_individuals,
                                                             intralinks_for_each_plant_sp) == FALSE)){
        
        data_for_intralinks_for_each_plant_sp <- 
          generate_intralinks_for_each_plant_sp(plant_pollinator_interlink_matrix,
                                                number_intralinks_per_pollinator, pol_i)
        
        plant_sps <- data_for_intralinks_for_each_plant_sp[[1]]
        intralinks_for_each_plant_sp <- data_for_intralinks_for_each_plant_sp[[2]]
        
        cont2 <- cont2 + 1
        # cat("cont2: ",cont2,", ")
      }
      
      if((cont2 < 100) & (test_intralinks_for_each_plant_sp(plant_sps,
                                                          list_plant_individuals,
                                                          intralinks_for_each_plant_sp))){
       
        for (plant_sp_i in 1:length(plant_sps)) {
          
          plant_sp_selected <- plant_sps[plant_sp_i]
          individuals_selected <- sample(list_plant_individuals[[plant_sp_selected]],
                                         intralinks_for_each_plant_sp[plant_sp_i])
          
          interaction_matrix[individuals_selected,pol_i] <- 1
          
        }
         
      }
      
    }
    
    ##############
    
  }

  return(list(interaction_matrix,number_plant_individuals_per_plant_sp))
}

test_intralinks_for_each_plant_sp <- function(plant_sps,list_plant_individuals,
                                              intralinks_for_each_plant_sp){
  
  if(all(intralinks_for_each_plant_sp>0)==FALSE){
    
    return(FALSE)
    
  }else{
    
    test_aux <- NULL
    
    for(plant_i in 1:length(plant_sps)){
      
      plant_sp_i <- plant_sps[plant_i]
      
      test_aux <- 
        c(test_aux,
          length(list_plant_individuals[[plant_sp_i]]) >=  intralinks_for_each_plant_sp[plant_i])
      
    }
    
    return(all(test_aux))
    
  }
  
}


generate_intralinks_for_each_plant_sp <- function(plant_pollinator_interlink_matrix,
                                                  number_intralinks_per_pollinator,
                                                  pol_i){
  
  if(sum(plant_pollinator_interlink_matrix[,pol_i])>0){
    plant_sps <- which(plant_pollinator_interlink_matrix[,pol_i]>0) %>% as.numeric()
  }else{
    plant_sps <- sample(1:5,size = 1)
  }
  
  
  intralinks_for_each_plant_sp <- rand_vect(length(plant_sps),
                                            number_intralinks_per_pollinator[pol_i])
  
  return(list(plant_sps,intralinks_for_each_plant_sp))
  
}


weighted_interaction_matrix_generator <- function(number_plant_individuals,
                                                  number_pollinator_sp,
                                                  number_intralinks,
                                                  average_intralink_strenght,
                                                  number_interlinks){
  
  data_interaction_matrix <- generate_interaction_matrix(number_pollinator_sp,
                                                         number_plant_individuals,
                                                         number_intralinks,
                                                         number_interlinks)
  
  interaction_matrix_links <- data_interaction_matrix[[1]]
  number_plant_individuals_per_plant_sp <- data_interaction_matrix[[2]]
  
  sanity <- all(rowSums(interaction_matrix_links)>0) &
    (sum(interaction_matrix_links>0) == number_intralinks)
  
  while(!sanity){
    
    data_interaction_matrix <- generate_interaction_matrix(number_pollinator_sp,
                                                           number_plant_individuals,
                                                           number_intralinks,
                                                           number_interlinks)
    
    interaction_matrix_links <- data_interaction_matrix[[1]]
    number_plant_individuals_per_plant_sp <- data_interaction_matrix[[2]]
    
    sanity <- all(rowSums(interaction_matrix_links)>0) &
      (sum(interaction_matrix_links>0) == number_intralinks)
  }
  
  
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
  
  rowSums(weighted_interaction_matrix)==instrength_plant_individuals
  
  #Sanity_checks
  
  # weighted_interlinks <- interlink_tester(weighted_interaction_matrix,
  #                                     number_plant_individuals_per_plant_sp)
  # 
  # cat(paste0("the number of interlinks in weighted interaction matrix is correct: ",
  #            weighted_interlinks == number_interlinks,"\n"))
  # 
  # cat(paste0("the number of intralinks in weighted interaction matrix is correct: ",
  #            sum(weighted_interaction_matrix>0) == number_intralinks,"\n"))
  # 
  # cat(paste0("the instrength of intralinks in weighted interaction matrix is correct: ",
  #            all(sum(rowSums(weighted_interaction_matrix))==
  #                  number_intralinks*average_intralink_strenght),"\n"))
  # 
  # 
  
  return(list(weighted_interaction_matrix,number_plant_individuals_per_plant_sp))
  
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
  
  poll_plant_sp_1 <- poll_names[colSums(as.matrix(poll_plant_sp_1_aux))>0]
  poll_plant_sp_2 <- poll_names[colSums(as.matrix(poll_plant_sp_2_aux))>0]
  poll_plant_sp_3 <- poll_names[colSums(as.matrix(poll_plant_sp_3_aux))>0]
  poll_plant_sp_4 <- poll_names[colSums(as.matrix(poll_plant_sp_4_aux))>0]
  poll_plant_sp_5 <- poll_names[colSums(as.matrix(poll_plant_sp_5_aux))>0]
  
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
                                                   number_pollinator_sp,
                                                   number_intralinks,
                                                   average_intralink_strenght,
                                                   number_interlinks){
  
  #number_pollinator_sp <- number_intralinks - number_interlinks
  
  data_weighted_interaction_matrix <- weighted_interaction_matrix_generator(number_plant_individuals,
                                                                            number_pollinator_sp,
                                                                            number_intralinks,
                                                                            average_intralink_strenght,
                                                                            number_interlinks)
  
  weighted_interaction_matrix <- data_weighted_interaction_matrix[[1]]
  number_plant_individuals_per_plant_sp <- data_weighted_interaction_matrix[[2]]
  
  
  real_interlinks <- interlink_tester(weighted_interaction_matrix,
                                      number_plant_individuals_per_plant_sp)
  # Sanity check
  test1 <- all(sum(rowSums(weighted_interaction_matrix))==number_intralinks*average_intralink_strenght)
  test2 <- sum(weighted_interaction_matrix > 0) == number_intralinks
  test3 <- number_interlinks == real_interlinks
  
  while(all(test1,test2,test3)==FALSE){
    
    data_weighted_interaction_matrix <- weighted_interaction_matrix_generator(number_plant_individuals,
                                                                              number_pollinator_sp,
                                                                              number_intralinks,
                                                                              average_intralink_strenght,
                                                                              number_interlinks)
    
    weighted_interaction_matrix <- data_weighted_interaction_matrix[[1]]
    number_plant_individuals_per_plant_sp <- data_weighted_interaction_matrix[[2]]
    
    
    real_interlinks <- interlink_tester(weighted_interaction_matrix,
                                        number_plant_individuals_per_plant_sp)
    
    # Sanity check
    test1 <- all(sum(rowSums(weighted_interaction_matrix))==number_intralinks*average_intralink_strenght)
    test2 <- sum(weighted_interaction_matrix > 0) == number_intralinks
    test3 <- number_interlinks == real_interlinks
    
  }
  
  
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
  
  return(list(list_interactions,number_plant_individuals_per_plant_sp))
  
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
                                                number_pollinator_sp, number_intralinks,
                                                average_intralink_strenght,
                                                number_interlinks,
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
  

  interaction_list_without_phen_data <- generate_interaction_list_without_phen(number_plant_individuals,
                                                                               number_pollinator_sp,
                                                                               number_intralinks,
                                                                               average_intralink_strenght,
                                                                               number_interlinks)
  
  
  interaction_list_without_phen <- interaction_list_without_phen_data[[1]]
  number_plant_individuals_per_plant_sp <- interaction_list_without_phen_data[[1]]
  
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


