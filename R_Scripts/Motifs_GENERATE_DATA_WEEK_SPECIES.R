# load libraries
library(tidyverse)
library(bipartite)
library(matlib)
library(igraph)

####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################


fitness_data2 <- read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID_RAPE.csv")

# Filtering & relabeling
fitness2 <- fitness_data2 %>% filter(Year==2019)
  
# Calculating week number
fitness2 <- fitness2 %>% select(Day,Month,Year,Plot,Subplot,Plant_Simple,ID_Simple,Visits) %>%
  mutate(date_raw=as.Date(paste(Day,Month,Year,sep="/"), "%d/%m/%Y"),
         Week=as.numeric(format(date_raw, "%V")))

fitness2 <- rename(fitness2,ID=ID_Simple)

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

  #print(names_plot)
  #print(length(list_incid_matrix_19))
  
  first_motif_encuounter <- F
  
  for (plot_i in 1:length(list_incid_matrix_19)){
    
    
    #print(names_plot[plot_i])
    
    # Incidence matrix for each plot network

    incid_matrix_i <- list_incid_matrix_19[[plot_i ]] 

    graph_i <- graph_from_incidence_matrix(incid_matrix_i, weighted = T, directed = F)
    
    #plot(graph_i)
    #E(graph_i)
    

    # motiv in igraph (triplet is equal to a two-path graph, i.e., a star graph with 2 nodes)
    
    pattern <- make_star(3, mode = "undirected")
    #E(pattern)
    #plot(pattern)
    
    iso <- subgraph_isomorphisms(pattern, graph_i)      # takes a while
    
    #WARNING subgraph_isomorphisms works with directed graphs and
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
      
      motif_3$plot_id <- as.numeric(names_plot[plot_i])
      
      #for (i in 1:length(motifs)){print(V(motifs[[i]])$name)}
      
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
      #write_csv(motif_3, paste("Triplets 3_plot",plot_i,".csv",sep=""))
      
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

  # Homo_Motifs: misma especie otros subplots
  
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

########################################################################
# ESTIMATING CARACOLES'S homo + hetero motifs
########################################################################

for (week_i in unique(fitness2$Week)){

  print("WEEK")
  print(week_i)
  
  
  fitness_week_i <- fitness2 %>% filter(Week==week_i)
  
  #fitness_week_i <- fitness_week_i %>% group_by(Plot,Subplot,Plant_Simple,ID,num.plantas,Fruit,Seed) %>%
  #  count(wt=Visits) %>% rename(Visits_tot = n)
  
  #fitness_week_i$Subplot_Plant_Label <- paste(fitness_week_i$Subplot,fitness_week_i$Plant_Simple,sep = " ")
  #fitness_week_i <- fitness_week_i %>% mutate(Seeds_tot = num.plantas*Seed)
  
  #Nuevos fitness maria
  fitness_week_i <- fitness_week_i %>% group_by(Plot,Subplot,Plant_Simple,ID) %>%
    count(wt=Visits) %>% rename(Visits_tot = n)
  
  fitness_week_i$Subplot_Plant_Label <- paste(fitness_week_i$Subplot,fitness_week_i$Plant_Simple,sep = " ")
  fitness_week_i <- fitness_week_i
  
  visit_list_week <- fitness_week_i %>% ungroup() %>% select(Plot,ID,Subplot_Plant_Label,Visits_tot)
  
  #x <- motifs_extraction(visit_list_week)
  
  visit_list_week <- homo_hete_motifs(visit_list_week)
  visit_list_week <- visit_list_week %>% mutate(Week=week_i)
  
  if (week_i==min(fitness2$Week)){
    visit_list <- visit_list_week
  }else{
    visit_list <- visit_list %>% bind_rows(visit_list_week) 
  }
}

write_csv(visit_list, "Processed_data/Motifs_WEEK/Caracoles_WEEK_SPECIES.csv")

