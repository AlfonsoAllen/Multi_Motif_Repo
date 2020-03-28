# load libraries
library(tidyverse)
library(bipartite)
library(matlib)
library(igraph)

####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################

fitness_data <- read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019.csv")

fitness_data %>% group_by(G_F) %>% count()


fitness <- fitness_data %>% group_by(Plot,Subplot,Plant_Simple,G_F,num.plantas,Fruit,Seed) %>%
  count(wt=Visits) %>% rename(Visits_tot = n)

#####################################
# Filtering & relabeling
#####################################

fitness <- fitness %>% filter(!Subplot == "OUT" & !is.na(G_F))

fitness$Subplot_Plant_Label <- paste(fitness$Subplot,fitness$Plant_Simple,sep = " ")
fitness <- fitness %>% mutate(Seeds_tot = num.plantas*Seed)


##################################################################
# FUNCTION (NAME): motifs_extraction
# INPUT (1) -> visit_list: it contains (Plot,G_F,Subplot_Plant_Label(e.g. A1 LEMA),Visits_tot)
# OUTPUT (1) -> list with motifs: it contains
# (plot_id,number_nodes,number_plants,plant_1,subplot_plant_1,
# plant_2,subplot_plant_2,poll_1,poll_2
##################################################################

motifs_extraction <- function(visit_list) {
  
  ##############################################################
  # GENERATE A LIST OF BIPARTITE NETWORKS (1 NETWORK PER PLOT)
  ##############################################################

  testdata_19 <-   data.frame(higher = visit_list$G_F,
                              lower = visit_list$Subplot_Plant_Label,
                              webID = visit_list$Plot,
                              freq = visit_list$Visits_tot)
  
  list_incid_matrix_19 <- frame2webs(testdata_19,type.out="list")

  for (plot_i in 1:9){
    
    # Incidence matrix for each plot network

    incid_matrix_i <- list_incid_matrix_19[[plot_i ]] 

    graph_i <- graph_from_incidence_matrix(incid_matrix_i, weighted = T)

    # motiv in igraph (triplet is equal to a two-path graph, i.e., a star graph with 2 nodes)
    
    pattern <- make_star(3, mode = "undirected")
    iso <- subgraph_isomorphisms(pattern, graph_i)      # takes a while
    motifs <- lapply(iso, function (x) { induced_subgraph(graph_i, x) })

    ##############################################################
    # TRIPLETS ANALYSIS
    ##############################################################
    
    # We collect triplets information in "motif_3" tibble
    
    tbl_colnames <- c("plot_id","number_nodes","number_plants","plant_1","subplot_plant_1",
                      "plant_2","subplot_plant_2","poll_1","poll_2")

    motif_3 <- as_tibble(data.frame(matrix(nrow=length(motifs),ncol=length(tbl_colnames))))
    colnames(motif_3) <- tbl_colnames
    
    motif_3$plot_id <- plot_i
    
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
    print(paste("Triplets_plot",plot_i,sep=""))
    #write_csv(motif_3, paste("Triplets 3_plot",plot_i,".csv",sep=""))
    
    if (plot_i==1){
      motif_3 <- unique(motif_3) #to avoid duplicities due to subgraph isomorphisms
      motif_3_list <- motif_3
    }else{
      motif_3 <- unique(motif_3) #to avoid duplicities due to subgraph isomorphisms
      motif_3_list <- motif_3_list %>% bind_rows(motif_3) 
    }
  }
  return(motif_3_list)
}

##################################################################
# FUNCTION (NAME): homo_hete_motifs
# INPUT (1) -> visit_list: it contains (Plot,G_F,Subplot_Plant_Label(e.g. A1 LEMA),Visits_tot)
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
    
    homo_motif <- motif_3_Carac %>% filter(Same_plant & poll_1 == output_funct$G_F[i] &
                                             (motif_3_Carac$descript_plant_1 == descrip|motif_3_Carac$descript_plant_2==descrip))
    num_homo_motif <- sum(homo_motif$Same_plant,na.rm = TRUE)
    
    hete_motif <- motif_3_Carac %>% filter(!Same_plant & number_plants==2 & poll_1 == output_funct$G_F[i] &
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

visit_list <- fitness %>% ungroup() %>% select(Plot,G_F,Subplot_Plant_Label,Visits_tot)
visit_list <- homo_hete_motifs(visit_list)

write_csv(visit_list, "Processed_data/Motifs_YEAR/Caracoles_YEAR.csv")

########################################################################
# NULL MODEL
########################################################################

abundances <- read_csv2("Raw_Data/abundances_2019.csv")
abundances_19 <- abundances %>% filter(year==2019) %>% 
  select(plot,subplot,species,individuals) %>% arrange(plot,subplot)

# Since there is no ME in abundances and pollination data shows that plant, we set 
# MEEL, and MESU to ME

abundances_19$species[abundances_19$species=="MEEL"] <- "ME"
abundances_19$species[abundances_19$species=="MESU"] <- "ME"

abundances_19 <- abundances_19 %>% group_by(plot, subplot,species) %>% count(wt=individuals)%>%
  rename(individuals=n)

abundances_19_flowers <- abundances_19 %>%
  filter(species %in% c("CHFU","CHMI","LEMA","ME","PUPA","RAPE")) 

total_indv <-  abundances_19_flowers %>% group_by(plot,subplot) %>% count(wt=individuals) %>%
  rename(total_indviduals=n)
abundances_19_flowers <- abundances_19_flowers %>% left_join(total_indv,by=c("plot","subplot"))
abundances_19_flowers  <- mutate(abundances_19_flowers, proportion=individuals/total_indviduals)

visits_GF <- fitness %>% ungroup() %>% group_by(Plot,Subplot,G_F)%>% count(wt=Visits_tot)%>%
  rename(Visits_tot=n)%>% ungroup()


##################################################################
##################################################################
# FUNCTION (NAME): random_visits_funct
# INPUT (1) -> visits_GF: it contains (Plot,Subplot,G_F,Visits_tot)
# OUTPUT (1) -> random_visits: (Plot,G_F,Subplot_Plant_Label(e.g. A1 LEMA),Visits_tot)
##################################################################

random_visits_funct <- function(visits_GF){

  for (i in 1:nrow(visits_GF)){
    
    ab_plot_subplot <- abundances_19_flowers %>% filter(plot==visits_GF$Plot[i],
                                                        subplot==visits_GF$Subplot[i])
    
    sample_i <- sample(ab_plot_subplot$species,
                       size = visits_GF$Visits_tot[i],
                       prob=ab_plot_subplot$proportion,replace = TRUE)
    
    sample_i <- as.data.frame(table(sample_i))
    
    #We create a tibble to load the random visitations
    
    tbl_colnames <- c("Plot","Subplot","Plant_Simple","G_F",
                      "Visits_tot","Subplot_Plant_Label")
    
    random_visits_i <- as_tibble(data.frame(matrix(nrow=nrow(sample_i),ncol=length(tbl_colnames))))
    colnames(random_visits_i) <- tbl_colnames
    
    for (j in 1: nrow(random_visits_i)){
  
      random_visits_i$Plot[j] <- visits_GF$Plot[i]
      random_visits_i$Subplot[j] <- visits_GF$Subplot[i]
      random_visits_i$Plant_Simple[j] <- levels(sample_i[j,1])
      random_visits_i$G_F[j] <- visits_GF$G_F[i]
      random_visits_i$Visits_tot[j] <- sample_i[j,2]
      random_visits_i$Subplot_Plant_Label[j] <- paste(visits_GF$Subplot[i],sample_i[j,1],sep = " ")
    }
    
    if(i==1){
      
      random_visits <- random_visits_i
        
    }else{
      
      random_visits <- bind_rows(random_visits,random_visits_i)
      }
    
  }
  
  return(random_visits)
  
}

set.seed(123)

for (count_i in 1:100){

  rd_visits <- random_visits_funct(visits_GF)
  
  rd_visits_adapt <- rd_visits %>% select(Plot,G_F,Subplot_Plant_Label,Visits_tot)
  
  rd_visit_list <- homo_hete_motifs(as_tibble(rd_visits_adapt))
  
  print(paste("random_",count_i,".csv",sep = ""))
  
  write_csv(rd_visit_list, paste("Processed_data/Motifs_YEAR/random_",count_i,"_YEAR.csv",sep = ""))

}



