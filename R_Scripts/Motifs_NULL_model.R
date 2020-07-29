# load libraries
library(tidyverse)
library(bipartite)
library(matlib)
library(igraph)
source("R_Scripts/random_visits_WEEKS_motifs.R")
####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################


fitness_data2 <- read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")

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

set.seed(123) # for reproducible results

number_repetitions = 500

visit_list_REAL <- read_csv("Processed_data/Motifs_WEEK/Caracoles_WEEK_SPECIES.csv")

pooled_data_REAL <- visit_list_REAL %>% group_by(Plot,Subplot_Plant_Label) %>% 
  summarise(homo_motif_real=sum(homo_motif),hete_motif_real=sum(hete_motif))

homo_table <- pooled_data_REAL %>% select(-hete_motif_real)
hete_table <- pooled_data_REAL %>% select(-homo_motif_real)


for (rep in 1:number_repetitions){

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
    
    visit_list_week_REAL <- fitness_week_i %>% ungroup() %>% select(Plot,ID,Subplot_Plant_Label,Visits_tot)
    
    #random visits_week
    
    visit_list_week <- random_visits_WEEKS(visit_list_week_REAL)
    
    # Motif estimation
    
    visit_list_week <- homo_hete_motifs(visit_list_week)
    visit_list_week <- visit_list_week %>% mutate(Week=week_i)
    
    if (week_i==min(fitness2$Week)){
      visit_list <- visit_list_week
    }else{
      visit_list <- visit_list %>% bind_rows(visit_list_week) 
    }
  }

  
  pooled_data_random <- visit_list %>% group_by(Plot,Subplot_Plant_Label) %>% 
    summarise(homo_motif=sum(homo_motif), hete_motif=sum(hete_motif))
  
  homo_random <- pooled_data_random %>% select(-hete_motif)
  hete_random <- pooled_data_random %>% select(-homo_motif)
  
  homo_table <- homo_table %>% left_join(homo_random,by=c("Plot","Subplot_Plant_Label"))
  
  colnames(homo_table)[4:length(colnames(homo_table))] <- as.character(1:(length(colnames(homo_table))-3))
  
  hete_table <- hete_table %>% left_join(hete_random,by=c("Plot","Subplot_Plant_Label"))
  colnames(hete_table)[4:length(colnames(hete_table))] <- as.character(1:(length(colnames(hete_table))-3))
}
 
# Commented for security reasons

# write_csv(homo_table, "Processed_data/Motifs_WEEK/HOMO_MOTIFS_WEEK_SPECIES.csv")
# write_csv(hete_table, "Processed_data/Motifs_WEEK/HETE_MOTIFS_WEEK_SPECIES.csv")

homo_table <- read_csv("Processed_data/Motifs_WEEK/HOMO_MOTIFS_WEEK_SPECIES.csv")
hete_table <- read_csv("Processed_data/Motifs_WEEK/HETE_MOTIFS_WEEK_SPECIES.csv")




homo_table[is.na(homo_table)] <- 0
hete_table[is.na(hete_table)] <- 0

# Same order for plants and plots: Sanity check
sum(homo_table[,c(1,2)] != hete_table[,c(1,2)])

sum_table <- homo_table

for(i in 1:nrow(sum_table)){
  for(j in 3:ncol(sum_table)){
    
    sum_table[i,j] <- homo_table[i,j] + hete_table[i,j]
    
  }
}


homo_table$significance <- NA
hete_table$significance <- NA
sum_table$significance <- NA

for (i in 1:nrow(homo_table)){
  
  # Homo motuf significance
  
  CI <- quantile(as.numeric(homo_table[i,4:(ncol(homo_table)-1)]), probs = c(0.025, 0.975))
  
  if (as.numeric(homo_table[i,1]) < CI[[1]] | as.numeric(homo_table[i,1]) > CI[[2]]){
    homo_table$significance[i] <- TRUE
  }else{
    homo_table$significance[i] <- FALSE
  }
  
  # Hetemotif significance
  
  CI <- quantile(as.numeric(hete_table[i,4:(ncol(hete_table)-1)]), probs = c(0.025, 0.975))
  
  if (as.numeric(hete_table[i,1]) < CI[[1]] | as.numeric(hete_table[i,1]) > CI[[2]]){
    hete_table$significance[i] <- TRUE
  }else{
    hete_table$significance[i] <- FALSE
  }
  
  # Sum significance
  
  CI <- quantile(as.numeric(sum_table[i,4:(ncol(sum_table)-1)]), probs = c(0.025, 0.975))
  
  if (as.numeric(sum_table[i,1]) < CI[[1]] | as.numeric(sum_table[i,1]) > CI[[2]]){
    sum_table$significance[i] <- TRUE
  }else{
    sum_table$significance[i] <- FALSE
  }
  
}

###########################
# PROCESSING SIGNIFICANCE
###########################

# HOMOMOTiFS-----

sum(homo_table$significance)/nrow(homo_table) #0.2249443

homo_table <- homo_table %>% mutate(aux_c=Subplot_Plant_Label) %>%
  separate(aux_c,c("Subplot","Plant_Label")," ") %>% select(-Subplot)

homo_table_Plants <- homo_table %>% select(Plot,Plant_Label) %>%  group_by(Plot,Plant_Label) %>% 
  count() %>% rename(total_focals=n)

homo_table_SIGNIF_Plants <- homo_table %>% filter(significance==TRUE) %>% 
  select(Plot,Plant_Label) %>%  group_by(Plot,Plant_Label) %>% count() %>%
  rename(significant_focals=n)

homo_significance_Plants <- homo_table_Plants %>% left_join(homo_table_SIGNIF_Plants,by=c("Plot","Plant_Label"))
homo_significance_Plants[is.na(homo_significance_Plants)] <- 0.0
homo_significance_Plants <- homo_significance_Plants %>% mutate(percent_signif=100*significant_focals/total_focals)

homo_significance_Figure <- homo_significance_Plants %>% 
  mutate(non_significant=total_focals-significant_focals) %>%
  gather(significance,amount_plants,c("significant_focals","non_significant")) %>%
  arrange(Plot,Plant_Label)

homo_significance_Figure$significance[homo_significance_Figure$significance=="non_significant"] <- 
  "Non significant"
homo_significance_Figure$significance[homo_significance_Figure$significance=="significant_focals"] <- 
  "Significant"


# HETEMOTiFS-----

sum(hete_table$significance)/nrow(hete_table) #0.8262806

hete_table <- hete_table %>% mutate(aux_c=Subplot_Plant_Label) %>%
  separate(aux_c,c("Subplot","Plant_Label")," ") %>% select(-Subplot)

hete_table_Plants <- hete_table %>% select(Plot,Plant_Label) %>%  group_by(Plot,Plant_Label) %>% 
  count() %>% rename(total_focals=n)

hete_table_SIGNIF_Plants <- hete_table %>% filter(significance==TRUE) %>% 
  select(Plot,Plant_Label) %>%  group_by(Plot,Plant_Label) %>% count() %>%
  rename(significant_focals=n)

hete_significance_Plants <- hete_table_Plants %>% left_join(hete_table_SIGNIF_Plants,by=c("Plot","Plant_Label"))
hete_significance_Plants[is.na(hete_significance_Plants)] <- 0.0
hete_significance_Plants <- hete_significance_Plants %>% mutate(percent_signif=100*significant_focals/total_focals)

hete_significance_Figure <- hete_significance_Plants %>% 
  mutate(non_significant=total_focals-significant_focals) %>%
  gather(significance,amount_plants,c("significant_focals","non_significant")) %>%
  arrange(Plot,Plant_Label)

hete_significance_Figure$significance[hete_significance_Figure$significance=="non_significant"] <- 
  "Non significant"
hete_significance_Figure$significance[hete_significance_Figure$significance=="significant_focals"] <- 
  "Significant"

# SUM MOTiFS-----

sum(sum_table$significance)/nrow(sum_table) #0.1804009

sum_table <- sum_table %>% mutate(aux_c=Subplot_Plant_Label) %>%
  separate(aux_c,c("Subplot","Plant_Label")," ") %>% select(-Subplot)

sum_table_Plants <- sum_table %>% select(Plot,Plant_Label) %>%  group_by(Plot,Plant_Label) %>% 
  count() %>% rename(total_focals=n)

sum_table_SIGNIF_Plants <- sum_table %>% filter(significance==TRUE) %>% 
  select(Plot,Plant_Label) %>%  group_by(Plot,Plant_Label) %>% count() %>%
  rename(significant_focals=n)

sum_significance_Plants <- sum_table_Plants %>% left_join(sum_table_SIGNIF_Plants,by=c("Plot","Plant_Label"))
sum_significance_Plants[is.na(sum_significance_Plants)] <- 0.0
sum_significance_Plants <- sum_significance_Plants %>% mutate(percent_signif=100*significant_focals/total_focals)

sum_significance_Figure <- sum_significance_Plants %>% 
  mutate(non_significant=total_focals-significant_focals) %>%
  gather(significance,amount_plants,c("significant_focals","non_significant")) %>%
  arrange(Plot,Plant_Label)

sum_significance_Figure$significance[sum_significance_Figure$significance=="non_significant"] <- 
  "Non significant"
sum_significance_Figure$significance[sum_significance_Figure$significance=="significant_focals"] <- 
  "Significant"



######################
# PLOTS

plot_labs <-c(
  "Plot 1",
  "Plot 2",
  "Plot 3",
  "Plot 4",
  "Plot 5",
  "Plot 6",
  "Plot 7",
  "Plot 8",
  "Plot 9"
)
names(plot_labs) <- c(
  '1',
  '2',
  '3',
  '4',
  '5',
  "6",
  '7',
  '8',
  "9"
)

# Homo----

ggplot(homo_significance_Figure, aes(fill=significance, y=amount_plants, x=Plant_Label)) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw()+labs(title="Homo-motifs",
                  x ="Plant Species", y = "Number of focal individuals",fill=NULL)
  

ggplot(homo_significance_Figure, aes(fill=significance, y=amount_plants, x=Plant_Label)) + 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  theme_bw()+labs(title="Homo-motifs",
                  x ="Plant Species", y = "Number of focal individuals",fill=NULL)


# Example of distribution

i=5

homo_example  <-  tibble(triplets=as.numeric(homo_table[i,4:(ncol(homo_table)-1)]))

ggplot(homo_example)+
  geom_histogram(aes(x=triplets),binwidth=1)+
  theme_bw()+
  geom_vline(aes(xintercept=as.numeric(homo_table[i,3])), colour="deepskyblue",linetype = "dashed",size=1)+
  labs(x="Total amount of homospecific networks", y = "Number of randomized networks",
       title = paste("Plot",homo_table$Plot[i], homo_table$Subplot_Plant_Label[i],sep=" "))



# Hete----

ggplot(hete_significance_Figure, aes(fill=significance, y=amount_plants, x=Plant_Label)) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw()+labs(title="Hete-motifs",
                  x ="Plant Species", y = "Number of focal individuals",fill=NULL)


ggplot(hete_significance_Figure, aes(fill=significance, y=amount_plants, x=Plant_Label)) + 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  theme_bw()+labs(title="Hete-motifs",
                  x ="Plant Species", y = "Number of focal individuals",fill=NULL)


# SUM----

ggplot(sum_significance_Figure, aes(fill=significance, y=amount_plants, x=Plant_Label)) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw()+labs(title="Homo-motifs + Hete-motifs",
                  x ="Plant Species", y = "Number of focal individuals",fill=NULL)


ggplot(sum_significance_Figure, aes(fill=significance, y=amount_plants, x=Plant_Label)) + 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  theme_bw()+labs(title="Homo-motifs + Hete-motifs",
                  x ="Plant Species", y = "Number of focal individuals",fill=NULL)


#################################
# SIGNIFICANT HETEMOTIFS: ARE LARGER OR SMALLER THAN THE OBSERVED VALUES?
#################################

hete_table$Comparison <- NA

for (i in 1:nrow(hete_table)){
  
  # confidence interval
  
  CI <- quantile(as.numeric(hete_table[i,4:(ncol(hete_table)-3)]), probs = c(0.025, 0.975))
  
  if(hete_table$significance[i] != TRUE){
    hete_table$Comparison[i] <- NA
  }else if(as.numeric(hete_table[i,1]) < CI[[1]]){
    hete_table$Comparison[i] <- "Smaller"
  }else if(as.numeric(hete_table[i,1]) > CI[[2]]){
    hete_table$Comparison[i] <- "Larger"
  }
 
  
}

hete_table %>% group_by(Comparison) %>% count()
