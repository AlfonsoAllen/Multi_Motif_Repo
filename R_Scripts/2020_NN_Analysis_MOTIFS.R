# load libraries
library(tidyverse)
library(bipartite)
library(matlib)
library(igraph)
source("R_Scripts/functions.R")
####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################

fitness_data2 <- read_csv2("Raw_Data/final_Pollinators_2020.csv")
# Remove points from ID names
fitness_data2$ID <- sub("\\.", "", fitness_data2$ID)
fitness_data2$ID_Simple <- sub("\\.", "", fitness_data2$ID_Simple)

# filter tabanidae
fitness_data2 <- fitness_data2 %>% filter(ID != "Tabanidae")

# Filtering & relabeling
fitness2 <- fitness_data2 %>% filter(!is.na(Plant),Plant!="0",Subplot!="OUT",Plant!="Ground")

# Calculating week number
fitness2 <- fitness2 %>% select(Day,Month,Year,Plot,Subplot,Plant,ID_Simple,Visits) %>%
  mutate(date_raw=as.Date(paste(Day,Month,Year,sep="/"), "%d/%m/%Y"),
         Week=as.numeric(format(date_raw, "%V")))

fitness2 <- rename(fitness2,ID=ID_Simple,Plant_Simple=Plant)


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

visit_list_REAL <- read_csv("Processed_data/Motifs_WEEK/2020_Caracoles_WEEK_SPECIES.csv")

pooled_data_REAL <- visit_list_REAL %>% group_by(Plot,Subplot_Plant_Label) %>% 
  summarise(homo_motif_real=sum(homo_motif),hete_motif_real=sum(hete_motif))

homo_table <- pooled_data_REAL %>% dplyr::select(-hete_motif_real)
hete_table <- pooled_data_REAL %>% dplyr::select(-homo_motif_real)


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
    
    visit_list_week_REAL <- fitness_week_i %>% ungroup() %>% dplyr::select(Plot,ID,Subplot_Plant_Label,Visits_tot)
    
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
    summarise(homo_motif=sum(homo_motif,na.rm = T), hete_motif=sum(hete_motif,na.rm = T))
  
  homo_random <- pooled_data_random %>% dplyr::select(-hete_motif)
  hete_random <- pooled_data_random %>% dplyr::select(-homo_motif)
  
  homo_table <- homo_table %>% left_join(homo_random,by=c("Plot","Subplot_Plant_Label"))
  
  colnames(homo_table)[4:length(colnames(homo_table))] <- as.character(1:(length(colnames(homo_table))-3))
  
  hete_table <- hete_table %>% left_join(hete_random,by=c("Plot","Subplot_Plant_Label"))
  colnames(hete_table)[4:length(colnames(hete_table))] <- as.character(1:(length(colnames(hete_table))-3))
}
 
# Commented for security reasons

# write_csv(homo_table, "Processed_data/Motifs_WEEK/2020_HOMO_MOTIFS_WEEK_SPECIES.csv")
# write_csv(hete_table, "Processed_data/Motifs_WEEK/2020_HETE_MOTIFS_WEEK_SPECIES.csv")

homo_table <- read_csv("Processed_data/Motifs_WEEK/2020_HOMO_MOTIFS_WEEK_SPECIES.csv")
hete_table <- read_csv("Processed_data/Motifs_WEEK/2020_HETE_MOTIFS_WEEK_SPECIES.csv")




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
# percentage de homo tripl. that are significant
100*sum(homo_table$significance)/nrow(homo_table) #0.2249443

homo_table <- homo_table %>% mutate(aux_c=Subplot_Plant_Label) %>%
  separate(aux_c,c("Subplot","Plant_Label")," ") %>% dplyr::select(-Subplot)

homo_table_Plants <- homo_table %>% dplyr::select(Plot,Plant_Label) %>%  group_by(Plot,Plant_Label) %>% 
  count() %>% rename(total_focals=n)

homo_table_SIGNIF_Plants <- homo_table %>% filter(significance==TRUE) %>% 
  dplyr::select(Plot,Plant_Label) %>%  group_by(Plot,Plant_Label) %>% count() %>%
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

100*sum(hete_table$significance)/nrow(hete_table) #0.8262806

hete_table <- hete_table %>% mutate(aux_c=Subplot_Plant_Label) %>%
  separate(aux_c,c("Subplot","Plant_Label")," ") %>% dplyr::select(-Subplot)

hete_table_Plants <- hete_table %>% dplyr::select(Plot,Plant_Label) %>%  group_by(Plot,Plant_Label) %>% 
  count() %>% rename(total_focals=n)

hete_table_SIGNIF_Plants <- hete_table %>% filter(significance==TRUE) %>% 
  dplyr::select(Plot,Plant_Label) %>%  group_by(Plot,Plant_Label) %>% count() %>%
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

sum(sum_table$significance)/nrow(sum_table) #0.44

sum_table <- sum_table %>% mutate(aux_c=Subplot_Plant_Label) %>%
  separate(aux_c,c("Subplot","Plant_Label")," ") %>% dplyr::select(-Subplot)

sum_table_Plants <- sum_table %>% dplyr::select(Plot,Plant_Label) %>%  group_by(Plot,Plant_Label) %>% 
  count() %>% rename(total_focals=n)

sum_table_SIGNIF_Plants <- sum_table %>% filter(significance==TRUE) %>% 
  dplyr::select(Plot,Plant_Label) %>%  group_by(Plot,Plant_Label) %>% count() %>%
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
                  x ="Plant Species", y = "Number of focal individuals",fill=NULL)+
  theme(axis.text.x = element_text(angle = 90))


# Example of distribution

i=5

homo_example  <-  tibble(triplets=as.numeric(homo_table[i,4:(ncol(homo_table)-1)]))

# Save 700 x 500
png("New_Figures/figA111.png", width=1961*2, height = 1961*2*362/517, res=300*2)
ggplot(homo_example)+
  geom_histogram(aes(x=triplets),binwidth=1)+
  theme_bw()+
  geom_vline(aes(xintercept=as.numeric(homo_table[i,3])), colour="deepskyblue",linetype = "dashed",size=1)+
  labs(x="Total amount of homospecific subgraphs", y = "Number of randomized networks",
       title = expression(paste("Plot ","1 ",italic("L. maroccanus"), sep=" ")))
dev.off()


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

homo_table$Comparison <- NA

which(homo_table$Plot==2 & homo_table$Subplot_Plant_Label=="F3 LEMA")

for (i in 1:nrow(homo_table)){
  
  # confidence interval
  
  CI <- quantile(as.numeric(homo_table[i,4:(ncol(homo_table)-3)]), probs = c(0.025, 0.975))
  
  if(as.numeric(homo_table[i,1]) < CI[[1]]){
    homo_table$Comparison[i] <- "Significantly smaller"
  }else if(as.numeric(homo_table[i,1]) > CI[[2]]){
    homo_table$Comparison[i] <- "Significantly larger"
  }else{
    homo_table$Comparison[i] <- "Non-significant"
  }
}

homo_table %>% group_by(Comparison,Plant_Label) %>% count()

homo_table$cant_focal <- 1

homo_table_exp <- homo_table
homo_table_exp$Plant_Label[homo_table_exp$Plant_Label == "BEMA"] <- "B. macrocarpa"
homo_table_exp$Plant_Label[homo_table_exp$Plant_Label == "CETE"] <- "C. tenuiflorum"
homo_table_exp$Plant_Label[homo_table_exp$Plant_Label == "CHFU"] <- "C. fuscatum"
homo_table_exp$Plant_Label[homo_table_exp$Plant_Label == "CHMI"] <- "C. mixtum"
homo_table_exp$Plant_Label[homo_table_exp$Plant_Label == "LEMA"] <- "L. maroccanus"
homo_table_exp$Plant_Label[homo_table_exp$Plant_Label == "MESU"] <- "M. sulcatus"
homo_table_exp$Plant_Label[homo_table_exp$Plant_Label == "PUPA"] <- "P. paludosa"
homo_table_exp$Plant_Label[homo_table_exp$Plant_Label == "SCLA"] <- "S. laciniata"
homo_table_exp$Plant_Label[homo_table_exp$Plant_Label == "SOAS"] <- "S. asper"
homo_table_exp$Plant_Label[homo_table_exp$Plant_Label == "SPRU"] <- "S. rubra"

png("New_Figures/figA112.png", width=1961*2, height = 1961*2*500/600, res=300*2)
ggplot(homo_table_exp, aes(fill=Comparison, y=cant_focal, x=Plant_Label)) + 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  theme_bw()+labs(title="Homospecific subgraphs",
                  x ="Plant Species", y = "Number of focal individuals",fill=NULL)+
  theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(angle = 90,vjust=0.5, hjust=1))+ 
  theme(axis.text.x = element_text(face = "italic"))
dev.off()
#save 600 x 500

ggplot(homo_table, aes(fill=Comparison, y=cant_focal, x=Plant_Label)) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw()+labs(title="Homo-motifs",
                  x ="Plant Species", y = "Number of focal individuals",fill=NULL)+
  theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(angle = 90))

#Hete----
hete_table$Comparison <- NA


for (i in 1:nrow(hete_table)){
  
  # confidence interval
  
  CI <- quantile(as.numeric(hete_table[i,4:(ncol(hete_table)-3)]), probs = c(0.025, 0.975))
  
  if(as.numeric(hete_table[i,1]) < CI[[1]]){
    hete_table$Comparison[i] <- "Significantly smaller"
  }else if(as.numeric(hete_table[i,1]) > CI[[2]]){
    hete_table$Comparison[i] <- "Significantly larger"
  }else{
    hete_table$Comparison[i] <- "Non-significant"
  }
}

hete_table %>% group_by(Comparison,Plant_Label) %>% count()

hete_table$cant_focal <- 1

hete_table_exp <- hete_table
hete_table_exp$Plant_Label[hete_table_exp$Plant_Label == "BEMA"] <- "B. macrocarpa"
hete_table_exp$Plant_Label[hete_table_exp$Plant_Label == "CETE"] <- "C. tenuiflorum"
hete_table_exp$Plant_Label[hete_table_exp$Plant_Label == "CHFU"] <- "C. fuscatum"
hete_table_exp$Plant_Label[hete_table_exp$Plant_Label == "CHMI"] <- "C. mixtum"
hete_table_exp$Plant_Label[hete_table_exp$Plant_Label == "LEMA"] <- "L. maroccanus"
hete_table_exp$Plant_Label[hete_table_exp$Plant_Label == "MESU"] <- "M. sulcatus"
hete_table_exp$Plant_Label[hete_table_exp$Plant_Label == "PUPA"] <- "P. paludosa"
hete_table_exp$Plant_Label[hete_table_exp$Plant_Label == "SCLA"] <- "S. laciniata"
hete_table_exp$Plant_Label[hete_table_exp$Plant_Label == "SOAS"] <- "S. asper"
hete_table_exp$Plant_Label[hete_table_exp$Plant_Label == "SPRU"] <- "S. rubra"

png("New_Figures/figA113.png", width=1961*2, height = 1961*2*500/600, res=300*2)
ggplot(hete_table_exp, aes(fill=Comparison, y=cant_focal, x=Plant_Label)) + 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  theme_bw()+labs(title="Heterospecific subgraphs",
                  x ="Plant Species", y = "Number of focal individuals",fill=NULL)+
  theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(angle = 90,vjust=0.5, hjust=1))+
  theme(axis.text.x = element_text(face = "italic"))
dev.off()

homo_table %>% separate(Subplot_Plant_Label,c("Subplot","Plant_Simple")," ") %>%
  group_by(Comparison,Plant_Simple) %>% count() %>% mutate(percen=100*n/nrow(hete_table))


homo_table %>% separate(Subplot_Plant_Label,c("Subplot","Plant_Simple")," ") %>%
  group_by(Plant_Simple) %>% count() %>% mutate(percen=100*n/nrow(hete_table))

homo_table %>% separate(Subplot_Plant_Label,c("Subplot","Plant_Simple")," ") %>%
  group_by(Plot) %>% count() %>% mutate(percen=100*n/nrow(hete_table))%>% arrange(desc(n))

homo_table %>% separate(Subplot_Plant_Label,c("Subplot","Plant_Simple")," ") %>%
  group_by(Plant_Simple,Plot) %>% count() %>% mutate(percen=100*n/nrow(hete_table)) %>% 
  filter(percen>1) %>% arrange(desc(n))

homo_table %>% group_by(Comparison) %>% count() %>% mutate(percen=100*n/nrow(hete_table))
hete_table %>% group_by(Comparison) %>% count() %>% mutate(percen=100*n/nrow(hete_table))


############
homo_data <- homo_table %>% dplyr::select(Plot,Subplot_Plant_Label,homo_motif_real) %>% 
  separate(Subplot_Plant_Label,c("Sub","Plant_Simple")," ")%>% rename(motif=homo_motif_real) %>%
  mutate(type="Homosp. triplet")

hete_data <- hete_table %>% dplyr::select(Plot,Subplot_Plant_Label,hete_motif_real) %>% 
  separate(Subplot_Plant_Label,c("Sub","Plant_Simple")," ") %>% rename(motif=hete_motif_real) %>%
  mutate(type="Heterosp. triplet")

data_total <- bind_rows(homo_data,hete_data)

means <- aggregate((motif) ~  Plant_Simple + Plot,
                   homo_data, mean)

ggplot(homo_data,aes(x=Plant_Simple,y=motif))+
  geom_boxplot()+
  geom_point(aes(color=Plant_Simple), position = "jitter", alpha=0.2)+
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = means, aes(label = round(`(motif)`,2), y = 38))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("# Homo-motifs")+
  theme_bw()+ theme(legend.position = "none")#+stat_compare_means()
#labs(Color

##########
# Test significance

# Test significance

# Plot1
i=1
homo_data_i <- homo_data %>% filter(Plot==i)#Significant
kruskal.test(motif ~ Plant_Simple, data = homo_data_i)
pairwise.wilcox.test(homo_data_i$motif, homo_data_i$Plant_Simple,
                     p.adjust.method = "BH")#No sig: ME-CHFU

# Plot2
i=2
homo_data_i <- homo_data %>% filter(Plot==i) #Significant
kruskal.test(motif ~ Plant_Simple, data = homo_data_i)
pairwise.wilcox.test(homo_data_i$motif, homo_data_i$Plant_Simple,
                     p.adjust.method = "BH") 
# Plot3
i=3
homo_data_i <- homo_data %>% filter(Plot==i)
kruskal.test(motif ~ Plant_Simple, data = homo_data_i)
pairwise.wilcox.test(homo_data_i$motif, homo_data_i$Plant_Simple,
                     p.adjust.method = "BH")
# Plot4
i=4
homo_data_i <- homo_data %>% filter(Plot==i) #No significant
kruskal.test(motif ~ Plant_Simple, data = homo_data_i)
pairwise.wilcox.test(homo_data_i$motif, homo_data_i$Plant_Simple,
                     p.adjust.method = "BH")
homo_data_i %>% group_by(Plant_Simple) %>% count()
# Plot5
i=5
homo_data_i <- homo_data %>% filter(Plot==i) #Significant
kruskal.test(motif ~ Plant_Simple, data = homo_data_i)
pairwise.wilcox.test(homo_data_i$motif, homo_data_i$Plant_Simple,
                     p.adjust.method = "BH")

homo_data_i %>% group_by(Plant_Simple) %>% count()
# Plot6
i=6
homo_data_i <- homo_data %>% filter(Plot==i) 
kruskal.test(motif ~ Plant_Simple, data = homo_data_i) #significant
pairwise.wilcox.test(homo_data_i$motif, homo_data_i$Plant_Simple,
                     p.adjust.method = "BH") #N0: ME-CHFU,ME-LEMA
# Plot7
i=7
homo_data_i <- homo_data %>% filter(Plot==i)
kruskal.test(motif ~ Plant_Simple, data = homo_data_i) #Significant
pairwise.wilcox.test(homo_data_i$motif, homo_data_i$Plant_Simple,
                     p.adjust.method = "BH") #ME-CHMI,PUPA-LEMA

# Plot8
i=8
homo_data_i <- homo_data %>% filter(Plot==i)
kruskal.test(motif ~ Plant_Simple, data = homo_data_i) #significant
pairwise.wilcox.test(homo_data_i$motif, homo_data_i$Plant_Simple,
                     p.adjust.method = "BH") #LEMA-CHFU

# Plot9
i=9
homo_data_i <- homo_data %>% filter(Plot==i)
kruskal.test(motif ~ Plant_Simple, data = homo_data_i) #significant
pairwise.wilcox.test(homo_data_i$motif, homo_data_i$Plant_Simple,
                     p.adjust.method = "BH") #CHFU-ME

#########
#HETEROMOTIFS

means_hete <- aggregate((motif) ~  Plant_Simple + Plot,
                   hete_data, mean)

ggplot(hete_data,aes(x=Plant_Simple,y=motif))+
  geom_boxplot()+
  geom_point(aes(color=Plant_Simple), position = "jitter", alpha=0.2)+
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = means_hete, aes(label = round(`(motif)`,2), y = 10))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("# Hete-motifs")+
  theme_bw()+ theme(legend.position = "none")#+stat_compare_means()
#labs(Color

##########
# Test significance

# Test significance

# Plot1
i=1
hete_data_i <- hete_data %>% filter(Plot==i)#Significant
kruskal.test(motif ~ Plant_Simple, data = hete_data_i)
pairwise.wilcox.test(hete_data_i$motif, hete_data_i$Plant_Simple,
                     p.adjust.method = "BH")#No sig: LEMA-CHFU, ME-CHFU

# Plot2
i=2
hete_data_i <- hete_data %>% filter(Plot==i) # No Significant
kruskal.test(motif ~ Plant_Simple, data = hete_data_i)
pairwise.wilcox.test(hete_data_i$motif, hete_data_i$Plant_Simple,
                     p.adjust.method = "BH") 
# Plot3
i=3
hete_data_i <- hete_data %>% filter(Plot==i)
kruskal.test(motif ~ Plant_Simple, data = hete_data_i)
pairwise.wilcox.test(hete_data_i$motif, hete_data_i$Plant_Simple,
                     p.adjust.method = "BH") #LEMA-CHFU
# Plot4
i=4
hete_data_i <- hete_data %>% filter(Plot==i) #No significant
kruskal.test(motif ~ Plant_Simple, data = hete_data_i)
pairwise.wilcox.test(hete_data_i$motif, hete_data_i$Plant_Simple,
                     p.adjust.method = "BH")
# Plot5
i=5
hete_data_i <- hete_data %>% filter(Plot==i) # NO Significant
kruskal.test(motif ~ Plant_Simple, data = hete_data_i)
pairwise.wilcox.test(hete_data_i$motif, hete_data_i$Plant_Simple,
                     p.adjust.method = "BH")
# Plot6
i=6
hete_data_i <- hete_data %>% filter(Plot==i) 
kruskal.test(motif ~ Plant_Simple, data = hete_data_i) # NO significant
pairwise.wilcox.test(hete_data_i$motif, hete_data_i$Plant_Simple,
                     p.adjust.method = "BH") 
# Plot7
i=7
hete_data_i <- hete_data %>% filter(Plot==i)
kruskal.test(motif ~ Plant_Simple, data = hete_data_i) #Significant
pairwise.wilcox.test(hete_data_i$motif, hete_data_i$Plant_Simple,
                     p.adjust.method = "BH") #Solo ME y PUPA LEMA son signif

# Plot8
i=8
hete_data_i <- hete_data %>% filter(Plot==i)
kruskal.test(motif ~ Plant_Simple, data = hete_data_i) #significant
pairwise.wilcox.test(hete_data_i$motif, hete_data_i$Plant_Simple,
                     p.adjust.method = "BH") #Solo pupa signific

# Plot9
i=9
hete_data_i <- hete_data %>% filter(Plot==i)
kruskal.test(motif ~ Plant_Simple, data = hete_data_i) #NO significant
pairwise.wilcox.test(hete_data_i$motif, hete_data_i$Plant_Simple,
                     p.adjust.method = "BH")


homo_data2 <- homo_data %>% rename(homo_motif = motif) %>% select(-type)
hete_data2 <- hete_data %>% rename(hete_motif = motif) %>% select(-type)

data_total2 <- homo_data2 %>% left_join(hete_data2,by=c("Plot","Sub","Plant_Simple"))
data_total2$Plant_Simple <- as.factor(data_total2$Plant_Simple)

# Sanity check
data_total2 %>% group_by(Plot,Sub,Plant_Simple) %>% count() %>% filter(n>1)

library(RColorBrewer)
ggplot(data_total2,aes(x=homo_motif,y = hete_motif, color = Plant_Simple, shape = Plant_Simple))+
  geom_point(aes(color=Plant_Simple), position = "jitter")+
  scale_shape_manual(values=1:nlevels(data_total2$Plant_Simple))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  geom_abline(aes(slope=1,intercept=0),linetype = "dashed")+
  scale_color_brewer(palette="Paired")+
  #ggtitle(paste0("Plot ",i)) +
  xlab("# Homospecific motifs") + ylab("# Heterospecific motifs")+
  theme_bw()+
  labs(color=NULL,shape=NULL)+
  theme(legend.position = "bottom")
#labs(Color

library(scales)

data_total2_exp <- data_total2
data_total2_exp$Plant_Simple <- as.character(data_total2_exp$Plant_Simple)
data_total2_exp$Plant_Simple[data_total2_exp$Plant_Simple == "BEMA"] <- "B. macrocarpa"
data_total2_exp$Plant_Simple[data_total2_exp$Plant_Simple == "CETE"] <- "C. tenuiflorum"
data_total2_exp$Plant_Simple[data_total2_exp$Plant_Simple == "CHFU"] <- "C. fuscatum"
data_total2_exp$Plant_Simple[data_total2_exp$Plant_Simple == "CHMI"] <- "C. mixtum"
data_total2_exp$Plant_Simple[data_total2_exp$Plant_Simple == "LEMA"] <- "L. maroccanus"
data_total2_exp$Plant_Simple[data_total2_exp$Plant_Simple == "MESU"] <- "M. sulcatus"
data_total2_exp$Plant_Simple[data_total2_exp$Plant_Simple == "PUPA"] <- "P. paludosa"
data_total2_exp$Plant_Simple[data_total2_exp$Plant_Simple == "SCLA"] <- "S. laciniata"
data_total2_exp$Plant_Simple[data_total2_exp$Plant_Simple == "SOAS"] <- "S. asper"
data_total2_exp$Plant_Simple[data_total2_exp$Plant_Simple == "SPRU"] <- "S. rubra"

data_total2_exp$Plant_Simple <- as.factor(data_total2_exp$Plant_Simple)

png("New_Figures/fig5.png", width=1961*2, height = 1961*2*300/600, res=300*2)
ggplot(data_total2_exp,aes(x=homo_motif+0.5,y = hete_motif+0.5))+
  geom_point(alpha=0.2, position = "jitter")+
  #scale_shape_manual(values=1:nlevels(data_total2$Plant_Simple))+
  #facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  facet_wrap(vars(Plant_Simple),nrow = 2,ncol = 5)+
  geom_abline(aes(slope=1,intercept=0),linetype = "dashed")+
  scale_color_brewer(palette="Paired")+
  #ggtitle(paste0("Plot ",i)) +
  xlab("# Homospecific subgraphs + 0.5") + ylab("# Heterospecific subgraphs + 0.5")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()+
  labs(color=NULL,shape=NULL)+
  theme(legend.position = "bottom")+ theme(strip.text = element_text(face = "italic"))
dev.off()

# save 600 x 300

plant <- "L. maroccanus"
homo_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(homo_motif) %>% pull()
hete_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(hete_motif) %>% pull()
wilcox.test(homo_plant, hete_plant, paired = TRUE) #p-value < 2.2e-16

plant <- "C. fuscatum"
homo_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(homo_motif) %>% pull()
hete_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(hete_motif) %>% pull()
wilcox.test(homo_plant, hete_plant, paired = TRUE) #p-value < 2.2e-16

plant <- "P. paludosa"
homo_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(homo_motif) %>% pull()
hete_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(hete_motif) %>% pull()
wilcox.test(homo_plant, hete_plant, paired = TRUE) #p-value = 0.1356

plant <- "B. macrocarpa"
homo_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(homo_motif) %>% pull()
hete_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(hete_motif) %>% pull()
wilcox.test(homo_plant, hete_plant, paired = TRUE) #p-value = 0.05791

plant <- "C. mixtum"
homo_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(homo_motif) %>% pull()
hete_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(hete_motif) %>% pull()
wilcox.test(homo_plant, hete_plant, paired = TRUE) #p-value = 0.07186

plant <- "C. tenuiflorum"
homo_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(homo_motif) %>% pull()
hete_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(hete_motif) %>% pull()
wilcox.test(homo_plant, hete_plant, paired = TRUE) #p-value = 0.0035

plant <- "M. sulcatus"
homo_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(homo_motif) %>% pull()
hete_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(hete_motif) %>% pull()
wilcox.test(homo_plant, hete_plant, paired = TRUE) #p-value = 0.5862

plant <- "S. asper"
homo_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(homo_motif) %>% pull()
hete_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(hete_motif) %>% pull()
wilcox.test(homo_plant, hete_plant, paired = TRUE) #p-value = 0.5

plant <- "S. laciniata"
homo_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(homo_motif) %>% pull()
hete_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(hete_motif) %>% pull()
wilcox.test(homo_plant, hete_plant, paired = TRUE) #p-value = 0.05906

plant <- "S. rubra"
homo_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(homo_motif) %>% pull()
hete_plant <- data_total2_exp %>% filter(Plant_Simple == plant) %>% 
  select(hete_motif) %>% pull()
wilcox.test(homo_plant, hete_plant, paired = TRUE) #p-value = 1

