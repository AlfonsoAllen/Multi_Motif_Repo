
library(tidyverse)


####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################

fitness_data2 <- read_csv("Raw_data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")

fitness2 <- fitness_data2 %>% filter(Year==2019)

fitness2 <- fitness2 %>% select(-Order,-Family,-Superfamily,-ID) %>% rename(ID=ID_Simple) %>%
  mutate(date_raw=as.Date(paste(Day,Month,Year,sep="/"), "%d/%m/%Y"),
         Week=as.numeric(format(date_raw, "%V")))


fitness2$Line <- NA

for (i in 1:nrow(fitness2 )){
  if(fitness2$Plot[i] %in% c(1,2,3)){fitness2$Line[i] <- 1}
  else if(fitness2$Plot[i] %in% c(4,5,6)){fitness2$Line[i] <- 2}
  else{fitness2$Line[i] <- 3}
} 

fitness3 <- fitness2 %>% group_by(Line,Plot,Subplot,Plant_Simple,Week,ID) %>%
  summarize(Seeds_GF = mean(Seed,na.rm = T),
            Fruit_GF = mean(Fruit,na.rm = T),
            visits_GF = sum(Visits,na.rm = T))

#####################################
# Uploading motifs data
#####################################

caracoles_motif <- read_csv("Processed_data/Motifs_WEEK/Caracoles_WEEK_SPECIES.csv")

caracoles_motif <- caracoles_motif %>% separate(Subplot_Plant_Label,c("Subplot",
                                                                      "Plant_Simple")," ")

#####################################
# Merging motifs data and fitness
#####################################

fitness_aux <- caracoles_motif %>% left_join(fitness3, by=c("Plot","Subplot","Plant_Simple","ID","Week"))

#Sanity check: visits from caracoles and fitness3 are equal
fitness_aux %>% filter(visits_GF!=Visits_tot)

#####################################
# COMPETITION
#####################################
competition <- read_csv2("Raw_Data/competition.txt")

competition$ME_iden <- NA
competition$ME_iden[competition$focal=="MEEL"] <- "MEEL" # we add this dummy variable to identify ME
competition$ME_iden[competition$focal=="MESU"] <- "MESU"
competition$focal[competition$focal=="MEEL"] <- "ME"
competition$focal[competition$focal=="MESU"] <- "ME"


#CHMI: 55 semillas planta
#MEEL (ME): 1 en compet. no hay en poll.
#MESU (ME): 2 en compet. 1 en poll.
#CHFU: 76 poll.
#LEMA: 64 poll.
#PUPA: 35 poll.

competition_fil <- competition %>% 
  filter(year==2019,focal %in% c("PUPA","LEMA","CHFU", "CHMI","ME")) %>%
  select(-neighbour,-number) %>% unique() %>% rename(Year=year,Plot=plot,Subplot=subplot,
                                                     Plant_Simple=focal,
                                                     Fruit_GF = fruit,
                                                     Seeds_GF = seed)

fitness <- competition_fil %>% full_join(fitness_aux,by=c("Plot","Subplot",
                                                          "Plant_Simple"))

#According to fitness dataframe, poll. dataset does not contain MEEL observations

for (i in 1:nrow(fitness)){
  
  if (is.na(fitness$Fruit_GF.y[i]) | is.nan(fitness$Fruit_GF.y[i])){
    fitness$Fruit_GF.y[i] <- fitness$Fruit_GF.x[i]
    if (fitness$Plant_Simple[i]=="CHMI"){fitness$Seeds_GF.y[i] <- 55*fitness$Fruit_GF.y[i]}
    else if(fitness$Plant_Simple[i]=="PUPA"){fitness$Seeds_GF.y[i] <- 35*fitness$Fruit_GF.y[i]}
    else if(fitness$Plant_Simple[i]=="LEMA"){fitness$Seeds_GF.y[i] <- 64*fitness$Fruit_GF.y[i]}
    else if(fitness$Plant_Simple[i]=="CHFU"){fitness$Seeds_GF.y[i] <- 76*fitness$Fruit_GF.y[i]}
    else{fitness$Seeds_GF.y[i] <- fitness$Seeds_GF.x[i]}
    #CHMI: 55 semillas planta
    #MEEL (ME): 1 en compet. no hay en poll.
    #MESU (ME): 2 en compet. 1 en poll.
    #CHFU: 76 poll.
    #LEMA: 64 poll.
    #PUPA: 35 poll.
  }
}

# Test number of fruits: There are 12 differences between compet. dataset and poll. dataset
fitness %>% select(Plot,Subplot,Plant_Simple, Fruit_GF.x,Fruit_GF.y)%>%
  filter(Fruit_GF.x!=Fruit_GF.y) %>% group_by(Plot,Subplot,Plant_Simple,Fruit_GF.x,Fruit_GF.y) %>%
  count()

# We use fruit and seed data from pollinator dataset

fitness_final <- fitness %>% select(Plot,Subplot,Plant_Simple,Seeds_GF.y,
                                    Fruit_GF.y,visits_GF,ID,homo_motif,
                                    hete_motif) %>%
  rename(Seeds_GF = Seeds_GF.y, Fruit_GF = Fruit_GF.y)

fitness_final$visits_GF[is.na(fitness_final$visits_GF)] <- 0
fitness_final$homo_motif[is.na(fitness_final$homo_motif)] <- 0
fitness_final$hete_motif[is.na(fitness_final$hete_motif)] <- 0
fitness_final$ID[is.na(fitness_final$ID)] <- "None"


str(fitness_final)

#9     E3      LEMA       Seeds_GF=NA     Fruit_GF=NA    visits_GF=3 Melyridae  0   0
fitness_final <- fitness_final %>% filter(!is.na(Seeds_GF))


#############################################
# ADDING CENTRALITY MEASSURES (PHENOLOGICAL OVERLAP)
#############################################

for (i in 1:9){
  
  file_i = paste0("C:/Users/USUARIO/Desktop/Multi_Motif_Repo/Processed_data/Muxviz_Pheno_Overlap/centrality_table_Plot",i)
  Centrality_i <- read_delim(file_i,";")
  
  # Remove data for eigenvectors
  Centrality_i <- Centrality_i %>% filter(Layer=="1-Multi") %>% mutate(Plot=i) %>%
    select(-Layer,-Node,-Eigenvector) %>% separate(Label,c("Subplot","Plant_Simple")," ")
  
  
  if(i==1){centrality <- Centrality_i}else{centrality <- bind_rows(centrality,Centrality_i)}
  
}

fitness_final <- fitness_final %>% left_join(centrality, by=c("Plot","Subplot","Plant_Simple"))

# Removing centrality NAs
fitness_final[is.na(fitness_final)] <- 0

fitness_final$Plot <- as.factor(fitness_final$Plot)
fitness_final$Subplot <- as.factor(fitness_final$Subplot)
fitness_final$ID <- as.factor(fitness_final$ID)
fitness_final$Plant_Simple <- as.factor(fitness_final$Plant_Simple)







#############################################
# EXPLORING THE STRUCTURAL DATA
##############################################


########################################################
# Scaling powerlaw between degree and strength
########################################################

# For a randomly distributed set of weights, it can be shown that a linear
# relation (exponent=1) exists. If the importance of a given node in the
# network is lower than predicted by its degree, we would observe
# exponent<1. In mutualistic webs, a superlinear behaviour is found (that is,
# exponent>1), indicating that species with many connections tend to display
# stronger interactions than

ggplot(centrality)+
  geom_point(aes(x=log10(DegreeIn),y=log10(StrengthIn),color=Plant_Simple))+
  geom_smooth(aes(x=log10(DegreeIn),y=log10(StrengthIn),color=Plant_Simple),method = "lm")+
  geom_path(aes(x=log10(DegreeIn),y=log10(DegreeIn)))


#######################################
# Asymmetry values (between A_ij and A_ji)
# the asymmetry of this pair is defined as
# phi(i, j) = abs(phi_ij-phi_ji)/max{phi_ij, phi_ji}.
#######################################

full_edge_list <- read_csv("Processed_data/Modularity_Pheno_Overlap/Edge_list_Phen_Over_PLOT.csv")

#Sanity check: There should be an even number of observations per plot
full_edge_list %>% group_by(Plot) %>% count() 

full_edge_list_mod <- full_edge_list %>%
  mutate(full_from = paste(Plot, layer_from, node_from, sep = " "),
         full_to = paste(Plot, layer_to, node_to, sep = " "),
         interlink = ifelse(layer_to!=layer_from,T,F)) %>%
  select(Plot,full_from,full_to,interlink,weight)

full_edge_list_mod$weight2 <- NA
full_edge_list_mod$phi <- NA

for(i in 1:nrow(full_edge_list_mod)){
  
  if(is.na(full_edge_list_mod$weight2[i])){

    index_reverse_i <- which(full_edge_list_mod$full_from==full_edge_list_mod$full_to[i]
          &full_edge_list_mod$full_to==full_edge_list_mod$full_from[i])
    
    full_edge_list_mod$weight2[i] <- full_edge_list_mod$weight[index_reverse_i]
    full_edge_list_mod$weight2[index_reverse_i] <- -1
    
    full_edge_list_mod$phi[i] <- 
      abs(full_edge_list_mod$weight[i]-full_edge_list_mod$weight2[i])/max(full_edge_list_mod$weight2[i],full_edge_list_mod$weight[i])
    
  }
  
}
full_edge_list_mod <- full_edge_list_mod %>% filter(weight2!=-1)

#Removing interlinks with weight and weight2 equal to zero, respectively.
full_edge_list_mod <- full_edge_list_mod %>% filter(!is.nan(phi))

# Frequency histogram for phi

ggplot(full_edge_list_mod, aes(x = phi)) +geom_histogram()

# Frequency histogram for phi, when only intralinks are considered
full_edge_list_mod_intra <- full_edge_list_mod %>% filter(!interlink)
ggplot(full_edge_list_mod_intra, aes(x = phi,fill=as.factor(Plot))) + 
  geom_histogram()

r <- hist(full_edge_list_mod_intra$phi)
log_hist <- tibble(breaks=r$breaks[-1],counts= r$counts)
ggplot(log_hist, aes(x = log10(breaks),y=log10(counts))) + geom_path()

       