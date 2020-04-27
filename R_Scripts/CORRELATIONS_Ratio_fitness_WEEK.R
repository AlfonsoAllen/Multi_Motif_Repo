# load libraries
library(tidyverse)

####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################

fitness_data2 <- read_csv("Raw_data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")

fitness2 <- fitness_data2 %>% filter(Year==2019)

fitness2 <- fitness2 %>% select(-Order,-Family,-Superfamily,-ID) %>% rename(ID=ID_Simple) %>%
  mutate(date_raw=as.Date(paste(Day,Month,Year,sep="/"), "%d/%m/%Y"),
         Week=as.numeric(format(date_raw, "%V")))


#####################################
# Uploading motifs data
#####################################

caracoles_motif <- read_csv("Processed_data/Motifs_WEEK/Caracoles_WEEK_SPECIES.csv")

#Adding Subplot and Plant_Simple
for (i in 1:nrow(caracoles_motif)){
  x <- strsplit(caracoles_motif$Subplot_Plant_Label[i], " ")
  caracoles_motif$Subplot[i] <- x[[1]][1]
  caracoles_motif$Plant_Simple[i] <- x[[1]][2]
}

caracoles_motif <- caracoles_motif %>% filter(!Plant_Simple%in%c("HOMA","Lysimachia_arvensis")) %>%
  select(-Visits_tot,-Subplot_Plant_Label)

#####################################
# Merging motifs data and fitness
#####################################

fitness <- caracoles_motif %>% left_join(fitness2, by=c("Plot","Subplot","Plant_Simple","ID","Week"))

# Adding GF contributions

fitness_SUM_Seed2 <- fitness %>% group_by(Plot,Subplot,Plant_Simple) %>%
  summarize(Seeds_GF = mean(Seed),
            Homo_Sum=sum(homo_motif),
            Hete_Sum=sum(hete_motif),
            Motif_Sum=sum(homo_motif+hete_motif)) %>% 
  mutate(Homo_ratio=100*Homo_Sum/Motif_Sum)

fitness_SUM_Seed2 <- filter(fitness_SUM_Seed2,Motif_Sum!=0)

fitness_SUM_Seed <- fitness_SUM_Seed2 #%>% filter(Week==13)

#####################################
# ADDING MODULARITY MEASSURES
#####################################

for (i in 1:9){
  
  
  file_i = paste0("Processed_data/Modularity_Pheno_Overlap/Modularity_Plot",i,".csv")
  modularity_i <- read_csv(file_i)
  
  modularity_i <- modularity_i %>% filter(type=="plant") %>%
    select(-node_id,-layer_id,-layer_name,-type) %>% separate(species,c("Subplot","Plant_Simple")," ")
  
  
  if(i==1){modularity <- modularity_i}else{modularity <- bind_rows(modularity,modularity_i)}
  
}


fitness_SUM_Seed <- fitness_SUM_Seed %>% left_join(modularity,by=c("Plot","Subplot","Plant_Simple"))


ggplot(fitness_SUM_Seed)+
  geom_point(aes(x = (module),y =(Seeds_GF),color = as.factor(Plant_Simple)))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  labs(title="Phen. Overlap Modules \n (each point represents the result for a given subplot and plant)",
       x ="Module", y = "Mean(# seeds)",color="Plant")

ggplot(fitness_SUM_Seed)+
  geom_boxplot(aes(x = as.factor(module),y =(Seeds_GF),color = as.factor(module)))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  labs(title="Phen. Overlap Modules \n (each point represents the result for a given subplot and plant)",
       x ="Module", y = "Mean(# seeds)",color="Plant")

########################################################################
# SOME MOTIF VS SEEDS GRAPHS
########################################################################

#-----
# pdf("Processed_data/Motifs_WEEK/Examples_motifs_seed_correlation_graphs/Homo_Hete_Seed_Plot.pdf",
#     width = 11.69, # The width of the plot in inches
#     height = 8.27)

ggplot(fitness_SUM_Seed2)+
  geom_point(aes(x = (Homo_ratio),y =(Seeds_GF)))+
  geom_smooth(aes(x=(Homo_ratio),y=(Seeds_GF)),method = "lm",se = F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  labs(title="Aggregation of triplets \n (each point represents the result for a given subplot and plant)",
       x ="% homotifs", y = "Mean(# seeds)",color="Plant")

# dev.off()

#-----
# pdf("Processed_data/Motifs_WEEK/Examples_motifs_seed_correlation_graphs/Homo_Seed_Plot.pdf",
#     width = 11.69, # The width of the plot in inches
#     height = 8.27)


############################
#-----
# pdf("Processed_data/Motifs_WEEK/Examples_motifs_seed_correlation_graphs/Homo_Hete_Seed_plant_Plot.pdf",
#     width = 11.69, # The width of the plot in inches
#     height = 8.27)

ggplot(fitness_SUM_Seed2)+
  geom_point(aes(x = (Homo_ratio),y =(Seeds_GF),color = as.factor(Plant_Simple)))+
  geom_smooth(aes(x=(Homo_ratio),y=(Seeds_GF),color = as.factor(Plant_Simple)),method = "lm",se = F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  labs(title="Aggregation of triplets \n (each point represents the result for a given subplot and plant)",
       x ="Sum(# homo + hetero triplets)", y = "Mean(# seeds)",color="Plant")


###############
###############



#-----
# pdf("Processed_data/Motifs_WEEK/Examples_motifs_seed_correlation_graphs/Homo_Hete_Seed.pdf",
#     width = 11.69, # The width of the plot in inches
#     height = 8.27)

ggplot(fitness_SUM_Seed2)+
  geom_point(aes(x = (Homo_ratio),y =(Seeds_GF)))+
  geom_smooth(aes(x=(Homo_ratio),y=(Seeds_GF)),method = "lm",se = F)+
  labs(title="Aggregation of triplets \n (each point represents the result for a given subplot and plant)",
       x ="Sum(# homo + hetero triplets)", y = "Mean(# seeds)",color="Plant")

# dev.off()

############################
#-----
# pdf("Processed_data/Motifs_WEEK/Examples_motifs_seed_correlation_graphs/Homo_Hete_Seed_plant.pdf",
#     width = 11.69, # The width of the plot in inches
#     height = 8.27)

ggplot(fitness_SUM_Seed2)+
  geom_point(aes(x = (Homo_ratio),y =(Seeds_GF),color = as.factor(Plant_Simple)))+
  geom_smooth(aes(x=(Homo_ratio),y=(Seeds_GF),color = as.factor(Plant_Simple)),method = "lm",se = F)+
  labs(title="Aggregation of triplets \n (each point represents the result for a given subplot and plant)",
       x ="Sum(# homo + hetero triplets)", y = "Mean(# seeds)",color="Plant")

# dev.off()
