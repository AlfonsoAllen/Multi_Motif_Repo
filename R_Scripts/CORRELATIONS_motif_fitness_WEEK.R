# load libraries
library(tidyverse)

####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################

fitness_data2 <- read_csv("Raw_data/Metadata_Pollinators_Abundances_Seeds_2019.csv")

fitness2 <- fitness_data2 %>% filter(Year==2019,!Subplot == "OUT" & !is.na(G_F))

fitness2 <- fitness2 %>% select(-Order,-Family,-Superfamily,-ID2) %>%
  mutate(date_raw=as.Date(paste(Day,Month,Year,sep="/"), "%d/%m/%Y"),
         Week=as.numeric(format(date_raw, "%V")))


#####################################
# Uploading motifs data
#####################################

caracoles_motif <- read_csv("Processed_data/Motifs_WEEK/Caracoles_WEEK.csv")

#Adding Subplot and Plant_Simple
for (i in 1:nrow(caracoles_motif)){
  x <- strsplit(caracoles_motif$Subplot_Plant_Label[i], " ")
  caracoles_motif$Subplot[i] <- x[[1]][1]
  caracoles_motif$Plant_Simple[i] <- x[[1]][2]
}

caracoles_motif <- caracoles_motif %>% select(-Visits_tot,-Subplot_Plant_Label)

#####################################
# Merging motifs data and fitness
#####################################

fitness <- caracoles_motif %>% left_join(fitness2, by=c("Plot","Subplot","Plant_Simple","G_F","Week"))

########################################################################
# SOME MOTIF VS SEEDS GRAPHS
########################################################################

# Adding GF contributions

fitness_SUM_Seed2 <- fitness %>% group_by(Plot,Subplot,Plant_Simple) %>%
  summarize(Seeds_GF = mean(Seed),
            Homo_Sum=sum(homo_motif),
            Hete_Sum=sum(hete_motif))

fitness_SUM_Seed <- fitness_SUM_Seed2 #%>% filter(Week==13)

#-----
pdf("Processed_data/Motifs_WEEK/Examples_motifs_seed_correlation_graphs/Homo_Hete_Seed.pdf",
    width = 11.69, # The width of the plot in inches
    height = 8.27)

ggplot(fitness_SUM_Seed2)+
  geom_point(aes(x = (Homo_Sum + Hete_Sum),y =(Seeds_GF)))+
  geom_smooth(aes(x=(Homo_Sum + Hete_Sum),y=(Seeds_GF)),method = "lm",se = F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  labs(title="Aggregation of triplets \n (each point represents the result for a given subplot and plant)",
       x ="Sum(# homo + hetero triplets)", y = "Mean(# seeds)",color="Plant")

dev.off()

#-----
pdf("Processed_data/Motifs_WEEK/Examples_motifs_seed_correlation_graphs/Homo_Seed.pdf",
    width = 11.69, # The width of the plot in inches
    height = 8.27)


ggplot(fitness_SUM_Seed2)+
  geom_point(aes(x = (Homo_Sum),y =(Seeds_GF)))+
  geom_smooth(aes(x=(Homo_Sum),y=(Seeds_GF)),method = "lm",se = F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  labs(title="Aggregation of homospecific triplets \n (each point represents the result for a given subplot and plant)",
       x ="Sum(# homo triplets)", y = "Mean(# seeds)")

dev.off()


#-----
pdf("Processed_data/Motifs_WEEK/Examples_motifs_seed_correlation_graphs/Hete_Seed.pdf",
    width = 11.69, # The width of the plot in inches
    height = 8.27)


ggplot(fitness_SUM_Seed2)+
  geom_point(aes(x = (Hete_Sum),y =(Seeds_GF)))+
  geom_smooth(aes(x=(Hete_Sum),y=(Seeds_GF)),method = "lm",se = F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  labs(title="Aggregation of  heterospecific triplets \n (each point represents the result for a given subplot and plant)",
       x ="Sum(# hetero triplets)", y = "Mean(# seeds)")

dev.off()


############################
#-----
pdf("Processed_data/Motifs_WEEK/Examples_motifs_seed_correlation_graphs/Homo_Hete_Seed_plant.pdf",
    width = 11.69, # The width of the plot in inches
    height = 8.27)

ggplot(fitness_SUM_Seed2)+
  geom_point(aes(x = (Homo_Sum + Hete_Sum),y =(Seeds_GF),color = as.factor(Plant_Simple)))+
  geom_smooth(aes(x=(Homo_Sum + Hete_Sum),y=(Seeds_GF),color = as.factor(Plant_Simple)),method = "lm",se = F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  labs(title="Aggregation of triplets \n (each point represents the result for a given subplot and plant)",
       x ="Sum(# homo + hetero triplets)", y = "Mean(# seeds)",color="Plant")

dev.off()

#-----
pdf("Processed_data/Motifs_WEEK/Examples_motifs_seed_correlation_graphs/Homo_Seed_plant.pdf",
    width = 11.69, # The width of the plot in inches
    height = 8.27)


ggplot(fitness_SUM_Seed2)+
  geom_point(aes(x = (Homo_Sum),y =(Seeds_GF),color = as.factor(Plant_Simple)))+
  geom_smooth(aes(x=(Homo_Sum),y=(Seeds_GF),color = as.factor(Plant_Simple)),method = "lm",se = F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  labs(title="Aggregation of homospecific triplets \n (each point represents the result for a given subplot and plant)",
       x ="Sum(# homo triplets)", y = "Mean(# seeds)",color="Plant")

dev.off()


#-----
pdf("Processed_data/Motifs_WEEK/Examples_motifs_seed_correlation_graphs/Hete_Seed_plant.pdf",
    width = 11.69, # The width of the plot in inches
    height = 8.27)


ggplot(fitness_SUM_Seed2)+
  geom_point(aes(x = (Hete_Sum),y =(Seeds_GF),color = as.factor(Plant_Simple)))+
  geom_smooth(aes(x=(Hete_Sum),y=(Seeds_GF),color = as.factor(Plant_Simple)),method = "lm",se = F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  labs(title="Aggregation of  heterospecific triplets \n (each point represents the result for a given subplot and plant)",
       x ="Sum(# hetero triplets)", y = "Mean(# seeds)",color="Plant")

dev.off()

