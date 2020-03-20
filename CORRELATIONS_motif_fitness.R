# load libraries
library(tidyverse)

####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################

fitness_data <- read_csv2("FV_2019_FINAL_visits_abundances_seeds.csv")

fitness_data %>% group_by(G_F) %>% count()


fitness <- fitness_data %>% group_by(Plot,Subplot,Plant_Simple,G_F,num.plantas,Fruit,Seed) %>%
  count(wt=Visits) %>% rename(Visits_tot = n)

#####################################
# Filtering & relabeling
#####################################

fitness <- fitness %>% filter(!Subplot == "OUT" & !is.na(G_F))

fitness$Subplot_Plant_Label <- paste(fitness$Subplot,fitness$Plant_Simple,sep = " ")
fitness <- fitness %>% mutate(Seeds_tot = num.plantas*Seed)

#####################################
# Uploading motifs data
#####################################

caracoles_motif <- read_csv("Motifs/Caracoles.csv")

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

fitness <- fitness %>% left_join(caracoles_motif, by=c("Plot","Subplot","Plant_Simple","G_F"))

########################################################################
# SOME MOTIF VS SEEDS GRAPHS
########################################################################

# Adding GF contributions

fitness_SUM_Seed <- fitness %>% group_by(Plot,Subplot,Plant_Simple) %>%
  summarize(Seeds_GF = mean(Seed),
            Homo_Sum=sum(homo_motif),
            Hete_Sum=sum(hete_motif))


ggplot(fitness_SUM_Seed)+
  geom_point(aes(x=(Homo_Sum+Hete_Sum),y=(Seeds_GF),position = "jitter"))+
  geom_smooth(aes(x=(Homo_Sum+Hete_Sum),y=(Seeds_GF)),method = "lm",se = F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)

ggplot(fitness_SUM_Seed)+
  geom_point(aes(x=(Homo_Sum),y=(Seeds_GF)),position = "jitter")+
  geom_smooth(aes(x=(Homo_Sum),y=(Seeds_GF)),method = "lm",se = F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)

ggplot(fitness_SUM_Seed)+
  geom_point(aes(x=(Hete_Sum),y=(Seeds_GF)),position = "jitter")+
  geom_smooth(aes(x=(Hete_Sum),y=(Seeds_GF)),method = "lm",se = F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)




ggplot(fitness_SUM_Seed)+
  geom_point(aes(x=(Homo_Sum+Hete_Sum),y=(Seeds_GF),
                 color=Plant_Simple),position = "jitter")+
  geom_smooth(aes(x=(Homo_Sum+Hete_Sum),y=(Seeds_GF),
                  color=Plant_Simple),method = "lm",se = F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)

ggplot(fitness_SUM_Seed)+
  geom_point(aes(x=(Homo_Sum),y=(Seeds_GF),
                 color=Plant_Simple),position = "jitter")+
  geom_smooth(aes(x=(Homo_Sum),y=(Seeds_GF),
                  color=Plant_Simple),method = "lm",se = F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)

ggplot(fitness_SUM_Seed)+
  geom_point(aes(x=(Hete_Sum),y=(Seeds_GF),
                 color=Plant_Simple),position = "jitter")+
  geom_smooth(aes(x=(Hete_Sum),y=(Seeds_GF),
                  color=Plant_Simple),method = "lm",se = F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)




