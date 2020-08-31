
library(tidyverse)


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


# read data ---------------------------------------------------------------

fitness_final_aux0 <- read.csv(file = "data_models_phenol_overlap.csv",
                         header = TRUE,
                         stringsAsFactors = FALSE)

PageRank_results <- read.csv(file = "PageRank_results.csv",
                              header = TRUE,
                              stringsAsFactors = FALSE)

PageRank_results1 <- PageRank_results %>% separate(Label,c("Subplot","Plant0")," ") %>% dplyr::select(-Plant0) 

# NOTE THAT PREVIOUS CENTRALITY MEASURES WERE WRONG BECAUSE THEY WEREN'T NORMALIZED

fitness_final_aux <- fitness_final_aux0 %>% 
  left_join(PageRank_results1,by=c("Plot","Subplot","Plant_Simple"))

fitness_final_aux[is.na(fitness_final_aux)] <- 0
fitness_final_aux$Ratio[fitness_final_aux$ID=="None"] <- 1
fitness_final_aux$pecentage_same_plants[fitness_final_aux$ID=="None"] <- 1

fitness_final_aux %>% filter(ID!="None") %>% count()
fitness_final_aux %>% filter(ID!="None") %>% count(wt=visits_GF)
fitness_final_aux %>% filter(Real_PR_Multi>0) %>% count()

# Add G_F

G_F_list <- read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv") %>%
  dplyr::select(G_F,ID_Simple) %>% rename(ID=ID_Simple) %>% unique()

G_F_list <- bind_rows(G_F_list,tibble(G_F="None",ID="None"))

# Fix "Odontomyia_sp."

G_F_list$G_F[G_F_list$ID=="Odontomyia_sp."] <- "Small_flies"
G_F_list <- unique(G_F_list)

# Sanity check
G_F_list %>% group_by(ID) %>% count() %>% filter(n>1)

fitness_orig <- fitness_final_aux %>% dplyr::left_join(G_F_list,by = "ID")

#write_csv(fitness_final,"data_models_phenol_overlap_Random_GF.csv")

# Turn ID, GF and Plot into factors
fitness_orig$Plot <- as.factor(fitness_orig$Plot)
fitness_orig$ID <- as.factor(fitness_orig$ID)
fitness_orig$G_F <- as.factor(fitness_orig$G_F)


# remove fitness = 0

fitness.data <- subset(fitness_orig,Seeds_GF > 0)
#fitness.data <- subset(fitness.data,ID!="None")

fitness.data %>% count(wt=visits_GF)

# class of every column
# fitness.data %>% map_chr(class)

fitness_LEMA <- subset(fitness.data,Plant_Simple == "LEMA")
fitness_PUPA <- subset(fitness.data,Plant_Simple == "PUPA")
fitness_CHFU <- subset(fitness.data,Plant_Simple == "CHFU")

#################################
# POLLINATORS
##################################
fitness_orig_pol <- fitness_orig %>% filter(ID!="None") %>% group_by(G_F) %>% count(wt=visits_GF) 
fitness_orig_pol %>% mutate(perc=100*n/sum(fitness_orig_pol$n))

fitness_orig_plant <- fitness_orig %>% filter(ID!="None") %>% group_by(Plant_Simple) %>% count(wt=visits_GF) 
fitness_orig_plant %>% mutate(perc=100*n/sum(fitness_orig_plant$n))

fitness_orig_plot <- fitness_orig %>% filter(ID!="None") %>% group_by(Plot) %>% count(wt=visits_GF) 
fitness_orig_plot %>% mutate(perc=100*n/sum(fitness_orig_plot$n))

fitness_LEMA_pol <- fitness_LEMA %>% filter(ID!="None") %>% group_by(G_F) %>% count() 
fitness_LEMA_pol %>% mutate(perc=100*n/sum(fitness_LEMA_pol$n))

fitness_CHFU_pol <- fitness_CHFU %>% filter(ID!="None") %>% group_by(G_F) %>% count() 
fitness_CHFU_pol %>% mutate(perc=100*n/sum(fitness_CHFU_pol$n))

fitness_PUPA_pol <- fitness_PUPA%>% filter(ID!="None")  %>% group_by(G_F) %>% count() 
fitness_PUPA_pol %>% mutate(perc=100*n/sum(fitness_PUPA_pol$n))

library(viridis)
library(RColorBrewer)

brewer.pal(11, name = "Paired")
display.brewer.pal(n = 11, name = "Paired")

# Change G_F labelling

fitness_orig_label <- fitness_orig
fitness_orig_label$G_F <- as.character(fitness_orig_label$G_F)

fitness_orig_label$G_F[fitness_orig_label$G_F=="Flower_beetles"] <- "Flower beetles"
fitness_orig_label$G_F[fitness_orig_label$G_F=="House_flies"] <- "House flies"
fitness_orig_label$G_F[fitness_orig_label$G_F=="Small_beetles"] <- "Small beetles"
fitness_orig_label$G_F[fitness_orig_label$G_F=="Small_flies"] <- "Small flies"

fitness_orig_label$G_F <- as.factor(fitness_orig_label$G_F)

ggplot(fitness_orig_label %>% filter(G_F!="None"), aes(fill=G_F, y=visits_GF, x=Plant_Simple)) + 
  geom_bar(position="stack", stat="identity")+ theme_bw()+
  scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Plant species", y = "Number of visits",fill=NULL)+ theme(legend.position="bottom")

ggplot(fitness_orig_label %>% filter(G_F!="None"), aes(fill=G_F, y=visits_GF, x=Plant_Simple)) + 
  geom_bar(position="stack", stat="identity")+ theme_bw()+
  scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Plant species", y = "Number of visits",fill=NULL)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  theme(legend.position="bottom")

# Stacked + percent
ggplot(fitness_orig %>% filter(G_F!="None"), aes(fill=G_F, y=visits_GF, x=Plant_Simple)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Plant species", y = "Percentage",fill=NULL)
  

###########################
# NUMBER OF PLANTS PER PLOT

plants_PLOT <- fitness_orig  %>% filter(visits_GF>0) %>% select(Plot,Subplot,Plant_Simple) %>% unique() %>%
  group_by(Plot,Plant_Simple) %>% count()



ggplot(plants_PLOT, aes(fill=Plant_Simple, y=n, x=Plot)) + 
  geom_bar(position="stack", stat="identity")+ theme_bw()+
  scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Plot", y = "Number of visited focal individuals",fill=NULL)+ theme(legend.position="bottom")








######################
ggplot(fitness_orig %>% filter(G_F!="None"),aes(x=(G_F),y=scale(Seeds_GF),fill=G_F))+
  geom_boxplot()+theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
  facet_wrap(vars(Plant_Simple),nrow = 3,ncol = 3)+
  scale_color_brewer(palette = 'Paired')+labs(x =NULL, y = "Seeds",fill=NULL)
  


ggplot(fitness_orig %>% filter(G_F!="None"),aes(x=(homo_motif),y=(Seeds_GF),color=G_F))+
  geom_point()+geom_smooth(method = 'lm',formula = y~x, se=F)+
  facet_wrap(vars(Plant_Simple),nrow = 3,ncol = 3)+
  scale_color_brewer(palette = 'Paired')



ggplot(fitness.data,aes(x=(homo_motif),y=entropy,color=Plant_Simple))+
  geom_point()+geom_smooth(method = 'lm',formula = y~x)

ggplot(fitness.data%>% filter(entropy>0),aes(x=(homo_motif),y=entropy,color=Plant_Simple))+
  geom_boxplot()+facet_wrap(vars(Plot),nrow = 3,ncol = 3)
ggplot(fitness.data,aes(x=as.factor(homo_motif),y=entropy,color=Plant_Simple))+
  geom_boxplot()+facet_wrap(vars(Plot),nrow = 3,ncol = 3)

ggplot(fitness.data,aes(x=as.factor(homo_motif),y=entropy,color=Plant_Simple))+
  geom_violin()

ggplot(fitness.data,aes(x=Real_PR_Multi,y=entropy,color=Plant_Simple))+
  geom_point()+geom_smooth(method = 'lm',formula = y~x)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)

ggplot(fitness.data,aes(x=Real_PR_Multi,y=entropy,color=Plant_Simple))+
  geom_point()+geom_smooth(method = 'lm',formula = y~x)
ggplot(fitness.data,aes(x=Real_PR_Multi,y=homo_motif,color=Plant_Simple))+
  geom_point()+geom_smooth(method = 'lm',formula = y~x)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)

ggplot(fitness.data,aes(x=Real_PR_Multi,y=homo_motif,color=Plant_Simple))+
  geom_point()+geom_smooth(method = 'lm',formula = y~x)

ggplot(fitness.data,aes(x=Real_PR_Multi,y=log(Seeds_GF),color=Plant_Simple))+
  geom_point()+geom_smooth(method = 'lm',formula = y~x)+
  facet_wrap(vars(Plant_Simple),nrow = 3,ncol = 3)

ggplot(fitness.data,aes(x=homo_motif,y=log(Seeds_GF),color=Plant_Simple))+
  geom_point()+geom_smooth(method = 'lm',formula = y~x)+
  facet_wrap(vars(G_F),nrow = 4,ncol = 3)


ggplot(fitness.data,aes(x=homo_motif,y=log(Seeds_GF),color=Plant_Simple))+
  geom_point()+geom_smooth(method = 'lm',formula = y~x)+
  facet_wrap(vars(ID),nrow = 4,ncol = 10)
