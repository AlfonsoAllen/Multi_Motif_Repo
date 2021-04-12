
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

fitness_final_aux <- read.csv(file = "Processed_data/2020_NN_NEW_data_models_phenol_overlap.csv",
                         header = TRUE,
                         stringsAsFactors = FALSE)


fitness_final_aux %>% filter(ID!="None") %>% count()
fitness_final_aux %>% filter(ID=="None") %>% count()
fitness_final_aux %>% filter(ID!="None") %>% count(wt=visits_GF)
fitness_final_aux %>% filter(Real_PR_Multi>0) %>% count()

# Add G_F

G_F_list <- read_csv2("Raw_Data/final_Pollinators_2020.csv") %>%
  filter(ID != "Tabanidae") %>%
  dplyr::select(G_F,ID_Simple) %>% unique() %>% rename(ID=ID_Simple)

# Remove points from ID names
G_F_list$ID <- sub("\\.", "", G_F_list$ID)

G_F_list <- bind_rows(G_F_list,tibble(G_F="None",ID="None"))

G_F_list <- unique(G_F_list)

G_F_list$G_F %>% unique() %>% sort()

# Sanity check
G_F_list %>% group_by(ID) %>% count() %>% filter(n>1)

fitness_orig <- fitness_final_aux %>% dplyr::left_join(G_F_list,by = "ID")

#write_csv(fitness_final,"data_models_phenol_overlap_Random_GF.csv")

# Turn ID, GF and Plot into factors
fitness_orig$Plot <- as.factor(fitness_orig$Plot)
fitness_orig$ID <- as.factor(fitness_orig$ID)
fitness_orig$G_F <- as.factor(fitness_orig$G_F)


# remove fitness = 0
fitness.data <- subset(fitness_orig,Seeds_GF >= 0)
#fitness.data <- subset(fitness.data,ID!="None")

fitness.data %>% count(wt=visits_GF)

# class of every column
# fitness.data %>% map_chr(class)

fitness_LEMA <- subset(fitness.data,Plant == "LEMA")
fitness_PUPA <- subset(fitness.data,Plant == "PUPA")
fitness_CHFU <- subset(fitness.data,Plant == "CHFU")

#################################
# POLLINATORS
##################################
fitness_orig_pol <- fitness_orig %>% filter(ID!="None") %>% group_by(G_F) %>% count(wt=visits_GF) 
fitness_orig_pol %>% mutate(perc=100*n/sum(fitness_orig_pol$n)) %>% 
  arrange(desc(perc))

fitness_orig_plant <- fitness_orig %>% filter(ID!="None") %>% group_by(Plant) %>% count(wt=visits_GF) 
fitness_orig_plant %>% mutate(perc=100*n/sum(fitness_orig_plant$n))%>% 
  arrange(desc(perc))

fitness_orig_plot <- fitness_orig %>% filter(ID!="None") %>% group_by(Plot) %>% count(wt=visits_GF) 
fitness_orig_plot %>% mutate(perc=100*n/sum(fitness_orig_plot$n)) %>% 
  arrange(desc(perc))

fitness_LEMA_pol <- fitness_LEMA %>% filter(ID!="None") %>% group_by(G_F) %>% count() 
fitness_LEMA_pol %>% mutate(perc=100*n/sum(fitness_LEMA_pol$n))%>% 
  arrange(desc(perc))

fitness_CHFU_pol <- fitness_CHFU %>% filter(ID!="None") %>% group_by(G_F) %>% count() 
fitness_CHFU_pol %>% mutate(perc=100*n/sum(fitness_CHFU_pol$n))%>% 
  arrange(desc(perc))

fitness_PUPA_pol <- fitness_PUPA%>% filter(ID!="None")  %>% group_by(G_F) %>% count() 
fitness_PUPA_pol %>% mutate(perc=100*n/sum(fitness_PUPA_pol$n))%>% 
  arrange(desc(perc))

library(viridis)
library(RColorBrewer)

# Change G_F labelling

fitness_orig_label <- fitness_orig
fitness_orig_label$G_F <- as.character(fitness_orig_label$G_F)

unique(fitness_orig_label$G_F)

fitness_orig_label$G_F[fitness_orig_label$G_F=="Flower_beetles"] <- "Flower beetles"
fitness_orig_label$G_F[fitness_orig_label$G_F=="House_flies"] <- "House flies"
fitness_orig_label$G_F[fitness_orig_label$G_F=="Small_beetles"] <- "Small beetles"
fitness_orig_label$G_F[fitness_orig_label$G_F=="Big_beetles"] <- "Big beetles"
fitness_orig_label$G_F[fitness_orig_label$G_F=="Small_flies"] <- "Small flies"
fitness_orig_label$G_F[fitness_orig_label$G_F=="Solitary_bees"] <- "Solitary bees"
fitness_orig_label$G_F[fitness_orig_label$G_F=="Wasp"] <- "Wasps"

fitness_orig_label$G_F <- as.factor(fitness_orig_label$G_F)

#A solution is to use the function colorRampPalette() which can extend any list of colors:
  
# Define the number of colors you want
nb.cols <- 11
mycolors <- colorRampPalette(brewer.pal(11, "Paired"))(nb.cols)

library("scales")


fitness_orig_label_expl <- fitness_orig_label
fitness_orig_label_expl$Plant[fitness_orig_label_expl$Plant == "BEMA"] <- "B. macrocarpa"
fitness_orig_label_expl$Plant[fitness_orig_label_expl$Plant == "CETE"] <- "C. tenuiflorum"
fitness_orig_label_expl$Plant[fitness_orig_label_expl$Plant == "CHFU"] <- "C. fuscatum"
fitness_orig_label_expl$Plant[fitness_orig_label_expl$Plant == "CHMI"] <- "C. mixtum"
fitness_orig_label_expl$Plant[fitness_orig_label_expl$Plant == "LEMA"] <- "L. maroccanus"
fitness_orig_label_expl$Plant[fitness_orig_label_expl$Plant == "MESU"] <- "M. sulcatus"
fitness_orig_label_expl$Plant[fitness_orig_label_expl$Plant == "PUPA"] <- "P. paludosa"
fitness_orig_label_expl$Plant[fitness_orig_label_expl$Plant == "SCLA"] <- "S. laciniata"
fitness_orig_label_expl$Plant[fitness_orig_label_expl$Plant == "SOAS"] <- "S. asper"
fitness_orig_label_expl$Plant[fitness_orig_label_expl$Plant == "SPRU"] <- "S. rubra"

ggplot(fitness_orig_label_expl %>% filter(G_F!="None"), aes(fill=G_F, y=visits_GF, x=Plant)) + 
  geom_bar(position="stack", stat="identity")+ theme_bw()+
  #scale_fill_brewer(palette = 'Paired')+ 
  scale_fill_manual(values = mycolors) +
  #scale_y_continuous(trans = scales::pseudo_log_trans(base = 10,sigma = 100))+
  labs(x ="Plant species", y = "Number of visits",fill=NULL)+
  theme(legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(fitness_orig_label_expl %>% filter(G_F!="None"), aes(fill=G_F, y=visits_GF, x=Plant)) + 
  geom_bar(position="stack", stat="identity")+ theme_bw()+
  #scale_fill_brewer(palette = 'Paired')+ 
  scale_fill_manual(values = mycolors) +
  labs(x ="Plant species", y = "Number of visits",fill=NULL)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.text.x = element_text(face = "italic"))

# Stacked + percent
ggplot(fitness_orig %>% filter(G_F!="None"), aes(fill=G_F, y=visits_GF, x=Plant)) + 
  geom_bar(position="fill", stat="identity")+
  #scale_fill_brewer(palette = 'Paired')+ 
  scale_fill_manual(values = mycolors) +
  labs(x ="Plant species", y = "Percentage of visits",fill=NULL)+
  theme(legend.position="bottom")
  


###########################
# NUMBER OF PLANTS PER PLOT

plants_PLOT <- fitness_orig  %>% filter(visits_GF>0) %>% select(Plot,Subplot,Plant) %>% unique() %>%
  group_by(Plot,Plant) %>% count()

plants_PLOT_expl <- plants_PLOT
plants_PLOT_expl$Plant[plants_PLOT_expl$Plant == "BEMA"] <- "B. macrocarpa"
plants_PLOT_expl$Plant[plants_PLOT_expl$Plant == "CETE"] <- "C. tenuiflorum"
plants_PLOT_expl$Plant[plants_PLOT_expl$Plant == "CHFU"] <- "C. fuscatum"
plants_PLOT_expl$Plant[plants_PLOT_expl$Plant == "CHMI"] <- "C. mixtum"
plants_PLOT_expl$Plant[plants_PLOT_expl$Plant == "LEMA"] <- "L. maroccanus"
plants_PLOT_expl$Plant[plants_PLOT_expl$Plant == "MESU"] <- "M. sulcatus"
plants_PLOT_expl$Plant[plants_PLOT_expl$Plant == "PUPA"] <- "P. paludosa"
plants_PLOT_expl$Plant[plants_PLOT_expl$Plant == "SCLA"] <- "S. laciniata"
plants_PLOT_expl$Plant[plants_PLOT_expl$Plant == "SOAS"] <- "S. asper"
plants_PLOT_expl$Plant[plants_PLOT_expl$Plant == "SPRU"] <- "S. rubra"


ggplot(plants_PLOT_expl, aes(fill=Plant, y=n, x=Plot)) + 
  geom_bar(position="stack", stat="identity")+ theme_bw()+
  scale_fill_brewer(palette = 'Paired')+ 
  labs(x ="Plot", y = "Number of visited focal individuals",fill=NULL)+ theme(legend.position="bottom")+
  theme(legend.text = element_text(face = "italic"))

