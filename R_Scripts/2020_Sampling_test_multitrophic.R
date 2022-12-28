
# Some consultants are unaware of the difference between rarefaction and species
# accumulation curves, and often calculate and present rarefaction curves as species
# accumulation curves in fauna reports. Rarefaction curves are useful for comparing
# species richness values for different sampling efforts. Rarefaction cannot be used for
# extrapolation as it does not provide an estimate of asymptotic richness.

library(tidyverse)
library(vegan)
library(iNEXT)

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


####################################################################
# Loadind Plant-pollinator dataset (Caracoles): visits, abundances, seeds
####################################################################
fitness_data2 <- read_csv2("Raw_Data/final_Pollinators_2020.csv")

# Inspect ID, ID_Simple

fitness_data2$ID %>% unique() %>% sort()
fitness_data2$ID_Simple %>% unique() %>% sort()

fitness_data2$ID[grep(" ",fitness_data2$ID)] # No labels contain spaces
fitness_data2$ID_Simple[grep(" ",fitness_data2$ID_Simple)]# No labels contain spaces

fitness_data2$ID[grep("\\.",fitness_data2$ID)] # Labels contain dots
fitness_data2$ID_Simple[grep("\\.",fitness_data2$ID_Simple)]# Labels contain dots

# Remove points from ID names
fitness_data2$ID <- sub("\\.", "", fitness_data2$ID)
fitness_data2$ID_Simple <- sub("\\.", "", fitness_data2$ID_Simple)

# filter tabanidae
fitness_data2 <- fitness_data2 %>% filter(ID != "Tabanidae")

# Filtering & relabeling
fitness2 <- fitness_data2 %>% filter(!is.na(Plant),
                                     Plant!="0",
                                     Subplot!="OUT",
                                     Plant!="Ground")

# Transform day month to day of the year

fitness2$Day_Year <- paste(fitness2$Day,
                           fitness2$Month,
                           fitness2$Year,
                           sep = "/") %>% as.Date(format='%d/%m/%Y') %>%
  
  lubridate::yday()



counts_visitors_subplot <- fitness2 %>% mutate(Position = paste(Plot,Subplot,sep = "_")) %>% 
  group_by(Plot,ID_Simple,Day_Year) %>%
  count(wt=Visits) %>% rename(Visits_tot = n,ID = ID_Simple) %>%
  spread(ID,Visits_tot)



######################
# INSECTS PER PLOT

# Species accumulation curve by adding sites in a random order, for 100 permutations.

data_acc_counts_visitors_subplot <- NULL

places <- counts_visitors_subplot$Plot %>% unique() #Position is no valid because vegan
# is not able to handle plots with just one round of sampling

for(i in 1:length(places)){
  
  new_rounds <- counts_visitors_subplot %>% ungroup() %>%
    filter(Plot == places[i]) %>% select(Day_Year) %>% unique() %>% pull()
    
  count_i <- counts_visitors_subplot %>% ungroup() %>%
    filter(Plot == places[i]) %>%
    select(-Plot,-Day_Year)
    
  count_i[is.na(count_i)] <- 0
    
  ind <- sapply(count_i, function(x) sum(x==0)) != nrow(count_i)
  count_i <- count_i[,ind]
    
  accurve <- specaccum(count_i, method="random", permutations=100)
  accurve_coll <- specaccum(count_i, method="collector")
    
  data_aux <- data.frame(Sites = accurve$sites,
                          Richness = accurve$richness,
                          Richness_coll = accurve_coll$richness,
                          SD = accurve$sd,
                          Plot = places[i],
                         Day_Year = new_rounds)
    
  data_acc_counts_visitors_subplot <- bind_rows(data_acc_counts_visitors_subplot,
                                                data_aux)
  
}

# Collector sampling
ggplot(data_acc_counts_visitors_subplot) +
  geom_point(aes(x = Sites, y = Richness_coll, color = as.factor(Plot))) +
  geom_line(aes(x=Sites, y=Richness_coll,
                color = as.factor(Plot), 
                group = as.factor(Plot))) +
  facet_wrap(~ Plot, labeller=labeller(Plot= plot_labs))+
  theme_bw()+
  labs(y="Plant richness",x="Round number",color=NULL,
       title = "Richness of insect visitors accumulated as sampling effort increases, pooling the data for all plant species and subplots")


ggplot(data_acc_counts_visitors_subplot) +
  geom_point(aes(x = Day_Year, y = Richness_coll, color = as.factor(Plot))) +
  # geom_line(aes(x=Day_Year, y=Richness_coll,
  #               color = as.factor(Plot), 
  #               group = as.factor(Plot))) +
  facet_wrap(~ Plot, labeller=labeller(Plot= plot_labs))+
  theme_bw()+
  theme(legend.position = "none")+
  labs(y="Richness of insect visitors accumulated",x="Day of the year",color=NULL,
       title = "Richness of insect visitors accumulated as sampling effort increases, pooling the data for all plant species and subplots")


ggplot(data_acc_counts_visitors_subplot) +
  geom_point(aes(x = Sites, y = Richness, color = as.factor(Plot))) +
  geom_line(aes(x= Sites, y=Richness, color = as.factor(Plot),
                group = as.factor(Plot))) +
  geom_ribbon(aes(x= Sites,
                  ymin=(Richness-2*SD),ymax=(Richness+2*SD),
                  color = as.factor(Plot),
                  fill = as.factor(Plot)),alpha=0.2)+
  facet_wrap(~ Plot)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(y="Richness of insect visitors accumulated",x="Round number",color=NULL,
       title = "Richness of insect visitors accumulated as sampling effort increases, pooling the data for all plant species and subplots")


ggplot(data_acc_counts_visitors_subplot) +
  geom_point(aes(x = Day_Year, y = Richness, color = as.factor(Plot))) +
  geom_line(aes(x=Day_Year, y=Richness, color = as.factor(Plot),
                group = as.factor(Plot))) +
  geom_ribbon(aes(x=Day_Year,
                  ymin=(Richness-2*SD),ymax=(Richness+2*SD),
                  color = as.factor(Plot),
                  fill = as.factor(Plot)),alpha=0.2)+
  facet_wrap(~ Plot)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(y="Richness of insect visitors accumulated",x="Day of the year",color=NULL,
       title = "Richness of insect visitors accumulated as sampling effort increases, pooling the data for all plant species and subplots")

# Chao 2#
counts_visitors_subplot[is.na(counts_visitors_subplot)] <- 0
richness_aux <- counts_visitors_subplot %>% group_by(Plot) %>% summarise_all(sum) %>%
  select(-Day_Year)

richness_aux[is.na(richness_aux)] <- 0
richness_aux$r_obser <-  0
richness_aux$r_chao <-  0

for (j in 1:nrow(richness_aux)) {
  x <- as.numeric(richness_aux[j,2:(ncol(richness_aux)-2)])
  chao  <-  ChaoRichness(x, datatype = "abundance", conf = 0.95)
  richness_aux$r_obser[j] <-  chao$Observed
  richness_aux$r_chao[j] <-  chao$Estimator
}

ratio_chao <- tibble(x = rep(65,9), 
                     y = rep(20,9), 
                     label = round(100*richness_aux$r_obser/richness_aux$r_chao,2),
                     Plot = 1:9)

ggplot(data_acc_counts_visitors_subplot) +
  geom_point(aes(x = Day_Year, y = Richness_coll, color = as.factor(Plot))) +
  # geom_line(aes(x=Day_Year, y=Richness_coll,
  #               color = as.factor(Plot), 
  #               group = as.factor(Plot))) +
  geom_text(data = ratio_chao, 
            aes(x = x, y = y, label = label))+
  facet_wrap(~ Plot, labeller=labeller(Plot= plot_labs))+
  theme_bw()+
  theme(legend.position = "none")+
  labs(y="Richness of insect visitors accumulated",x="Day of the year",color=NULL,
       title = "Richness of insect visitors accumulated as sampling effort increases, pooling the data for all plant species and subplots")


##################################
# INTERACTIONS
##################################


counts_interactions_subplot <- fitness2 %>% 
  mutate(ID = paste(ID_Simple,Plant,sep = "_")) %>% 
  group_by(Plot,ID,Day_Year) %>%
  count(wt=Visits) %>% rename(Visits_tot = n) %>%
  spread(ID,Visits_tot)




# Interactions accumulation curve by adding sites in a random order, for 100 permutations.

data_acc_counts_interactions_subplot <- NULL

places <- counts_interactions_subplot$Plot %>% unique() #Position is no valid because vegan
# is not able to handle plots with just one round of sampling

for(i in 1:length(places)){
  
  new_rounds <- counts_interactions_subplot %>% ungroup() %>%
    filter(Plot == places[i]) %>% select(Day_Year) %>% unique() %>% pull()
  
  count_i <- counts_interactions_subplot %>% ungroup() %>%
    filter(Plot == places[i]) %>%
    select(-Plot,-Day_Year)
  
  count_i[is.na(count_i)] <- 0
  
  ind <- sapply(count_i, function(x) sum(x==0)) != nrow(count_i)
  count_i <- count_i[,ind]
  
  accurve <- specaccum(count_i, method="random", permutations=100)
  accurve_coll <- specaccum(count_i, method="collector")
  
  data_aux <- data.frame(Sites = accurve$sites,
                         Richness = accurve$richness,
                         Richness_coll = accurve_coll$richness,
                         SD = accurve$sd,
                         Plot = places[i],
                         Day_Year = new_rounds)
  
  data_acc_counts_interactions_subplot <- bind_rows(data_acc_counts_interactions_subplot,
                                                data_aux)
  
}

# Collector sampling
ggplot(data_acc_counts_interactions_subplot) +
  geom_point(aes(x = Sites, y = Richness_coll, color = as.factor(Plot))) +
  geom_line(aes(x=Sites, y=Richness_coll,
                color = as.factor(Plot), 
                group = as.factor(Plot))) +
  facet_wrap(~ Plot, labeller=labeller(Plot= plot_labs))+
  theme_bw()+
  labs(y="Interaction richness",x="Round number",color=NULL,
       title = "Richness of interactions accumulated as sampling effort increases, pooling the data for all plant species and subplots")


ggplot(data_acc_counts_interactions_subplot) +
  geom_point(aes(x = Day_Year, y = Richness_coll, color = as.factor(Plot))) +
  # geom_line(aes(x=Day_Year, y=Richness_coll,
  #               color = as.factor(Plot), 
  #               group = as.factor(Plot))) +
  facet_wrap(~ Plot, labeller=labeller(Plot= plot_labs))+
  theme_bw()+
  theme(legend.position = "none")+
  labs(y="Richness of interactions accumulated",x="Day of the year",color=NULL,
       title = "Richness of interactions accumulated as sampling effort increases, pooling the data for all plant species and subplots")


ggplot(data_acc_counts_interactions_subplot) +
  geom_point(aes(x = Sites, y = Richness, color = as.factor(Plot))) +
  geom_line(aes(x= Sites, y=Richness, color = as.factor(Plot),
                group = as.factor(Plot))) +
  geom_ribbon(aes(x= Sites,
                  ymin=(Richness-2*SD),ymax=(Richness+2*SD),
                  color = as.factor(Plot),
                  fill = as.factor(Plot)),alpha=0.2)+
  facet_wrap(~ Plot)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(y="Richness of interactions accumulated",x="Round number",color=NULL,
       title = "Richness of interactions accumulated as sampling effort increases, pooling the data for all plant species and subplots")


ggplot(data_acc_counts_interactions_subplot) +
  geom_point(aes(x = Day_Year, y = Richness, color = as.factor(Plot))) +
  geom_line(aes(x=Day_Year, y=Richness, color = as.factor(Plot),
                group = as.factor(Plot))) +
  geom_ribbon(aes(x=Day_Year,
                  ymin=(Richness-2*SD),ymax=(Richness+2*SD),
                  color = as.factor(Plot),
                  fill = as.factor(Plot)),alpha=0.2)+
  facet_wrap(~ Plot)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(y="Richness of interactions accumulated",x="Day of the year",color=NULL,
       title = "Richness of interactions accumulated as sampling effort increases, pooling the data for all plant species and subplots")

# Chao 2#
counts_interactions_subplot[is.na(counts_interactions_subplot)] <- 0
richness_int_aux <- counts_interactions_subplot %>% group_by(Plot) %>% summarise_all(sum) %>%
  select(-Day_Year)

richness_int_aux[is.na(richness_aux)] <- 0
richness_int_aux$r_obser <-  0
richness_int_aux$r_chao <-  0

for (j in 1:nrow(richness_int_aux)) {
  x <- as.numeric(richness_int_aux[j,2:(ncol(richness_int_aux)-2)])
  chao  <-  ChaoRichness(x, datatype = "abundance", conf = 0.95)
  richness_int_aux$r_obser[j] <-  chao$Observed
  richness_int_aux$r_chao[j] <-  chao$Estimator
}

ratio_int_chao <- tibble(x = rep(65,9), 
                     y = rep(30,9), 
                     label = round(100*richness_int_aux$r_obser/richness_int_aux$r_chao,2),
                     Plot = 1:9)

ggplot(data_acc_counts_interactions_subplot) +
  geom_point(aes(x = Day_Year, y = Richness_coll, color = as.factor(Plot))) +
  # geom_line(aes(x=Day_Year, y=Richness_coll,
  #               color = as.factor(Plot), 
  #               group = as.factor(Plot))) +
  geom_text(data = ratio_int_chao, 
            aes(x = x, y = y, label = label))+
  facet_wrap(~ Plot, labeller=labeller(Plot= plot_labs))+
  theme_bw()+
  theme(legend.position = "none")+
  labs(y="Richness of interactions accumulated",x="Day of the year",color=NULL,
       title = "Richness of interactions accumulated as sampling effort increases, pooling the data for all plant species and subplots")

#################
# RAREFACTION CURVES

x <- counts_interactions_subplot %>% ungroup() %>% select(-Day_Year) %>% 
  group_by(Plot) %>% summarise_all(sum)
col_names_x <- paste0("Plot ",x$Plot)

x <- x %>% ungroup() %>% select(-Plot)

rownames(x) <- col_names_x

x.list <- setNames(split(x, seq(nrow(x))), rownames(x))
str(x.list)
iNEXT(x.list, q=0, datatype="abundance")

iNEXT(unlist(x.list), q=0, datatype="abundance")

out <- iNEXT(unlist(x.list[[2]]), q=c(0), datatype="abundance", endpoint=400)
# Sample‐size‐based R/E curves, separating by "site""

ggiNEXT(out, type=2)+
  scale_shape_manual(values = rep(19,9))+
  theme_bw(base_size = 18) + 
  theme(legend.position="none")+
  labs(x="Number of interactions")

gb4 <- ggplot_build(g)
gb4$data[[1]]$size <- 2
gb4$data[[2]]$size <- 1
gt4 <- ggplot_gtable(gb4)

library(grid)
grid.draw(gt4)

png("New_Figures/figA31.png",
    width = 11.69*0.5, # The width of the plot in inches
    height = 8.27*0.65, units = "in", res=300*2)
grid.draw(gt4)
dev.off()
