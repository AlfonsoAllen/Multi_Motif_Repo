
library(tidyverse)

decode_plot_label <- function(data_frame_info){
  
  data_frame_info$rep <- NA
  data_frame_info$number_intra <- NA
  data_frame_info$number_inter <- NA
  data_frame_info$number_individuals <- NA
  data_frame_info$number_pollinators <- NA
  
  for (i in 1:nrow(data_frame_info)) {
    info_bits <- strsplit(data_frame_info$Plot[i],"_")
    info_bits <- gsub("[^0-9.-]","", info_bits[[1]]) %>% as.numeric()
    data_frame_info$rep[i] <- info_bits[1]
    data_frame_info$number_intra[i] <- info_bits[2]
    data_frame_info$number_inter[i] <- info_bits[3]
    data_frame_info$number_individuals[i] <- info_bits[4]
    data_frame_info$number_pollinators[i] <- info_bits[5]
  }
  
  return(data_frame_info)
  
}

################################################################################
# Obtain intralayer degree from real data
################################################################################

# Load data
pollination <- read_csv2("Raw_Data/final_Pollinators_2020.csv")
# Remove points from ID names
pollination$ID <- sub("\\.", "", pollination$ID)
pollination$ID_Simple <- sub("\\.", "", pollination$ID_Simple)

# filter tabanidae
pollination <- pollination %>% filter(ID != "Tabanidae")


unique(pollination$ID_Simple[grep(" ",pollination$ID_Simple,ignore.case = T)])
#No labels with spaces -> Good!

data_visits_observed <- pollination %>% 
  dplyr::select(Plot,Subplot,Plant,ID_Simple,Visits) %>%
  group_by(Plot,Subplot,Plant) %>% count(wt=Visits) %>%
  ungroup() %>% dplyr::select(Plot,n)

################################################################################
# Obtain intralayer degree from simulated data
################################################################################

# Load data obtained from simulations-------------
data_simulation_raw_interlinks <- read_csv("Processed_data/Data_simulation/data_changing_interlinks_V2.csv") %>% mutate(type="Changing interlinks")
data_simulation_raw_intralinks <- read_csv("Processed_data/Data_simulation/data_changing_intralinks_V2.csv") %>% mutate(type="Changing intralinks")
data_simulation_raw_pollinators <- read_csv("Processed_data/Data_simulation/data_changing_pollinators_V2.csv") %>% mutate(type="Changing pollinators")
data_simulation_raw_individuals <- read_csv("Processed_data/Data_simulation/data_changing_individuals_V2.csv") %>% mutate(type="Changing plant individuals")

data_visits_simulation_intralinks <- data_simulation_raw_intralinks %>% 
  dplyr::select(Plot,Subplot,Plant,ID,type,Visits) %>% 
  group_by(Plot,Subplot,Plant,type) %>% count(wt=Visits) %>%
  ungroup() %>% dplyr::select(Plot,type,n)

unique(data_visits_simulation_intralinks$Plot) # 300 random networks


data_visits_simulation_interlinks <- data_simulation_raw_interlinks %>% 
  dplyr::select(Plot,Subplot,Plant,ID,type,Visits) %>% 
  group_by(Plot,Subplot,Plant,type) %>% count(wt=Visits) %>%
  ungroup() %>% dplyr::select(Plot,type,n)

data_visits_simulation_pollinators <- data_simulation_raw_pollinators %>% 
  dplyr::select(Plot,Subplot,Plant,ID,type,Visits) %>% 
  group_by(Plot,Subplot,Plant,type) %>% count(wt=Visits) %>%
  ungroup() %>% dplyr::select(Plot,type,n)

data_visits_simulation_individuals <- data_simulation_raw_individuals %>% 
  dplyr::select(Plot,Subplot,Plant,ID,type,Visits) %>% 
  group_by(Plot,Subplot,Plant,type) %>% count(wt=Visits) %>%
  ungroup() %>% dplyr::select(Plot,type,n)



data_visits_simulation_intralinks$Plot <- paste0(data_visits_simulation_intralinks$Plot,
                                            25,"_pollinators")
data_visits_simulation_interlinks$Plot <- paste0(data_visits_simulation_interlinks$Plot,
                                            25,"_pollinators")
data_visits_simulation_individuals$Plot <- paste0(data_visits_simulation_individuals$Plot,
                                             25,"_pollinators")

data_visits_simulation_intralinks_deco <- decode_plot_label(data_visits_simulation_intralinks) %>% 
  dplyr::select(Plot,type,n,number_intra)
data_visits_simulation_interlinks_deco <- decode_plot_label(data_visits_simulation_interlinks) %>% 
  dplyr::select(Plot,type,n,number_inter)
data_visits_simulation_individuals_deco <- decode_plot_label(data_visits_simulation_individuals) %>% 
  dplyr::select(Plot,type,n,number_individuals)
data_visits_simulation_pollinators_deco <- decode_plot_label(data_visits_simulation_pollinators) %>% 
  dplyr::select(Plot,type,n,number_pollinators)


max_visits <- max(max(data_visits_simulation_intralinks_deco$n),
                  max(data_visits_simulation_interlinks_deco$n),
                  max(data_visits_simulation_individuals_deco$n),
                  max(data_visits_simulation_pollinators_deco$n))


# Plot observed VS simulated degree-------------
intralink_plot <- ggplot()+
  geom_line(data = data_visits_observed %>% filter(Plot==1), aes(x=n,color=Plot), stat="density", alpha=0.4, color = "blue",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==2), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "red",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==3), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "green",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==4), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "pink",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==5), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "orange",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==6), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "cyan",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==7), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "brown",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==8), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "magenta",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==9), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "tomato",linewidth = 2)+
  geom_line(data = data_visits_simulation_intralinks_deco %>%
              mutate(number_intra = paste0(number_intra," intralinks")), aes(x=n,color=Plot), stat="density", alpha=0.25, show.legend = FALSE)+
  scale_color_manual(values =  rep("gray24",1000))+
  facet_wrap(~number_intra)+
  xlim(0,max_visits )+
  labs(x="Intralayer degree", y = "Density", title = unique(data_visits_simulation_intralinks$type))+
  theme_bw()

intralink_plot

interlink_plot <- ggplot()+
  geom_line(data = data_visits_observed %>% filter(Plot==1), aes(x=n,color=Plot), stat="density", alpha=0.4, color = "blue",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==2), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "red",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==3), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "green",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==4), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "pink",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==5), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "orange",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==6), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "cyan",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==7), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "brown",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==8), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "magenta",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==9), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "tomato",linewidth = 2)+
  geom_line(data = data_visits_simulation_interlinks_deco %>%
              mutate(number_inter = paste0(number_inter," interlinks")), aes(x=n,color=Plot), stat="density", alpha=0.25, show.legend = FALSE)+
  scale_color_manual(values =  rep("gray24",1000))+
  facet_wrap(~number_inter)+
  xlim(0,max_visits )+
  labs(x="Intralayer degree", y = "Density", title = unique(data_visits_simulation_interlinks$type))+
  theme_bw()

interlink_plot


individuals_plot <- ggplot()+
  geom_line(data = data_visits_observed %>% filter(Plot==1), aes(x=n,color=Plot), stat="density", alpha=0.4, color = "blue",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==2), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "red",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==3), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "green",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==4), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "pink",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==5), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "orange",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==6), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "cyan",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==7), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "brown",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==8), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "magenta",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==9), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "tomato",linewidth = 2)+
  geom_line(data = data_visits_simulation_individuals_deco %>%
              mutate(number_individuals = paste0(number_individuals," plant individuals")), aes(x=n,color=Plot), stat="density", alpha=0.25, show.legend = FALSE)+
  scale_color_manual(values =  rep("gray24",1000))+
  facet_wrap(~number_individuals)+
  xlim(0,max_visits )+
  labs(x="Intralayer degree", y = "Density", title = unique(data_visits_simulation_individuals$type))+
  theme_bw()

individuals_plot


pollinators_plot <- ggplot()+
  geom_line(data = data_visits_observed %>% filter(Plot==1), aes(x=n,color=Plot), stat="density", alpha=0.4, color = "blue",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==2), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "red",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==3), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "green",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==4), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "pink",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==5), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "orange",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==6), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "cyan",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==7), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "brown",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==8), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "magenta",linewidth = 2)+
  geom_line(data = data_visits_observed %>% filter(Plot==9), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "tomato",linewidth = 2)+
  geom_line(data = data_visits_simulation_pollinators_deco %>%
              mutate(number_pollinators = paste0(number_pollinators," pollinator sp.")), aes(x=n,color=Plot), stat="density", alpha=0.25, show.legend = FALSE)+
  scale_color_manual(values =  rep("gray24",1000))+
  facet_wrap(~number_pollinators)+
  xlim(0,max_visits )+
  labs(x="Intralayer degree", y = "Density", title = unique(data_visits_simulation_pollinators$type))+
  theme_bw()

pollinators_plot

library(patchwork)

intralink_plot / interlink_plot / individuals_plot / pollinators_plot
