
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

data_degree_observed <- pollination %>% dplyr::select(Plot,Subplot,Plant,ID_Simple) %>% unique() %>% group_by(Plot,Subplot,Plant) %>% count() %>%
  ungroup() %>% dplyr::select(Plot,n)

################################################################################
# Obtain intralayer degree from simulated data
################################################################################

# Load data obtained from simulations-------------
data_simulation_raw_interlinks <- read_csv("Processed_data/Data_simulation/data_changing_interlinks_V2.csv") %>% mutate(type="Changing interlinks")
data_simulation_raw_intralinks <- read_csv("Processed_data/Data_simulation/data_changing_intralinks_V2.csv") %>% mutate(type="Changing intralinks")
data_simulation_raw_pollinators <- read_csv("Processed_data/Data_simulation/data_changing_pollinators_V2.csv") %>% mutate(type="Changing pollinators")
data_simulation_raw_individuals <- read_csv("Processed_data/Data_simulation/data_changing_individuals_V2.csv") %>% mutate(type="Changing plant individuals")
data_simulation_raw_increasing_mix <- read_csv("Processed_data/Data_simulation/data_increasing_mix_V2.csv") %>% mutate(type="Increasing plant individuals and intralinks")
data_simulation_raw_decreasing_mix <- read_csv("Processed_data/Data_simulation/data_decreasing_mix_V2.csv") %>% mutate(type="Decreasing plant individuals and intralinks")

data_degree_simulation_intralinks <- data_simulation_raw_intralinks %>% 
  dplyr::select(Plot,Subplot,Plant,ID,type) %>% unique() %>% 
  group_by(Plot,Subplot,Plant,type) %>% count() %>%
  ungroup() %>% dplyr::select(Plot,type,n)

unique(data_degree_simulation_intralinks$Plot) # 300 random networks


data_degree_simulation_interlinks <- data_simulation_raw_interlinks %>% 
  dplyr::select(Plot,Subplot,Plant,ID,type) %>% unique() %>% 
  group_by(Plot,Subplot,Plant,type) %>% count() %>%
  ungroup() %>% dplyr::select(Plot,type,n)

data_degree_simulation_pollinators <- data_simulation_raw_pollinators %>% 
  dplyr::select(Plot,Subplot,Plant,ID,type) %>% unique() %>% 
  group_by(Plot,Subplot,Plant,type) %>% count() %>%
  ungroup() %>% dplyr::select(Plot,type,n)

data_degree_simulation_individuals <- data_simulation_raw_individuals %>% 
  dplyr::select(Plot,Subplot,Plant,ID,type) %>% unique() %>% 
  group_by(Plot,Subplot,Plant,type) %>% count() %>%
  ungroup() %>% dplyr::select(Plot,type,n)

data_degree_simulation_increasing_mix <- data_simulation_raw_increasing_mix %>% 
  dplyr::select(Plot,Subplot,Plant,ID,type) %>% unique() %>% 
  group_by(Plot,Subplot,Plant,type) %>% count() %>%
  ungroup() %>% dplyr::select(Plot,type,n)

data_degree_simulation_decreasing_mix <- data_simulation_raw_decreasing_mix %>% 
  dplyr::select(Plot,Subplot,Plant,ID,type) %>% unique() %>% 
  group_by(Plot,Subplot,Plant,type) %>% count() %>%
  ungroup() %>% dplyr::select(Plot,type,n)



data_degree_simulation_intralinks$Plot <- paste0(data_degree_simulation_intralinks$Plot,
                                            25,"_pollinators")
data_degree_simulation_interlinks$Plot <- paste0(data_degree_simulation_interlinks$Plot,
                                            25,"_pollinators")
data_degree_simulation_individuals$Plot <- paste0(data_degree_simulation_individuals$Plot,
                                             25,"_pollinators")
data_degree_simulation_increasing_mix$Plot <- paste0(data_degree_simulation_increasing_mix$Plot,
                                                 25,"_pollinators")
data_degree_simulation_decreasing_mix$Plot <- paste0(data_degree_simulation_decreasing_mix$Plot,
                                                  25,"_pollinators")

data_degree_simulation_intralinks_deco <- decode_plot_label(data_degree_simulation_intralinks) %>% 
  dplyr::select(Plot,type,n,number_intra,number_individuals)
data_degree_simulation_interlinks_deco <- decode_plot_label(data_degree_simulation_interlinks) %>% 
  dplyr::select(Plot,type,n,number_inter)
data_degree_simulation_individuals_deco <- decode_plot_label(data_degree_simulation_individuals) %>% 
  dplyr::select(Plot,type,n,number_intra,number_individuals)
data_degree_simulation_pollinators_deco <- decode_plot_label(data_degree_simulation_pollinators) %>% 
  dplyr::select(Plot,type,n,number_pollinators)
data_degree_simulation_increasing_mix_deco <- decode_plot_label(data_degree_simulation_increasing_mix) %>% 
  dplyr::select(Plot,type,n,number_intra,number_individuals)
data_degree_simulation_decreasing_mix_deco <- decode_plot_label(data_degree_simulation_decreasing_mix) %>% 
  dplyr::select(Plot,type,n,number_intra,number_individuals)

max_degree <- max(max(data_degree_simulation_intralinks_deco$n),
                  max(data_degree_simulation_interlinks_deco$n),
                  max(data_degree_simulation_individuals_deco$n),
                  max(data_degree_simulation_pollinators_deco$n),
                  max(data_degree_simulation_increasing_mix_deco$n),
                  max(data_degree_simulation_decreasing_mix_deco$n))


# Plot observed VS simulated degree-------------
intralink_plot <- ggplot()+
  geom_line(data = data_degree_observed %>% filter(Plot==1), aes(x=n,color=Plot), stat="density", alpha=0.4, color = "blue",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==2), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "red",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==3), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "green",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==4), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "pink",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==5), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "orange",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==6), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "cyan",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==7), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "brown",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==8), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "magenta",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==9), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "tomato",linewidth = 1.5)+
  geom_line(data = data_degree_simulation_intralinks_deco %>%
              mutate(number_intra = paste0(number_intra," intralinks, ",number_individuals, " plant individuals")), aes(x=n,color=Plot), stat="density", alpha=0.25, show.legend = FALSE)+
  scale_color_manual(values =  rep("gray24",1000))+
  facet_wrap(~number_intra)+
  xlim(0,max_degree )+
  labs(x="Intralayer degree", y = "Density", title = unique(data_degree_simulation_intralinks$type))+
  theme_bw()

intralink_plot

intralink_plot_filtered <- ggplot()+
  geom_line(data = data_degree_observed %>% filter(Plot==1), aes(x=n,color=Plot), stat="density", alpha=0.4, color = "blue",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==2), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "red",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==3), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "green",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==4), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "pink",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==5), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "orange",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==6), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "cyan",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==7), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "brown",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==8), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "magenta",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==9), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "tomato",linewidth = 1.5)+
  geom_line(data = data_degree_simulation_intralinks_deco %>% filter(number_intra == 100) %>%
              mutate(number_intra = paste0(number_intra," intralinks, ",number_individuals, " plant individuals")), aes(x=n,color=Plot), stat="density", alpha=0.25, show.legend = FALSE)+
  scale_color_manual(values =  rep("gray24",1000))+
  facet_wrap(~number_intra)+
  xlim(0,10)+
  labs(x="Intralayer degree", y = "Density")+
  theme_bw()

intralink_plot_filtered

interlink_plot <- ggplot()+
  geom_line(data = data_degree_observed %>% filter(Plot==1), aes(x=n,color=Plot), stat="density", alpha=0.4, color = "blue",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==2), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "red",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==3), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "green",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==4), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "pink",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==5), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "orange",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==6), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "cyan",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==7), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "brown",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==8), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "magenta",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==9), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "tomato",linewidth = 1.5)+
  geom_line(data = data_degree_simulation_interlinks_deco %>%
              mutate(number_inter = paste0(number_inter," interlinks")), aes(x=n,color=Plot), stat="density", alpha=0.25, show.legend = FALSE)+
  scale_color_manual(values =  rep("gray24",1000))+
  facet_wrap(~number_inter)+
  xlim(0,max_degree )+
  labs(x="Intralayer degree", y = "Density", title = unique(data_degree_simulation_interlinks$type))+
  theme_bw()

interlink_plot


individuals_plot <- ggplot()+
  geom_line(data = data_degree_observed %>% filter(Plot==1), aes(x=n,color=Plot), stat="density", alpha=0.4, color = "blue",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==2), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "red",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==3), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "green",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==4), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "pink",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==5), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "orange",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==6), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "cyan",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==7), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "brown",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==8), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "magenta",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==9), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "tomato",linewidth = 1.5)+
  geom_line(data = data_degree_simulation_individuals_deco %>%
              mutate(number_individuals = paste0(number_intra," intralinks, ",number_individuals, " plant individuals")), aes(x=n,color=Plot), stat="density", alpha=0.25, show.legend = FALSE)+
  scale_color_manual(values =  rep("gray24",1000))+
  facet_wrap(~number_individuals)+
  xlim(0,max_degree )+
  labs(x="Intralayer degree", y = "Density", title = unique(data_degree_simulation_individuals$type))+
  theme_bw()

individuals_plot

individuals_plot_filtered <- ggplot()+
  geom_line(data = data_degree_observed %>% filter(Plot==1), aes(x=n,color=Plot), stat="density", alpha=0.4, color = "blue",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==2), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "red",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==3), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "green",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==4), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "pink",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==5), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "orange",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==6), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "cyan",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==7), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "brown",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==8), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "magenta",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==9), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "tomato",linewidth = 1.5)+
  geom_line(data = data_degree_simulation_individuals_deco %>% filter(number_individuals == 100) %>%
              mutate(number_individuals = paste0(number_intra," intralinks, ",number_individuals, " plant individuals")), aes(x=n,color=Plot), stat="density", alpha=0.25, show.legend = FALSE)+
  scale_color_manual(values =  rep("gray24",1000))+
  facet_wrap(~number_individuals)+
  xlim(0,10)+
  labs(x="Intralayer degree", y = "Density")+
  theme_bw()

individuals_plot_filtered


pollinators_plot <- ggplot()+
  geom_line(data = data_degree_observed %>% filter(Plot==1), aes(x=n,color=Plot), stat="density", alpha=0.4, color = "blue",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==2), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "red",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==3), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "green",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==4), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "pink",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==5), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "orange",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==6), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "cyan",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==7), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "brown",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==8), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "magenta",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==9), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "tomato",linewidth = 1.5)+
  geom_line(data = data_degree_simulation_pollinators_deco %>%
              mutate(number_pollinators = paste0(number_pollinators," pollinator sp.")), aes(x=n,color=Plot), stat="density", alpha=0.25, show.legend = FALSE)+
  scale_color_manual(values =  rep("gray24",1000))+
  facet_wrap(~number_pollinators)+
  xlim(0,max_degree )+
  labs(x="Intralayer degree", y = "Density", title = unique(data_degree_simulation_pollinators$type))+
  theme_bw()

pollinators_plot


increasing_mix_plot <- ggplot()+
  geom_line(data = data_degree_observed %>% filter(Plot==1), aes(x=n,color=Plot), stat="density", alpha=0.4, color = "blue",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==2), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "red",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==3), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "green",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==4), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "pink",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==5), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "orange",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==6), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "cyan",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==7), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "brown",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==8), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "magenta",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==9), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "tomato",linewidth = 1.5)+
  geom_line(data = data_degree_simulation_increasing_mix_deco %>%
              mutate(number_individuals = paste0(number_intra," intralinks, ",number_individuals, " plant individuals")), aes(x=n,color=Plot), stat="density", alpha=0.25, show.legend = FALSE)+
  scale_color_manual(values =  rep("gray24",1000))+
  facet_wrap(~number_individuals)+
  xlim(0,10 )+
  labs(x="Intralayer degree", y = "Density")+
  theme_bw()

increasing_mix_plot


decreasing_mix_plot <- ggplot()+
  geom_line(data = data_degree_observed %>% filter(Plot==1), aes(x=n,color=Plot), stat="density", alpha=0.4, color = "blue",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==2), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "red",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==3), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "green",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==4), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "pink",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==5), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "orange",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==6), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "cyan",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==7), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "brown",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==8), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "magenta",linewidth = 1.5)+
  geom_line(data = data_degree_observed %>% filter(Plot==9), aes(x=n,color=Plot), stat="density", alpha=0.4,  color = "tomato",linewidth = 1.5)+
  geom_line(data = data_degree_simulation_decreasing_mix_deco %>%
              mutate(number_individuals = paste0(number_intra," intralinks, ",number_individuals, " plant individuals")), aes(x=n,color=Plot), stat="density", alpha=0.25, show.legend = FALSE)+
  scale_color_manual(values =  rep("gray24",1000))+
  facet_wrap(~number_individuals)+
  xlim(0,10 )+
  labs(x="Intralayer degree", y = "Density")+
  theme_bw()

decreasing_mix_plot
library(patchwork)

intralink_plot / interlink_plot / individuals_plot / pollinators_plot

png("New_Figures/simulations_degree_distribution.png",
    width = 11.69*0.5, # The width of the plot in inches
    height = 11.69*0.4, units = "in", res=300*2)
(intralink_plot_filtered & decreasing_mix_plot) / (increasing_mix_plot & individuals_plot_filtered)
dev.off()
