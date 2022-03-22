
library(tidyverse)
library(RColorBrewer)
library(patchwork)

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


# Load data
raw_intralinks_CENTRALITY <- read_csv("Processed_data/Data_simulation/data_changing_intralinks_CENTRALITY_V2.csv")
raw_interlinks_CENTRALITY <- read_csv("Processed_data/Data_simulation/data_changing_interlinks_CENTRALITY_V2.csv")
raw_individuals_CENTRALITY <- read_csv("Processed_data/Data_simulation/data_changing_individuals_CENTRALITY_V2.csv")
raw_pollinators_CENTRALITY <- read_csv("Processed_data/Data_simulation/data_changing_pollinators_CENTRALITY_V2.csv")



total_PR_intralinks_CENTRALITY <- raw_intralinks_CENTRALITY[grep("ind",
                                                                   raw_intralinks_CENTRALITY$species,
                                                                   ignore.case = T),] %>%
  separate(species,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,Real_PR_Multi)

total_PR_interlinks_CENTRALITY <- raw_interlinks_CENTRALITY[grep("ind",
                                                                 raw_interlinks_CENTRALITY$species,
                                                                   ignore.case = T),] %>%
  separate(species,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,Real_PR_Multi)

total_PR_individuals_CENTRALITY <- raw_individuals_CENTRALITY[grep("ind",
                                                                   raw_individuals_CENTRALITY$species,
                                                                    ignore.case = T),] %>%
  separate(species,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,Real_PR_Multi)


total_PR_pollinators_CENTRALITY <- raw_pollinators_CENTRALITY[grep("ind",
                                                                   raw_pollinators_CENTRALITY$species,
                                                                    ignore.case = T),] %>%
  separate(species,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,Real_PR_Multi)


total_PR_intralinks_CENTRALITY$Plot <- paste0(total_PR_intralinks_CENTRALITY$Plot,
                                          10,"_pollinators")
total_PR_interlinks_CENTRALITY$Plot <- paste0(total_PR_interlinks_CENTRALITY$Plot,
                                             10,"_pollinators")
total_PR_individuals_CENTRALITY$Plot <- paste0(total_PR_individuals_CENTRALITY$Plot,
                                             10,"_pollinators")

total_PR_intralinks_CENTRALITY_deco <- decode_plot_label(total_PR_intralinks_CENTRALITY)
total_PR_interlinks_CENTRALITY_deco <- decode_plot_label(total_PR_interlinks_CENTRALITY)
total_PR_individuals_CENTRALITY_deco <- decode_plot_label(total_PR_individuals_CENTRALITY)
total_PR_pollinators_CENTRALITY_deco <- decode_plot_label(total_PR_pollinators_CENTRALITY)

summary_PR_intralinks <- total_PR_intralinks_CENTRALITY_deco %>%
  group_by(Plant,number_intra,number_inter,number_individuals,number_pollinators) %>% 
  summarise(Mean = mean(Real_PR_Multi, na.rm = T), SD = sd(Real_PR_Multi, na.rm = T))

summary_PR_intralinks$Plant[summary_PR_intralinks$Plant=="plant_sp_1"] <- "Sp. 1"
summary_PR_intralinks$Plant[summary_PR_intralinks$Plant=="plant_sp_2"] <- "Sp. 2"
summary_PR_intralinks$Plant[summary_PR_intralinks$Plant=="plant_sp_3"] <- "Sp. 3"
summary_PR_intralinks$Plant[summary_PR_intralinks$Plant=="plant_sp_4"] <- "Sp. 4"
summary_PR_intralinks$Plant[summary_PR_intralinks$Plant=="plant_sp_5"] <- "Sp. 5"

summary_PR_interlinks <- total_PR_interlinks_CENTRALITY_deco %>%
  group_by(Plant,number_intra,number_inter,number_individuals,number_pollinators) %>% 
  summarise(Mean = mean(Real_PR_Multi, na.rm = T), SD = sd(Real_PR_Multi, na.rm = T))

summary_PR_interlinks$Plant[summary_PR_interlinks$Plant=="plant_sp_1"] <- "Sp. 1"
summary_PR_interlinks$Plant[summary_PR_interlinks$Plant=="plant_sp_2"] <- "Sp. 2"
summary_PR_interlinks$Plant[summary_PR_interlinks$Plant=="plant_sp_3"] <- "Sp. 3"
summary_PR_interlinks$Plant[summary_PR_interlinks$Plant=="plant_sp_4"] <- "Sp. 4"
summary_PR_interlinks$Plant[summary_PR_interlinks$Plant=="plant_sp_5"] <- "Sp. 5"

summary_PR_individuals <- total_PR_individuals_CENTRALITY_deco %>%
  group_by(Plant,number_intra,number_inter,number_individuals,number_pollinators) %>% 
  summarise(Mean = mean(Real_PR_Multi, na.rm = T), SD = sd(Real_PR_Multi, na.rm = T))

summary_PR_individuals$Plant[summary_PR_individuals$Plant=="plant_sp_1"] <- "Sp. 1"
summary_PR_individuals$Plant[summary_PR_individuals$Plant=="plant_sp_2"] <- "Sp. 2"
summary_PR_individuals$Plant[summary_PR_individuals$Plant=="plant_sp_3"] <- "Sp. 3"
summary_PR_individuals$Plant[summary_PR_individuals$Plant=="plant_sp_4"] <- "Sp. 4"
summary_PR_individuals$Plant[summary_PR_individuals$Plant=="plant_sp_5"] <- "Sp. 5"

summary_PR_pollinators <- total_PR_pollinators_CENTRALITY_deco %>%
  group_by(Plant,number_intra,number_inter,number_individuals,number_pollinators) %>% 
  summarise(Mean = mean(Real_PR_Multi, na.rm = T), SD = sd(Real_PR_Multi, na.rm = T))

summary_PR_pollinators$Plant[summary_PR_pollinators$Plant=="plant_sp_1"] <- "Sp. 1"
summary_PR_pollinators$Plant[summary_PR_pollinators$Plant=="plant_sp_2"] <- "Sp. 2"
summary_PR_pollinators$Plant[summary_PR_pollinators$Plant=="plant_sp_3"] <- "Sp. 3"
summary_PR_pollinators$Plant[summary_PR_pollinators$Plant=="plant_sp_4"] <- "Sp. 4"
summary_PR_pollinators$Plant[summary_PR_pollinators$Plant=="plant_sp_5"] <- "Sp. 5"

p_intra <- ggplot(summary_PR_intralinks, aes(x=Plant, y=Mean, fill=Plant))+
  facet_wrap(~number_intra)+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  geom_errorbar(aes(ymin = Mean, ymax = Mean+SD),width = 0.2,
                position = position_dodge(0.9))+
  scale_fill_brewer(palette = "Paired")+theme_bw()+
  labs(title="Dependence of the number of PageRank on the number of intralinks\n(30 interlinks, 50 plant individuals, 10 pollinator sps.)")
  
p_inter <- ggplot(summary_PR_interlinks, aes(x=Plant, y=Mean, fill=Plant))+facet_wrap(~number_inter)+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  geom_errorbar(aes(ymin = Mean, ymax = Mean+SD),width = 0.2,
                position = position_dodge(0.9))+
  scale_fill_brewer(palette = "Paired")+theme_bw()+
  labs(title="Dependence of the number of PageRank on the number of interlinks\n(125 intralinks, 50 plant individuals, 10 pollinator sps.)")

p_indiv <- ggplot(summary_PR_individuals, aes(x=Plant, y=Mean, fill=Plant))+facet_wrap(~number_individuals)+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  geom_errorbar(aes(ymin = Mean, ymax = Mean+SD),width = 0.2,
                position = position_dodge(0.9))+
  scale_fill_brewer(palette = "Paired")+theme_bw()+
  labs(title="Dependence of the number of PageRank on the total number of plant individuals\n(125 intralinks, 30 interlinks, 10 pollinator sps.)")

p_polli <- ggplot(summary_PR_pollinators, aes(x=Plant, y=Mean, fill=Plant))+
  facet_wrap(~number_pollinators)+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  geom_errorbar(aes(ymin = Mean, ymax = Mean+SD),width = 0.2,
                position = position_dodge(0.9))+
  scale_fill_brewer(palette = "Paired")+theme_bw()+
  labs(title="Dependence of the number of PageRank on the total number of pollinator sp.\n(125 intralinks, 30 interlinks, 50 plant individuals)")



png("New_Figures/simulations_PAGERANK.png",
    width = 11.69*0.75, # The width of the plot in inches
    height = 11.69, units = "in", res=300*2)
p_indiv/p_intra/p_inter/p_polli
dev.off()

