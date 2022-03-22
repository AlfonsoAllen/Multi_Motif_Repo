
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
raw_intralinks_MOTIFS <- read_csv("Processed_data/Data_simulation/data_changing_intralinks_MOTIFS_V2.csv")
raw_interlinks_MOTIFS <- read_csv("Processed_data/Data_simulation/data_changing_interlinks_MOTIFS_V2.csv")
raw_individuals_MOTIFS <- read_csv("Processed_data/Data_simulation/data_changing_individuals_MOTIFS_V2.csv")
raw_pollinators_MOTIFS <- read_csv("Processed_data/Data_simulation/data_changing_pollinators_MOTIFS_V2.csv")


total_hete_intralinks_MOTIFS <- raw_intralinks_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=hete_motif) %>% rename(hete_motif=n)

total_hete_interlinks_MOTIFS <- raw_interlinks_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=hete_motif) %>% rename(hete_motif=n)

total_hete_individuals_MOTIFS <- raw_individuals_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=hete_motif) %>% rename(hete_motif=n)

total_hete_pollinators_MOTIFS <- raw_pollinators_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=hete_motif) %>% rename(hete_motif=n)


total_hete_intralinks_MOTIFS$Plot <- paste0(total_hete_intralinks_MOTIFS$Plot,
                                          10,"_pollinators")
total_hete_interlinks_MOTIFS$Plot <- paste0(total_hete_interlinks_MOTIFS$Plot,
                                             10,"_pollinators")
total_hete_individuals_MOTIFS$Plot <- paste0(total_hete_individuals_MOTIFS$Plot,
                                             10,"_pollinators")

total_hete_intralinks_MOTIFS_deco <- decode_plot_label(total_hete_intralinks_MOTIFS)
total_hete_interlinks_MOTIFS_deco <- decode_plot_label(total_hete_interlinks_MOTIFS)
total_hete_individuals_MOTIFS_deco <- decode_plot_label(total_hete_individuals_MOTIFS)
total_hete_pollinators_MOTIFS_deco <- decode_plot_label(total_hete_pollinators_MOTIFS)

summary_hete_intralinks <- total_hete_intralinks_MOTIFS_deco %>%
  group_by(Plant,number_intra,number_inter,number_individuals,number_pollinators) %>% 
  summarise(Mean = mean(hete_motif, na.rm = T), SD = sd(hete_motif, na.rm = T))

summary_hete_intralinks$Plant[summary_hete_intralinks$Plant=="plant_sp_1"] <- "Sp. 1"
summary_hete_intralinks$Plant[summary_hete_intralinks$Plant=="plant_sp_2"] <- "Sp. 2"
summary_hete_intralinks$Plant[summary_hete_intralinks$Plant=="plant_sp_3"] <- "Sp. 3"
summary_hete_intralinks$Plant[summary_hete_intralinks$Plant=="plant_sp_4"] <- "Sp. 4"
summary_hete_intralinks$Plant[summary_hete_intralinks$Plant=="plant_sp_5"] <- "Sp. 5"

summary_hete_interlinks <- total_hete_interlinks_MOTIFS_deco %>%
  group_by(Plant,number_intra,number_inter,number_individuals,number_pollinators) %>% 
  summarise(Mean = mean(hete_motif, na.rm = T), SD = sd(hete_motif, na.rm = T))

summary_hete_interlinks$Plant[summary_hete_interlinks$Plant=="plant_sp_1"] <- "Sp. 1"
summary_hete_interlinks$Plant[summary_hete_interlinks$Plant=="plant_sp_2"] <- "Sp. 2"
summary_hete_interlinks$Plant[summary_hete_interlinks$Plant=="plant_sp_3"] <- "Sp. 3"
summary_hete_interlinks$Plant[summary_hete_interlinks$Plant=="plant_sp_4"] <- "Sp. 4"
summary_hete_interlinks$Plant[summary_hete_interlinks$Plant=="plant_sp_5"] <- "Sp. 5"

summary_hete_individuals <- total_hete_individuals_MOTIFS_deco %>%
  group_by(Plant,number_intra,number_inter,number_individuals,number_pollinators) %>% 
  summarise(Mean = mean(hete_motif, na.rm = T), SD = sd(hete_motif, na.rm = T))

summary_hete_individuals$Plant[summary_hete_individuals$Plant=="plant_sp_1"] <- "Sp. 1"
summary_hete_individuals$Plant[summary_hete_individuals$Plant=="plant_sp_2"] <- "Sp. 2"
summary_hete_individuals$Plant[summary_hete_individuals$Plant=="plant_sp_3"] <- "Sp. 3"
summary_hete_individuals$Plant[summary_hete_individuals$Plant=="plant_sp_4"] <- "Sp. 4"
summary_hete_individuals$Plant[summary_hete_individuals$Plant=="plant_sp_5"] <- "Sp. 5"

summary_hete_pollinators <- total_hete_pollinators_MOTIFS_deco %>%
  group_by(Plant,number_intra,number_inter,number_individuals,number_pollinators) %>% 
  summarise(Mean = mean(hete_motif, na.rm = T), SD = sd(hete_motif, na.rm = T))

summary_hete_pollinators$Plant[summary_hete_pollinators$Plant=="plant_sp_1"] <- "Sp. 1"
summary_hete_pollinators$Plant[summary_hete_pollinators$Plant=="plant_sp_2"] <- "Sp. 2"
summary_hete_pollinators$Plant[summary_hete_pollinators$Plant=="plant_sp_3"] <- "Sp. 3"
summary_hete_pollinators$Plant[summary_hete_pollinators$Plant=="plant_sp_4"] <- "Sp. 4"
summary_hete_pollinators$Plant[summary_hete_pollinators$Plant=="plant_sp_5"] <- "Sp. 5"

p_intra <- ggplot(summary_hete_intralinks, aes(x=Plant, y=Mean, fill=Plant))+
  facet_wrap(~number_intra)+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  geom_errorbar(aes(ymin = Mean, ymax = Mean+SD),width = 0.2,
                position = position_dodge(0.9))+
  scale_fill_brewer(palette = "Paired")+theme_bw()+
  labs(title="Dependence of the number of heterospecific motifs on the number of intralinks\n(30 interlinks, 50 plant individuals , 10 pollinator sps.)")
  
p_inter <- ggplot(summary_hete_interlinks, aes(x=Plant, y=Mean, fill=Plant))+facet_wrap(~number_inter)+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  geom_errorbar(aes(ymin = Mean, ymax = Mean+SD),width = 0.2,
                position = position_dodge(0.9))+
  scale_fill_brewer(palette = "Paired")+theme_bw()+
  labs(title="Dependence of the number of heterospecific motifs on the number of interlinks\n(125 intralinks, 50 plant individuals , 10 pollinator sps.)")

p_indiv <- ggplot(summary_hete_individuals, aes(x=Plant, y=Mean, fill=Plant))+facet_wrap(~number_individuals)+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  geom_errorbar(aes(ymin = Mean, ymax = Mean+SD),width = 0.2,
                position = position_dodge(0.9))+
  scale_fill_brewer(palette = "Paired")+theme_bw()+
  labs(title="Dependence of the number of heterospecific motifs on the total number of plant individuals\n(125 intralinks, 30 interlinks, 10 pollinator sps.)")

p_polli <- ggplot(summary_hete_pollinators, aes(x=Plant, y=Mean, fill=Plant))+
  facet_wrap(~number_pollinators)+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  geom_errorbar(aes(ymin = Mean, ymax = Mean+SD),width = 0.2,
                position = position_dodge(0.9))+
  scale_fill_brewer(palette = "Paired")+theme_bw()+
  labs(title="Dependence of the number of heterospecific motifs on the total number of pollinator sp.\n(125 intralinks, 30 interlinks, 50 plant individuals)")



png("New_Figures/simulations_HETE.png",
    width = 11.69*0.75, # The width of the plot in inches
    height = 11.69, units = "in", res=300*2)
p_indiv/p_intra/p_inter/p_polli
dev.off()
