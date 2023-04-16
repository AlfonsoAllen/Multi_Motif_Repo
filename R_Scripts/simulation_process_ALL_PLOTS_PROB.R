library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(scales)

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


# Load MOTIF data
raw_intralinks_MOTIFS <- read_csv("Processed_data/Data_simulation/data_changing_intralinks_MOTIFS_V2.csv")
raw_interlinks_MOTIFS <- read_csv("Processed_data/Data_simulation/data_changing_interlinks_MOTIFS_V2.csv")
raw_individuals_MOTIFS <- read_csv("Processed_data/Data_simulation/data_changing_individuals_MOTIFS_V2.csv")
raw_pollinators_MOTIFS <- read_csv("Processed_data/Data_simulation/data_changing_pollinators_MOTIFS_V2.csv")


total_homo_intralinks_MOTIFS <- raw_intralinks_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=homo_motif) %>% rename(homo_motif=n)

total_homo_interlinks_MOTIFS <- raw_interlinks_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=homo_motif) %>% rename(homo_motif=n)

total_homo_individuals_MOTIFS <- raw_individuals_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=homo_motif) %>% rename(homo_motif=n)

total_homo_pollinators_MOTIFS <- raw_pollinators_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=homo_motif) %>% rename(homo_motif=n)

total_hete_intralinks_MOTIFS <- raw_intralinks_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=hete_motif) %>% rename(hete_motif=n)

total_hete_interlinks_MOTIFS <- raw_interlinks_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=hete_motif) %>% rename(hete_motif=n)

total_hete_individuals_MOTIFS <- raw_individuals_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=hete_motif) %>% rename(hete_motif=n)

total_hete_pollinators_MOTIFS <- raw_pollinators_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=hete_motif) %>% rename(hete_motif=n)


total_homo_intralinks_MOTIFS$Plot <- paste0(total_homo_intralinks_MOTIFS$Plot,
                                            25,"_pollinators")
total_homo_interlinks_MOTIFS$Plot <- paste0(total_homo_interlinks_MOTIFS$Plot,
                                            25,"_pollinators")
total_homo_individuals_MOTIFS$Plot <- paste0(total_homo_individuals_MOTIFS$Plot,
                                             25,"_pollinators")
total_hete_intralinks_MOTIFS$Plot <- paste0(total_hete_intralinks_MOTIFS$Plot,
                                            25,"_pollinators")
total_hete_interlinks_MOTIFS$Plot <- paste0(total_hete_interlinks_MOTIFS$Plot,
                                            25,"_pollinators")
total_hete_individuals_MOTIFS$Plot <- paste0(total_hete_individuals_MOTIFS$Plot,
                                             25,"_pollinators")

total_homo_intralinks_MOTIFS_deco <- decode_plot_label(total_homo_intralinks_MOTIFS)
total_homo_interlinks_MOTIFS_deco <- decode_plot_label(total_homo_interlinks_MOTIFS)
total_homo_individuals_MOTIFS_deco <- decode_plot_label(total_homo_individuals_MOTIFS)
total_homo_pollinators_MOTIFS_deco <- decode_plot_label(total_homo_pollinators_MOTIFS)
total_hete_intralinks_MOTIFS_deco <- decode_plot_label(total_hete_intralinks_MOTIFS)
total_hete_interlinks_MOTIFS_deco <- decode_plot_label(total_hete_interlinks_MOTIFS)
total_hete_individuals_MOTIFS_deco <- decode_plot_label(total_hete_individuals_MOTIFS)
total_hete_pollinators_MOTIFS_deco <- decode_plot_label(total_hete_pollinators_MOTIFS)

total_homo_intralinks_MOTIFS_deco$Plant[total_homo_intralinks_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_homo_intralinks_MOTIFS_deco$Plant[total_homo_intralinks_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_homo_intralinks_MOTIFS_deco$Plant[total_homo_intralinks_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_homo_intralinks_MOTIFS_deco$Plant[total_homo_intralinks_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_homo_intralinks_MOTIFS_deco$Plant[total_homo_intralinks_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_homo_interlinks_MOTIFS_deco$Plant[total_homo_interlinks_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_homo_interlinks_MOTIFS_deco$Plant[total_homo_interlinks_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_homo_interlinks_MOTIFS_deco$Plant[total_homo_interlinks_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_homo_interlinks_MOTIFS_deco$Plant[total_homo_interlinks_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_homo_interlinks_MOTIFS_deco$Plant[total_homo_interlinks_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_homo_individuals_MOTIFS_deco$Plant[total_homo_individuals_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_homo_individuals_MOTIFS_deco$Plant[total_homo_individuals_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_homo_individuals_MOTIFS_deco$Plant[total_homo_individuals_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_homo_individuals_MOTIFS_deco$Plant[total_homo_individuals_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_homo_individuals_MOTIFS_deco$Plant[total_homo_individuals_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_homo_pollinators_MOTIFS_deco$Plant[total_homo_pollinators_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_homo_pollinators_MOTIFS_deco$Plant[total_homo_pollinators_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_homo_pollinators_MOTIFS_deco$Plant[total_homo_pollinators_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_homo_pollinators_MOTIFS_deco$Plant[total_homo_pollinators_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_homo_pollinators_MOTIFS_deco$Plant[total_homo_pollinators_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_hete_intralinks_MOTIFS_deco$Plant[total_hete_intralinks_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_hete_intralinks_MOTIFS_deco$Plant[total_hete_intralinks_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_hete_intralinks_MOTIFS_deco$Plant[total_hete_intralinks_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_hete_intralinks_MOTIFS_deco$Plant[total_hete_intralinks_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_hete_intralinks_MOTIFS_deco$Plant[total_hete_intralinks_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_hete_interlinks_MOTIFS_deco$Plant[total_hete_interlinks_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_hete_interlinks_MOTIFS_deco$Plant[total_hete_interlinks_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_hete_interlinks_MOTIFS_deco$Plant[total_hete_interlinks_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_hete_interlinks_MOTIFS_deco$Plant[total_hete_interlinks_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_hete_interlinks_MOTIFS_deco$Plant[total_hete_interlinks_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_hete_individuals_MOTIFS_deco$Plant[total_hete_individuals_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_hete_individuals_MOTIFS_deco$Plant[total_hete_individuals_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_hete_individuals_MOTIFS_deco$Plant[total_hete_individuals_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_hete_individuals_MOTIFS_deco$Plant[total_hete_individuals_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_hete_individuals_MOTIFS_deco$Plant[total_hete_individuals_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_hete_pollinators_MOTIFS_deco$Plant[total_hete_pollinators_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_hete_pollinators_MOTIFS_deco$Plant[total_hete_pollinators_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_hete_pollinators_MOTIFS_deco$Plant[total_hete_pollinators_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_hete_pollinators_MOTIFS_deco$Plant[total_hete_pollinators_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_hete_pollinators_MOTIFS_deco$Plant[total_hete_pollinators_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"

# Load PROBABILITY data -------
raw_intralinks_PROBABILITY <- read_csv("Processed_data/Data_simulation/data_changing_intralinks_PROBABILITIES.csv")
raw_interlinks_PROBABILITY <- read_csv("Processed_data/Data_simulation/data_changing_interlinks_PROBABILITIES.csv")
raw_individuals_PROBABILITY <- read_csv("Processed_data/Data_simulation/data_changing_individuals_PROBABILITIES.csv")
raw_pollinators_PROBABILITY <- read_csv("Processed_data/Data_simulation/data_changing_pollinators_PROBABILITIES.csv")

raw_intralinks_PROBABILITY_UNCOUPLED <- read_csv("Processed_data/Data_simulation/data_changing_intralinks_PROBABILITIES_UNCOUPLED.csv")
raw_interlinks_PROBABILITY_UNCOUPLED <- read_csv("Processed_data/Data_simulation/data_changing_interlinks_PROBABILITIES_UNCOUPLED.csv")
raw_individuals_PROBABILITY_UNCOUPLED <- read_csv("Processed_data/Data_simulation/data_changing_individuals_PROBABILITIES_UNCOUPLED.csv")
raw_pollinators_PROBABILITY_UNCOUPLED <- read_csv("Processed_data/Data_simulation/data_changing_pollinators_PROBABILITIES_UNCOUPLED.csv")


total_heter_prob_intralinks_PROBABILITY <- raw_intralinks_PROBABILITY[grep("ind",
                                                                    raw_intralinks_PROBABILITY$name,
                                                                    ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,heter_prob)

total_heter_prob_interlinks_PROBABILITY <- raw_interlinks_PROBABILITY[grep("ind",
                                                                    raw_interlinks_PROBABILITY$name,
                                                                    ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,heter_prob)

total_heter_prob_individuals_PROBABILITY <- raw_individuals_PROBABILITY[grep("ind",
                                                                      raw_individuals_PROBABILITY$name,
                                                                      ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,heter_prob)


total_heter_prob_pollinators_PROBABILITY <- raw_pollinators_PROBABILITY[grep("ind",
                                                                      raw_pollinators_PROBABILITY$name,
                                                                      ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,heter_prob)


total_consp_prob_UNCOUPLED_intralinks_PROBABILITY <- raw_intralinks_PROBABILITY_UNCOUPLED[grep("ind",
                                                                 raw_intralinks_PROBABILITY$name,
                                                                 ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,consp_prob_UNCOUPLED)

total_consp_prob_UNCOUPLED_interlinks_PROBABILITY <- raw_interlinks_PROBABILITY_UNCOUPLED[grep("ind",
                                                                 raw_interlinks_PROBABILITY$name,
                                                                 ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,consp_prob_UNCOUPLED)

total_consp_prob_UNCOUPLED_individuals_PROBABILITY <- raw_individuals_PROBABILITY_UNCOUPLED[grep("ind",
                                                                   raw_individuals_PROBABILITY$name,
                                                                   ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,consp_prob_UNCOUPLED)


total_consp_prob_UNCOUPLED_pollinators_PROBABILITY <- raw_pollinators_PROBABILITY_UNCOUPLED[grep("ind",
                                                                   raw_pollinators_PROBABILITY$name,
                                                                   ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,consp_prob_UNCOUPLED)


total_heter_prob_intralinks_PROBABILITY$Plot <- paste0(total_heter_prob_intralinks_PROBABILITY$Plot,
                                                 25,"_pollinators")
total_heter_prob_interlinks_PROBABILITY$Plot <- paste0(total_heter_prob_interlinks_PROBABILITY$Plot,
                                                 25,"_pollinators")
total_heter_prob_individuals_PROBABILITY$Plot <- paste0(total_heter_prob_individuals_PROBABILITY$Plot,
                                                  25,"_pollinators")

total_consp_prob_UNCOUPLED_intralinks_PROBABILITY$Plot <- paste0(total_consp_prob_UNCOUPLED_intralinks_PROBABILITY$Plot,
                                              25,"_pollinators")
total_consp_prob_UNCOUPLED_interlinks_PROBABILITY$Plot <- paste0(total_consp_prob_UNCOUPLED_interlinks_PROBABILITY$Plot,
                                              25,"_pollinators")
total_consp_prob_UNCOUPLED_individuals_PROBABILITY$Plot <- paste0(total_consp_prob_UNCOUPLED_individuals_PROBABILITY$Plot,
                                               25,"_pollinators")


total_heter_prob_intralinks_PROBABILITY_deco <- decode_plot_label(total_heter_prob_intralinks_PROBABILITY)
total_heter_prob_interlinks_PROBABILITY_deco <- decode_plot_label(total_heter_prob_interlinks_PROBABILITY)
total_heter_prob_individuals_PROBABILITY_deco <- decode_plot_label(total_heter_prob_individuals_PROBABILITY)
total_heter_prob_pollinators_PROBABILITY_deco <- decode_plot_label(total_heter_prob_pollinators_PROBABILITY)

total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco <- decode_plot_label(total_consp_prob_UNCOUPLED_intralinks_PROBABILITY)
total_consp_prob_UNCOUPLED_interlinks_PROBABILITY_deco <- decode_plot_label(total_consp_prob_UNCOUPLED_interlinks_PROBABILITY)
total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco <- decode_plot_label(total_consp_prob_UNCOUPLED_individuals_PROBABILITY)
total_consp_prob_UNCOUPLED_pollinators_PROBABILITY_deco <- decode_plot_label(total_consp_prob_UNCOUPLED_pollinators_PROBABILITY)

total_heter_prob_intralinks_PROBABILITY_deco$Plant[total_heter_prob_intralinks_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_heter_prob_intralinks_PROBABILITY_deco$Plant[total_heter_prob_intralinks_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_heter_prob_intralinks_PROBABILITY_deco$Plant[total_heter_prob_intralinks_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_heter_prob_intralinks_PROBABILITY_deco$Plant[total_heter_prob_intralinks_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_heter_prob_intralinks_PROBABILITY_deco$Plant[total_heter_prob_intralinks_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_heter_prob_interlinks_PROBABILITY_deco$Plant[total_heter_prob_interlinks_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_heter_prob_interlinks_PROBABILITY_deco$Plant[total_heter_prob_interlinks_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_heter_prob_interlinks_PROBABILITY_deco$Plant[total_heter_prob_interlinks_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_heter_prob_interlinks_PROBABILITY_deco$Plant[total_heter_prob_interlinks_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_heter_prob_interlinks_PROBABILITY_deco$Plant[total_heter_prob_interlinks_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_heter_prob_individuals_PROBABILITY_deco$Plant[total_heter_prob_individuals_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_heter_prob_individuals_PROBABILITY_deco$Plant[total_heter_prob_individuals_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_heter_prob_individuals_PROBABILITY_deco$Plant[total_heter_prob_individuals_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_heter_prob_individuals_PROBABILITY_deco$Plant[total_heter_prob_individuals_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_heter_prob_individuals_PROBABILITY_deco$Plant[total_heter_prob_individuals_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_heter_prob_pollinators_PROBABILITY_deco$Plant[total_heter_prob_pollinators_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_heter_prob_pollinators_PROBABILITY_deco$Plant[total_heter_prob_pollinators_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_heter_prob_pollinators_PROBABILITY_deco$Plant[total_heter_prob_pollinators_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_heter_prob_pollinators_PROBABILITY_deco$Plant[total_heter_prob_pollinators_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_heter_prob_pollinators_PROBABILITY_deco$Plant[total_heter_prob_pollinators_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_consp_prob_UNCOUPLED_interlinks_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_interlinks_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_consp_prob_UNCOUPLED_interlinks_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_interlinks_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_consp_prob_UNCOUPLED_interlinks_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_interlinks_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_consp_prob_UNCOUPLED_interlinks_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_interlinks_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_consp_prob_UNCOUPLED_interlinks_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_interlinks_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_consp_prob_UNCOUPLED_pollinators_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_pollinators_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_consp_prob_UNCOUPLED_pollinators_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_pollinators_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_consp_prob_UNCOUPLED_pollinators_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_pollinators_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_consp_prob_UNCOUPLED_pollinators_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_pollinators_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_consp_prob_UNCOUPLED_pollinators_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_pollinators_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

#----------------------------------

p_intra_box3_homo <- ggplot(total_homo_intralinks_MOTIFS_deco,
                       aes(x=number_intra, y=homo_motif, group=number_intra))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun =mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of intralinks")

p_inter_box3_homo <- ggplot(total_homo_interlinks_MOTIFS_deco,
                       aes(x=number_inter, y=homo_motif, group=number_inter))+
  # geom_jitter(aes(color=as.factor(number_inter)),shape=16, position=position_jitter(5.5), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of interlinks")

p_indiv_box3_homo <- ggplot(total_homo_individuals_MOTIFS_deco, 
                       aes(x=number_individuals, y=homo_motif, group=number_individuals))+
  # geom_jitter(aes(color=as.factor(number_individuals)),shape=16, position=position_jitter(3.8), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="Number of HoT",x="Number of individual plants", title = "Homospecific triplets (HoT)")

p_polli_box3_homo <- ggplot(total_homo_pollinators_MOTIFS_deco, 
                       aes(x=number_pollinators, y=homo_motif, group=number_pollinators))+
  # geom_jitter(aes(color=as.factor(number_pollinators)),shape=16, position=position_jitter(1.9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of pollinator species")

patchwork_homo <- p_indiv_box3_homo+p_polli_box3_homo+p_intra_box3_homo+p_inter_box3_homo+ 
  plot_layout(ncol = 4) + plot_layout(tag_level = "new")
  # plot_annotation(
  #   title = 'Homospecific motifs',
  #   subtitle = 'Parameters for central box plots: 125 intralinks, 30 interlinks, 50 plant individuals, 25 pollinator sps.'
  # )


#----------------------------------

p_intra_box3_hete <- ggplot(total_hete_intralinks_MOTIFS_deco,
                            aes(x=number_intra, y=hete_motif, group=number_intra))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of intralinks")

p_inter_box3_hete <- ggplot(total_hete_interlinks_MOTIFS_deco,
                            aes(x=number_inter, y=hete_motif, group=number_inter))+
  # geom_jitter(aes(color=as.factor(number_inter)),shape=16, position=position_jitter(5.5), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of interlinks")

p_indiv_box3_hete <- ggplot(total_hete_individuals_MOTIFS_deco, 
                            aes(x=number_individuals, y=hete_motif, group=number_individuals))+
  # geom_jitter(aes(color=as.factor(number_individuals)),shape=16, position=position_jitter(3.8), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="Number of HeT",x="Number of individual plants", title = "Heterospecific triplets (HeT)")

p_polli_box3_hete <- ggplot(total_hete_pollinators_MOTIFS_deco, 
                            aes(x=number_pollinators, y=hete_motif, group=number_pollinators))+
  # geom_jitter(aes(color=as.factor(number_pollinators)),shape=16, position=position_jitter(1.9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of pollinator species")

patchwork_hete <- p_indiv_box3_hete+p_polli_box3_hete+p_intra_box3_hete+p_inter_box3_hete+ 
  plot_layout(ncol = 4) + plot_layout(tag_level = "new")
  # plot_annotation(
  #   title = 'Heterospecific motifs',
  #   subtitle = 'Parameters for central box plots: 125 intralinks, 30 interlinks, 50 plant individuals, 25 pollinator sps.'
  # )

#----------------------------------

p_intra_box3_heter_prob_log <- ggplot(total_heter_prob_intralinks_PROBABILITY_deco,
                             aes(x=number_intra, y=heter_prob, group=number_intra))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of intralinks")

p_inter_box3_heter_prob_log <- ggplot(total_heter_prob_interlinks_PROBABILITY_deco,
                             aes(x=number_inter, y=heter_prob, group=number_inter))+
  # geom_jitter(aes(color=as.factor(number_inter)),shape=16, position=position_jitter(5.5), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of interlinks")

p_indiv_box3_heter_prob_log <- ggplot(total_heter_prob_individuals_PROBABILITY_deco, 
                             aes(x=number_individuals, y=heter_prob, group=number_individuals))+
  # geom_jitter(aes(color=as.factor(number_individuals)),shape=16, position=position_jitter(3.8), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="Prob. of receiving information\nfrom heterosp.",x="Number of individual plants", title = "Prob. of receiving information\nfrom heterosp.")

p_polli_box3_heter_prob_log <- ggplot(total_heter_prob_pollinators_PROBABILITY_deco, 
                             aes(x=number_pollinators, y=heter_prob, group=number_pollinators))+
  # geom_jitter(aes(color=as.factor(number_pollinators)),shape=16, position=position_jitter(1.9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of pollinator species")

patchwork_heter_prob_log <- p_indiv_box3_heter_prob_log + p_polli_box3_heter_prob_log + p_intra_box3_heter_prob_log + p_inter_box3_heter_prob_log + 
  plot_layout(ncol = 4) + plot_layout(tag_level = "new")
  # plot_annotation(
  #   title = 'In-Strength',
  #   subtitle = 'Parameters for central box plots: 125 intralinks, 30 interlinks, 50 plant individuals, 25 pollinator sps.'
  # )

p_intra_box3_heter_prob <- ggplot(total_heter_prob_intralinks_PROBABILITY_deco,
                                      aes(x=number_intra, y=heter_prob, group=number_intra))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of intralinks")

p_inter_box3_heter_prob <- ggplot(total_heter_prob_interlinks_PROBABILITY_deco,
                                      aes(x=number_inter, y=heter_prob, group=number_inter))+
  # geom_jitter(aes(color=as.factor(number_inter)),shape=16, position=position_jitter(5.5), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of interlinks")

p_indiv_box3_heter_prob <- ggplot(total_heter_prob_individuals_PROBABILITY_deco, 
                                      aes(x=number_individuals, y=heter_prob, group=number_individuals))+
  # geom_jitter(aes(color=as.factor(number_individuals)),shape=16, position=position_jitter(3.8), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="Prob. of receiving information\nfrom heterosp.",x="Number of individual plants", title = "Prob. of receiving information\nfrom heterosp.")

p_polli_box3_heter_prob <- ggplot(total_heter_prob_pollinators_PROBABILITY_deco, 
                                      aes(x=number_pollinators, y=heter_prob, group=number_pollinators))+
  # geom_jitter(aes(color=as.factor(number_pollinators)),shape=16, position=position_jitter(1.9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of pollinator species")

patchwork_heter_prob <- p_indiv_box3_heter_prob + p_polli_box3_heter_prob + p_intra_box3_heter_prob + p_inter_box3_heter_prob + 
  plot_layout(ncol = 4) + plot_layout(tag_level = "new")
# plot_annotation(
#   title = 'In-Strength',
#   subtitle = 'Parameters for central box plots: 125 intralinks, 30 interlinks, 50 plant individuals, 25 pollinator sps.'
# )

#----------------------------------

p_intra_box3_consp_prob_UNCOUPLED <- ggplot(total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco,
                               aes(x=number_intra, y=consp_prob_UNCOUPLED, group=number_intra))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of intralinks")

p_inter_box3_consp_prob_UNCOUPLED <- ggplot(total_consp_prob_UNCOUPLED_interlinks_PROBABILITY_deco,
                               aes(x=number_inter, y=consp_prob_UNCOUPLED, group=number_inter))+
  # geom_jitter(aes(color=as.factor(number_inter)),shape=16, position=position_jitter(5.5), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of interlinks")

p_indiv_box3_consp_prob_UNCOUPLED <- ggplot(total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco, 
                               aes(x=number_individuals, y=consp_prob_UNCOUPLED, group=number_individuals))+
  # geom_jitter(aes(color=as.factor(number_individuals)),shape=16, position=position_jitter(3.8), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="Prob. of receiving information\nfrom consp. ",x="Number of individual plants", title = "Prob. of receiving information\nfrom consp. (uncoupled)")

p_polli_box3_consp_prob_UNCOUPLED <- ggplot(total_consp_prob_UNCOUPLED_pollinators_PROBABILITY_deco, 
                               aes(x=number_pollinators, y=consp_prob_UNCOUPLED, group=number_pollinators))+
  # geom_jitter(aes(color=as.factor(number_pollinators)),shape=16, position=position_jitter(1.9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of pollinator species")

patchwork_consp_prob_UNCOUPLED <- p_indiv_box3_consp_prob_UNCOUPLED + p_polli_box3_consp_prob_UNCOUPLED + p_intra_box3_consp_prob_UNCOUPLED + p_inter_box3_consp_prob_UNCOUPLED + 
  plot_layout(ncol = 4) + plot_layout(tag_level = "new")
  # plot_annotation(
  #   title = 'PageRank',
  #   subtitle = 'Parameters for central box plots: 125 intralinks, 30 interlinks, 50 plant individuals, 25 pollinator sps.'
  # )


p_intra_box3_consp_prob_UNCOUPLED_log <- ggplot(total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco,
                                            aes(x=number_intra, y=consp_prob_UNCOUPLED, group=number_intra))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of intralinks")

p_inter_box3_consp_prob_UNCOUPLED_log <- ggplot(total_consp_prob_UNCOUPLED_interlinks_PROBABILITY_deco,
                                            aes(x=number_inter, y=consp_prob_UNCOUPLED, group=number_inter))+
  # geom_jitter(aes(color=as.factor(number_inter)),shape=16, position=position_jitter(5.5), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of interlinks")

p_indiv_box3_consp_prob_UNCOUPLED_log <- ggplot(total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco, 
                                            aes(x=number_individuals, y=consp_prob_UNCOUPLED, group=number_individuals))+
  # geom_jitter(aes(color=as.factor(number_individuals)),shape=16, position=position_jitter(3.8), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="Prob. of receiving information\nfrom consp.",x="Number of individual plants", title = "Prob. of receiving information\nfrom consp. (uncoupled)")

p_polli_box3_consp_prob_UNCOUPLED_log <- ggplot(total_consp_prob_UNCOUPLED_pollinators_PROBABILITY_deco, 
                                            aes(x=number_pollinators, y=consp_prob_UNCOUPLED, group=number_pollinators))+
  # geom_jitter(aes(color=as.factor(number_pollinators)),shape=16, position=position_jitter(1.9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun=mean, geom="point", shape=23, size=3,fill="red")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of pollinator species")

patchwork_consp_prob_UNCOUPLED_log <- p_indiv_box3_consp_prob_UNCOUPLED_log + p_polli_box3_consp_prob_UNCOUPLED_log + p_intra_box3_consp_prob_UNCOUPLED_log + p_inter_box3_consp_prob_UNCOUPLED_log + 
  plot_layout(ncol = 4) + plot_layout(tag_level = "new")

#----------------------------------

png("New_Figures/simulations_boxplot4_log.png",
    width = 11.69*0.9, # The width of the plot in inches
    height = 11.69, units = "in", res=300*2)
patchwork_consp_prob_UNCOUPLED_log/patchwork_heter_prob_log/patchwork_homo/patchwork_hete
  # plot_annotation(tag_levels = c("I", "a"))
dev.off()

png("New_Figures/simulations_boxplot4.png",
    width = 11.69*0.9, # The width of the plot in inches
    height = 11.69, units = "in", res=300*2)
patchwork_consp_prob_UNCOUPLED/patchwork_heter_prob/patchwork_homo/patchwork_hete
# plot_annotation(tag_levels = c("I", "a"))
dev.off()

