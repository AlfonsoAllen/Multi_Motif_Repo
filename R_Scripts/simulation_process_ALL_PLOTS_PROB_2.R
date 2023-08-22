library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(scales)

decode_plot_label <- function(data_frame_info){
  
  data_frame_info$rep <- NA
  data_frame_info$number_intra <- NA
  data_frame_info$number_inter <- NA
  data_frame_info$number_individuals <- NA
  data_frame_info$number_increasing_mix <- NA
  
  for (i in 1:nrow(data_frame_info)) {
    info_bits <- strsplit(data_frame_info$Plot[i],"_")
    info_bits <- gsub("[^0-9.-]","", info_bits[[1]]) %>% as.numeric()
    data_frame_info$rep[i] <- info_bits[1]
    data_frame_info$number_intra[i] <- info_bits[2]
    data_frame_info$number_inter[i] <- info_bits[3]
    data_frame_info$number_individuals[i] <- info_bits[4]
    data_frame_info$number_increasing_mix[i] <- info_bits[5]
  }
  
  return(data_frame_info)
  
}


# Load MOTIF data
raw_intralinks_MOTIFS <- read_csv("Processed_data/Data_simulation/data_changing_intralinks_MOTIFS_V2.csv")
raw_decreasing_mix_MOTIFS <- read_csv("Processed_data/Data_simulation/data_changing_decreasing_mix_MOTIFS_V2.csv")
raw_individuals_MOTIFS <- read_csv("Processed_data/Data_simulation/data_changing_individuals_MOTIFS_V2.csv")
raw_increasing_mix_MOTIFS <- read_csv("Processed_data/Data_simulation/data_changing_increasing_mix_MOTIFS_V2.csv")


total_homo_intralinks_MOTIFS <- raw_intralinks_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=homo_motif) %>% rename(homo_motif=n)

total_homo_decreasing_mix_MOTIFS <- raw_decreasing_mix_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=homo_motif) %>% rename(homo_motif=n)

total_homo_individuals_MOTIFS <- raw_individuals_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=homo_motif) %>% rename(homo_motif=n)

total_homo_increasing_mix_MOTIFS <- raw_increasing_mix_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=homo_motif) %>% rename(homo_motif=n)

total_hete_intralinks_MOTIFS <- raw_intralinks_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=hete_motif) %>% rename(hete_motif=n)

total_hete_decreasing_mix_MOTIFS <- raw_decreasing_mix_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=hete_motif) %>% rename(hete_motif=n)

total_hete_individuals_MOTIFS <- raw_individuals_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=hete_motif) %>% rename(hete_motif=n)

total_hete_increasing_mix_MOTIFS <- raw_increasing_mix_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=hete_motif) %>% rename(hete_motif=n)


total_homo_intralinks_MOTIFS$Plot <- paste0(total_homo_intralinks_MOTIFS$Plot,
                                            25,"_increasing_mix")
total_homo_decreasing_mix_MOTIFS$Plot <- paste0(total_homo_decreasing_mix_MOTIFS$Plot,
                                            25,"_increasing_mix")
total_homo_individuals_MOTIFS$Plot <- paste0(total_homo_individuals_MOTIFS$Plot,
                                             25,"_increasing_mix")
total_hete_intralinks_MOTIFS$Plot <- paste0(total_hete_intralinks_MOTIFS$Plot,
                                            25,"_increasing_mix")
total_hete_decreasing_mix_MOTIFS$Plot <- paste0(total_hete_decreasing_mix_MOTIFS$Plot,
                                            25,"_increasing_mix")
total_hete_individuals_MOTIFS$Plot <- paste0(total_hete_individuals_MOTIFS$Plot,
                                             25,"_increasing_mix")

total_homo_intralinks_MOTIFS_deco <- decode_plot_label(total_homo_intralinks_MOTIFS)
total_homo_decreasing_mix_MOTIFS_deco <- decode_plot_label(total_homo_decreasing_mix_MOTIFS)
total_homo_individuals_MOTIFS_deco <- decode_plot_label(total_homo_individuals_MOTIFS)
total_homo_increasing_mix_MOTIFS_deco <- decode_plot_label(total_homo_increasing_mix_MOTIFS)
total_hete_intralinks_MOTIFS_deco <- decode_plot_label(total_hete_intralinks_MOTIFS)
total_hete_decreasing_mix_MOTIFS_deco <- decode_plot_label(total_hete_decreasing_mix_MOTIFS)
total_hete_individuals_MOTIFS_deco <- decode_plot_label(total_hete_individuals_MOTIFS)
total_hete_increasing_mix_MOTIFS_deco <- decode_plot_label(total_hete_increasing_mix_MOTIFS)

total_homo_intralinks_MOTIFS_deco$Plant[total_homo_intralinks_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_homo_intralinks_MOTIFS_deco$Plant[total_homo_intralinks_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_homo_intralinks_MOTIFS_deco$Plant[total_homo_intralinks_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_homo_intralinks_MOTIFS_deco$Plant[total_homo_intralinks_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_homo_intralinks_MOTIFS_deco$Plant[total_homo_intralinks_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_homo_decreasing_mix_MOTIFS_deco$Plant[total_homo_decreasing_mix_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_homo_decreasing_mix_MOTIFS_deco$Plant[total_homo_decreasing_mix_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_homo_decreasing_mix_MOTIFS_deco$Plant[total_homo_decreasing_mix_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_homo_decreasing_mix_MOTIFS_deco$Plant[total_homo_decreasing_mix_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_homo_decreasing_mix_MOTIFS_deco$Plant[total_homo_decreasing_mix_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_homo_individuals_MOTIFS_deco$Plant[total_homo_individuals_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_homo_individuals_MOTIFS_deco$Plant[total_homo_individuals_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_homo_individuals_MOTIFS_deco$Plant[total_homo_individuals_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_homo_individuals_MOTIFS_deco$Plant[total_homo_individuals_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_homo_individuals_MOTIFS_deco$Plant[total_homo_individuals_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_homo_increasing_mix_MOTIFS_deco$Plant[total_homo_increasing_mix_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_homo_increasing_mix_MOTIFS_deco$Plant[total_homo_increasing_mix_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_homo_increasing_mix_MOTIFS_deco$Plant[total_homo_increasing_mix_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_homo_increasing_mix_MOTIFS_deco$Plant[total_homo_increasing_mix_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_homo_increasing_mix_MOTIFS_deco$Plant[total_homo_increasing_mix_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_hete_intralinks_MOTIFS_deco$Plant[total_hete_intralinks_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_hete_intralinks_MOTIFS_deco$Plant[total_hete_intralinks_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_hete_intralinks_MOTIFS_deco$Plant[total_hete_intralinks_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_hete_intralinks_MOTIFS_deco$Plant[total_hete_intralinks_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_hete_intralinks_MOTIFS_deco$Plant[total_hete_intralinks_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_hete_decreasing_mix_MOTIFS_deco$Plant[total_hete_decreasing_mix_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_hete_decreasing_mix_MOTIFS_deco$Plant[total_hete_decreasing_mix_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_hete_decreasing_mix_MOTIFS_deco$Plant[total_hete_decreasing_mix_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_hete_decreasing_mix_MOTIFS_deco$Plant[total_hete_decreasing_mix_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_hete_decreasing_mix_MOTIFS_deco$Plant[total_hete_decreasing_mix_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_hete_individuals_MOTIFS_deco$Plant[total_hete_individuals_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_hete_individuals_MOTIFS_deco$Plant[total_hete_individuals_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_hete_individuals_MOTIFS_deco$Plant[total_hete_individuals_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_hete_individuals_MOTIFS_deco$Plant[total_hete_individuals_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_hete_individuals_MOTIFS_deco$Plant[total_hete_individuals_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_hete_increasing_mix_MOTIFS_deco$Plant[total_hete_increasing_mix_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_hete_increasing_mix_MOTIFS_deco$Plant[total_hete_increasing_mix_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_hete_increasing_mix_MOTIFS_deco$Plant[total_hete_increasing_mix_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_hete_increasing_mix_MOTIFS_deco$Plant[total_hete_increasing_mix_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_hete_increasing_mix_MOTIFS_deco$Plant[total_hete_increasing_mix_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"

# Load PROBABILITY data -------
raw_intralinks_PROBABILITY <- read_csv("Processed_data/Data_simulation/data_changing_intralinks_PROBABILITIES.csv")
raw_decreasing_mix_PROBABILITY <- read_csv("Processed_data/Data_simulation/data_changing_decreasing_mix_PROBABILITIES.csv")
raw_individuals_PROBABILITY <- read_csv("Processed_data/Data_simulation/data_changing_individuals_PROBABILITIES.csv")
raw_increasing_mix_PROBABILITY <- read_csv("Processed_data/Data_simulation/data_changing_increasing_mix_PROBABILITIES.csv")

raw_intralinks_PROBABILITY_UNCOUPLED <- read_csv("Processed_data/Data_simulation/data_changing_intralinks_PROBABILITIES_UNCOUPLED.csv")
raw_decreasing_mix_PROBABILITY_UNCOUPLED <- read_csv("Processed_data/Data_simulation/data_changing_decreasing_mix_PROBABILITIES_UNCOUPLED.csv")
raw_individuals_PROBABILITY_UNCOUPLED <- read_csv("Processed_data/Data_simulation/data_changing_individuals_PROBABILITIES_UNCOUPLED.csv")
raw_increasing_mix_PROBABILITY_UNCOUPLED <- read_csv("Processed_data/Data_simulation/data_changing_increasing_mix_PROBABILITIES_UNCOUPLED.csv")


total_heter_prob_intralinks_PROBABILITY <- raw_intralinks_PROBABILITY[grep("ind",
                                                                    raw_intralinks_PROBABILITY$name,
                                                                    ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,heter_prob)

total_heter_prob_decreasing_mix_PROBABILITY <- raw_decreasing_mix_PROBABILITY[grep("ind",
                                                                    raw_decreasing_mix_PROBABILITY$name,
                                                                    ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,heter_prob)

total_heter_prob_individuals_PROBABILITY <- raw_individuals_PROBABILITY[grep("ind",
                                                                      raw_individuals_PROBABILITY$name,
                                                                      ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,heter_prob)


total_heter_prob_increasing_mix_PROBABILITY <- raw_increasing_mix_PROBABILITY[grep("ind",
                                                                      raw_increasing_mix_PROBABILITY$name,
                                                                      ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,heter_prob)


total_consp_prob_UNCOUPLED_intralinks_PROBABILITY <- raw_intralinks_PROBABILITY_UNCOUPLED[grep("ind",
                                                                 raw_intralinks_PROBABILITY$name,
                                                                 ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,consp_prob_UNCOUPLED)

total_consp_prob_UNCOUPLED_decreasing_mix_PROBABILITY <- raw_decreasing_mix_PROBABILITY_UNCOUPLED[grep("ind",
                                                                 raw_decreasing_mix_PROBABILITY$name,
                                                                 ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,consp_prob_UNCOUPLED)

total_consp_prob_UNCOUPLED_individuals_PROBABILITY <- raw_individuals_PROBABILITY_UNCOUPLED[grep("ind",
                                                                   raw_individuals_PROBABILITY$name,
                                                                   ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,consp_prob_UNCOUPLED)


total_consp_prob_UNCOUPLED_increasing_mix_PROBABILITY <- raw_increasing_mix_PROBABILITY_UNCOUPLED[grep("ind",
                                                                   raw_increasing_mix_PROBABILITY$name,
                                                                   ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,consp_prob_UNCOUPLED)


total_heter_prob_intralinks_PROBABILITY$Plot <- paste0(total_heter_prob_intralinks_PROBABILITY$Plot,
                                                 25,"_increasing_mix")
total_heter_prob_decreasing_mix_PROBABILITY$Plot <- paste0(total_heter_prob_decreasing_mix_PROBABILITY$Plot,
                                                 25,"_increasing_mix")
total_heter_prob_individuals_PROBABILITY$Plot <- paste0(total_heter_prob_individuals_PROBABILITY$Plot,
                                                  25,"_increasing_mix")

total_consp_prob_UNCOUPLED_intralinks_PROBABILITY$Plot <- paste0(total_consp_prob_UNCOUPLED_intralinks_PROBABILITY$Plot,
                                              25,"_increasing_mix")
total_consp_prob_UNCOUPLED_decreasing_mix_PROBABILITY$Plot <- paste0(total_consp_prob_UNCOUPLED_decreasing_mix_PROBABILITY$Plot,
                                              25,"_increasing_mix")
total_consp_prob_UNCOUPLED_individuals_PROBABILITY$Plot <- paste0(total_consp_prob_UNCOUPLED_individuals_PROBABILITY$Plot,
                                               25,"_increasing_mix")


total_heter_prob_intralinks_PROBABILITY_deco <- decode_plot_label(total_heter_prob_intralinks_PROBABILITY)
total_heter_prob_decreasing_mix_PROBABILITY_deco <- decode_plot_label(total_heter_prob_decreasing_mix_PROBABILITY)
total_heter_prob_individuals_PROBABILITY_deco <- decode_plot_label(total_heter_prob_individuals_PROBABILITY)
total_heter_prob_increasing_mix_PROBABILITY_deco <- decode_plot_label(total_heter_prob_increasing_mix_PROBABILITY)

total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco <- decode_plot_label(total_consp_prob_UNCOUPLED_intralinks_PROBABILITY)
total_consp_prob_UNCOUPLED_decreasing_mix_PROBABILITY_deco <- decode_plot_label(total_consp_prob_UNCOUPLED_decreasing_mix_PROBABILITY)
total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco <- decode_plot_label(total_consp_prob_UNCOUPLED_individuals_PROBABILITY)
total_consp_prob_UNCOUPLED_increasing_mix_PROBABILITY_deco <- decode_plot_label(total_consp_prob_UNCOUPLED_increasing_mix_PROBABILITY)

total_heter_prob_intralinks_PROBABILITY_deco$Plant[total_heter_prob_intralinks_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_heter_prob_intralinks_PROBABILITY_deco$Plant[total_heter_prob_intralinks_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_heter_prob_intralinks_PROBABILITY_deco$Plant[total_heter_prob_intralinks_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_heter_prob_intralinks_PROBABILITY_deco$Plant[total_heter_prob_intralinks_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_heter_prob_intralinks_PROBABILITY_deco$Plant[total_heter_prob_intralinks_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_heter_prob_decreasing_mix_PROBABILITY_deco$Plant[total_heter_prob_decreasing_mix_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_heter_prob_decreasing_mix_PROBABILITY_deco$Plant[total_heter_prob_decreasing_mix_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_heter_prob_decreasing_mix_PROBABILITY_deco$Plant[total_heter_prob_decreasing_mix_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_heter_prob_decreasing_mix_PROBABILITY_deco$Plant[total_heter_prob_decreasing_mix_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_heter_prob_decreasing_mix_PROBABILITY_deco$Plant[total_heter_prob_decreasing_mix_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_heter_prob_individuals_PROBABILITY_deco$Plant[total_heter_prob_individuals_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_heter_prob_individuals_PROBABILITY_deco$Plant[total_heter_prob_individuals_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_heter_prob_individuals_PROBABILITY_deco$Plant[total_heter_prob_individuals_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_heter_prob_individuals_PROBABILITY_deco$Plant[total_heter_prob_individuals_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_heter_prob_individuals_PROBABILITY_deco$Plant[total_heter_prob_individuals_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_heter_prob_increasing_mix_PROBABILITY_deco$Plant[total_heter_prob_increasing_mix_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_heter_prob_increasing_mix_PROBABILITY_deco$Plant[total_heter_prob_increasing_mix_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_heter_prob_increasing_mix_PROBABILITY_deco$Plant[total_heter_prob_increasing_mix_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_heter_prob_increasing_mix_PROBABILITY_deco$Plant[total_heter_prob_increasing_mix_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_heter_prob_increasing_mix_PROBABILITY_deco$Plant[total_heter_prob_increasing_mix_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_consp_prob_UNCOUPLED_decreasing_mix_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_decreasing_mix_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_consp_prob_UNCOUPLED_decreasing_mix_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_decreasing_mix_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_consp_prob_UNCOUPLED_decreasing_mix_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_decreasing_mix_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_consp_prob_UNCOUPLED_decreasing_mix_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_decreasing_mix_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_consp_prob_UNCOUPLED_decreasing_mix_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_decreasing_mix_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_consp_prob_UNCOUPLED_increasing_mix_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_increasing_mix_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_consp_prob_UNCOUPLED_increasing_mix_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_increasing_mix_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_consp_prob_UNCOUPLED_increasing_mix_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_increasing_mix_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_consp_prob_UNCOUPLED_increasing_mix_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_increasing_mix_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_consp_prob_UNCOUPLED_increasing_mix_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_increasing_mix_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

#----------------------------------

total_homo_MOTIFS_deco <- bind_rows(total_homo_intralinks_MOTIFS_deco %>% filter(number_intra == 100),
                                    total_homo_individuals_MOTIFS_deco %>% filter(number_individuals == 50),
                                    total_homo_decreasing_mix_MOTIFS_deco,
                                    total_homo_increasing_mix_MOTIFS_deco
)

plot_homo <- ggplot(total_homo_MOTIFS_deco,
                       aes(x=number_intra, y=homo_motif, group=number_intra))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun =mean, geom="point", shape=23, size=3,fill="red")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1,1e2))+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="Number of homospecific subgraphs",x="Number of intralinks")+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))


total_hete_MOTIFS_deco <- bind_rows(total_hete_intralinks_MOTIFS_deco %>% filter(number_intra == 100),
                                    total_hete_individuals_MOTIFS_deco %>% filter(number_individuals == 50),
                                    total_hete_decreasing_mix_MOTIFS_deco,
                                    total_hete_increasing_mix_MOTIFS_deco
)

plot_hete <- ggplot(total_hete_MOTIFS_deco,
                    aes(x=number_intra, y=hete_motif, group=number_intra))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun =mean, geom="point", shape=23, size=3,fill="red")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1,1e2))+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="Number of heterospecific subgraphs",x="Number of intralinks")+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))


total_heter_prob_PROBABILITY_deco <- bind_rows(total_heter_prob_intralinks_PROBABILITY_deco %>% filter(number_intra == 100),
                                    total_heter_prob_individuals_PROBABILITY_deco %>% filter(number_individuals == 50),
                                    total_heter_prob_decreasing_mix_PROBABILITY_deco,
                                    total_heter_prob_increasing_mix_PROBABILITY_deco
)

plot_heter_prob <- ggplot(total_heter_prob_PROBABILITY_deco,
                    aes(x=number_intra, y=heter_prob, group=number_intra))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun =mean, geom="point", shape=23, size=3,fill="red")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e-5,1e-1))+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="Prob. of receiving pollen\nfrom heterosp.",x="Number of intralinks")+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))


total_consp_prob_UNCOUPLED_PROBABILITY_deco <- bind_rows(total_consp_prob_UNCOUPLED_intralinks_PROBABILITY_deco %>% filter(number_intra == 100),
                                               total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco %>% filter(number_individuals == 50),
                                               total_consp_prob_UNCOUPLED_decreasing_mix_PROBABILITY_deco,
                                               total_consp_prob_UNCOUPLED_increasing_mix_PROBABILITY_deco
)

plot_consp_prob_UNCOUPLED <- ggplot(total_consp_prob_UNCOUPLED_PROBABILITY_deco,
                          aes(x=number_intra, y=consp_prob_UNCOUPLED, group=number_intra))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun =mean, geom="point", shape=23, size=3,fill="red")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e-5,1e-1))+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="Prob. of receiving pollen\nfrom consp. (uncopled layers)",x="Number of intralinks")+
  theme(legend.text = element_text(size=15))+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))



#--------
library(patchwork)
png("New_Figures/simulations_boxplot.png",
    width = 11.69, # The width of the plot in inches
    height = 11.69, units = "in", res=300*2)

(plot_homo & plot_hete)/(plot_consp_prob_UNCOUPLED & plot_heter_prob)


  # plot_annotation(tag_levels = c("I", "a"))
dev.off()


#####################################
# Plot modularity results

modularity_INCREASING_MIX <- read_csv("Processed_data/Data_simulation/data_changing_increasing_mix_MODULARITY.csv")
modularity_DECREASING_MIX <- read_csv("Processed_data/Data_simulation/data_changing_decreasing_mix_MODULARITY.csv")
modularity_INTRALINKS <- read_csv("Processed_data/Data_simulation/data_changing_intralinks_MODULARITY.csv")
modularity_INDIVIDUALS <- read_csv("Processed_data/Data_simulation/data_changing_individuals_MODULARITY.csv")

modularity_INTRALINKS$Plot <- paste0(modularity_INTRALINKS$Plot,
                                            25,"_intralinks")
modularity_DECREASING_MIX$Plot <- paste0(modularity_DECREASING_MIX$Plot,
                                                25,"_decreasing_mix")

modularity_INCREASING_MIX$Plot <- paste0(modularity_INCREASING_MIX$Plot,
                                         25,"_increasing_mix")

modularity_INDIVIDUALS$Plot <- paste0(modularity_INDIVIDUALS_MIX$Plot,
                                             25,"_individuals")

modularity_INTRALINKS_deco <- decode_plot_label(modularity_INTRALINKS) %>% filter(number_intra == 100)
modularity_DECREASING_MIX_deco <- decode_plot_label(modularity_DECREASING_MIX)
modularity_INCREASING_MIX_deco <- decode_plot_label(modularity_INCREASING_MIX)
modularity_INDIVIDUALS_deco <- decode_plot_label(modularity_INDIVIDUALS) %>% filter(number_individuals == 100)


modularity_simulations <- bind_rows(modularity_INCREASING_MIX_deco,
                                    modularity_DECREASING_MIX_deco,
                                    modularity_INTRALINKS_deco,
                                    modularity_INDIVIDUALS_deco)

# Plot observed VS simulated degree-------------
code_length_plot <- ggplot(data = modularity_simulations, aes(x=as.factor(number_intra),y = L))+
  geom_hline(yintercept = 3.709, alpha=0.4, color = "blue",linewidth = 1.5)+
  geom_hline(yintercept = 3.971, alpha=0.4,  color = "red",linewidth = 1.5)+
  geom_hline(yintercept = 3.696, alpha=0.4,  color = "green",linewidth = 1.5)+
  geom_hline(yintercept = 3.259, alpha=0.4,  color = "black",linewidth = 1.5)+
  geom_hline(yintercept = 3.126, alpha=0.4,  color = "orange",linewidth = 1.5)+
  geom_hline(yintercept = 2.997, alpha=0.4,  color = "tomato",linewidth = 1.5)+
  geom_hline(yintercept = 3.680, alpha=0.4,  color = "brown",linewidth = 1.5)+
  geom_hline(yintercept = 3.539, alpha=0.4,  color = "magenta",linewidth = 1.5)+
  geom_hline(yintercept = 3.701, alpha=0.4,  color = "cyan",linewidth = 1.5)+
  geom_boxplot()+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun =mean, geom="point", shape=23, size=3,fill="red")+
  labs(x="Number of intralinks", y = "Code length, L (bits)")+
  theme_bw()

code_length_plot


number_modules_plot <- ggplot(data = modularity_simulations, aes(x=as.factor(number_intra),y = m))+
  geom_hline(yintercept = 10, alpha=0.4, color = "blue",linewidth = 1.5)+
  geom_hline(yintercept = 11, alpha=0.4,  color = "red",linewidth = 1.5)+
  geom_hline(yintercept = 14, alpha=0.4,  color = "green",linewidth = 1.5)+
  geom_hline(yintercept = 7, alpha=0.4,  color = "black",linewidth = 1.5)+
  geom_hline(yintercept = 6, alpha=0.4,  color = "orange",linewidth = 1.5)+
  geom_hline(yintercept = 3, alpha=0.4,  color = "tomato",linewidth = 1.5)+
  geom_hline(yintercept = 14, alpha=0.4, linetype='dashed', color = "brown",linewidth = 1.5)+
  geom_hline(yintercept = 16, alpha=0.4,  color = "magenta",linewidth = 1.5)+
  geom_hline(yintercept = 16, alpha=0.4,  linetype='dashed', color = "cyan",linewidth = 1.5)+
  geom_boxplot()+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun =mean, geom="point", shape=23, size=3,fill="red")+
  labs(x="Number of intralinks", y = "Number of modules, m")+
  theme_bw()

number_modules_plot

library(patchwork)
png("New_Figures/simulations_modularity.png",
    width = 11.69*0.9, # The width of the plot in inches
    height = 11.69*0.4, units = "in", res=300*2)

(code_length_plot & number_modules_plot)
dev.off()
