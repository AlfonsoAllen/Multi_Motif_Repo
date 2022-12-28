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

# Load centrality data -------
raw_intralinks_CENTRALITY <- read_csv("Processed_data/Data_simulation/data_changing_intralinks_CENTRALITY_V2.csv")
raw_interlinks_CENTRALITY <- read_csv("Processed_data/Data_simulation/data_changing_interlinks_CENTRALITY_V2.csv")
raw_individuals_CENTRALITY <- read_csv("Processed_data/Data_simulation/data_changing_individuals_CENTRALITY_V2.csv")
raw_pollinators_CENTRALITY <- read_csv("Processed_data/Data_simulation/data_changing_pollinators_CENTRALITY_V2.csv")

total_InStr_intralinks_CENTRALITY <- raw_intralinks_CENTRALITY[grep("ind",
                                                                    raw_intralinks_CENTRALITY$species,
                                                                    ignore.case = T),] %>%
  separate(species,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,StrengthIn)

total_InStr_interlinks_CENTRALITY <- raw_interlinks_CENTRALITY[grep("ind",
                                                                    raw_interlinks_CENTRALITY$species,
                                                                    ignore.case = T),] %>%
  separate(species,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,StrengthIn)

total_InStr_individuals_CENTRALITY <- raw_individuals_CENTRALITY[grep("ind",
                                                                      raw_individuals_CENTRALITY$species,
                                                                      ignore.case = T),] %>%
  separate(species,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,StrengthIn)


total_InStr_pollinators_CENTRALITY <- raw_pollinators_CENTRALITY[grep("ind",
                                                                      raw_pollinators_CENTRALITY$species,
                                                                      ignore.case = T),] %>%
  separate(species,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,StrengthIn)


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


total_ACRatio_intralinks_CENTRALITY <- raw_intralinks_CENTRALITY[grep("ind",
                                                                      raw_intralinks_CENTRALITY$species,
                                                                      ignore.case = T),] %>%
  separate(species,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,Ratio)

total_ACRatio_interlinks_CENTRALITY <- raw_interlinks_CENTRALITY[grep("ind",
                                                                      raw_interlinks_CENTRALITY$species,
                                                                      ignore.case = T),] %>%
  separate(species,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,Ratio)

total_ACRatio_individuals_CENTRALITY <- raw_individuals_CENTRALITY[grep("ind",
                                                                        raw_individuals_CENTRALITY$species,
                                                                        ignore.case = T),] %>%
  separate(species,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,Ratio)


total_ACRatio_pollinators_CENTRALITY <- raw_pollinators_CENTRALITY[grep("ind",
                                                                        raw_pollinators_CENTRALITY$species,
                                                                        ignore.case = T),] %>%
  separate(species,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,Ratio)


total_InStr_intralinks_CENTRALITY$Plot <- paste0(total_InStr_intralinks_CENTRALITY$Plot,
                                                 25,"_pollinators")
total_InStr_interlinks_CENTRALITY$Plot <- paste0(total_InStr_interlinks_CENTRALITY$Plot,
                                                 25,"_pollinators")
total_InStr_individuals_CENTRALITY$Plot <- paste0(total_InStr_individuals_CENTRALITY$Plot,
                                                  25,"_pollinators")

total_PR_intralinks_CENTRALITY$Plot <- paste0(total_PR_intralinks_CENTRALITY$Plot,
                                              25,"_pollinators")
total_PR_interlinks_CENTRALITY$Plot <- paste0(total_PR_interlinks_CENTRALITY$Plot,
                                              25,"_pollinators")
total_PR_individuals_CENTRALITY$Plot <- paste0(total_PR_individuals_CENTRALITY$Plot,
                                               25,"_pollinators")

total_ACRatio_intralinks_CENTRALITY$Plot <- paste0(total_ACRatio_intralinks_CENTRALITY$Plot,
                                                   25,"_pollinators")
total_ACRatio_interlinks_CENTRALITY$Plot <- paste0(total_ACRatio_interlinks_CENTRALITY$Plot,
                                                   25,"_pollinators")
total_ACRatio_individuals_CENTRALITY$Plot <- paste0(total_ACRatio_individuals_CENTRALITY$Plot,
                                                    25,"_pollinators")

total_InStr_intralinks_CENTRALITY_deco <- decode_plot_label(total_InStr_intralinks_CENTRALITY)
total_InStr_interlinks_CENTRALITY_deco <- decode_plot_label(total_InStr_interlinks_CENTRALITY)
total_InStr_individuals_CENTRALITY_deco <- decode_plot_label(total_InStr_individuals_CENTRALITY)
total_InStr_pollinators_CENTRALITY_deco <- decode_plot_label(total_InStr_pollinators_CENTRALITY)

total_PR_intralinks_CENTRALITY_deco <- decode_plot_label(total_PR_intralinks_CENTRALITY)
total_PR_interlinks_CENTRALITY_deco <- decode_plot_label(total_PR_interlinks_CENTRALITY)
total_PR_individuals_CENTRALITY_deco <- decode_plot_label(total_PR_individuals_CENTRALITY)
total_PR_pollinators_CENTRALITY_deco <- decode_plot_label(total_PR_pollinators_CENTRALITY)

total_ACRatio_intralinks_CENTRALITY_deco <- decode_plot_label(total_ACRatio_intralinks_CENTRALITY)
total_ACRatio_interlinks_CENTRALITY_deco <- decode_plot_label(total_ACRatio_interlinks_CENTRALITY)
total_ACRatio_individuals_CENTRALITY_deco <- decode_plot_label(total_ACRatio_individuals_CENTRALITY)
total_ACRatio_pollinators_CENTRALITY_deco <- decode_plot_label(total_ACRatio_pollinators_CENTRALITY)

total_InStr_intralinks_CENTRALITY_deco$Plant[total_InStr_intralinks_CENTRALITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_InStr_intralinks_CENTRALITY_deco$Plant[total_InStr_intralinks_CENTRALITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_InStr_intralinks_CENTRALITY_deco$Plant[total_InStr_intralinks_CENTRALITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_InStr_intralinks_CENTRALITY_deco$Plant[total_InStr_intralinks_CENTRALITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_InStr_intralinks_CENTRALITY_deco$Plant[total_InStr_intralinks_CENTRALITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_InStr_interlinks_CENTRALITY_deco$Plant[total_InStr_interlinks_CENTRALITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_InStr_interlinks_CENTRALITY_deco$Plant[total_InStr_interlinks_CENTRALITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_InStr_interlinks_CENTRALITY_deco$Plant[total_InStr_interlinks_CENTRALITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_InStr_interlinks_CENTRALITY_deco$Plant[total_InStr_interlinks_CENTRALITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_InStr_interlinks_CENTRALITY_deco$Plant[total_InStr_interlinks_CENTRALITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_InStr_individuals_CENTRALITY_deco$Plant[total_InStr_individuals_CENTRALITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_InStr_individuals_CENTRALITY_deco$Plant[total_InStr_individuals_CENTRALITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_InStr_individuals_CENTRALITY_deco$Plant[total_InStr_individuals_CENTRALITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_InStr_individuals_CENTRALITY_deco$Plant[total_InStr_individuals_CENTRALITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_InStr_individuals_CENTRALITY_deco$Plant[total_InStr_individuals_CENTRALITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_InStr_pollinators_CENTRALITY_deco$Plant[total_InStr_pollinators_CENTRALITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_InStr_pollinators_CENTRALITY_deco$Plant[total_InStr_pollinators_CENTRALITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_InStr_pollinators_CENTRALITY_deco$Plant[total_InStr_pollinators_CENTRALITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_InStr_pollinators_CENTRALITY_deco$Plant[total_InStr_pollinators_CENTRALITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_InStr_pollinators_CENTRALITY_deco$Plant[total_InStr_pollinators_CENTRALITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_PR_intralinks_CENTRALITY_deco$Plant[total_PR_intralinks_CENTRALITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_PR_intralinks_CENTRALITY_deco$Plant[total_PR_intralinks_CENTRALITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_PR_intralinks_CENTRALITY_deco$Plant[total_PR_intralinks_CENTRALITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_PR_intralinks_CENTRALITY_deco$Plant[total_PR_intralinks_CENTRALITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_PR_intralinks_CENTRALITY_deco$Plant[total_PR_intralinks_CENTRALITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_PR_interlinks_CENTRALITY_deco$Plant[total_PR_interlinks_CENTRALITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_PR_interlinks_CENTRALITY_deco$Plant[total_PR_interlinks_CENTRALITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_PR_interlinks_CENTRALITY_deco$Plant[total_PR_interlinks_CENTRALITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_PR_interlinks_CENTRALITY_deco$Plant[total_PR_interlinks_CENTRALITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_PR_interlinks_CENTRALITY_deco$Plant[total_PR_interlinks_CENTRALITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_PR_individuals_CENTRALITY_deco$Plant[total_PR_individuals_CENTRALITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_PR_individuals_CENTRALITY_deco$Plant[total_PR_individuals_CENTRALITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_PR_individuals_CENTRALITY_deco$Plant[total_PR_individuals_CENTRALITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_PR_individuals_CENTRALITY_deco$Plant[total_PR_individuals_CENTRALITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_PR_individuals_CENTRALITY_deco$Plant[total_PR_individuals_CENTRALITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_PR_pollinators_CENTRALITY_deco$Plant[total_PR_pollinators_CENTRALITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_PR_pollinators_CENTRALITY_deco$Plant[total_PR_pollinators_CENTRALITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_PR_pollinators_CENTRALITY_deco$Plant[total_PR_pollinators_CENTRALITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_PR_pollinators_CENTRALITY_deco$Plant[total_PR_pollinators_CENTRALITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_PR_pollinators_CENTRALITY_deco$Plant[total_PR_pollinators_CENTRALITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_ACRatio_intralinks_CENTRALITY_deco$Plant[total_ACRatio_intralinks_CENTRALITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_ACRatio_intralinks_CENTRALITY_deco$Plant[total_ACRatio_intralinks_CENTRALITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_ACRatio_intralinks_CENTRALITY_deco$Plant[total_ACRatio_intralinks_CENTRALITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_ACRatio_intralinks_CENTRALITY_deco$Plant[total_ACRatio_intralinks_CENTRALITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_ACRatio_intralinks_CENTRALITY_deco$Plant[total_ACRatio_intralinks_CENTRALITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_ACRatio_interlinks_CENTRALITY_deco$Plant[total_ACRatio_interlinks_CENTRALITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_ACRatio_interlinks_CENTRALITY_deco$Plant[total_ACRatio_interlinks_CENTRALITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_ACRatio_interlinks_CENTRALITY_deco$Plant[total_ACRatio_interlinks_CENTRALITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_ACRatio_interlinks_CENTRALITY_deco$Plant[total_ACRatio_interlinks_CENTRALITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_ACRatio_interlinks_CENTRALITY_deco$Plant[total_ACRatio_interlinks_CENTRALITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_ACRatio_individuals_CENTRALITY_deco$Plant[total_ACRatio_individuals_CENTRALITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_ACRatio_individuals_CENTRALITY_deco$Plant[total_ACRatio_individuals_CENTRALITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_ACRatio_individuals_CENTRALITY_deco$Plant[total_ACRatio_individuals_CENTRALITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_ACRatio_individuals_CENTRALITY_deco$Plant[total_ACRatio_individuals_CENTRALITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_ACRatio_individuals_CENTRALITY_deco$Plant[total_ACRatio_individuals_CENTRALITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_ACRatio_pollinators_CENTRALITY_deco$Plant[total_ACRatio_pollinators_CENTRALITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_ACRatio_pollinators_CENTRALITY_deco$Plant[total_ACRatio_pollinators_CENTRALITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_ACRatio_pollinators_CENTRALITY_deco$Plant[total_ACRatio_pollinators_CENTRALITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_ACRatio_pollinators_CENTRALITY_deco$Plant[total_ACRatio_pollinators_CENTRALITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_ACRatio_pollinators_CENTRALITY_deco$Plant[total_ACRatio_pollinators_CENTRALITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

#----------------------------------

p_intra_box3_homo <- ggplot(total_homo_intralinks_MOTIFS_deco,
                       aes(x=number_intra, y=homo_motif, group=number_intra))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of intralinks")

p_inter_box3_homo <- ggplot(total_homo_interlinks_MOTIFS_deco,
                       aes(x=number_inter, y=homo_motif, group=number_inter))+
  # geom_jitter(aes(color=as.factor(number_inter)),shape=16, position=position_jitter(5.5), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of interlinks")

p_indiv_box3_homo <- ggplot(total_homo_individuals_MOTIFS_deco, 
                       aes(x=number_individuals, y=homo_motif, group=number_individuals))+
  # geom_jitter(aes(color=as.factor(number_individuals)),shape=16, position=position_jitter(3.8), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="Number of HoM",x="Number of individual plants", title = "Homospecific motifs (HoM)")

p_polli_box3_homo <- ggplot(total_homo_pollinators_MOTIFS_deco, 
                       aes(x=number_pollinators, y=homo_motif, group=number_pollinators))+
  # geom_jitter(aes(color=as.factor(number_pollinators)),shape=16, position=position_jitter(1.9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
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
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of intralinks")

p_inter_box3_hete <- ggplot(total_hete_interlinks_MOTIFS_deco,
                            aes(x=number_inter, y=hete_motif, group=number_inter))+
  # geom_jitter(aes(color=as.factor(number_inter)),shape=16, position=position_jitter(5.5), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of interlinks")

p_indiv_box3_hete <- ggplot(total_hete_individuals_MOTIFS_deco, 
                            aes(x=number_individuals, y=hete_motif, group=number_individuals))+
  # geom_jitter(aes(color=as.factor(number_individuals)),shape=16, position=position_jitter(3.8), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="Number of HeM",x="Number of individual plants", title = "Heterospecific motifs (HeM)")

p_polli_box3_hete <- ggplot(total_hete_pollinators_MOTIFS_deco, 
                            aes(x=number_pollinators, y=hete_motif, group=number_pollinators))+
  # geom_jitter(aes(color=as.factor(number_pollinators)),shape=16, position=position_jitter(1.9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
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

p_intra_box3_InStr <- ggplot(total_InStr_intralinks_CENTRALITY_deco,
                             aes(x=number_intra, y=StrengthIn, group=number_intra))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of intralinks")

p_inter_box3_InStr <- ggplot(total_InStr_interlinks_CENTRALITY_deco,
                             aes(x=number_inter, y=StrengthIn, group=number_inter))+
  # geom_jitter(aes(color=as.factor(number_inter)),shape=16, position=position_jitter(5.5), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of interlinks")

p_indiv_box3_InStr <- ggplot(total_InStr_individuals_CENTRALITY_deco, 
                             aes(x=number_individuals, y=StrengthIn, group=number_individuals))+
  # geom_jitter(aes(color=as.factor(number_individuals)),shape=16, position=position_jitter(3.8), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="In-strength",x="Number of individual plants", title = "In-strength")

p_polli_box3_InStr <- ggplot(total_InStr_pollinators_CENTRALITY_deco, 
                             aes(x=number_pollinators, y=StrengthIn, group=number_pollinators))+
  # geom_jitter(aes(color=as.factor(number_pollinators)),shape=16, position=position_jitter(1.9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of pollinator species")

patchwork_InStr <- p_indiv_box3_InStr + p_polli_box3_InStr + p_intra_box3_InStr + p_inter_box3_InStr + 
  plot_layout(ncol = 4) + plot_layout(tag_level = "new")
  # plot_annotation(
  #   title = 'In-Strength',
  #   subtitle = 'Parameters for central box plots: 125 intralinks, 30 interlinks, 50 plant individuals, 25 pollinator sps.'
  # )

#----------------------------------

p_intra_box3_PR <- ggplot(total_PR_intralinks_CENTRALITY_deco,
                               aes(x=number_intra, y=Real_PR_Multi, group=number_intra))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of intralinks")

p_inter_box3_PR <- ggplot(total_PR_interlinks_CENTRALITY_deco,
                               aes(x=number_inter, y=Real_PR_Multi, group=number_inter))+
  # geom_jitter(aes(color=as.factor(number_inter)),shape=16, position=position_jitter(5.5), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of interlinks")

p_indiv_box3_PR <- ggplot(total_PR_individuals_CENTRALITY_deco, 
                               aes(x=number_individuals, y=Real_PR_Multi, group=number_individuals))+
  # geom_jitter(aes(color=as.factor(number_individuals)),shape=16, position=position_jitter(3.8), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="PageRank",x="Number of individual plants", title = "PageRank")

p_polli_box3_PR <- ggplot(total_PR_pollinators_CENTRALITY_deco, 
                               aes(x=number_pollinators, y=Real_PR_Multi, group=number_pollinators))+
  # geom_jitter(aes(color=as.factor(number_pollinators)),shape=16, position=position_jitter(1.9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of pollinator species")

patchwork_PR <- p_indiv_box3_PR + p_polli_box3_PR + p_intra_box3_PR + p_inter_box3_PR + 
  plot_layout(ncol = 4) + plot_layout(tag_level = "new")
  # plot_annotation(
  #   title = 'PageRank',
  #   subtitle = 'Parameters for central box plots: 125 intralinks, 30 interlinks, 50 plant individuals, 25 pollinator sps.'
  # )

#----------------------------------

p_intra_box3_ACRatio <- ggplot(total_ACRatio_intralinks_CENTRALITY_deco,
                             aes(x=number_intra, y=Ratio, group=number_intra))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of intralinks")

p_inter_box3_ACRatio <- ggplot(total_ACRatio_interlinks_CENTRALITY_deco,
                             aes(x=number_inter, y=Ratio, group=number_inter))+
  # geom_jitter(aes(color=as.factor(number_inter)),shape=16, position=position_jitter(5.5), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of interlinks")

p_indiv_box3_ACRatio <- ggplot(total_ACRatio_individuals_CENTRALITY_deco, 
                             aes(x=number_individuals, y=Ratio, group=number_individuals))+
  # geom_jitter(aes(color=as.factor(number_individuals)),shape=16, position=position_jitter(3.8), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="ALC-ratio",x="Number of individual plants", title = "Among-layer centrality")

p_polli_box3_ACRatio <- ggplot(total_ACRatio_pollinators_CENTRALITY_deco, 
                             aes(x=number_pollinators, y=Ratio, group=number_pollinators))+
  # geom_jitter(aes(color=as.factor(number_pollinators)),shape=16, position=position_jitter(1.9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=.7)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y=NULL,x="Number of pollinator species")

patchwork_ACRatio <- p_indiv_box3_ACRatio + p_polli_box3_ACRatio + p_intra_box3_ACRatio + p_inter_box3_ACRatio + 
  plot_layout(ncol = 4) + plot_layout(tag_level = "new")
  # plot_annotation(
  #   title = 'Among-layer centrality ratio',
  #   subtitle = 'Parameters for central box plots: 125 intralinks, 30 interlinks, 50 plant individuals, 25 pollinator sps.'
  # )


png("New_Figures/simulations_boxplot3.png",
    width = 11.69*0.85, # The width of the plot in inches
    height = 11.69, units = "in", res=300*2)
patchwork_homo/patchwork_hete/patchwork_InStr/patchwork_PR/patchwork_ACRatio
  # plot_annotation(tag_levels = c("I", "a"))
dev.off()
