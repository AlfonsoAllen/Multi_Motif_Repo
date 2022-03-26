

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



total_intralinks_CENTRALITY <- raw_intralinks_CENTRALITY[grep("ind",
                                                                    raw_intralinks_CENTRALITY$species,
                                                                    ignore.case = T),] %>%
  separate(species,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,Ratio,StrengthIn,Real_PR_Multi)

total_interlinks_CENTRALITY <- raw_interlinks_CENTRALITY[grep("ind",
                                                                    raw_interlinks_CENTRALITY$species,
                                                                    ignore.case = T),] %>%
  separate(species,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,Ratio,StrengthIn,Real_PR_Multi)

total_individuals_CENTRALITY <- raw_individuals_CENTRALITY[grep("ind",
                                                                      raw_individuals_CENTRALITY$species,
                                                                      ignore.case = T),] %>%
  separate(species,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,Ratio,StrengthIn,Real_PR_Multi)


total_pollinators_CENTRALITY <- raw_pollinators_CENTRALITY[grep("ind",
                                                                      raw_pollinators_CENTRALITY$species,
                                                                      ignore.case = T),] %>%
  separate(species,c("Ind","Plant"), " ") %>%
  select(Plot,Plant,Ratio,StrengthIn,Real_PR_Multi)


total_intralinks_CENTRALITY$Plot <- paste0(total_intralinks_CENTRALITY$Plot,
                                                 10,"_pollinators")
total_interlinks_CENTRALITY$Plot <- paste0(total_interlinks_CENTRALITY$Plot,
                                                 10,"_pollinators")
total_individuals_CENTRALITY$Plot <- paste0(total_individuals_CENTRALITY$Plot,
                                                  10,"_pollinators")

total_intralinks_CENTRALITY_deco <- decode_plot_label(total_intralinks_CENTRALITY)
total_interlinks_CENTRALITY_deco <- decode_plot_label(total_interlinks_CENTRALITY)
total_individuals_CENTRALITY_deco <- decode_plot_label(total_individuals_CENTRALITY)
total_pollinators_CENTRALITY_deco <- decode_plot_label(total_pollinators_CENTRALITY)



total_centrality <- bind_rows(total_intralinks_CENTRALITY_deco,
                         total_interlinks_CENTRALITY_deco,
                         total_individuals_CENTRALITY_deco,
                         total_pollinators_CENTRALITY_deco) %>%
  mutate(log10_PR=log10(Real_PR_Multi),log10_ratio=log10(Ratio))


library(jtools)
library(lme4)
library(ggstance)
library(broom.mixed)
fm_InStr <- lm(StrengthIn ~ number_intra + number_inter + number_individuals + number_pollinators+Plant,
               total_centrality)
fm_Real_PR_Multi <- lm(log10_PR ~ number_intra + number_inter + number_individuals + number_pollinators+Plant,
               total_centrality)
fm_Ratio <- lm(log10_ratio ~ number_intra + number_inter + number_individuals + number_pollinators+Plant,
               total_centrality)

summ(fm_InStr)
summ(fm_Real_PR_Multi)
summ(fm_Ratio)

plot_summs(fm_InStr,fm_Real_PR_Multi,fm_Ratio, plot.distributions = TRUE,
           scale = TRUE,model.names = c("In-strength","log10(PageRank)","log10(ACRatio)")) 
