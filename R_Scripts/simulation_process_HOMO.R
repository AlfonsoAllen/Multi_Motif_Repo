
library(tidyverse)
library(RColorBrewer)
library(patchwork)

decode_plot_label <- function(data_frame_info){
  
  data_frame_info$rep <- NA
  data_frame_info$number_intra <- NA
  data_frame_info$number_inter <- NA
  data_frame_info$number_individuals <- NA
  
  for (i in 1:nrow(data_frame_info)) {
    info_bits <- strsplit(data_frame_info$Plot[i],"_")
    info_bits <- gsub("[^0-9.-]","", info_bits[[1]]) %>% as.numeric()
    data_frame_info$rep[i] <- info_bits[1]
    data_frame_info$number_intra[i] <- info_bits[2]
    data_frame_info$number_inter[i] <- info_bits[3]
    data_frame_info$number_individuals[i] <- info_bits[4]
    
  }
  
  return(data_frame_info)
  
}


# Load data
raw_intralinks_MOTIFS <- read_csv("Processed_data/Data_simulation/data_changing_intralinks_MOTIFS.csv")
raw_interlinks_MOTIFS <- read_csv("Processed_data/Data_simulation/data_changing_interlinks_MOTIFS.csv")
raw_individuals_MOTIFS <- read_csv("Processed_data/Data_simulation/data_changing_individuals_MOTIFS.csv")

total_homo_intralinks_MOTIFS <- raw_intralinks_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=homo_motif) %>% rename(homo_motif=n)

total_homo_interlinks_MOTIFS <- raw_interlinks_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=homo_motif) %>% rename(homo_motif=n)

total_homo_individuals_MOTIFS <- raw_individuals_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=homo_motif) %>% rename(homo_motif=n)

total_homo_intralinks_MOTIFS_deco <- decode_plot_label(total_homo_intralinks_MOTIFS)
total_homo_interlinks_MOTIFS_deco <- decode_plot_label(total_homo_interlinks_MOTIFS)
total_homo_individuals_MOTIFS_deco <- decode_plot_label(total_homo_individuals_MOTIFS)

summary_homo_intralinks <- total_homo_intralinks_MOTIFS_deco %>%
  group_by(Plant,number_intra,number_inter,number_individuals) %>% 
  summarise(Mean = mean(homo_motif, na.rm = T), SD = sd(homo_motif, na.rm = T))

summary_homo_intralinks$Plant[summary_homo_intralinks$Plant=="plant_sp_1"] <- "Sp. 1"
summary_homo_intralinks$Plant[summary_homo_intralinks$Plant=="plant_sp_2"] <- "Sp. 2"
summary_homo_intralinks$Plant[summary_homo_intralinks$Plant=="plant_sp_3"] <- "Sp. 3"
summary_homo_intralinks$Plant[summary_homo_intralinks$Plant=="plant_sp_4"] <- "Sp. 4"
summary_homo_intralinks$Plant[summary_homo_intralinks$Plant=="plant_sp_5"] <- "Sp. 5"

summary_homo_interlinks <- total_homo_interlinks_MOTIFS_deco %>%
  group_by(Plant,number_intra,number_inter,number_individuals) %>% 
  summarise(Mean = mean(homo_motif, na.rm = T), SD = sd(homo_motif, na.rm = T))

summary_homo_interlinks$Plant[summary_homo_interlinks$Plant=="plant_sp_1"] <- "Sp. 1"
summary_homo_interlinks$Plant[summary_homo_interlinks$Plant=="plant_sp_2"] <- "Sp. 2"
summary_homo_interlinks$Plant[summary_homo_interlinks$Plant=="plant_sp_3"] <- "Sp. 3"
summary_homo_interlinks$Plant[summary_homo_interlinks$Plant=="plant_sp_4"] <- "Sp. 4"
summary_homo_interlinks$Plant[summary_homo_interlinks$Plant=="plant_sp_5"] <- "Sp. 5"

summary_homo_individuals <- total_homo_individuals_MOTIFS_deco %>%
  group_by(Plant,number_intra,number_inter,number_individuals) %>% 
  summarise(Mean = mean(homo_motif, na.rm = T), SD = sd(homo_motif, na.rm = T))

summary_homo_individuals$Plant[summary_homo_individuals$Plant=="plant_sp_1"] <- "Sp. 1"
summary_homo_individuals$Plant[summary_homo_individuals$Plant=="plant_sp_2"] <- "Sp. 2"
summary_homo_individuals$Plant[summary_homo_individuals$Plant=="plant_sp_3"] <- "Sp. 3"
summary_homo_individuals$Plant[summary_homo_individuals$Plant=="plant_sp_4"] <- "Sp. 4"
summary_homo_individuals$Plant[summary_homo_individuals$Plant=="plant_sp_5"] <- "Sp. 5"

p_intra <- ggplot(summary_homo_intralinks, aes(x=Plant, y=Mean, fill=Plant))+facet_wrap(~number_intra)+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  geom_errorbar(aes(ymin = Mean, ymax = Mean+SD),width = 0.2,
                position = position_dodge(0.9))+
  scale_fill_brewer(palette = "Paired")+theme_bw()+
  labs(title="Dependence of the number of homospecific motifs on the number of intralinks\n(number of interlinks = 75, total number of plant individuals = 85)")
  
p_inter <- ggplot(summary_homo_interlinks, aes(x=Plant, y=Mean, fill=Plant))+facet_wrap(~number_inter)+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  geom_errorbar(aes(ymin = Mean, ymax = Mean+SD),width = 0.2,
                position = position_dodge(0.9))+
  scale_fill_brewer(palette = "Paired")+theme_bw()+
  labs(title="Dependence of the number of homospecific motifs on the number of interlinks\n(number of intralinks = 125, total number of plant individuals = 85)")

p_indiv <- ggplot(summary_homo_individuals, aes(x=Plant, y=Mean, fill=Plant))+facet_wrap(~number_individuals)+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  geom_errorbar(aes(ymin = Mean, ymax = Mean+SD),width = 0.2,
                position = position_dodge(0.9))+
  scale_fill_brewer(palette = "Paired")+theme_bw()+
  labs(title="Dependence of the number of homospecific motifs on the total number of plant individuals\n(number of intralinks = 125, number of interlinks = 75)")




png("New_Figures/simulations_HOMO.png",
    width = 11.69*0.75, # The width of the plot in inches
    height = 11.69*0.75, units = "in", res=300*2)
p_indiv/p_intra/p_inter
dev.off()
