
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


total_homo_intralinks_MOTIFS <- raw_intralinks_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=homo_motif) %>% rename(homo_motif=n)

total_homo_interlinks_MOTIFS <- raw_interlinks_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=homo_motif) %>% rename(homo_motif=n)

total_homo_individuals_MOTIFS <- raw_individuals_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=homo_motif) %>% rename(homo_motif=n)

total_homo_pollinators_MOTIFS <- raw_pollinators_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=homo_motif) %>% rename(homo_motif=n)


total_homo_intralinks_MOTIFS$Plot <- paste0(total_homo_intralinks_MOTIFS$Plot,
                                          10,"_pollinators")
total_homo_interlinks_MOTIFS$Plot <- paste0(total_homo_interlinks_MOTIFS$Plot,
                                             10,"_pollinators")
total_homo_individuals_MOTIFS$Plot <- paste0(total_homo_individuals_MOTIFS$Plot,
                                             10,"_pollinators")

total_homo_intralinks_MOTIFS_deco <- decode_plot_label(total_homo_intralinks_MOTIFS)
total_homo_interlinks_MOTIFS_deco <- decode_plot_label(total_homo_interlinks_MOTIFS)
total_homo_individuals_MOTIFS_deco <- decode_plot_label(total_homo_individuals_MOTIFS)
total_homo_pollinators_MOTIFS_deco <- decode_plot_label(total_homo_pollinators_MOTIFS)

summary_homo_intralinks <- total_homo_intralinks_MOTIFS_deco %>%
  group_by(Plant,number_intra,number_inter,number_individuals,number_pollinators) %>% 
  summarise(Mean = mean(homo_motif, na.rm = T), SD = sd(homo_motif, na.rm = T))

summary_homo_intralinks$Plant[summary_homo_intralinks$Plant=="plant_sp_1"] <- "Sp. 1"
summary_homo_intralinks$Plant[summary_homo_intralinks$Plant=="plant_sp_2"] <- "Sp. 2"
summary_homo_intralinks$Plant[summary_homo_intralinks$Plant=="plant_sp_3"] <- "Sp. 3"
summary_homo_intralinks$Plant[summary_homo_intralinks$Plant=="plant_sp_4"] <- "Sp. 4"
summary_homo_intralinks$Plant[summary_homo_intralinks$Plant=="plant_sp_5"] <- "Sp. 5"

summary_homo_interlinks <- total_homo_interlinks_MOTIFS_deco %>%
  group_by(Plant,number_intra,number_inter,number_individuals,number_pollinators) %>% 
  summarise(Mean = mean(homo_motif, na.rm = T), SD = sd(homo_motif, na.rm = T))

summary_homo_interlinks$Plant[summary_homo_interlinks$Plant=="plant_sp_1"] <- "Sp. 1"
summary_homo_interlinks$Plant[summary_homo_interlinks$Plant=="plant_sp_2"] <- "Sp. 2"
summary_homo_interlinks$Plant[summary_homo_interlinks$Plant=="plant_sp_3"] <- "Sp. 3"
summary_homo_interlinks$Plant[summary_homo_interlinks$Plant=="plant_sp_4"] <- "Sp. 4"
summary_homo_interlinks$Plant[summary_homo_interlinks$Plant=="plant_sp_5"] <- "Sp. 5"

summary_homo_individuals <- total_homo_individuals_MOTIFS_deco %>%
  group_by(Plant,number_intra,number_inter,number_individuals,number_pollinators) %>% 
  summarise(Mean = mean(homo_motif, na.rm = T), SD = sd(homo_motif, na.rm = T))

summary_homo_individuals$Plant[summary_homo_individuals$Plant=="plant_sp_1"] <- "Sp. 1"
summary_homo_individuals$Plant[summary_homo_individuals$Plant=="plant_sp_2"] <- "Sp. 2"
summary_homo_individuals$Plant[summary_homo_individuals$Plant=="plant_sp_3"] <- "Sp. 3"
summary_homo_individuals$Plant[summary_homo_individuals$Plant=="plant_sp_4"] <- "Sp. 4"
summary_homo_individuals$Plant[summary_homo_individuals$Plant=="plant_sp_5"] <- "Sp. 5"

summary_homo_pollinators <- total_homo_pollinators_MOTIFS_deco %>%
  group_by(Plant,number_intra,number_inter,number_individuals,number_pollinators) %>% 
  summarise(Mean = mean(homo_motif, na.rm = T), SD = sd(homo_motif, na.rm = T))

summary_homo_pollinators$Plant[summary_homo_pollinators$Plant=="plant_sp_1"] <- "Sp. 1"
summary_homo_pollinators$Plant[summary_homo_pollinators$Plant=="plant_sp_2"] <- "Sp. 2"
summary_homo_pollinators$Plant[summary_homo_pollinators$Plant=="plant_sp_3"] <- "Sp. 3"
summary_homo_pollinators$Plant[summary_homo_pollinators$Plant=="plant_sp_4"] <- "Sp. 4"
summary_homo_pollinators$Plant[summary_homo_pollinators$Plant=="plant_sp_5"] <- "Sp. 5"

p_intra <- ggplot(summary_homo_intralinks, aes(x=Plant, y=Mean, fill=Plant))+
  facet_wrap(~number_intra)+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  geom_errorbar(aes(ymin = Mean, ymax = Mean+SD),width = 0.2,
                position = position_dodge(0.9))+
  scale_fill_brewer(palette = "Paired")+theme_bw()+
  labs(title="Dependence of the number of homospecific motifs on the number of intralinks\n(30 interlinks, 50 plant individuals , 10 pollinator sps.)")
  
p_inter <- ggplot(summary_homo_interlinks, aes(x=Plant, y=Mean, fill=Plant))+facet_wrap(~number_inter)+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  geom_errorbar(aes(ymin = Mean, ymax = Mean+SD),width = 0.2,
                position = position_dodge(0.9))+
  scale_fill_brewer(palette = "Paired")+theme_bw()+
  labs(title="Dependence of the number of homospecific motifs on the number of interlinks\n(125 intralinks, 50 plant individuals , 10 pollinator sps.)")

p_indiv <- ggplot(summary_homo_individuals, aes(x=Plant, y=Mean, fill=Plant))+facet_wrap(~number_individuals)+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  geom_errorbar(aes(ymin = Mean, ymax = Mean+SD),width = 0.2,
                position = position_dodge(0.9))+
  scale_fill_brewer(palette = "Paired")+theme_bw()+
  labs(title="Dependence of the number of homospecific motifs on the total number of plant individuals\n(125 intralinks, 30 interlinks, 10 pollinator sps.)")

p_polli <- ggplot(summary_homo_pollinators, aes(x=Plant, y=Mean, fill=Plant))+
  facet_wrap(~number_pollinators)+
  geom_bar(stat = "identity", color = "black", position = position_dodge())+
  geom_errorbar(aes(ymin = Mean, ymax = Mean+SD),width = 0.2,
                position = position_dodge(0.9))+
  scale_fill_brewer(palette = "Paired")+theme_bw()+
  labs(title="Dependence of the number of homospecific motifs on the total number of pollinator sp.\n(125 intralinks, 30 interlinks, 50 plant individuals)")



png("New_Figures/simulations_HOMO.png",
    width = 11.69*0.75, # The width of the plot in inches
    height = 11.69, units = "in", res=300*2)
p_indiv/p_intra/p_inter/p_polli
dev.off()


#---------

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

p_intra_box <- ggplot(total_homo_intralinks_MOTIFS_deco, aes(x=Plant, y=homo_motif))+
  facet_wrap(~number_intra)+
  geom_jitter(aes(color=Plant),shape=16, position=position_jitter(0.2), alpha=0.3)+
  geom_boxplot(alpha = 0.0)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Paired")+theme_bw()+theme(legend.position="none")+
  labs(y=NULL,title="Dependence of the number of homospecific motifs on the number of intralinks\n(30 interlinks, 50 plant individuals , 10 pollinator sps.)")

p_inter_box <- ggplot(total_homo_interlinks_MOTIFS_deco, aes(x=Plant, y=homo_motif))+
  facet_wrap(~number_inter)+
  geom_jitter(aes(color=Plant),shape=16, position=position_jitter(0.2), alpha=0.3)+
  geom_boxplot(alpha = 0.0)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Paired")+theme_bw()+theme(legend.position="none")+
  labs(y=NULL,title="Dependence of the number of homospecific motifs on the number of interlinks\n(125 intralinks, 50 plant individuals , 10 pollinator sps.)")

p_indiv_box <- ggplot(total_homo_individuals_MOTIFS_deco, aes(x=Plant, y=homo_motif))+
  facet_wrap(~number_individuals)+
  geom_jitter(aes(color=Plant),shape=16, position=position_jitter(0.2), alpha=0.3)+
  geom_boxplot(alpha = 0.0)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Paired")+theme_bw()+theme(legend.position="none")+
  labs(y=NULL,title="Dependence of the number of homospecific motifs on the total number of plant individuals\n(125 intralinks, 30 interlinks, 10 pollinator sps.)")

p_polli_box <- ggplot(total_homo_pollinators_MOTIFS_deco, aes(x=Plant, y=homo_motif))+
  facet_wrap(~number_pollinators)+
  geom_jitter(aes(color=Plant),shape=16, position=position_jitter(0.2), alpha=0.3)+
  geom_boxplot(alpha = 0.0)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Paired")+theme_bw()+theme(legend.position="none")+
  labs(y=NULL,title="Dependence of the number of homospecific motifs on the total number of pollinator sp.\n(125 intralinks, 30 interlinks, 50 plant individuals)")



png("New_Figures/simulations_HOMO_boxplot.png",
    width = 11.69*0.75, # The width of the plot in inches
    height = 11.69, units = "in", res=300*2)
p_indiv_box/p_intra_box/p_inter_box/p_polli_box
dev.off()


#----------------------------------

p_intra_box2 <- ggplot(total_homo_intralinks_MOTIFS_deco,
                       aes(x=number_intra, y=homo_motif, group=number_intra))+
  geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 0.0,lwd=1.3)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+theme(legend.position="none")+
  labs(y=NULL,x="Number of intralinks",
       title="Dependence of the number of homospecific motifs on the number of intralinks\n(30 interlinks, 50 plant individuals , 10 pollinator sps.)")

p_inter_box2 <- ggplot(total_homo_interlinks_MOTIFS_deco,
                       aes(x=number_inter, y=homo_motif, group=number_inter))+
  geom_jitter(aes(color=as.factor(number_inter)),shape=16, position=position_jitter(5.5), alpha=0.3)+
  geom_boxplot(alpha = 0.0,lwd=1.3)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+theme(legend.position="none")+
  labs(y=NULL,x="Number of interlinks",
       title="Dependence of the number of homospecific motifs on the number of interlinks\n(125 intralinks, 50 plant individuals , 10 pollinator sps.)")

p_indiv_box2 <- ggplot(total_homo_individuals_MOTIFS_deco, 
                       aes(x=number_individuals, y=homo_motif, group=number_individuals))+
  geom_jitter(aes(color=as.factor(number_individuals)),shape=16, position=position_jitter(3.8), alpha=0.3)+
  geom_boxplot(alpha = 0.0,lwd=1.3)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+theme(legend.position="none")+
  labs(y=NULL,x="Number of individual plants",
       title="Dependence of the number of homospecific motifs on the total number of plant individuals\n(125 intralinks, 30 interlinks, 10 pollinator sps.)")

p_polli_box2 <- ggplot(total_homo_pollinators_MOTIFS_deco, 
                       aes(x=number_pollinators, y=homo_motif, group=number_pollinators))+
  geom_jitter(aes(color=as.factor(number_pollinators)),shape=16, position=position_jitter(1.9), alpha=0.3)+
  geom_boxplot(alpha = 0.0,lwd=1.3)+
  stat_summary(fun.y = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=3,fill="red")+
  scale_color_brewer(palette = "Set2")+theme_bw()+theme(legend.position="none")+
  labs(y=NULL,x="Number of pollinator species",
       title="Dependence of the number of homospecific motifs on the total number of pollinator sp.\n(125 intralinks, 30 interlinks, 50 plant individuals)")



png("New_Figures/simulations_HOMO_boxplot2.png",
    width = 11.69*0.75, # The width of the plot in inches
    height = 11.69, units = "in", res=300*2)
p_indiv_box2/p_intra_box2/p_inter_box2/p_polli_box2
dev.off()
