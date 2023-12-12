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
raw_individuals_MOTIFS <- read_csv("Processed_data/Data_simulation/data_changing_INDIVIDUALS_MOTIFS_V3.csv")


total_homo_individuals_MOTIFS <- raw_individuals_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=homo_motif) %>% rename(homo_motif=n)

total_hete_individuals_MOTIFS <- raw_individuals_MOTIFS %>% group_by(Plot,Plant,Subplot) %>%
  count(wt=hete_motif) %>% rename(hete_motif=n)



total_homo_individuals_MOTIFS$Plot <- paste0(total_homo_individuals_MOTIFS$Plot,
                                             25,"_increasing_mix")
total_hete_individuals_MOTIFS$Plot <- paste0(total_hete_individuals_MOTIFS$Plot,
                                             25,"_increasing_mix")


total_homo_individuals_MOTIFS_deco <- decode_plot_label(total_homo_individuals_MOTIFS)
total_hete_individuals_MOTIFS_deco <- decode_plot_label(total_hete_individuals_MOTIFS)

total_homo_individuals_MOTIFS_deco$Plant[total_homo_individuals_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_homo_individuals_MOTIFS_deco$Plant[total_homo_individuals_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_homo_individuals_MOTIFS_deco$Plant[total_homo_individuals_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_homo_individuals_MOTIFS_deco$Plant[total_homo_individuals_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_homo_individuals_MOTIFS_deco$Plant[total_homo_individuals_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"


total_hete_individuals_MOTIFS_deco$Plant[total_hete_individuals_MOTIFS_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_hete_individuals_MOTIFS_deco$Plant[total_hete_individuals_MOTIFS_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_hete_individuals_MOTIFS_deco$Plant[total_hete_individuals_MOTIFS_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_hete_individuals_MOTIFS_deco$Plant[total_hete_individuals_MOTIFS_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_hete_individuals_MOTIFS_deco$Plant[total_hete_individuals_MOTIFS_deco$Plant=="plant_sp_5"] <- "Sp. 5"

# Load PROBABILITY data -------

raw_individuals_PROBABILITY <- read_csv("Processed_data/Data_simulation/data_changing_INDIVIDUALS_PROBABILITIES_V3.csv")
raw_individuals_PROBABILITY_UNCOUPLED <- read_csv("Processed_data/Data_simulation/data_changing_INDIVIDUALS_PROBABILITIES_UNCOUPLED_V3.csv")


total_heter_prob_individuals_PROBABILITY <- raw_individuals_PROBABILITY[grep("ind",
                                                                      raw_individuals_PROBABILITY$name,
                                                                      ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  dplyr::select(Plot,Plant,heter_prob)



total_consp_prob_UNCOUPLED_individuals_PROBABILITY <- raw_individuals_PROBABILITY_UNCOUPLED[grep("ind",
                                                                   raw_individuals_PROBABILITY$name,
                                                                   ignore.case = T),] %>%
  separate(name,c("Ind","Plant"), " ") %>%
  dplyr::select(Plot,Plant,consp_prob_UNCOUPLED)


total_heter_prob_individuals_PROBABILITY$Plot <- paste0(total_heter_prob_individuals_PROBABILITY$Plot,
                                                  25,"_increasing_mix")

total_consp_prob_UNCOUPLED_individuals_PROBABILITY$Plot <- paste0(total_consp_prob_UNCOUPLED_individuals_PROBABILITY$Plot,
                                               25,"_increasing_mix")

total_heter_prob_individuals_PROBABILITY_deco <- decode_plot_label(total_heter_prob_individuals_PROBABILITY)
total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco <- decode_plot_label(total_consp_prob_UNCOUPLED_individuals_PROBABILITY)

total_heter_prob_individuals_PROBABILITY_deco$Plant[total_heter_prob_individuals_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_heter_prob_individuals_PROBABILITY_deco$Plant[total_heter_prob_individuals_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_heter_prob_individuals_PROBABILITY_deco$Plant[total_heter_prob_individuals_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_heter_prob_individuals_PROBABILITY_deco$Plant[total_heter_prob_individuals_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_heter_prob_individuals_PROBABILITY_deco$Plant[total_heter_prob_individuals_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant=="plant_sp_1"] <- "Sp. 1"
total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant=="plant_sp_2"] <- "Sp. 2"
total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant=="plant_sp_3"] <- "Sp. 3"
total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant=="plant_sp_4"] <- "Sp. 4"
total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant[total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco$Plant=="plant_sp_5"] <- "Sp. 5"

#----------------------------------

total_homo_MOTIFS_deco <- bind_rows(total_homo_individuals_MOTIFS_deco)

plot_homo <- ggplot(total_homo_MOTIFS_deco,
                       aes(x=number_individuals, y=homo_motif, group=number_individuals))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun =mean, geom="point", shape=23, size=3,fill="red")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1,5e1))+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="Number of homospecific subgraphs",x="Number of plant individuals")+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))


total_hete_MOTIFS_deco <- bind_rows(total_hete_individuals_MOTIFS_deco)

plot_hete <- ggplot(total_hete_MOTIFS_deco,
                    aes(x=number_individuals, y=hete_motif, group=number_individuals))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun =mean, geom="point", shape=23, size=3,fill="red")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1,5e1))+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="Number of heterospecific subgraphs",x="Number of plant individuals")+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))


total_heter_prob_PROBABILITY_deco <- bind_rows(
                                    total_heter_prob_individuals_PROBABILITY_deco
)

plot_heter_prob <- ggplot(total_heter_prob_PROBABILITY_deco,
                    aes(x=number_individuals, y=heter_prob, group=number_individuals))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun =mean, geom="point", shape=23, size=3,fill="red")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e-4,1e-1))+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="Prob. of receiving pollen\nfrom heterosp.",x="Number of plant individuals")+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))


total_consp_prob_UNCOUPLED_PROBABILITY_deco <- bind_rows(
                                               total_consp_prob_UNCOUPLED_individuals_PROBABILITY_deco
)

plot_consp_prob_UNCOUPLED <- ggplot(total_consp_prob_UNCOUPLED_PROBABILITY_deco,
                          aes(x=number_individuals, y=consp_prob_UNCOUPLED, group=number_individuals))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun =mean, geom="point", shape=23, size=3,fill="red")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e-4,1e-1))+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="Prob. of receiving pollen\nfrom consp. (uncopled layers)",x="Number of plant individuals")+
  theme(legend.text = element_text(size=15))+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))


ggplot(total_consp_prob_UNCOUPLED_PROBABILITY_deco,
       aes(x=number_individuals, y=consp_prob_UNCOUPLED, group=number_individuals))+
  # geom_jitter(aes(color=as.factor(number_intra)),shape=16, position=position_jitter(9), alpha=0.3)+
  geom_boxplot(alpha = 1.0,lwd=0.7)+
  stat_summary(fun = mean, color = "red", geom = "line", aes(group = 1),lwd=1.5)+
  stat_summary(fun =mean, geom="point", shape=23, size=3,fill="red")+
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x)),
  #               limits = c(1e-4,1e-1))+
  scale_color_brewer(palette = "Set2")+theme_bw()+
  theme(legend.position="none", plot.title = element_text(size = 11))+
  labs(y="Prob. of receiving pollen\nfrom consp. (uncopled layers)",x="Number of plant individuals")+
  theme(legend.text = element_text(size=15))+
  theme(axis.text=element_text(size=15),  axis.title=element_text(size=17,face="bold"))+                                                                # Change font size
  theme(strip.text.x = element_text(size = 18))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE),colour = guide_legend(override.aes = list(size=5)))

#--------
library(patchwork)
png("New_Figures/simulations_INDIVIDUALS_boxplot_V3.png",
    width = 11.69, # The width of the plot in inches
    height = 11.69, units = "in", res=300*2)

(plot_homo & plot_hete)/(plot_consp_prob_UNCOUPLED & plot_heter_prob)


  # plot_annotation(tag_levels = c("I", "a"))
dev.off()



##############################################################################
# Significant pairwise differences in probabilities

library(sandwich)
library(multcomp)

amod <- aov(hete_motif ~ condition, total_hete_individuals_MOTIFS_deco %>%
              mutate(condition=as.factor(number_individuals)))
amod_glht <- glht(amod, mcp(condition="Tukey"), vcov=vcovHC)
summary(amod_glht)
confint(amod_glht)
confint_report <- confint(amod_glht)

# Configurar la relación de aspecto y los márgenes
par(mar = c(5, 5, 2, 2))  # c(bottom, left, top, right)
plot(confint(amod_glht), main="95% family-wise confidence level for number of heterosp. motifs")

results_hete <- data.frame(confint_report[["confint"]])
results_hete$comparisson <- rownames(results_hete)

png("New_Figures/pairwise_heter_motifs_individuals_V3.png",
    width = 11.69*0.55, # The width of the plot in inches
    height = 11.69*0.35, units = "in", res=300*2)

ggplot(results_hete,aes(y = as.factor(comparisson)))+
  geom_point(aes(x=Estimate,size=2))+
  geom_errorbar(aes(xmin=lwr, xmax=upr), width=.2)+
  geom_vline(xintercept = 0,linetype = "dashed")+
  theme_bw()+
  guides(size = "none")+
  labs(title="Pairwise comparisons of the mean number of heterosp. motifs", x=NULL, y = NULL)+
  theme(axis.text=element_text(size=14),  plot.title=element_text(size=19))+
  theme(strip.text = element_text(size=15,face = "italic"))

dev.off()




amod <- aov(consp_prob_UNCOUPLED ~ condition, total_consp_prob_UNCOUPLED_PROBABILITY_deco %>%
              mutate(condition=as.factor(number_individuals)))
amod_glht <- glht(amod, mcp(condition="Tukey"), vcov=vcovHC)
summary(amod_glht)
confint(amod_glht)
confint_report <- confint(amod_glht)

# Configurar la relación de aspecto y los márgenes
par(mar = c(5, 5, 2, 2))  # c(bottom, left, top, right)
plot(confint(amod_glht), main="95% family-wise confidence level for the prob. of receiving pollen from consp. (uncopled layers)")

results_consp <- data.frame(confint_report[["confint"]])
results_consp$comparisson <- rownames(results_consp)

png("New_Figures/pairwise_comp_consp_prob_UNCOUPLED_individuals_V3.png",
    width = 11.69*0.55, # The width of the plot in inches
    height = 11.69*0.35, units = "in", res=300*2)

ggplot(results_consp,aes(y = as.factor(comparisson)))+
  geom_point(aes(x=Estimate,size=2))+
  geom_errorbar(aes(xmin=lwr, xmax=upr), width=.2)+
  geom_vline(xintercept = 0,linetype = "dashed")+
  theme_bw()+
  guides(size = "none")+
  labs(title="Pairwise comparisons of the mean prob. of\nreceiving pollen from consp. (uncopled layers)", x=NULL, y = NULL)+
  theme(axis.text=element_text(size=14),  plot.title=element_text(size=19))+
  theme(strip.text = element_text(size=15,face = "italic"))

dev.off()


amod <- aov(heter_prob ~ condition, total_heter_prob_PROBABILITY_deco %>%
              mutate(condition=as.factor(number_individuals)))
amod_glht <- glht(amod, mcp(condition="Tukey"), vcov=vcovHC)
summary(amod_glht)
confint(amod_glht)
plot(confint(amod_glht), main="95% family-wise confidence level for the prob. of receiving pollen from heterosp.")
confint_report <- confint(amod_glht)


results_heter <- data.frame(confint_report[["confint"]])
results_heter$comparisson <- rownames(results_heter)

png("New_Figures/pairwise_comp_heter_prob_individuals_V3.png",
    width = 11.69*0.55, # The width of the plot in inches
    height = 11.69*0.35, units = "in", res=300*2)

ggplot(results_heter,aes(y = as.factor(comparisson)))+
  geom_point(aes(x=Estimate,size=2))+
  geom_errorbar(aes(xmin=lwr, xmax=upr), width=.2)+
  geom_vline(xintercept = 0,linetype = "dashed")+
  theme_bw()+
  guides(size = "none")+
  labs(title="Pairwise comparisons of the mean prob. of\nreceiving pollen from heterosp.", x=NULL, y = NULL)+
  theme(axis.text=element_text(size=14),  plot.title=element_text(size=19))+
  theme(strip.text = element_text(size=15,face = "italic"))

dev.off()
