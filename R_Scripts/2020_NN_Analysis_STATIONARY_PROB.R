
library(tidyverse)

Prob_results <- read_csv("Processed_data/2020_NN_plant_stationary_prob_results.csv")

Prob_results$consp_prob <- round(Prob_results$consp_prob,10)
Prob_results$heter_prob <- round(Prob_results$heter_prob,10)

plot_labs <-c(
  "Plot 1",
  "Plot 2",
  "Plot 3",
  "Plot 4",
  "Plot 5",
  "Plot 6",
  "Plot 7",
  "Plot 8",
  "Plot 9"
)
names(plot_labs) <- c(
  '1',
  '2',
  '3',
  '4',
  '5',
  "6",
  '7',
  '8',
  "9"
)



library(RColorBrewer)

ggplot(Prob_results,aes(x=layer,y=log10(consp_prob)))+
  #geom_violin()+ 
  
  geom_boxplot()+geom_point(aes(color=layer),
                            position = "jitter",alpha=0.3)+
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  xlab("Plant Species") + ylab("log10(conspecific probability)") +
  labs(color = NULL)+ theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



means <- aggregate((consp_prob) ~  layer + Plot,
                   Prob_results, mean)

#Test significance of differences
library(ggpubr)
Prob_results_exp <- Prob_results_ini
Prob_results_exp$layer[Prob_results_exp$layer == "BEMA"] <- "B. macrocarpa"
Prob_results_exp$layer[Prob_results_exp$layer == "CETE"] <- "C. tenuiflorum"
Prob_results_exp$layer[Prob_results_exp$layer == "CHFU"] <- "C. fuscatum"
Prob_results_exp$layer[Prob_results_exp$layer == "CHMI"] <- "C. mixtum"
Prob_results_exp$layer[Prob_results_exp$layer == "LEMA"] <- "L. maroccanus"
Prob_results_exp$layer[Prob_results_exp$layer == "MESU"] <- "M. sulcatus"
Prob_results_exp$layer[Prob_results_exp$layer == "PUPA"] <- "P. paludosa"
Prob_results_exp$layer[Prob_results_exp$layer == "SCLA"] <- "S. laciniata"
Prob_results_exp$layer[Prob_results_exp$layer == "SOAS"] <- "S. asper"
Prob_results_exp$layer[Prob_results_exp$layer == "SPRU"] <- "S. rubra"

# We add equiprobability line


nodes <- tibble(Plot = c(1:9), nodes = c(49,66,66,26,19,14,43,84,74)+
                            0*c(23,30,30,11,8,7,22,41,34))

Prob_results_exp_nodes <- Prob_results_exp %>% left_join(nodes, by = "Plot") 

ggplot(Prob_results_exp_nodes,
       aes(x=layer,y=(consp_prob)))+
  geom_boxplot()+
  geom_point(aes(color=layer),position = "jitter",alpha=0.3)+
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+scale_y_continuous(trans='log10')+
  stat_summary(fun=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  #geom_text(data = means, aes(label = round(`(consp_prob)`,2), y = 0.1))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #geom_hline(data = Prob_results_exp_nodes, aes(yintercept = 1/nodes),linetype ="dashed")+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("Arrival probability of conspecific pollen")+
  labs(fill = NULL)+ theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))+ 
  theme(axis.text.x = element_text(face = "italic"))
#+stat_compare_means()

#Save 600 x 400

png("New_Figures/fig4_cons.png", width=1961*2, height = 1961*2*400/600, res=300*2)
ggplot(Prob_results_exp_nodes,
       aes(x=layer,y=(consp_prob)))+
  geom_boxplot()+
  geom_point(aes(color=layer),position = "jitter",alpha=0.3)+
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+scale_y_continuous(trans='log10')+
  stat_summary(fun=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  #geom_text(data = means, aes(label = round(`(consp_prob)`,2), y = 0.1))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #geom_hline(data = Prob_results_exp_nodes, aes(yintercept = 1/nodes),linetype ="dashed")+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("Arrival probability of conspecific pollen")+
  labs(fill = NULL)+ theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))+ 
  theme(axis.text.x = element_text(face = "italic"))
dev.off()

means_plot <- aggregate((consp_prob) ~  Plot,
          Prob_results_exp, mean)

cor(means_plot$`(consp_prob)`,nodes$nodes,method = "spearman")
cor.test(means_plot$`(consp_prob)`,nodes$nodes,method = "spearman")

means <- aggregate((consp_prob) ~  layer,
                   Prob_results_exp, mean)

#Test significance of differences
library(ggpubr)

ggplot(Prob_results_exp,aes(x=layer,y=(consp_prob)))+
  geom_boxplot()+
  geom_point(aes(color=layer),position = "jitter",alpha=0.3)+
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+
  stat_summary(fun=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = means, aes(label = round(`(consp_prob)`,3), y = 0.135))+
  #facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("Arrival probability of conspecific pollen")+ 
  labs(fill = NULL)+ theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#+stat_compare_means()


##########################################

png("New_Figures/fig4_hter.png", width=1961*2, height = 1961*2*400/600, res=300*2)
ggplot(Prob_results_exp_nodes %>% mutate(heter_prob=heter_prob+1e-20),
       aes(x=layer,y=(heter_prob)))+
  geom_boxplot()+
  geom_point(aes(color=layer),position = "jitter",alpha=0.3)+
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+scale_y_continuous(trans='log10')+
  stat_summary(fun=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  #geom_text(data = means, aes(label = round(`(consp_prob)`,2), y = 0.1))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #geom_hline(data = Prob_results_exp_nodes, aes(yintercept = 1/nodes),linetype ="dashed")+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("Arrival probability of heterospecific pollen")+
  labs(fill = NULL)+ theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))+ 
  theme(axis.text.x = element_text(face = "italic"))
dev.off()

means_plot <- aggregate((heter_prob) ~  Plot,
                        Prob_results_exp, mean)

cor(means_plot$`(heter_prob)`,nodes$nodes,method = "spearman")
cor.test(means_plot$`(heter_prob)`,nodes$nodes,method = "spearman")

means <- aggregate((heter_prob) ~  layer,
                   Prob_results_exp, mean)

#Test significance of differences
library(ggpubr)

ggplot(Prob_results_exp,aes(x=layer,y=(heter_prob)))+
  geom_boxplot()+
  geom_point(aes(color=layer),position = "jitter",alpha=0.3)+
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+
  stat_summary(fun=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = means, aes(label = round(`(heter_prob)`,3), y = 0.135))+
  #facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("Arrival probability of heterospecific pollen")+ 
  labs(fill = NULL)+ theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#+stat_compare_means()

#################

png("New_Figures/fig4_cons_VS_heter.png", width=1961*2, height = 1961*2*400/600, res=300*2)
ggplot(Prob_results_exp_nodes,# %>% filter(heter_prob>1e-20),
       aes(x=(consp_prob),y=(heter_prob)))+
  geom_abline(aes(slope= 1, intercept = 0),linetype ="dashed",color="deepskyblue", size=1.1)+
  geom_point(aes(color=layer),position="jitter",size=2,alpha=0.3)+
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+#scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')+
  # stat_summary(fun=mean, colour="darkred", geom="point", 
  #              shape=18, size=3,show.legend = FALSE) + 
  #geom_text(data = means, aes(label = round(`(consp_prob)`,2), y = 0.1))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Arrival probability of conspecific pollen") + ylab("Arrival probability of heterospecific pollen")+
  labs(fill = NULL)+ theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))+ 
  theme(axis.text.x = element_text(face = "italic"))+
  labs(color = NULL)+ theme(legend.position="bottom")+
  theme(legend.text = element_text(face = "italic"))#+
  # theme(axis.text=element_text(size=15),  axis.title=element_text(size=15,face="bold"))+                                                                # Change font size
  # theme(strip.text.x = element_text(size = 15))+
  # theme(legend.text = element_text(size=15),legend.title=element_text(size=15))+ 
  # guides(fill = guide_legend(override.aes = list(size=2)))
dev.off()
