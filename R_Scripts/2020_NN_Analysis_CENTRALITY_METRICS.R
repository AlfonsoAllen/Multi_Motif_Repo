
library(tidyverse)

PageRank_results_ini <- read_csv("Processed_data/2020_NN_NEW_data_models_phenol_overlap.csv")

# sanity check
PageRank_results_ini %>% filter(ID!="None") %>% select(Plot,Plant,Real_PR_Multi) %>%
  unique() %>% group_by(Plot) %>% count(wt=Real_PR_Multi)
# Page rank values are smaller than 1 because we are not considering the insects


# Extract pagerank results for plants
PageRank_results <- PageRank_results_ini %>% 
  filter(ID!="None") %>% rename(Plant_Simple=Plant) %>%
  select(Plot,Subplot,Plant_Simple,Real_PR_Multi,Real_PR_Layer,StrengthIn,
         Ratio,visits_GF) %>%
  unique()


PageRank_results %>% group_by(Plant_Simple) %>% count()

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

ggplot(PageRank_results,aes(x=Plant_Simple,y=log10(Real_PR_Multi)))+
  #geom_violin()+ 
  
  geom_boxplot()+geom_point(aes(color=Plant_Simple),
                            position = "jitter",alpha=0.3)+
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  xlab("Plant Species") + ylab("log10(Page Rank)") +
  labs(color = NULL)+ theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



means <- aggregate((Real_PR_Multi) ~  Plant_Simple + Plot,
                   PageRank_results, mean)

#Test significance of differences
library(ggpubr)

ggplot(PageRank_results,
       aes(x=Plant_Simple,y=(Real_PR_Multi)))+
  geom_boxplot()+
  geom_point(aes(color=Plant_Simple),position = "jitter",alpha=0.3)+
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+scale_y_continuous(trans='log10')+
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  #geom_text(data = means, aes(label = round(`(Real_PR_Multi)`,2), y = 0.1))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("PageRank")+
  labs(fill = NULL)+ theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#+stat_compare_means()

aggregate((Real_PR_Multi) ~  Plot,
          PageRank_results, mean)

means <- aggregate((Real_PR_Multi) ~  Plant_Simple,
                   PageRank_results, mean)

#Test significance of differences
library(ggpubr)

ggplot(PageRank_results,aes(x=Plant_Simple,y=(Real_PR_Multi)))+
  geom_boxplot()+
  geom_point(aes(color=Plant_Simple),position = "jitter",alpha=0.3)+
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = means, aes(label = round(`(Real_PR_Multi)`,3), y = 0.135))+
  #facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("PageRank")+ 
  labs(fill = NULL)+ theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#+stat_compare_means()



#Test significance of differences
library(ggpubr)

PageRank_results1 <- PageRank_results %>% filter(Plot==1)
PageRank_results2 <- PageRank_results %>% filter(Plot==2)
PageRank_results3 <- PageRank_results %>% filter(Plot==3)
PageRank_results4 <- PageRank_results %>% filter(Plot==4)
PageRank_results5 <- PageRank_results %>% filter(Plot==5)
PageRank_results6 <- PageRank_results %>% filter(Plot==6)
PageRank_results7 <- PageRank_results %>% filter(Plot==7)
PageRank_results8 <- PageRank_results %>% filter(Plot==8)
PageRank_results9 <- PageRank_results %>% filter(Plot==9)

# If the p-value is less than the significance level 0.05, 
# we can conclude that there are significant differences between
# the plant species.

# Plot1
kruskal.test(Real_PR_Multi ~ Plant_Simple, data = PageRank_results1)
pairwise.wilcox.test(PageRank_results1$Real_PR_Multi, PageRank_results1$Plant_Simple,
                     p.adjust.method = "BH")

# Plot2
kruskal.test(Real_PR_Multi ~ Plant_Simple, data = PageRank_results2)
pairwise.wilcox.test(PageRank_results2$Real_PR_Multi, PageRank_results2$Plant_Simple,
                     p.adjust.method = "BH")
# Plot3
kruskal.test(Real_PR_Multi ~ Plant_Simple, data = PageRank_results3)
pairwise.wilcox.test(PageRank_results3$Real_PR_Multi, PageRank_results3$Plant_Simple,
                     p.adjust.method = "BH")
# Plot4
kruskal.test(Real_PR_Multi ~ Plant_Simple, data = PageRank_results4)
pairwise.wilcox.test(PageRank_results4$Real_PR_Multi, PageRank_results4$Plant_Simple,
                     p.adjust.method = "BH")
# Plot5
kruskal.test(Real_PR_Multi ~ Plant_Simple, data = PageRank_results5)
pairwise.wilcox.test(PageRank_results5$Real_PR_Multi, PageRank_results5$Plant_Simple,
                     p.adjust.method = "BH")
# Plot6
kruskal.test(Real_PR_Multi ~ Plant_Simple, data = PageRank_results6)
pairwise.wilcox.test(PageRank_results6$Real_PR_Multi, PageRank_results6$Plant_Simple,
                     p.adjust.method = "BH")
# Plot7
kruskal.test(Real_PR_Multi ~ Plant_Simple, data = PageRank_results7)
pairwise.wilcox.test(PageRank_results7$Real_PR_Multi, PageRank_results7$Plant_Simple,
                     p.adjust.method = "BH")

# Plot8
kruskal.test(Real_PR_Multi ~ Plant_Simple, data = PageRank_results8)
pairwise.wilcox.test(PageRank_results8$Real_PR_Multi, PageRank_results8$Plant_Simple,
                     p.adjust.method = "BH")

# Plot9
kruskal.test(Real_PR_Multi ~ Plant_Simple, data = PageRank_results9)
pairwise.wilcox.test(PageRank_results9$Real_PR_Multi, PageRank_results9$Plant_Simple,
                     p.adjust.method = "BH")

#############################################

means_caracoles <- aggregate((Real_PR_Multi) ~  Plant_Simple,
                   PageRank_results, mean)

ggplot(PageRank_results,aes(x=Plant_Simple,y=(Real_PR_Multi)))+
  geom_boxplot()+
  geom_point(aes(color=Plant_Simple),position = "jitter",alpha=0.3)+
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = means_caracoles, aes(label = round(`(Real_PR_Multi)`,3), y = 0.075))+
  #facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("PageRank")+ 
  labs(fill = NULL)+ theme(legend.position="none")+stat_compare_means()

library(ggpubr)

pairwise.wilcox.test(PageRank_results$Real_PR_Multi, PageRank_results$Plot,
                     p.adjust.method = "BH")

aggregate((Real_PR_Multi) ~ Plot,
          PageRank_results, mean)


library(ggpmisc)

ggplot(PageRank_results,aes(x=(StrengthIn),y=(Real_PR_Multi)))+
  geom_point(aes(color=as.factor(Plant_Simple)),position = "jitter",alpha=0.5)+
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_cor(method = "pearson", label.x = 0.01, label.y = 0.04)+
  #coord_trans(y = 'log10')+#scale_y_continuous(trans='log10')+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  xlab("In-Strength") + ylab("Multilayer Page Rank") +
  labs(color = NULL)+ theme(legend.position="bottom")+ylim(0,0.05)
  



ggplot(PageRank_results)+
  geom_point(aes(x=log10(Real_PR_Layer),y=log10(Real_PR_Multi),color=Plant_Simple,
                 shape=Plant_Simple),
             position = "jitter")+
  geom_smooth(aes(x=log10(Real_PR_Layer),y=log10(Real_PR_Multi),color=Plant_Simple),method="lm", se=F)+
  geom_abline(intercept = 0, slope = 1,color='steelblue', 
            size=1, alpha=0.4,linetype = "dashed")+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  theme_bw()+
  #ggtitle(paste0("Plot ",i)) +
  xlab("log10(LC)") + ylab("log10(MC)") +
  labs(color = "Plant Species", shape = "Plant Species")+ theme(legend.position="bottom")

means <- aggregate((Ratio) ~  Plant_Simple + Plot,
                   PageRank_results, mean)


x <- PageRank_results %>% filter( Ratio > 1)

means <- aggregate((Ratio) ~  Plant_Simple+Plot,
                   PageRank_results, mean)

ggplot(PageRank_results,aes(x=Plant_Simple,y=(Ratio),fill=Plant_Simple))+
  geom_boxplot()+
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = means, aes(label = round(`(Ratio)`,2), y = `(Ratio)` + 1.1))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  scale_y_continuous(trans='log10')+
  xlab("Plant Species") + ylab("R=MC/LC")+
  theme_bw()+ theme(legend.position = "none")+#+stat_compare_means()
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#labs(Color = "Plant Species")

means_ratio <- aggregate(Ratio ~  Plant_Simple+Plot,
                   PageRank_results, mean)

number_plants <- PageRank_results %>% select(Plot,Subplot,Plant_Simple) %>% unique() %>%
  group_by(Plot,Plant_Simple) %>%
  count()

total_number_plants <- PageRank_results %>% select(Plot,Subplot,Plant_Simple) %>% unique() %>%
  group_by(Plot) %>%
  count() %>% rename(total=n)

ratio_plants <- means_ratio %>% 
  left_join(number_plants,by=c("Plot","Plant_Simple")) %>%
  left_join(total_number_plants,by=c("Plot")) %>%
  mutate(perc=100*n/total)

ggplot(ratio_plants,aes(x=(perc),y=(Ratio)))+
  geom_point(aes(color=as.factor(Plant_Simple)),position = "jitter",alpha=0.5)+
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_cor(method = "pearson", label.x = 0.01, label.y = 0.04)+
  #coord_trans(y = 'log10')+#scale_y_continuous(trans='log10')+
  #facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  xlab("number of focal plants per plot") + ylab("Ratio") +
  labs(color = NULL)+ theme(legend.position="bottom")

means <- aggregate((Ratio) ~  Plant_Simple,
                   PageRank_results, mean)

ggplot(PageRank_results,aes(x=Plant_Simple,y=(Ratio),fill=Plant_Simple))+
  geom_boxplot()+
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = means, aes(label = round(`(Ratio)`,2), y = `(Ratio)` + 1.1))+
  #facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  scale_y_continuous(trans='log10')+
  xlab("Plant Species") + ylab("R=MC/LC")+
  theme_bw()+ theme(legend.position = "none")+#+stat_compare_means()
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#labs(Color = "Plant Species")

means <- aggregate((Ratio) ~  Plant_Simple,
                   PageRank_results, mean)


ggplot(PageRank_results,aes(x=Plant_Simple,y=(Ratio),fill=Plant_Simple))+
  geom_boxplot()+
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = means, aes(label = round(`(Ratio)`,3), y = `(Ratio)` + 0.3))+
  #facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("R=MC/LC")+
  theme_bw()+ theme(legend.position = "none")#+stat_compare_means()
#labs(Color = "Plant Species")



# Test significance

# Plot1
kruskal.test(Ratio ~ Plant_Simple, data = PageRank_results1)
pairwise.wilcox.test(PageRank_results1$Ratio, PageRank_results1$Plant_Simple,
                     p.adjust.method = "BH")

# Plot2
kruskal.test(Ratio ~ Plant_Simple, data = PageRank_results2)
pairwise.wilcox.test(PageRank_results2$Ratio, PageRank_results2$Plant_Simple,
                     p.adjust.method = "BH")
# Plot3
kruskal.test(Ratio ~ Plant_Simple, data = PageRank_results3)
pairwise.wilcox.test(PageRank_results3$Ratio, PageRank_results3$Plant_Simple,
                     p.adjust.method = "BH")
# Plot4
kruskal.test(Ratio ~ Plant_Simple, data = PageRank_results4)
pairwise.wilcox.test(PageRank_results4$Ratio, PageRank_results4$Plant_Simple,
                     p.adjust.method = "BH")
# Plot5
kruskal.test(Ratio ~ Plant_Simple, data = PageRank_results5)
pairwise.wilcox.test(PageRank_results5$Ratio, PageRank_results5$Plant_Simple,
                     p.adjust.method = "BH")
# Plot6
kruskal.test(Ratio ~ Plant_Simple, data = PageRank_results6)
pairwise.wilcox.test(PageRank_results6$Ratio, PageRank_results6$Plant_Simple,
                     p.adjust.method = "BH")
# Plot7
kruskal.test(Ratio ~ Plant_Simple, data = PageRank_results7)
pairwise.wilcox.test(PageRank_results7$Ratio, PageRank_results7$Plant_Simple,
                     p.adjust.method = "BH")

# Plot8
kruskal.test(Ratio ~ Plant_Simple, data = PageRank_results8)
pairwise.wilcox.test(PageRank_results8$Ratio, PageRank_results8$Plant_Simple,
                     p.adjust.method = "BH")

# Plot9
kruskal.test(Ratio ~ Plant_Simple, data = PageRank_results9)
pairwise.wilcox.test(PageRank_results9$Ratio, PageRank_results9$Plant_Simple,
                     p.adjust.method = "BH")



