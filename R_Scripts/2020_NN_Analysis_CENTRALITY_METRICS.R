
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
PageRank_results_exp <- PageRank_results
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "BEMA"] <- "B. macrocarpa"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "CETE"] <- "C. tenuiflorum"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "CHFU"] <- "C. fuscatum"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "CHMI"] <- "C. mixtum"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "LEMA"] <- "L. maroccanus"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "MESU"] <- "M. sulcatus"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "PUPA"] <- "P. paludosa"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "SCLA"] <- "S. laciniata"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "SOAS"] <- "S. asper"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "SPRU"] <- "S. rubra"

# We add equiprobability line


nodes <- tibble(Plot = c(1:9), nodes = c(49,66,66,26,19,14,43,84,74)+
                            c(23,30,30,11,8,7,22,41,34))

PageRank_results_exp_nodes <- PageRank_results_exp %>% left_join(nodes, by = "Plot") 

ggplot(PageRank_results_exp_nodes,
       aes(x=Plant_Simple,y=(Real_PR_Multi)))+
  geom_boxplot()+
  geom_point(aes(color=Plant_Simple),position = "jitter",alpha=0.3)+
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+scale_y_continuous(trans='log10')+
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  #geom_text(data = means, aes(label = round(`(Real_PR_Multi)`,2), y = 0.1))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  geom_hline(data = PageRank_results_exp_nodes, aes(yintercept = 1/nodes),linetype ="dashed")+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("Pollen arrival probability")+
  labs(fill = NULL)+ theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))+ 
  theme(axis.text.x = element_text(face = "italic"))
#+stat_compare_means()

#Save 600 x 400

png("New_Figures/fig4.png", width=1961*2, height = 1961*2*400/600, res=300*2)
ggplot(PageRank_results_exp_nodes,
       aes(x=Plant_Simple,y=(Real_PR_Multi)))+
  geom_boxplot()+
  geom_point(aes(color=Plant_Simple),position = "jitter",alpha=0.3)+
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+scale_y_continuous(trans='log10')+
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  #geom_text(data = means, aes(label = round(`(Real_PR_Multi)`,2), y = 0.1))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  geom_hline(data = PageRank_results_exp_nodes, aes(yintercept = 1/nodes),linetype ="dashed")+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("Pollen arrival probability")+
  labs(fill = NULL)+ theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))+ 
  theme(axis.text.x = element_text(face = "italic"))
dev.off()

means_plot <- aggregate((Real_PR_Multi) ~  Plot,
          PageRank_results_exp, mean)

cor(means_plot$`(Real_PR_Multi)`,nodes$nodes,method = "spearman")
cor.test(means_plot$`(Real_PR_Multi)`,nodes$nodes,method = "spearman")

means <- aggregate((Real_PR_Multi) ~  Plant_Simple,
                   PageRank_results_exp, mean)

#Test significance of differences
library(ggpubr)

ggplot(PageRank_results_exp,aes(x=Plant_Simple,y=(Real_PR_Multi)))+
  geom_boxplot()+
  geom_point(aes(color=Plant_Simple),position = "jitter",alpha=0.3)+
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = means, aes(label = round(`(Real_PR_Multi)`,3), y = 0.135))+
  #facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("Pollen arrival probability")+ 
  labs(fill = NULL)+ theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#+stat_compare_means()

# Calculate within group variances and between group variances
Plant_Simple_labels_exp <- PageRank_results_exp$Plant_Simple %>% unique()

run1 = PageRank_results_exp$Real_PR_Multi[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[1]]
group1 <- rep(Plant_Simple_labels_exp[1],length(run1))
plot1 = PageRank_results_exp$Plot[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[1]]

run2 = PageRank_results_exp$Real_PR_Multi[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[2]]
group2 <- rep(Plant_Simple_labels_exp[2],length(run2))
plot2 = PageRank_results_exp$Plot[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[2]]

run3 = PageRank_results_exp$Real_PR_Multi[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[3]]
group3 <- rep(Plant_Simple_labels_exp[3],length(run3))
plot3 = PageRank_results_exp$Plot[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[3]]

run4 = PageRank_results_exp$Real_PR_Multi[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[4]]
group4 <- rep(Plant_Simple_labels_exp[4],length(run4))
plot4 = PageRank_results_exp$Plot[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[4]]

run5 = PageRank_results_exp$Real_PR_Multi[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[5]]
group5 <- rep(Plant_Simple_labels_exp[5],length(run5))
plot5 = PageRank_results_exp$Plot[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[5]]

run6 = PageRank_results_exp$Real_PR_Multi[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[6]]
group6 <- rep(Plant_Simple_labels_exp[6],length(run6))
plot6 = PageRank_results_exp$Plot[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[6]]

run7 = PageRank_results_exp$Real_PR_Multi[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[7]]
group7 <- rep(Plant_Simple_labels_exp[7],length(run7))
plot7 = PageRank_results_exp$Plot[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[7]]

run8 = PageRank_results_exp$Real_PR_Multi[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[8]]
group8 <- rep(Plant_Simple_labels_exp[8],length(run8))
plot8 = PageRank_results_exp$Plot[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[8]]

run9 = PageRank_results_exp$Real_PR_Multi[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[9]]
group9 <- rep(Plant_Simple_labels_exp[9],length(run9))
plot9 = PageRank_results_exp$Plot[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[9]]

run10 = PageRank_results_exp$Real_PR_Multi[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[10]]
group10 <- rep(Plant_Simple_labels_exp[10],length(run10))
plot10 = PageRank_results_exp$Plot[
  PageRank_results_exp$Plant_Simple==Plant_Simple_labels_exp[10]]

runs = c(run1, run2, run3, run4, run5,
         run6, run7, run8, run9, run10)
group = c(group1, group2, group3, group4, group5,
         group6, group7, group8, group9, group10)
plot = c(plot1, plot2, plot3, plot4, plot5,
          plot6, plot7, plot8, plot9, plot10)

withinRunStats = function(x) c(sum = sum(x), mean = mean(x), var = var(x), n = length(x))
tapply(runs, group, withinRunStats)

data = data.frame(logPageRank = log(runs), Plant = factor(group),plot = factor(plot))
fit = lm(logPageRank ~ Plant, data)
fit
anova(fit)

# Error or within-group variance:
anova(fit)["Residuals", "Mean Sq"]

#Treatment or between-group variance:
anova(fit)["Plant", "Mean Sq"]

aovModel <- aov(logPageRank ~ Plant,data)
aovModel
summary(aovModel)



group_by(data , Plant) %>%
  summarise(
    count = n(),
    mean = mean(logPageRank, na.rm = TRUE),
    sd = sd(logPageRank, na.rm = TRUE)
  )

ggboxplot(data, x = "Plant", y = "logPageRank", 
          color = "Plant",
          ylab = "log PageRank", xlab = "Plant")
ggline(data, x = "Plant", y = "logPageRank", 
       add = c("mean_se", "jitter"),
       ylab = "log PageRank", xlab = "Plant")

res.aov <- aov(logPageRank ~ Plant, data = data)
# Summary of the analysis
summary(res.aov)
# As the p-value is less than the significance level 0.05, we can conclude that there are 
# significant differences between the groups highlighted with “*" in the model summary.

TukeyHSD(res.aov)
# It can be seen from the output, that the difference between
# L. maroccanus and C. fuscatum is significant with an adjusted p-value of 0.0000003,
# and that of P. paludosa-L. maroccanus with an adjusted p-value of 0.0016098

# 1. Homogeneity of variances
plot(res.aov, 1) # There are some outlayer that may affect the result
library(car)
leveneTest(logPageRank ~ Plant, data = data)
# From the output above we can see that the p-value is less than the significance level 
# of 0.05. This means that there is evidence to suggest that the variance across groups 
# is statistically significantly different. Therefore, we can not assume the homogeneity 
# of variances in the different treatment groups.

# ANOVA test with no assumption of equal variances
oneway.test(logPageRank ~ Plant, data = data)

pairwise.t.test(data$logPageRank, data$Plant,
                p.adjust.method = "BH", pool.sd = FALSE)

# 2. Normality
plot(res.aov, 2) # No normality
ks.test(res.aov)
plot(res.aov)

png("New_Figures/figA121.png", width=1461*2, height = 1461*2*361/515, res=300*2)
plot(res.aov,2)
dev.off()
#
kruskal.test(exp(y) ~ group, data = data)
# As the p-value is less than the significance level 0.05, we can conclude that there are 
# significant differences between the treatment groups.

pairwise.wilcox.test(data$y, data$group,
                     p.adjust.method = "BH")

pairwise.wilcox.test(exp(data$y), data$group,
                     p.adjust.method = "BH")
########
aggregate((Ratio) ~  Plot,
          PageRank_results_exp, mean)

means <- aggregate((Ratio) ~  Plant_Simple,
                   PageRank_results_exp, mean)


ggplot(PageRank_results_exp,aes(x=Plant_Simple,y=(Ratio)))+
  geom_boxplot()+
  geom_point(aes(color=Plant_Simple),position = "jitter",alpha=0.3)+
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = means, aes(label = round(`(Ratio)`,3), y = 2.5))+
  #facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("Among layer centrality ratio")+ 
  labs(fill = NULL)+ theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#+stat_compare_means()

#Test significance of differences
library(ggpubr)

kruskal.test(Real_PR_Multi ~ Plant_Simple, data = PageRank_results)

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
PageRank_results_exp <- PageRank_results
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "BEMA"] <- "B. macrocarpa"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "CETE"] <- "C. tenuiflorum"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "CHFU"] <- "C. fuscatum"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "CHMI"] <- "C. mixtum"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "LEMA"] <- "L. maroccanus"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "MESU"] <- "M. sulcatus"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "PUPA"] <- "P. paludosa"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "SCLA"] <- "S. laciniata"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "SOAS"] <- "S. asper"
PageRank_results_exp$Plant_Simple[PageRank_results_exp$Plant_Simple == "SPRU"] <- "S. rubra"

png("New_Figures/figA101.png", width=1961*2, height = 1961*2*361/595, res=300*2)
ggplot(PageRank_results_exp,aes(x=(StrengthIn),y=(Real_PR_Multi)))+
  geom_point(aes(color=as.factor(Plant_Simple)),position = "jitter",alpha=0.5)+
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_cor(method = "pearson", label.x = 0.01, label.y = 0.04)+
  #coord_trans(y = 'log10')+#scale_y_continuous(trans='log10')+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  xlab("Within layer centrality") + ylab("Pollen arrival probability") +
  labs(color = NULL)+ theme(legend.position="bottom")+ylim(0,0.05)+
  theme(legend.text = element_text(face = "italic"))
dev.off()
  



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

png("New_Figures/figA102.png", width=1961*2, height = 1961*2*441/515, res=300*2)
ggplot(PageRank_results_exp,aes(x=Plant_Simple,y=(Ratio),fill=Plant_Simple))+
  geom_boxplot()+
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  #geom_text(data = means, aes(label = round(`(Ratio)`,2), y = `(Ratio)` + 1.1))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  scale_y_continuous(trans='log10')+
  xlab("Plant Species") + ylab("Among-layer centrality")+
  theme_bw()+ theme(legend.position = "none")+#+stat_compare_means()
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  theme(axis.text.x = element_text(face = "italic"))
dev.off()

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



