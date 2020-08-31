
library(tidyverse)

#########################
# Layer Centrality (LC)
########################

for (i in 1:9){

# Layer Centrality
file_CL_i <- paste0("Processed_data/Muxviz_Pheno_Overlap/2020_centrality_perLayer_table_Plot",i)
centrality_layer <- read_delim(file_CL_i,";")

# Page rank results from Muxviz are not normalized
# We get the layer normalization constant = sum pargerank results
layer_norm_const <- centrality_layer %>% group_by(Layer) %>% 
  count(wt=PageRank) %>% rename(sum_weight=n)

centrality_layer2 <- centrality_layer %>% 
  left_join(layer_norm_const,by="Layer") %>% mutate(Real_PR_Layer=PageRank/sum_weight) %>%
  select(-Node,-Eigenvector,-Hub,-Authority,-Katz,-Multiplexity,-Kcore,-sum_weight,-PageRank)

#Sanity check
centrality_layer2 %>% group_by(Layer) %>% count(wt=Real_PR_Layer)

file_i <- paste0("Processed_data/Muxviz_Pheno_Overlap/Plot_",i,"/general_multilayer_layers_Plot",i,".txt")
layer_ID <- read.delim(file_i,header = T,sep =" ",stringsAsFactors =F) %>% 
  rename(Layer=layerID)
layer_ID$Layer <- as.character(layer_ID$Layer)
  
centrality_layer3 <- centrality_layer2 %>% left_join(layer_ID,by="Layer") %>%
  select(-Layer) %>% rename(Layer=layerLabel)
centrality_layer3$Layer[is.na(centrality_layer3$Layer)] <- "Aggr"

centrality_layer3 %>% group_by(Layer) %>% count(wt=Real_PR_Layer)
centrality_layer3$Plot <- i

centrality_layer3 %>% group_by(Plot,Layer) %>% summarise(Real_PR_Layer=sum(Real_PR_Layer))

if(i==1){
  LC_list <- centrality_layer3
}else{
  LC_list <- bind_rows(LC_list, centrality_layer3)
}

}

LC_list_fil <- LC_list %>% separate(Label,c("Ind","Plant_Simple")," ") %>% filter(!is.na(Plant_Simple)) %>% 
  filter(Layer==Plant_Simple) %>% mutate(Label=paste0(Ind," ",Plant_Simple)) %>% select(-Ind)


########################################
# Multilayer Centrality
########################################

for (i in 1:9){
  
file_CML_i <- paste0("Processed_data/Muxviz_Pheno_Overlap/2020_centrality_table_Plot",i)
centrality_multilayer <- read_delim(file_CML_i,";")
#layer_ID_Multi <- layer_ID
#layer_ID_Multi$Layer <- paste0(layer_ID_Multi$Layer,"-Multi")
#centrality_multilayer <- centrality_multilayer %>% 
#  left_join(layer_ID_Multi,by="Layer") %>% filter(Degree!=0)

#Sanity check
centrality_multilayer %>% filter(is.na(PageRank))

centrality_multilayer$Layer[centrality_multilayer$Layer!="Aggr"] <- "Multi"

centrality_multilayer_fil <- unique(centrality_multilayer)

centrality_multilayer_fil %>% group_by(Layer) %>% count(wt=PageRank)

# Page rank results from Muxviz are not normalized
# We get the layer normalization constant = sum pargerank results
Mlayer_norm_const <- centrality_multilayer_fil %>% group_by(Layer) %>% 
  count(wt=PageRank) %>% rename(sum_weight=n)

centrality_multilayer_fil2 <- centrality_multilayer_fil %>% 
  left_join(Mlayer_norm_const,by="Layer") %>% mutate(Real_PR_Multi=PageRank/sum_weight) %>%
  select(-Node,-Eigenvector,-Hub,-Authority,-Katz,-Multiplexity,-Kcore,-sum_weight,-PageRank)

#Sanity check
centrality_multilayer_fil2 %>% group_by(Layer) %>% count(wt=Real_PR_Multi)
centrality_multilayer_fil2 %>% group_by(Label) %>% count() %>% filter(n>2)
centrality_multilayer_fil2 %>% group_by(Label) %>% count() %>% filter(n==2)
centrality_multilayer_fil2 %>% filter(is.na(Real_PR_Multi))
centrality_multilayer_fil2 %>% filter(Real_PR_Multi==0)


centrality_multilayer_fil2$Plot <- i

if(i==1){
  MC_list <- centrality_multilayer_fil2
}else{
  MC_list <- bind_rows(MC_list, centrality_multilayer_fil2)
}

}

sum(MC_list$Real_PR_Multi[MC_list$Plot==1 & MC_list$Layer=="Multi"])

MC_list_fil <- MC_list %>% separate(Label,c("Ind","Plant_Simple")," ") %>% filter(!is.na(Plant_Simple)) %>% 
  mutate(Label=paste0(Ind," ",Plant_Simple)) %>% select(-Ind)

sum(MC_list_fil$Real_PR_Multi[MC_list_fil$Plot==1 & MC_list_fil$Layer=="Multi"])

##################


# Merge LC and MC lists for comparison

MC_list_final <- MC_list_fil %>% filter(Layer!="Aggr") %>% select(Real_PR_Multi, Plot, Label, Plant_Simple)
LC_list_final <- LC_list_fil %>% filter(Layer!="Aggr") %>% select(Real_PR_Layer, Plot, Label, Plant_Simple)

PageRank_results <- left_join(MC_list_final,LC_list_final,by=c("Plot","Label","Plant_Simple")) %>%
  mutate(Delta=Real_PR_Multi-Real_PR_Layer,Ratio=Real_PR_Multi/Real_PR_Layer)

sum(PageRank_results$Real_PR_Multi[PageRank_results$Plot==1])

#write_csv(PageRank_results,"2020_PageRank_results.csv")


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

ggplot(PageRank_results,aes(x=Plant_Simple,y=(Real_PR_Multi)))+
  geom_violin()+ 
  geom_point(aes(color=Plant_Simple),
             position = "jitter",alpha=0.2)+
  geom_boxplot(width=0.1)+
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  xlab("Plant Species") + ylab("Multilayer Page Rank") +
  labs(color = NULL)+ theme(legend.position="bottom")

means <- aggregate((Real_PR_Multi) ~  Plant_Simple + Plot,
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
  geom_text(data = means, aes(label = round(`(Real_PR_Multi)`,3), y = 0.04))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("PageRank")+ 
  labs(fill = NULL)+ theme(legend.position="none")+stat_compare_means()
  


PageRank_results1 <- PageRank_results %>% filter(Plot==1)
PageRank_results2 <- PageRank_results %>% filter(Plot==2)
PageRank_results3 <- PageRank_results %>% filter(Plot==3)
PageRank_results4 <- PageRank_results %>% filter(Plot==4)
PageRank_results5 <- PageRank_results %>% filter(Plot==5)
PageRank_results6 <- PageRank_results %>% filter(Plot==6)
PageRank_results7 <- PageRank_results %>% filter(Plot==7)
PageRank_results8 <- PageRank_results %>% filter(Plot==8)
PageRank_results9 <- PageRank_results %>% filter(Plot==9)

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
  geom_text(data = means_caracoles, aes(label = round(`(Real_PR_Multi)`,3), y = 0.04))+
  #facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("PageRank")+ 
  labs(fill = NULL)+ theme(legend.position="none")+stat_compare_means()

library(ggpubr)

pairwise.wilcox.test(PageRank_results$Real_PR_Multi, PageRank_results$Plot,
                     p.adjust.method = "BH")

aggregate((Real_PR_Multi) ~ Plot,
          PageRank_results, mean)



x <- MC_list %>% separate(Label,c("Sub","Spec")," ")
x$Spec[is.na(x$Spec)] <- "Visitor"

library(ggpmisc)

ggplot(x %>% filter(Spec!="Visitor"),aes(x=(StrengthIn),y=(Real_PR_Multi)))+
  geom_point(aes(color=as.factor(Spec)),position = "jitter",alpha=0.5)+
  geom_smooth(method = "lm")+
  theme_bw()+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  xlab("In-Strength") + ylab("Multilayer Page Rank") +
  labs(color = NULL)+ theme(legend.position="bottom")+
  stat_cor(method = "pearson", label.x = 0.01, label.y = 0.05)



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


PageRank_results %>% filter(Real_PR_Layer<Real_PR_Multi)

ggplot(PageRank_results,aes(x=Plant_Simple,y=(Ratio),fill=Plant_Simple))+
  geom_boxplot()+
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  geom_text(data = means, aes(label = round(`(Ratio)`,3), y = `(Ratio)` + 0.3))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
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


PageRank_results %>% filter(Plot==1) %>% summarise(Real_PR_Multi=sum(Real_PR_Multi))

###########################
# MODULAR STRUCTURE
###########################

for (Plot_i in 1:9){

  module_i <- read_csv(paste0("Processed_data/Modularity_Pheno_Overlap/Modularity_Plot",Plot_i,".csv")  )
  
  if (Plot_i==1){
    modules_final=module_i
  }else{
    modules_final=bind_rows(modules_final,module_i)
  }

}

plant_m <- modules_final %>% filter(type=="plant") %>% dplyr::select(Plot,module,layer_name) %>% 
  unique() %>% arrange(Plot,module,layer_name)%>% group_by(Plot,module) %>% count()

all_m <- modules_final %>%  dplyr::select(Plot,module,layer_name) %>% 
  unique() %>% arrange(Plot,module,layer_name)%>% group_by(Plot,module) %>% count()

modules_final %>% filter(module==2,Plot==8)
#Los módulos sólo de plantas tienen menos especies que
#los módulos con insectos, ya que éstos están presentes en otras capas
#module==2,Plot==8 tiene a pollinator Dilophus_sp. PUPA sin plantas pupa




modules_final %>% filter(type=="plant") %>% dplyr::select(Plot,module,layer_name) %>% 
  unique() %>% arrange(Plot,module,layer_name) %>% group_by(Plot,module) %>% count() %>%
  ggplot()+
  geom_histogram(aes(x=n))+theme_bw()+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("# Layers per module") + ylab("# Modules") + theme(legend.position = "none")
#Recoge a través de cuantas capas se extiende el módulo


# Plant modules data

modules_final %>% filter(type=="plant") %>% dplyr::select(Plot,module,layer_name) %>% 
  unique() %>% arrange(Plot,module,layer_name) %>% group_by(Plot,module) %>% count() %>%
  group_by(Plot,n) %>% count()
################################################

modules <- modules_final %>% 
  dplyr::select(Plot,module) %>% unique()

modules_spread <- modules %>% group_by(Plot,module) %>% 
  count() %>% spread(module,n)
modules_spread[is.na(modules_spread)] <- 0
modules_spread$total <- rowSums(modules_spread[,c(2:ncol(modules_spread))])

mean(modules_spread$total)
sd(modules_spread$total)

library(boot)
b <- boot(modules_spread$total, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(number_plant_species_per_module_i$n, probs = c(0.025, 0.975), na.rm = FALSE)

descrip_plants <- modules_final %>% filter(type=="plant") %>% dplyr::select(Plot,module,layer_name) %>% 
  unique() %>% arrange(Plot,module,layer_name) %>% group_by(Plot,module) %>% count() %>%
  group_by(Plot,n) %>% count()

sum(descrip_plants$nn)
sum(descrip_plants$nn[descrip_plants$n==1])
sum(descrip_plants$nn[descrip_plants$n==2])
sum(descrip_plants$nn[descrip_plants$n==3])
sum(descrip_plants$nn[descrip_plants$n==4])

#total amount of focal plants

modules_plants <- modules_final %>% filter(type=="plant") %>% 
  dplyr::select(Plot,module,layer_name,species)

species_plants_module <- modules_plants %>% group_by(Plot,module,layer_name) %>% 
  count() %>% spread(layer_name,n)
species_plants_module[is.na(species_plants_module)] <- 0
species_plants_module$total <- rowSums(species_plants_module[,c(3:ncol(species_plants_module))])

mean(species_plants_module$total)
sd(species_plants_module$total)

library(boot)
b <- boot(species_plants_module$total, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(species_plants_module$total, probs = c(0.025, 0.975), na.rm = FALSE)

# Total Pollinators

descrip_poll <- modules_final %>% filter(type!="plant") %>% dplyr::select(Plot,module,layer_name) %>% 
  arrange(Plot,module,layer_name) %>% group_by(Plot,module) %>% count() %>%
  group_by(Plot,n) %>% count()

sum(descrip_poll$nn)
sum(descrip_poll$nn[descrip_poll$n==1])
sum(descrip_poll$nn[descrip_poll$n==2])
sum(descrip_poll$nn[descrip_poll$n==3])
sum(descrip_poll$nn[descrip_poll$n==4])
sum(descrip_poll$nn[descrip_poll$n==5])

modules_poll <- modules_final %>% filter(type!="plant") %>% 
  dplyr::select(Plot,module,layer_name,species)

species_poll_module <- modules_poll %>% group_by(Plot,module,species) %>% 
  count() %>% spread(species,n)
species_poll_module[is.na(species_poll_module)] <- 0
species_poll_module$total <- rowSums(species_poll_module[,c(3:ncol(species_poll_module))])

mean(species_poll_module$total)
sd(species_poll_module$total)

library(boot)
b <- boot(species_poll_module$total, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(number_plant_species_per_module_i$n, probs = c(0.025, 0.975), na.rm = FALSE)



##################################################################3

perc_sp_plants_mod <- species_plants_module %>% rename(abs_CHFU=CHFU,
                                                       abs_CHMI=CHMI,
                                                       abs_LEMA=LEMA,
                                                       abs_ME=ME,
                                                       abs_PUPA=PUPA) %>%
  mutate(CHFU=abs_CHFU/total,
         CHMI=abs_CHMI/total,
         LEMA=abs_LEMA/total,
         ME=abs_ME/total,
         PUPA=abs_PUPA/total) %>%
  gather(layer_name,pecentage_same_plants,c(CHFU, CHMI, LEMA, ME, PUPA)) %>%
  dplyr::select(Plot, module,layer_name,pecentage_same_plants,total)

modules_plants2 <- modules_plants %>% left_join(perc_sp_plants_mod,by=c("Plot", "module", "layer_name"))


means_per <- aggregate(pecentage_same_plants ~  layer_name + Plot,
                       modules_plants2, mean)

ggplot(modules_plants2,aes(x=layer_name,y=pecentage_same_plants,fill=layer_name))+
  geom_boxplot()+
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3, show.legend = FALSE) + 
  geom_text(data = means_per, aes(label = round(pecentage_same_plants,3), y = 1.1))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller = labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("%  of individuals of same plant species in the module")+ theme(legend.position = "none")
#labs(Color = "Plant Species")

means_per2 <- aggregate(pecentage_same_plants ~  layer_name,
                       modules_plants2, mean)

ggplot(modules_plants2,aes(x=layer_name,y=pecentage_same_plants,fill=layer_name))+
  geom_boxplot()+
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3, show.legend = FALSE) + 
  geom_text(data = means_per2, aes(label = round(pecentage_same_plants,3), y = 1.1))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant Species") + ylab("%  of individuals of same plant species in the module")+ theme(legend.position = "none")
#labs(Color = "Plant Species")

ggplot(modules_plants2,aes(x=total,y=pecentage_same_plants,fill=layer_name))+
  geom_point(aes(color=layer_name,
                 shape=layer_name),
             position = "jitter")+
  geom_smooth(aes(x=total,y=pecentage_same_plants,color=layer_name),method="lm",
              formula=y~x,se=F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Total individual per module") + ylab("%  of individuals of same plant species in the module")+ theme(legend.position = "none")
#labs(Color = "Plant Species")

######
# SHANON DIVERSITY INDEX
#######

shannon <- species_plants_module %>% rename(abs_CHFU=CHFU,
                                                       abs_CHMI=CHMI,
                                                       abs_LEMA=LEMA,
                                                       abs_ME=ME,
                                                       abs_PUPA=PUPA) %>%
  mutate(CHFU=abs_CHFU/total,
         CHMI=abs_CHMI/total,
         LEMA=abs_LEMA/total,
         ME=abs_ME/total,
         PUPA=abs_PUPA/total) 


shannon$entropy <- NA

for (i in 1:nrow(shannon)){
  
  shannon$entropy[i] <- 0
  if (shannon$CHFU[i]>0){shannon$entropy[i] <-shannon$entropy[i] -shannon$CHFU[i]*log(shannon$CHFU[i]) }
  if (shannon$CHMI[i]>0){shannon$entropy[i] <-shannon$entropy[i] -shannon$CHMI[i]*log(shannon$CHMI[i]) }
  if (shannon$LEMA[i]>0){shannon$entropy[i] <-shannon$entropy[i] -shannon$LEMA[i]*log(shannon$LEMA[i]) }
  if (shannon$ME[i]>0){shannon$entropy[i] <-shannon$entropy[i] -shannon$ME[i]*log(shannon$ME[i]) }
  if (shannon$PUPA[i]>0){shannon$entropy[i] <-shannon$entropy[i] -shannon$PUPA[i]*log(shannon$PUPA[i]) }
}

ggplot(shannon %>% filter(entropy>=0),aes(x=total,y=entropy))+
  geom_point(position = "jitter")+
  geom_smooth(method="lm",
              formula=y~x,se=F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  xlab("Total individual per module") + ylab("entropy")+ theme(legend.position = "none")
#labs(Color = "Plant Species")


ggplot(shannon %>% filter(entropy>=0),aes(x=total,y=entropy))+
  geom_point(position = "jitter")+
  geom_smooth(method="lm",
              formula=y~x,se=F)+
  xlab("Total individual per module") + ylab("entropy")+ theme(legend.position = "none")
#labs(Color = "Plant Species")

#############
#
#######################

modulos <- shannon
modulos$Plant <- NA

for (i in 1:nrow(modulos)){
  
  if(modulos$abs_CHFU[i]>0 & modulos$entropy[i]==0){
    modulos$Plant[i] <- "Only CHFU"
  }
  else if(modulos$abs_PUPA[i]>0 & modulos$entropy[i]==0){
    modulos$Plant[i] <- "Only PUPA"
  }else if(modulos$abs_LEMA[i]>0 & modulos$entropy[i]==0){
    modulos$Plant[i] <- "Only LEMA"
  }else if(modulos$abs_CHMI[i]>0 & modulos$entropy[i]==0){
    modulos$Plant[i] <- "Only CHMI"
  }else if (modulos$abs_ME[i]>0 & modulos$entropy[i]==0){
    modulos$Plant[i] <- "Only ME"
  }else{modulos$Plant[i] <- "Mixed species"}
  
}


modulos_group <- modulos %>% group_by(Plot,Plant) %>% count() 

library(RColorBrewer)
ggplot(modulos_group, aes(fill=Plant, y=n, x=as.factor(Plot))) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw()+
  scale_fill_brewer(palette = 'Paired')+
  labs(x = "Plot", y = "Number of modules",fill=NULL)


modulos2 <- shannon %>% gather(Plant,Amount,c("abs_CHFU", "abs_CHMI", "abs_LEMA", "abs_ME", "abs_PUPA"))
modulos2$Plant[modulos2$Plant=="abs_CHFU"] <- "CHFU"
modulos2$Plant[modulos2$Plant=="abs_CHMI"] <- "CHMI"
modulos2$Plant[modulos2$Plant=="abs_LEMA"] <- "LEMA"
modulos2$Plant[modulos2$Plant== "abs_ME"] <-  "ME"
modulos2$Plant[modulos2$Plant=="abs_PUPA"] <- "PUPA"
  
  
ggplot(modulos2, aes(fill=Plant, y=Amount, x=as.factor(module))) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw()+
  scale_fill_brewer(palette = 'Paired')+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  labs(x = "Module", y = "Number of focal plants",fill=NULL)+
  theme(legend.position="bottom")
  

# Modules pollinators----
fitness_data2 <- read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")
fitness2 <- fitness_data2 %>% filter(Year==2019) %>% rename(species=ID_Simple) %>%
  select(Plot,species,G_F) %>% unique()

modules_final$species[modules_final$species=="Lassioglosum_.immunitum"] <- "Lassioglosum_ immunitum"
modules_pollinators <- modules_final %>% filter(type=="pollinator") %>% 
  left_join(fitness2,by=c("Plot","species")) %>% #rename(Plot=Plot.x) %>%
  group_by(module,Plot,G_F) %>% count()

ggplot(modules_pollinators, aes(fill=G_F, y=n, x=as.factor(module))) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw()+
  scale_fill_brewer(palette = 'Paired')+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  labs(x = "Module", y = "Number of pollinators",fill=NULL)+theme(legend.position="bottom")

#################################
# TERNARY DIAGRAMS
# #################################
# library(ggtern)
# library(ggalt) 
# 
# shannon %>% filter(abs_PUPA==0,abs_CHMI==0)
# 
# base = ggtern(shannon %>% filter(abs_PUPA==0,abs_CHMI==0,entropy!=0),aes(abs_LEMA,abs_ME,abs_CHFU,size=total,color=as.factor(Plot))) + 
#   theme_bw() + 
#   #theme_legend_position('bottommiddle') + 
#   #geom_encircle(alpha=0.5,size=1) + 
#   geom_point(alpha=0.4) +
#   labs(title    = "Example Plot")#+
#   #facet_wrap(vars(Plot),nrow = 3,ncol = 3)
# print(base)


#################################
# CENTRALITY VERSUS PERCENTAGE
#################################

# Incluyendo los outlayer de los bigotes de gato de porcentaje de individuos iguales

PageRank_results %>% group_by(Plot) %>% summarise(Real_PR_Multi=sum(Real_PR_Multi))

unique(modules_plants2)

new_PR_mod_info <- modules_plants2 %>% rename(Plant_Simple=layer_name,Label=species) %>%
  left_join(PageRank_results,by=c("Plot", "Label","Plant_Simple"))

new_PR_mod_info %>% group_by(Plot) %>% summarise(Real_PR_Multi=sum(Real_PR_Multi))


  ggplot(new_PR_mod_info)+
  geom_point(aes(x=pecentage_same_plants,y=log10(Ratio),color=Plant_Simple,
                 shape=Plant_Simple),
             position = "jitter")+
  geom_smooth(aes(x=pecentage_same_plants,y=log10(Ratio),color=Plant_Simple),method="lm", se=F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  #ggtitle(paste0("Plot ",i)) +
  xlab("% same plant individuals per module") + ylab("log10(Ratio)") +
  labs(color = "Plant Species", shape = "Plant Species")+ theme(legend.position="bottom")
  
  ggplot(new_PR_mod_info)+
    geom_point(aes(x=pecentage_same_plants,y=log10(Ratio),color=Plant_Simple,
                   shape=Plant_Simple),
               position = "jitter")+
    geom_smooth(aes(x=pecentage_same_plants,y=log10(Ratio),color=Plant_Simple),method="lm", se=F)+
    xlab("% same plant individuals per module") + ylab("log10(Ratio)") +
    labs(color = "Plant Species", shape = "Plant Species")+ theme(legend.position="bottom")

ggplot(new_PR_mod_info)+
  geom_point(aes(x=100*pecentage_same_plants,y=log10(Real_PR_Multi),color=Plant_Simple,
                 shape=Plant_Simple),
             position = "jitter")+
  geom_smooth(aes(x=100*pecentage_same_plants,y=log10(Real_PR_Multi),color=Plant_Simple),method="lm", se=F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  #ggtitle(paste0("Plot ",i)) +
  xlab("% same plant individuals per module") + ylab("log10(PageRank)") +
  labs(color = "Plant Species", shape = "Plant Species")+ theme(legend.position="bottom")

ggplot(new_PR_mod_info)+
  geom_point(aes(x=100*pecentage_same_plants,y=log10(Real_PR_Multi),color=Plant_Simple,
                 shape=Plant_Simple),
             position = "jitter")+
  geom_smooth(aes(x=100*pecentage_same_plants,y=log10(Real_PR_Multi),color=Plant_Simple),
              method="lm", formula=y~x,se=F)+
  xlab("% same plant individuals per module") + ylab("log10(PageRank)") +
  labs(color = "Plant Species", shape = "Plant Species")+ theme(legend.position="bottom")


####
new_PR_mod_info %>% group_by(Plot) %>% summarise(Real_PR_Multi=sum(Real_PR_Multi))

new_PR_mod_info$module <- as.factor(new_PR_mod_info$module)
means_mod <- aggregate(Real_PR_Multi ~  module + Plot,
                       new_PR_mod_info, mean)

ggplot(new_PR_mod_info,aes(x=as.factor(module),y=Real_PR_Multi))+
  geom_boxplot()+
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=3, show.legend = FALSE) + 
  geom_text(data = means_mod, aes(label = round(Real_PR_Multi,3), y = 0.03))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  xlab("Module") + ylab("log10(Real_PR_Multi)")+ theme(legend.position = "none")


mean_mod_info <- new_PR_mod_info %>% dplyr::select(-Plant_Simple,-Label) %>% 
  group_by(Plot, module) %>% summarise_all(mean)

ggplot(mean_mod_info)+
  geom_point(aes(x=log10(total),y=log10(Real_PR_Multi)),size=1.5,alpha=0.3,position = "jitter")+
  geom_smooth(aes(x=log10(total),y=log10(Real_PR_Multi)),method="lm", se=F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3, labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("log10(Plant individuals per module)") + ylab("log10(Average[PageRank])") +
  labs(color = "Plant Species", shape = "Plant Species")+ theme(legend.position="bottom")+
  theme_bw()

library(lme4)
model <- lm(log10(Real_PR_Multi) ~ log10(total)+as.factor(Plot), data = mean_mod_info)
summary(model)

sum_mod_info <- new_PR_mod_info %>% select(-Plant_Simple,-Label) %>% 
  group_by(Plot, module) %>% summarise(Real_PR_Multi=sum(Real_PR_Multi),total=mean(total))

new_PR_mod_info %>%
  group_by(Plot) %>% summarise(Real_PR_Multi=sum(Real_PR_Multi))

ggplot(sum_mod_info)+
  geom_point(aes(x=log10(total),y=log10(Real_PR_Multi)),position = "jitter")+
  geom_smooth(aes(x=log10(total),y=log10(Real_PR_Multi)),method="lm", se=F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  #ggtitle(paste0("Plot ",i)) +
  xlab("log10(Plant individuals per module)") + ylab("log10(Sum[PageRank])") +
  labs(color = "Plant Species", shape = "Plant Species")+ theme(legend.position="bottom")


# Extract slopes and interceps of this regression lines
sum_mod_info %>% 
  group_by(Plot) %>% 
  do({
    mod = lm(log10(Real_PR_Multi) ~ log10(total), data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })

###############################################
#TABLE MODULES
library(forcats)
library(ggridges)


modules_polli <- modules_final %>% filter(type!="plant") %>% 
  select(Plot,module,species) 
modules_polli[modules_polli %>% duplicated(),] 

modules_final %>% filter(type!="plant") %>% 
  select(Plot,species) %>% unique() %>% group_by(Plot) %>% count()

# Different Pollinators per plot
modules_final %>% filter(type!="plant") %>% 
  select(Plot,species) %>% unique() %>% group_by(Plot) %>% count() 

# Number of pollinator per module and plot
modules_final %>% filter(type!="plant") %>% 
  select(Plot,module,species) %>% unique() %>% group_by(Plot,module) %>% count() %>%
  ungroup() %>% group_by(Plot) %>% summarise(mean_number_pollinator_per_mod=mean(n),
                                             sd_number_pollinator_per_mod=sd(n))

# Number of pollinator per module
modules_final %>% filter(type!="plant") %>% 
  select(Plot,module,species) %>% unique() %>% group_by(Plot,module) %>% count() %>%
  ungroup() %>% summarise(mean_number_pollinator_per_mod=mean(n),
                                             sd_number_pollinator_per_mod=sd(n))

group_pollinators <- modules_final %>% filter(type!="plant") %>% 
  select(Plot,module,species) %>% unique() %>% arrange(Plot, module)
  

table_mod1 <- modules_plants2 %>% group_by(Plot) %>% count() %>% rename(number_plant_individuals=n)


table_mod2 <- modules_plants2 %>% select(Plot, module) %>% unique() %>% 
  group_by(Plot) %>% count() %>% rename(number_modules=n)

table_mod3 <- plant_m %>% group_by(Plot) %>% summarize(mean_number_plant_species_per_module=mean(n),
                                         sd_number_plant_species_per_module=sd(n))

library(boot)
# Useful info: https://blog.methodsconsultants.com/posts/understanding-bootstrap-confidence-interval-output-from-the-r-boot-package/
# mean_number_plant_species_per_module
i=1
number_plant_species_per_module_i <- plant_m %>% ungroup() %>% filter(Plot==i)
b <- boot(number_plant_species_per_module_i$n, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(number_plant_species_per_module_i$n, probs = c(0.025, 0.975), na.rm = FALSE)
i=2
number_plant_species_per_module_i <- plant_m %>% ungroup() %>% filter(Plot==i)
b <- boot(number_plant_species_per_module_i$n, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(number_plant_species_per_module_i$n, probs = c(0.025, 0.975), na.rm = FALSE)
i=3
number_plant_species_per_module_i <- plant_m %>% ungroup() %>% filter(Plot==i)
b <- boot(number_plant_species_per_module_i$n, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(number_plant_species_per_module_i$n, probs = c(0.025, 0.975), na.rm = FALSE)
i=4
number_plant_species_per_module_i <- plant_m %>% ungroup() %>% filter(Plot==i)
b <- boot(number_plant_species_per_module_i$n, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(number_plant_species_per_module_i$n, probs = c(0.025, 0.975), na.rm = FALSE)
i=5
number_plant_species_per_module_i <- plant_m %>% ungroup() %>% filter(Plot==i)
b <- boot(number_plant_species_per_module_i$n, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(number_plant_species_per_module_i$n, probs = c(0.025, 0.975), na.rm = FALSE)
i=6
number_plant_species_per_module_i <- plant_m %>% ungroup() %>% filter(Plot==i)
b <- boot(number_plant_species_per_module_i$n, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(number_plant_species_per_module_i$n, probs = c(0.025, 0.975), na.rm = FALSE)

i=7
number_plant_species_per_module_i <- plant_m %>% ungroup() %>% filter(Plot==i)
b <- boot(number_plant_species_per_module_i$n, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(number_plant_species_per_module_i$n, probs = c(0.025, 0.975), na.rm = FALSE)

i=8
number_plant_species_per_module_i <- plant_m %>% ungroup() %>% filter(Plot==i)
b <- boot(number_plant_species_per_module_i$n, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(number_plant_species_per_module_i$n, probs = c(0.025, 0.975), na.rm = FALSE)
i=9
number_plant_species_per_module_i <- plant_m %>% ungroup() %>% filter(Plot==i)
b <- boot(number_plant_species_per_module_i$n, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(number_plant_species_per_module_i$n, probs = c(0.025, 0.975), na.rm = FALSE)

table_mod4 <- shannon %>% group_by(Plot) %>% summarize(mean_total=mean(total),
                                         sd_total=sd(total),
                                         mean_entropy=mean(entropy),
                                         sd_entropy=sd(entropy),
                                         mean_CHFU=mean(CHFU),
                                         sd_CHFU=sd(CHFU),
                                         mean_CHMI=mean(CHMI),
                                         sd_CHMI=sd(CHMI),
                                         mean_LEMA=mean(LEMA),
                                         sd_LEMA=sd(LEMA),
                                         mean_ME=mean(ME),
                                         sd_ME=sd(ME),
                                         mean_PUPA=mean(PUPA),
                                         sd_PUPA=sd(PUPA),
)


i=1
entropy_i <- shannon %>% ungroup() %>% filter(Plot==i)
b <- boot(entropy_i$entropy, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

quantile(entropy_i$entropy, probs = c(0.025, 0.975), na.rm = FALSE)

i=2
entropy_i <- shannon %>% ungroup() %>% filter(Plot==i)
b <- boot(entropy_i$entropy, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(entropy_i$entropy, probs = c(0.025, 0.975), na.rm = FALSE)

i=3
entropy_i <- shannon %>% ungroup() %>% filter(Plot==i)
b <- boot(entropy_i$entropy, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(entropy_i$entropy, probs = c(0.025, 0.975), na.rm = FALSE)

i=4
entropy_i <- shannon %>% ungroup() %>% filter(Plot==i)
b <- boot(entropy_i$entropy, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(entropy_i$entropy, probs = c(0.025, 0.975), na.rm = FALSE)

i=5
entropy_i <- shannon %>% ungroup() %>% filter(Plot==i)
b <- boot(entropy_i$entropy, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(entropy_i$entropy, probs = c(0.025, 0.975), na.rm = FALSE)

i=6
entropy_i <- shannon %>% ungroup() %>% filter(Plot==i)
b <- boot(entropy_i$entropy, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(entropy_i$entropy, probs = c(0.025, 0.975), na.rm = FALSE)

i=7
entropy_i <- shannon %>% ungroup() %>% filter(Plot==i)
b <- boot(entropy_i$entropy, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(entropy_i$entropy, probs = c(0.025, 0.975), na.rm = FALSE)
i=8
entropy_i <- shannon %>% ungroup() %>% filter(Plot==i)
b <- boot(entropy_i$entropy, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(entropy_i$entropy, probs = c(0.025, 0.975), na.rm = FALSE)

i=9
entropy_i <- shannon %>% ungroup() %>% filter(Plot==i)
b <- boot(entropy_i$entropy, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
quantile(entropy_i$entropy, probs = c(0.025, 0.975), na.rm = FALSE)


table_mod5 <- mean_mod_info %>% group_by(Plot) %>% summarize(mean_mean_PR=mean(Real_PR_Multi),
                                               sd_mean_PR=sd(Real_PR_Multi))
table_mod6 <- sum_mod_info %>% group_by(Plot) %>% summarize(mean_sum_PR=mean(Real_PR_Multi),
                                              sd_sum_PR=sd(Real_PR_Multi))


table_mod <- table_mod1 %>% left_join(table_mod2,by="Plot") %>%
  left_join(table_mod3,by="Plot") %>%
  left_join(table_mod4,by="Plot") %>%
  left_join(table_mod5,by="Plot") %>%
  left_join(table_mod6,by="Plot")

write.csv(table_mod,"resume_table_data_modularity.csv")

###################################################

#####################################################
# NEW RESUME
###############################################
#TABLE MODULES
library(forcats)
library(ggridges)


table_mod1_2 <- modules_final %>% filter(type=="plant") %>% 
  select(Plot,module,layer_name) %>% unique() %>% group_by(Plot,module) %>% count() %>%
  rename(number_plants_species_per_mod=n)

table_mod2_2 <- modules_final %>% filter(type=="plant") %>% 
  select(Plot,module,species) %>% unique() %>% group_by(Plot,module) %>% count() %>%
  rename(number_focal_plants_per_mod=n)

i=1
number_focal_plants_per_mod_i <- table_mod2_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_focal_plants_per_mod_i$number_focal_plants_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=2
number_focal_plants_per_mod_i <- table_mod2_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_focal_plants_per_mod_i$number_focal_plants_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=3
number_focal_plants_per_mod_i <- table_mod2_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_focal_plants_per_mod_i$number_focal_plants_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=4
number_focal_plants_per_mod_i <- table_mod2_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_focal_plants_per_mod_i$number_focal_plants_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=5
number_focal_plants_per_mod_i <- table_mod2_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_focal_plants_per_mod_i$number_focal_plants_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=6
number_focal_plants_per_mod_i <- table_mod2_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_focal_plants_per_mod_i$number_focal_plants_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=7
number_focal_plants_per_mod_i <- table_mod2_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_focal_plants_per_mod_i$number_focal_plants_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=8
number_focal_plants_per_mod_i <- table_mod2_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_focal_plants_per_mod_i$number_focal_plants_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=9
number_focal_plants_per_mod_i <- table_mod2_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_focal_plants_per_mod_i$number_focal_plants_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)


table_mod3_2 <- modules_final %>% filter(type!="plant") %>% 
  select(Plot,module,species) %>% unique() %>% group_by(Plot,module) %>% count() %>%
  rename(number_pollinator_per_mod=n)

i=1
number_pollinator_per_mod_i <- table_mod3_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_pollinator_per_mod_i$number_pollinator_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=2
number_pollinator_per_mod_i <- table_mod3_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_pollinator_per_mod_i$number_pollinator_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=3
number_pollinator_per_mod_i <- table_mod3_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_pollinator_per_mod_i$number_pollinator_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=4
number_pollinator_per_mod_i <- table_mod3_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_pollinator_per_mod_i$number_pollinator_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=5
number_pollinator_per_mod_i <- table_mod3_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_pollinator_per_mod_i$number_pollinator_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=6
number_pollinator_per_mod_i <- table_mod3_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_pollinator_per_mod_i$number_pollinator_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=7
number_pollinator_per_mod_i <- table_mod3_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_pollinator_per_mod_i$number_pollinator_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=8
number_pollinator_per_mod_i <- table_mod3_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_pollinator_per_mod_i$number_pollinator_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=9
number_pollinator_per_mod_i <- table_mod3_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(number_pollinator_per_mod_i$number_pollinator_per_mod, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)


table_mod4_2 <- shannon %>% select(Plot,module,entropy)
table_mod5_2 <- mean_mod_info %>% select(Plot,module,Real_PR_Multi) %>% rename(av_PR_per_plant=Real_PR_Multi)

i=1
av_PR_per_plant_i <- table_mod5_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(av_PR_per_plant_i$av_PR_per_plant, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=2
av_PR_per_plant_i <- table_mod5_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(av_PR_per_plant_i$av_PR_per_plant, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=3
av_PR_per_plant_i <- table_mod5_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(av_PR_per_plant_i$av_PR_per_plant, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=4
av_PR_per_plant_i <- table_mod5_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(av_PR_per_plant_i$av_PR_per_plant, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=5
av_PR_per_plant_i <- table_mod5_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(av_PR_per_plant_i$av_PR_per_plant, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=6
av_PR_per_plant_i <- table_mod5_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(av_PR_per_plant_i$av_PR_per_plant, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=7
av_PR_per_plant_i <- table_mod5_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(av_PR_per_plant_i$av_PR_per_plant, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=8
av_PR_per_plant_i <- table_mod5_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(av_PR_per_plant_i$av_PR_per_plant, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=9
av_PR_per_plant_i <- table_mod5_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(av_PR_per_plant_i$av_PR_per_plant, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)



table_mod6_2 <- sum_mod_info %>% select(Plot,module,Real_PR_Multi) %>% rename(sum_PR_per_module=Real_PR_Multi)
table_mod5_2$module <- as.numeric(table_mod5_2$module)
table_mod6_2$module <- as.numeric(table_mod6_2$module)

i=1
sum_PR_per_module_i <- table_mod6_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(sum_PR_per_module_i$sum_PR_per_module, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=2
sum_PR_per_module_i <- table_mod6_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(sum_PR_per_module_i$sum_PR_per_module, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=3
sum_PR_per_module_i <- table_mod6_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(sum_PR_per_module_i$sum_PR_per_module, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=4
sum_PR_per_module_i <- table_mod6_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(sum_PR_per_module_i$sum_PR_per_module, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=5
sum_PR_per_module_i <- table_mod6_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(sum_PR_per_module_i$sum_PR_per_module, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=6
sum_PR_per_module_i <- table_mod6_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(sum_PR_per_module_i$sum_PR_per_module, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=7
sum_PR_per_module_i <- table_mod6_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(sum_PR_per_module_i$sum_PR_per_module, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=8
sum_PR_per_module_i <- table_mod6_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(sum_PR_per_module_i$sum_PR_per_module, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)
i=9
sum_PR_per_module_i <- table_mod6_2 %>% ungroup() %>% filter(Plot==i)
b <- boot(sum_PR_per_module_i$sum_PR_per_module, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)


table_mod_2 <- table_mod1_2 %>% left_join(table_mod2_2,by=c("Plot","module")) %>%
  left_join(table_mod3_2,by=c("Plot","module")) %>%
  left_join(table_mod4_2,by=c("Plot","module")) %>%
  left_join(table_mod5_2,by=c("Plot","module")) %>%
  left_join(table_mod6_2,by=c("Plot","module"))

ggplot(data=table_mod_2)+
  geom_point(aes(x=entropy,y=av_PR_per_plant),position = "jitter")+theme_bw()+
  geom_smooth(aes(x=entropy,y=av_PR_per_plant),method = "lm")+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  ylab(NULL)+xlab("Entropy")

library(gridExtra)
library(gtable)
library(grid)

names(table_mod_2)

plot1 <- ggplot(data=table_mod_2)+
  geom_histogram(aes(x=number_plants_species_per_mod))+theme_bw()+
  facet_wrap(vars(Plot),nrow = 1,ncol = 9,labeller=labeller(Plot= plot_labs))+
  ylab(NULL)+xlab("# Plant Species (Layers)")
plot2 <- ggplot(data=table_mod_2)+
  geom_histogram(aes(x=number_focal_plants_per_mod))+theme_bw()+
  facet_wrap(vars(Plot),nrow = 1,ncol = 9,labeller=labeller(Plot= plot_labs))+
  ylab(NULL)+xlab("# Focal plants per module")
plot3 <- ggplot(data=table_mod_2)+
  geom_histogram(aes(x=number_pollinator_per_mod))+theme_bw()+
  facet_wrap(vars(Plot),nrow = 1,ncol = 9,labeller=labeller(Plot= plot_labs))+
  ylab(NULL)+xlab("# Animal visitors per module")
plot4 <- ggplot(data=table_mod_2)+
  geom_histogram(aes(x=entropy))+theme_bw()+
  facet_wrap(vars(Plot),nrow = 1,ncol = 9,labeller=labeller(Plot= plot_labs))+
  ylab(NULL)+xlab("Shannon diversity")
plot5 <- ggplot(data=table_mod_2)+
  geom_histogram(aes(x=av_PR_per_plant))+theme_bw()+
  facet_wrap(vars(Plot),nrow = 1,ncol = 9,labeller=labeller(Plot= plot_labs))+
  ylab(NULL)+xlab("Average PageRank per focal plant")
plot6 <- ggplot(data=table_mod_2)+
  geom_histogram(aes(x=sum_PR_per_module))+theme_bw()+
  facet_wrap(vars(Plot),nrow = 1,ncol = 9,labeller=labeller(Plot= plot_labs))+
  ylab(NULL)+xlab("Sum PageRank per module")

grid.arrange(
  grobs = list(plot1, plot2, plot3,plot4, plot5,plot6),
  widths = c(1),
  layout_matrix = rbind(c(1),
                        c(2),
                        c(3),
                        c(4),
                        c(5),
                        c(6)),
  left = textGrob("Counts", rot = 90, vjust = 1)
)


plot12 <- ggplot(data=table_mod_2)+
  geom_histogram(aes(x=number_plants_species_per_mod))+theme_bw()+
  ylab(NULL)+xlab("# Plant Species (Layers)")
plot22 <- ggplot(data=table_mod_2)+
  geom_histogram(aes(x=number_focal_plants_per_mod))+theme_bw()+
  ylab(NULL)+xlab("# Focal plants per module")
plot32 <- ggplot(data=table_mod_2)+
  geom_histogram(aes(x=number_pollinator_per_mod))+theme_bw()+
  ylab(NULL)+xlab("# Animal visitors per module")
plot42 <- ggplot(data=table_mod_2)+
  geom_histogram(aes(x=entropy))+theme_bw()+
  ylab(NULL)+xlab("Shannon diversity")
plot52 <- ggplot(data=table_mod_2)+
  geom_histogram(aes(x=av_PR_per_plant))+theme_bw()+
  ylab(NULL)+xlab("Average PageRank per focal plant")
plot62 <- ggplot(data=table_mod_2)+
  geom_histogram(aes(x=sum_PR_per_module))+theme_bw()+
  ylab(NULL)+xlab("Sum PageRank per module")

grid.arrange(
  grobs = list(plot12, plot22, plot32,plot42, plot52,plot62),
  widths = c(1,1,1),
  layout_matrix = rbind(c(1,2,3),
                        c(4,5,6)),
  left = textGrob("Counts", rot = 90, vjust = 1)
)



################################################
# NEW (MODULE) FLOW PREDICTORS
#################################################

PageRank_results2 <- PageRank_results %>% left_join(rename(modules_plants2,Label=species), 
                               by=c("Plot","Label")) %>% select(-total,-layer_name) %>%
  left_join(shannon,by=c("Plot","module"))

write_csv(PageRank_results2,"PageRank_results.csv")

##################################################
####################################################
######
ggplot(new_PR_mod_info)+
  geom_point(aes(x=total,y=log10(Real_PR_Multi),color=Plant_Simple,
                 shape=Plant_Simple),
             position = "jitter")+
  geom_smooth(aes(x=total,y=log10(Real_PR_Multi),color=Plant_Simple),method="lm", se=F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Plant individuals per module") + ylab("log10(PageRank)") +
  labs(color = "Plant Species", shape = "Plant Species")+ theme(legend.position="bottom")

ggplot(new_PR_mod_info)+
  geom_point(aes(x=total,y=log10(Real_PR_Multi),color=Plant_Simple,
                 shape=Plant_Simple),
             position = "jitter")+
  geom_smooth(aes(x=total,y=log10(Real_PR_Multi),color=Plant_Simple),
              method="lm", formula=y~x,se=F)+
  xlab("Plant individuals per module") + ylab("log10(PageRank)") +
  labs(color = "Plant Species", shape = "Plant Species")+ theme(legend.position="bottom")


# Excluyendo los outlayer de los bigotes de gato de porcentaje de individuos iguales


new_PR_mod_info_outliers <- new_PR_mod_info %>%
  group_by(Plant_Simple,Plot) %>%
  summarise(Min = min(pecentage_same_plants),
            Max = max(pecentage_same_plants),
            Median = median(pecentage_same_plants),
            IQRange = IQR(pecentage_same_plants),
            Q1=quantile(pecentage_same_plants,0.25),
            Q3=quantile(pecentage_same_plants,0.25)
            ) %>%
  mutate(Lower_Limit=Q1-1.5*IQRange,
         Upper_Limit=Q3+1.5*IQRange)

new_PR_mod_info_outliers2 <-new_PR_mod_info_outliers %>% 
  left_join(new_PR_mod_info,by=c("Plant_Simple","Plot")) %>%
  filter(pecentage_same_plants>=Lower_Limit,pecentage_same_plants<=Upper_Limit)

ggplot(new_PR_mod_info_outliers2)+
  geom_point(aes(x=100*pecentage_same_plants,y=log10(Real_PR_Multi),color=Plant_Simple,
                 shape=Plant_Simple),
             position = "jitter")+
  geom_smooth(aes(x=100*pecentage_same_plants,y=log10(Real_PR_Multi),color=Plant_Simple),method="lm", se=F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  #ggtitle(paste0("Plot ",i)) +
  xlab("% same plant individuals per module") + ylab("log10(PageRank)") +
  labs(color = "Plant Species", shape = "Plant Species")+ theme(legend.position="bottom")

####################################################
# COMPARISON BETWEEN % SAME INDIVIDUALS AND FITNESS
#######################################################
fitness_final_aux0 <- read.csv(file = "data_models_phenol_overlap.csv",
                               header = TRUE,
                               stringsAsFactors = FALSE)
new_PR_mod_info1 <- new_PR_mod_info %>% separate(Label,c("Subplot","Plant0")," ") %>% 
  dplyr::select(-Plant0) 


new_PR_mod_info2 <- new_PR_mod_info1 %>% 
  left_join(fitness_final_aux0,by=c("Plot","Subplot","Plant_Simple"))


ggplot(new_PR_mod_info2%>%filter(Seeds_GF>0))+
  geom_point(aes(x=100*pecentage_same_plants,y=log10(Seeds_GF),color=Plant_Simple,
                 shape=Plant_Simple),
             position = "jitter")+
  geom_smooth(aes(x=100*pecentage_same_plants,y=log10(Seeds_GF),color=Plant_Simple),method="lm", se=F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  #ggtitle(paste0("Plot ",i)) +
  xlab("% same plant individuals per module") + ylab("log10(Seeds_GF)") +
  labs(color = "Plant Species", shape = "Plant Species")+ theme(legend.position="bottom")

ggplot(new_PR_mod_info2%>%filter(Seeds_GF>0))+
  geom_point(aes(x=log10(Real_PR_Multi),y=log10(Seeds_GF),color=Plant_Simple,
                 shape=Plant_Simple),
             position = "jitter")+
  geom_smooth(aes(x=log10(Real_PR_Multi),y=log10(Seeds_GF),color=Plant_Simple),method="lm", se=F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  #ggtitle(paste0("Plot ",i)) +
  xlab("log10(PageRank_Multi)") + ylab("log10(Seeds_GF)") +
  labs(color = "Plant Species", shape = "Plant Species")+ theme(legend.position="bottom")

ggplot(new_PR_mod_info2%>%filter(Seeds_GF>0))+
  geom_point(aes(x=log10(Real_PR_Layer),y=log10(Seeds_GF),color=Plant_Simple,
                 shape=Plant_Simple),
             position = "jitter")+
  geom_smooth(aes(x=log10(Real_PR_Layer),y=log10(Seeds_GF),color=Plant_Simple),method="lm", se=F)+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  #ggtitle(paste0("Plot ",i)) +
  xlab("log10(PageRank_Layer)") + ylab("log10(Seeds_GF)") +
  labs(color = "Plant Species", shape = "Plant Species")+ theme(legend.position="bottom")
