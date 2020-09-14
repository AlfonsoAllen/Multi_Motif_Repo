library(tidyverse)
library(boot)

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


###########################
# MODULAR STRUCTURE
###########################

for (Plot_i in 1:9){
  
  module_i <- read_csv(paste0("Processed_data/Modularity_Pheno_Overlap/NN_Modularity_Plot",Plot_i,".csv")  )
  
  if (Plot_i==1){
    modules_final=module_i
  }else{
    modules_final=bind_rows(modules_final,module_i)
  }
  
}

modules_final %>% filter(type=="plant") %>%  group_by(layer_name) %>% count()

modules_final %>% filter(type!="plant") %>%  group_by(layer_name) %>% count()

################################
# PROCESSING DATA
################################

# number of modules per plot

modules_final %>% dplyr::select(Plot,module) %>% unique() %>% 
  group_by(Plot) %>% count()


# Number of plants focal plants (with visitors) per plot

modules_final %>% filter(type=="plant") %>% dplyr::select(Plot,module,layer_name) %>% 
  group_by(Plot) %>% count()


# Number of pollinator species per plot 
# (we make no differences between Pol1_LEMA and Plo1_CHFU)

modules_final %>% filter(type!="plant") %>% dplyr::select(Plot,species) %>%
  unique() %>%
  group_by(Plot) %>% count()


# Number of plant species per module

plant_sp_m <- modules_final %>% filter(type=="plant") %>% dplyr::select(Plot,module,layer_name) %>% 
  unique() %>% arrange(Plot,module,layer_name)%>% group_by(Plot,module) %>% count()

plant_sp_m %>%  ggplot()+
  geom_histogram(aes(x=n))+theme_bw()+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("# Plant species per module") + ylab("# Modules") + theme(legend.position = "none")
#Recoge a través de cuantas capas se extiende el módulo

sum(plant_sp_m$n==1)
sum(plant_sp_m$n==2)
sum(plant_sp_m$n==3)
sum(plant_sp_m$n==4)

mean(plant_sp_m$n)
sd(plant_sp_m$n)

b <- boot(plant_sp_m$n, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)


# Number of poll species per module

poll_sp_m <- modules_final %>%  filter(type !="plant") %>% dplyr::select(Plot,module,species) %>% 
  unique() %>% arrange(Plot,module,species) %>% group_by(Plot,module) %>% count()

poll_sp_m %>%  ggplot()+
  geom_histogram(aes(x=n))+theme_bw()+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("# Animal species and morphospecies per module") + ylab("# Modules") + theme(legend.position = "none")

sum(poll_sp_m$n==0)
sum(poll_sp_m$n==1)
sum(poll_sp_m$n==2)
sum(poll_sp_m$n==3)
sum(poll_sp_m$n==4)
sum(poll_sp_m$n==5)
sum(poll_sp_m$n==6)
sum(poll_sp_m$n==9)

mean(poll_sp_m$n)
sd(poll_sp_m$n)

b <- boot(poll_sp_m$n, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)


# Number of plant ind per module

plant_ind_m <- modules_final %>% filter(type=="plant") %>% dplyr::select(Plot,module,layer_name) %>% 
  arrange(Plot,module,layer_name)%>% group_by(Plot,module) %>% count()

plant_ind_m %>%  ggplot()+
  geom_histogram(aes(x=n))+theme_bw()+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("# Plant individuals per module") + ylab("# Modules") + theme(legend.position = "none")
#Recoge a través de cuantas capas se extiende el módulo

mean(plant_ind_m$n)
sd(plant_ind_m$n)

b <- boot(plant_ind_m$n, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)


# Number of pollinator per module

poll_ind_m <- modules_final %>%  filter(type !="plant") %>% dplyr::select(Plot,module,species) %>%
  arrange(Plot,module,species)%>% group_by(Plot,module) %>% count()

poll_ind_m %>%  ggplot()+
  geom_histogram(aes(x=n))+theme_bw()+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("# Animal individuals per module") + ylab("# Modules") + theme(legend.position = "none")

mean(poll_ind_m$n)

b <- boot(poll_ind_m$n, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)


# Number of individuals per module

ind_m <- modules_final %>% dplyr::select(Plot,module) %>%
  group_by(Plot,module) %>% count()

ind_m %>%  ggplot()+
  geom_histogram(aes(x=n))+theme_bw()+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("# Individuals per module") + ylab("# Modules") + theme(legend.position = "none")

sum(ind_m$n==2)
sum(ind_m$n==3)
sum(ind_m$n==4)
sum(ind_m$n==5)
sum(ind_m$n==6)
sum(ind_m$n==29)


mean(ind_m$n)
sd(ind_m$n)

b <- boot(ind_m$n, function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)


# Modules plants----

plant_ind_m_species <- modules_final %>% filter(type=="plant") %>% dplyr::select(Plot,module,layer_name) %>% 
  arrange(Plot,module,layer_name)%>% group_by(Plot,module,layer_name) %>% count()

plant_ind_m_species %>% group_by(layer_name) %>% count()

library(RColorBrewer)
ggplot(plant_ind_m_species, 
       aes(fill=as.factor(layer_name), y=n, x=as.factor(module))) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw()+
  scale_fill_brewer(palette = 'Paired')+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  labs(x = "Module", y = "Number of plant individuals",fill=NULL)+theme(legend.position="bottom")


# Modules pollinators----
fitness_data2 <- read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")
fitness2 <- fitness_data2 %>% filter(Year==2019) %>% rename(species=ID_Simple) %>%
  dplyr:: select(Plot,species,G_F) %>% unique()

modules_final$species[modules_final$species=="Lassioglosum_.immunitum"] <- "Lassioglosum_ immunitum"
modules_pollinators <- modules_final %>% filter(type=="pollinator") %>% 
  left_join(fitness2,by=c("Plot","species")) %>% #rename(Plot=Plot.x) %>%
  group_by(module,Plot,G_F) %>% count()


modules_pollinators_label <- modules_pollinators
modules_pollinators_label$G_F <- as.character(modules_pollinators_label$G_F)

modules_pollinators_label$G_F[modules_pollinators_label$G_F=="Flower_beetles"] <- "Flower beetles"
modules_pollinators_label$G_F[modules_pollinators_label$G_F=="House_flies"] <- "House flies"
modules_pollinators_label$G_F[modules_pollinators_label$G_F=="Small_beetles"] <- "Small beetles"
modules_pollinators_label$G_F[modules_pollinators_label$G_F=="Small_flies"] <- "Small flies"

modules_pollinators_label$G_F <- as.factor(modules_pollinators_label$G_F)

ggplot(modules_pollinators_label, aes(fill=G_F, y=n, x=as.factor(module))) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw()+
  scale_fill_brewer(palette = 'Paired')+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  labs(x = "Module", y = "Number of pollinators",fill=NULL)+theme(legend.position="bottom")


