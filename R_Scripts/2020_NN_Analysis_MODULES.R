
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

  module_i <- read_csv(paste0("Processed_data/Modularity_Pheno_Overlap/2020_NN_Modularity_Plot",Plot_i,".csv")  )
  
  if (Plot_i==1){
    modules_final=module_i
  }else{
    modules_final=bind_rows(modules_final,module_i)
  }

}

################################
# PROCESSING DATA
################################

# Number of plants per plot

modules_final %>% filter(type=="plant") %>% dplyr::select(Plot,module,layer_name) %>% 
  group_by(Plot) %>% count()


modules_final %>% filter(type!="plant") %>% dplyr::select(Plot,module,layer_name) %>% 
  group_by(Plot) %>% count()

# Number of plant species

plant_sp_m <- modules_final %>% filter(type=="plant") %>% dplyr::select(Plot,module,layer_name) %>% 
  unique() %>% arrange(Plot,module,layer_name)%>% group_by(Plot,module) %>% count()

plant_sp_m %>%  ggplot()+
  geom_histogram(aes(x=n))+theme_bw()+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("# Plant species per module") + ylab("# Modules") + theme(legend.position = "none")
#Recoge a través de cuantas capas se extiende el módulo

plant_sp_m %>% filter(Plot==8) %>%  ggplot()+
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

i=1
plant_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=2
plant_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=3
plant_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=4
plant_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=5
plant_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=6
plant_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=7
plant_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=8
plant_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=9
plant_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)



# Number of poll species

poll_sp_m <- modules_final %>%  filter(type !="plant") %>% dplyr::select(Plot,module,species) %>% 
  unique() %>% arrange(Plot,module,species)%>% group_by(Plot,module) %>% count()

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


i=1
poll_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(poll_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=2
poll_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(poll_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=3
poll_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(poll_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=4
poll_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(poll_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=5
poll_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(poll_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=6
poll_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(poll_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=7
poll_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(poll_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=8
poll_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(poll_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=9
poll_sp_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(poll_sp_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)


# Number of plant ind

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

i=1
plant_ind_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_ind_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=2
plant_ind_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_ind_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=3
plant_ind_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_ind_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=4
plant_ind_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_ind_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=5
plant_ind_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_ind_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=6
plant_ind_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_ind_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=7
plant_ind_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_ind_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=8
plant_ind_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_ind_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)

i=9
plant_ind_m %>% filter(Plot==i) %>% group_by(Plot) %>% summarise(n=mean(n))
b <- boot(plant_ind_m %>% ungroup() %>% filter(Plot==i) %>% select(n) %>% pull(), 
          function(u,i) u[i], R = 1000)
boot.ci(b, type = c("norm", "basic", "perc"),conf = .95)


# Graph combined for suppl. mat.

p1 <- plant_sp_m %>%  ggplot()+
  geom_histogram(aes(x=n))+theme_bw()+
  #facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("# Plant species per module") + ylab("# Modules") + theme(legend.position = "none")


p2 <- poll_sp_m %>%  ggplot()+
  geom_histogram(aes(x=n))+theme_bw()+
  #facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("# Pollinator species per module") + ylab("# Modules") + theme(legend.position = "none")

p3 <- plant_ind_m %>%  ggplot()+
  geom_histogram(aes(x=n))+theme_bw()+
  #facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  #ggtitle(paste0("Plot ",i)) +
  xlab("# Plant individuals per module") + ylab("# Modules") + theme(legend.position = "none")

library(ggpubr)
ggarrange(p1,p2,p3,
          #labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)

# Number of poll individuals per module

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

# Number of ind

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





# Number of poll




# Plant modules data

modules_final %>% filter(type=="plant") %>% dplyr::select(Plot,module,layer_name) %>% 
  unique() %>% arrange(Plot,module,layer_name) %>% group_by(Plot,module) %>% count() %>%
  group_by(Plot,n) %>% count()
################################################

modules <- modules_final %>% 
  dplyr::select(Plot,module) %>% unique()


number_modul <- modules %>% group_by(Plot) %>% count()
mean(number_modul$n)

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
                                                       abs_MESU=MESU,
                                                       abs_PUPA=PUPA,
                                                       abs_BEMA=BEMA,  
                                                       abs_CETE=CETE,
                                                       abs_SCLA=SCLA,  
                                                       abs_SOAS=SOAS,  
                                                       abs_SPRU=SPRU) %>%
  mutate(CHFU=abs_CHFU/total,
         CHMI=abs_CHMI/total,
         LEMA=abs_LEMA/total,
         MESU=abs_MESU/total,
         PUPA=abs_PUPA/total,
         BEMA=abs_BEMA/total,  
         CETE=abs_CETE/total,
         SCLA=abs_SCLA/total,  
         SOAS=abs_SOAS/total,  
         SPRU=abs_SPRU/total) %>%
  gather(layer_name,pecentage_same_plants,c(CHFU, CHMI, LEMA, MESU, PUPA,BEMA,
                                            CETE,SCLA,SOAS,SPRU)) %>%
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
                                            abs_MESU=MESU,
                                            abs_PUPA=PUPA,
                                            abs_BEMA=BEMA,  
                                            abs_CETE=CETE,
                                            abs_SCLA=SCLA,  
                                            abs_SOAS=SOAS,  
                                            abs_SPRU=SPRU) %>%
  mutate(CHFU=abs_CHFU/total,
         CHMI=abs_CHMI/total,
         LEMA=abs_LEMA/total,
         MESU=abs_MESU/total,
         PUPA=abs_PUPA/total,
         BEMA=abs_BEMA/total,  
         CETE=abs_CETE/total,
         SCLA=abs_SCLA/total,  
         SOAS=abs_SOAS/total,  
         SPRU=abs_SPRU/total) 


shannon$entropy <- NA

for (i in 1:nrow(shannon)){
  
  shannon$entropy[i] <- 0
  if (shannon$CHFU[i]>0){shannon$entropy[i] <-shannon$entropy[i] -shannon$CHFU[i]*log(shannon$CHFU[i]) }
  if (shannon$CHMI[i]>0){shannon$entropy[i] <-shannon$entropy[i] -shannon$CHMI[i]*log(shannon$CHMI[i]) }
  if (shannon$LEMA[i]>0){shannon$entropy[i] <-shannon$entropy[i] -shannon$LEMA[i]*log(shannon$LEMA[i]) }
  if (shannon$MESU[i]>0){shannon$entropy[i] <-shannon$entropy[i] -shannon$MESU[i]*log(shannon$MESU[i]) }
  if (shannon$PUPA[i]>0){shannon$entropy[i] <-shannon$entropy[i] -shannon$PUPA[i]*log(shannon$PUPA[i]) }
  if (shannon$BEMA[i]>0){shannon$entropy[i] <-shannon$entropy[i] -shannon$BEMA[i]*log(shannon$BEMA[i]) }
  if (shannon$CETE[i]>0){shannon$entropy[i] <-shannon$entropy[i] -shannon$CETE[i]*log(shannon$CETE[i]) }
  if (shannon$SCLA[i]>0){shannon$entropy[i] <-shannon$entropy[i] -shannon$SCLA[i]*log(shannon$SCLA[i]) }
  if (shannon$SOAS[i]>0){shannon$entropy[i] <-shannon$entropy[i] -shannon$SOAS[i]*log(shannon$SOAS[i]) }
  if (shannon$SCLA[i]>0){shannon$entropy[i] <-shannon$entropy[i] -shannon$SCLA[i]*log(shannon$SCLA[i]) }
  
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




modulos2 <- shannon %>% gather(Plant,Amount,c("abs_CHFU", "abs_CHMI", "abs_LEMA", "abs_MESU", "abs_PUPA",
                                              "abs_BEMA","abs_CETE","abs_SCLA","abs_SOAS","abs_SPRU"))
modulos2$Plant[modulos2$Plant=="abs_CHFU"] <- "CHFU"
modulos2$Plant[modulos2$Plant=="abs_CHMI"] <- "CHMI"
modulos2$Plant[modulos2$Plant=="abs_LEMA"] <- "LEMA"
modulos2$Plant[modulos2$Plant== "abs_MESU"] <-  "MESU"
modulos2$Plant[modulos2$Plant=="abs_PUPA"] <- "PUPA"

modulos2$Plant[modulos2$Plant=="abs_BEMA"] <- "BEMA"
modulos2$Plant[modulos2$Plant=="abs_CETE"] <- "CETE"
modulos2$Plant[modulos2$Plant=="abs_SCLA"] <- "SCLA" 
modulos2$Plant[modulos2$Plant=="abs_SOAS"] <- "SOAS"
modulos2$Plant[modulos2$Plant=="abs_SPRU"] <- "SPRU"
  

ggplot(modulos2, aes(fill=Plant, y=Amount, x=as.factor(module))) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw()+
  scale_fill_brewer(palette = 'Paired')+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  labs(x = "Module", y = "Number of focal plants",fill=NULL)+
  theme(legend.position="bottom")

# No codes
modulos2_real <- modulos2
modulos2_real$Plant[modulos2_real$Plant == "BEMA"] <- "B. macrocarpa"
modulos2_real$Plant[modulos2_real$Plant == "CETE"] <- "C. tenuiflorum"
modulos2_real$Plant[modulos2_real$Plant == "CHFU"] <- "C. fuscatum"
modulos2_real$Plant[modulos2_real$Plant == "CHMI"] <- "C. mixtum"
modulos2_real$Plant[modulos2_real$Plant == "LEMA"] <- "L. maroccanus"
modulos2_real$Plant[modulos2_real$Plant == "MESU"] <- "M. sulcatus"
modulos2_real$Plant[modulos2_real$Plant == "PUPA"] <- "P. paludosa"
modulos2_real$Plant[modulos2_real$Plant == "SCLA"] <- "S. laciniata"
modulos2_real$Plant[modulos2_real$Plant == "SOAS"] <- "S. asper"
modulos2_real$Plant[modulos2_real$Plant == "SPRU"] <- "S. rubra"


ggplot(modulos2_real, aes(fill=Plant, y=Amount, x=as.factor(module))) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw()+
  scale_fill_brewer(palette = 'Paired')+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  labs(x = "Module", y = "Number of focal plants",fill=NULL)+
  theme(legend.position="bottom")+
  theme(legend.text = element_text(face = "italic"))

# Save 700 x 500

# Percentage of modules with PUPA and LEMA
mod_PUPA_LEMA <- modulos2 %>% filter(LEMA>0 & PUPA >0.5 & Plant == "PUPA") %>% unique() %>% nrow()
mod_PUPA <- modulos2 %>% filter(PUPA >0 & Plant == "PUPA") %>% unique() %>% nrow()
100*mod_PUPA_LEMA/mod_PUPA

# Percentage of modules in which at least 50% is LEMA
mod_LEMA_dom <- modulos2 %>% filter(LEMA>=0.5 & Plant == "LEMA") %>% unique() %>% nrow()
mod_LEMA <- modulos2 %>% filter(LEMA >0 & Plant == "LEMA") %>% unique() %>% nrow()
100*mod_LEMA_dom/mod_LEMA

# Percentage of modules in which at least 50% is CHFU
mod_CHFU_dom <- modulos2 %>% filter(CHFU>=0.5 & Plant == "CHFU") %>% unique() %>% nrow()
mod_CHFU <- modulos2 %>% filter(CHFU >0 & Plant == "CHFU") %>% unique() %>% nrow()
100*mod_CHFU_dom/mod_CHFU

ggplot(modulos2_real %>% filter(Plot==8), aes(fill=Plant, y=Amount, x=as.factor(module))) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw()+
  scale_fill_brewer(palette = 'Paired')+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  labs(x = "Module", y = "Number of focal plants",fill=NULL)+
  theme(legend.position="bottom")
  

# Modules pollinators----

fitness_data2 <- read_csv2("Raw_Data/final_Pollinators_2020.csv")
# Inspect ID, ID_Simple

fitness_data2$ID %>% unique() %>% sort()
fitness_data2$ID_Simple %>% unique() %>% sort()

fitness_data2$ID[grep(" ",fitness_data2$ID)] # No labels contain spaces
fitness_data2$ID_Simple[grep(" ",fitness_data2$ID_Simple)]# No labels contain spaces

fitness_data2$ID[grep("\\.",fitness_data2$ID)] # Labels contain dots
fitness_data2$ID_Simple[grep("\\.",fitness_data2$ID_Simple)]# Labels contain dots

# Remove points from ID names
fitness_data2$ID <- sub("\\.", "", fitness_data2$ID)
fitness_data2$ID_Simple <- sub("\\.", "", fitness_data2$ID_Simple)

# filter tabanidae
fitness_data2 <- fitness_data2 %>% filter(ID != "Tabanidae")
fitness2 <- fitness_data2 %>% filter(!is.na(Plant),Plant!="0",Subplot!="OUT",Plant!="Ground")

fitness2 <- fitness2 %>% dplyr::select(-Order,-Family,-Genus,-ID) %>% rename(ID=ID_Simple) %>%
  mutate(date_raw=as.Date(paste(Day,Month,Year,sep="/"), "%d/%m/%Y"),
         Week=as.numeric(format(date_raw, "%V")))%>% rename(species=ID) %>%
  dplyr:: select(Plot,species,G_F) %>% unique()

modules_pollinators <- modules_final %>% filter(type=="pollinator") %>% 
  left_join(fitness2,by=c("Plot","species")) %>% #rename(Plot=Plot.x) %>%
  group_by(module,Plot,G_F) %>% count()


modules_pollinators_label <- modules_pollinators
modules_pollinators_label$G_F <- as.character(modules_pollinators_label$G_F)

modules_pollinators_label$G_F[modules_pollinators_label$G_F=="Flower_beetles"] <- "Flower beetles"
modules_pollinators_label$G_F[modules_pollinators_label$G_F=="House_flies"] <- "House flies"
modules_pollinators_label$G_F[modules_pollinators_label$G_F=="Small_beetles"] <- "Small beetles"
modules_pollinators_label$G_F[modules_pollinators_label$G_F=="Small_flies"] <- "Small flies"
modules_pollinators_label$G_F[modules_pollinators_label$G_F=="Solitary_bees"] <- "Solitary bees"
modules_pollinators_label$G_F[modules_pollinators_label$G_F=="Wasp"] <- "Wasps"
modules_pollinators_label$G_F[modules_pollinators_label$G_F=="Big_beetles"] <- "Big beetles"

modules_pollinators_label$G_F <- as.factor(modules_pollinators_label$G_F)

ggplot(modules_pollinators_label, aes(fill=G_F, y=n, x=as.factor(module))) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw()+
  scale_fill_brewer(palette = 'Paired')+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  labs(x = "Module", y = "Number of pollinators",fill=NULL)+theme(legend.position="bottom")

# Save 700 x 500
