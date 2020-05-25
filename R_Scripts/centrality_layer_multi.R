
library(tidyverse)

#################
# Centrality layer
#################

i <- 8
# Layer Centrality
file_CL_i <- paste0("Processed_data/Muxviz_Pheno_Overlap/centrality_perLayer_table_Plot",i)
centrality_layer <- read_delim(file_CL_i,";")

# Page rank results from Muxviz are not normalized
# We get the layer normalization constant = sum pargerank results
layer_norm_const <- centrality_layer %>% group_by(Layer) %>% 
  count(wt=PageRank) %>% rename(sum_weight=n)

centrality_layer2 <- centrality_layer %>% 
  left_join(layer_norm_const,by="Layer") %>% mutate(Real_PR=PageRank/sum_weight) %>%
  select(-Node,-Eigenvector,-Hub,-Authority,-Katz,-Multiplexity,-Kcore,-sum_weight,-PageRank)

#Sanity check
centrality_layer2 %>% group_by(Layer) %>% count(wt=Real_PR)

file_i <- paste0("Processed_data/Muxviz_Pheno_Overlap/Plot_",i,"/general_multilayer_layers_Plot",i,".txt")
layer_ID <- read.delim(file_i,header = T,sep =" ",stringsAsFactors =F) %>% 
  rename(Layer=layerID)
layer_ID$Layer <- as.character(layer_ID$Layer)
  
centrality_layer3 <- centrality_layer2 %>% left_join(layer_ID,by="Layer") %>%
  select(-Layer) %>% rename(Layer=layerLabel)
centrality_layer3$Layer[is.na(centrality_layer3$Layer)] <- "Aggr"

centrality_layer3 %>% group_by(Layer) %>% count(wt=Real_PR)

########################################
# Multilayer Centrality
########################################

for (i in 1:9){

file_CML_i <- paste0("Processed_data/Muxviz_Pheno_Overlap/centrality_table_Plot",i)
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
  left_join(Mlayer_norm_const,by="Layer") %>% mutate(Real_PR=PageRank/sum_weight) %>%
  select(-Node,-Eigenvector,-Hub,-Authority,-Katz,-Multiplexity,-Kcore,-sum_weight,-PageRank)

#Sanity check
centrality_multilayer_fil2 %>% group_by(Layer) %>% count(wt=Real_PR)
centrality_multilayer_fil2 %>% group_by(Label) %>% count() %>% filter(n>2)
centrality_multilayer_fil2 %>% group_by(Label) %>% count() %>% filter(n==2)
centrality_multilayer_fil2 %>% filter(is.na(Real_PR))
centrality_multilayer_fil2 %>% filter(Real_PR==0)


# Nodes with interlinks have different strengths. Consequently, spreading them generates NAs
Links <- centrality_multilayer_fil2 %>% select(Layer,Label,Strength)
Links <- spread(Links,Layer,Strength) %>% mutate(Delta_Strength=Multi-Aggr) %>% 
  select(Label,Aggr,Delta_Strength) %>% rename(Strength_Aggr=Aggr)



Centr_ML <- centrality_multilayer_fil2 %>% left_join(Links,by="Label") %>%
  select(Layer,Label,Strength_Aggr,Delta_Strength,Real_PR) %>%
  spread(Layer,Real_PR)
Centr_ML$ID <- Centr_ML$Label
Centr_ML_ID <- Centr_ML %>%separate(ID,c("Subplot","ID")," ") %>% select(-Subplot)
Centr_ML_ID$ID[is.na(Centr_ML_ID$ID)] <- "Visitor"
Centr_ML_ID$ID <- as.factor(Centr_ML_ID$ID)

Centr_ML_ID$Plot <- i
Centr_ML_ID$Plot <- as.factor(Centr_ML_ID$Plot)

if (i==1){
  Centr_ML_ID_final <- Centr_ML_ID
}else{
  Centr_ML_ID_final <- bind_rows(Centr_ML_ID_final,Centr_ML_ID)
}

}

ggplot(Centr_ML_ID_final)+
  geom_point(aes(x=log10(Multi),y=log10(Aggr),color=ID, size=(1+Delta_Strength),
                 shape=ID),
             position = "jitter")+
  geom_line(aes(x=log10(Aggr),y=log10(Aggr)),color='steelblue', 
            size=1, alpha=0.4,linetype = "dashed")+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Multilayer: log10(PageRank)") + ylab("Bipartite: log10(PageRank)")
  #labs(Color = "ID",)


Centr_ML_ID_Plot <- Centr_ML_ID_final %>% 
  group_by(ID,Plot) %>% 
  summarise(Aggr=sum(Aggr), Multi=sum(Multi),
            Delta_Strength=sum(Delta_Strength),
            Strength_Aggr=sum(Strength_Aggr))

Centr_ML_ID_Plot <- Centr_ML_ID_Plot %>% mutate(new_ID=paste0(ID," ",Plot))
my_seq <- seq(0, 0.4, length = nrow(Centr_ML_ID_Plot %>% filter(ID!="Visitor")))

ggplot(Centr_ML_ID_Plot %>% filter(ID!="Visitor"))+
  geom_point(aes(x=(Multi),y=(Aggr),color=ID, fill = ID, size=Strength_Aggr,
                 shape=ID),
             position = "jitter")+
  geom_abline(intercept = 0, slope = 1,color='steelblue', 
              size=1, alpha=0.4,linetype = "dashed")+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  #scale_shape_manual(values=15:19)+
  #ggtitle(paste0("Plot ",i)) +
  xlab("Multilayer: PageRank") + ylab("Bipartite: PageRank")
#labs(Color = "ID",)
