library(tidyverse)

data <- read_csv("Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID_RAPE.csv")

plant_species_plot <- data %>% select(Plot,Plant_Simple) %>% unique() %>% group_by(Plot) %>%
  count()
mean(plant_species_plot$n)
sd(plant_species_plot$n)

plant_ind_plot <- data %>% select(Plot,Subplot,Plant_Simple) %>% unique() %>% group_by(Plot) %>%
  count()
mean(plant_ind_plot$n)
sd(plant_ind_plot$n)

poll_species_plot <- data %>% select(Plot,ID_Simple) %>% unique() %>% group_by(Plot) %>%
  count()
mean(poll_species_plot$n)
sd(poll_species_plot$n)

edge_list <- read_csv("Processed_data/Modularity_Pheno_Overlap/Edge_list_Phen_Over_PLOT.csv")

nodes <- edge_list %>% select(Plot,node_from,layer_from) %>% unique() %>% group_by(Plot) %>%
  count()
mean(nodes$n)
sd(nodes$n)

links <- edge_list %>% filter(layer_from==layer_to) %>% select(Plot,node_from) %>% unique() %>% group_by(Plot) %>%
  count()
mean(links$n)
sd(links$n)

links_inter <- edge_list %>% filter(layer_from!=layer_to) %>% select(Plot,node_from) %>% unique() %>% group_by(Plot) %>%
  count()
mean(links_inter$n)
sd(links_inter$n)

modules <- NULL

for (i in 1:9){
  plot_modules_NN_i <- read_csv(paste0("Processed_data/Modularity_Pheno_Overlap/NN_Modularity_Plot",
                                     i,".csv"))
  modules <- bind_rows(modules,plot_modules_NN_i)
                                
}

modules_l <- modules %>% select(Plot,module) %>% unique() %>% group_by(Plot) %>%
  count()
mean(modules_l$n)
sd(modules_l$n)
