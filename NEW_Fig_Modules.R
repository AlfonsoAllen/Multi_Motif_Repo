library(tidyverse)

for (Plot_i in 1:9){
  mod_i <- read_csv(paste0("Processed_data/Modularity_Pheno_Overlap/Modularity_Plot",Plot_i,".csv")  )
  
  if (Plot_i==1){
    modularity <- mod_i
  }else{
    modularity <- bind_rows(modularity, mod_i)
  }
  
}

mod_plot <- modularity %>% select(module, Plot) %>% unique() %>% group_by(Plot) %>% count()
mean(mod_plot$n)
sd(mod_plot$n)

mod_plot_layers <- modularity %>% select(module, Plot, layer_name) %>% unique() %>%
  group_by(Plot,module) %>% count() #%>% group_by(Plot, n) %>% count() 

ggplot(mod_plot_layers, aes(x=n)) + geom_histogram()+facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  theme(legend.position = 'none', axis.text.x = element_text(angle = 0))+
  xlab("# Species (Layer) per module")+ylab("# Modules")


modularity %>% filter(type=="plant") %>% group_by(Plot,module) %>% count() 


x <- modularity %>% group_by(module, layer_name, Plot) %>% count()

# Number of species in a combination of layer-module
modularity %>% group_by(module, layer_name, Plot) %>% count() %>% 
  ggplot()+
  geom_point(aes(layer_name, module, size=n, color=as.factor(module)))+
  geom_line(aes(layer_name, module, color=as.factor(module)))+
  theme_bw()+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  theme(legend.position = 'none', axis.text.x = element_text(angle = 0))+
  xlab("Species (Layer)")+ylab("Module ID")


