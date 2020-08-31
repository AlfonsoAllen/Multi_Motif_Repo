# c = 1 - sum( (k.it/k.i)^2) # among-module connectivity = participation coefficient P in Guimerà & Amaral
# 
# z = (k.is - ks.bar) / SD.ks # within-module degree
# 
# k.is = number of links of i to other species in its own module s
# ks.bar = average k.is of all species in module s
# SD.ks = standard deviation of k.is of all species in module s
# k.it = number of links of species i to module t
# k.i = degree of species i
# 
# Note that for any species alone (in its level) in a module the z-value will be NaN,
# since then SD.ks is 0. This is a limitation of the way the z-value is defined
# (in multiples of degree/strength standard deviations).
# 
# Olesen et al. (2006) give critical c and z values of 0.62 and 2.6, respectively.
# Species exceeding these values are deemed connectors or hubs of a network.
# The justification of these thresholds remains unclear to me. They may also not apply
# for the quantitative version.

cz_values_function_OUT <- function(edge_module_info_i){

  k.is <- edge_module_info_i %>% filter(module_to==module_from) %>% 
    group_by(node_from,layer_from,module_to,Plot) %>% count(wt=weight) %>% rename(k.is = n)
  ks.bar<- k.is %>% group_by(module_to,Plot) %>% summarise(ks.bar=mean(k.is,na.rm = T))
  SD.ks <- k.is %>% group_by(module_to,Plot) %>% summarise(SD.ks=sd(k.is,na.rm = T))
  
  k.it <- edge_module_info_i %>% filter(module_to!=module_from) %>% 
    group_by(node_from,layer_from,module_to,Plot) %>% count(wt=weight) %>% rename(k.it = n)
  
  k.i <- edge_module_info_i %>% group_by(node_from,layer_from,Plot) %>% count(wt=weight) %>% rename(k.i = n)
  
  z_edge_module_info_i <- k.is %>% 
    left_join(ks.bar,by=c("module_to","Plot")) %>% 
    left_join(SD.ks,by=c("module_to","Plot"))
  
  c_edge_module_info_i <- k.it %>% 
    left_join(k.i,by=c("node_from","layer_from","Plot")) %>%
    mutate(aux_c =  (k.it/k.i)^2)
  
  c <- c_edge_module_info_i %>% group_by(node_from,layer_from,Plot) %>% 
    count(wt=aux_c) %>% rename(sum_aux_c = n) %>% mutate(c=1-sum_aux_c)
  
  
  z_edge_module_info_i <- z_edge_module_info_i %>%
    mutate(z = (k.is - ks.bar) / SD.ks)
  
  cz_edge_module_info_i <- z_edge_module_info_i %>% 
    left_join(k.i,by=c("node_from","layer_from","Plot")) %>%
    left_join(c,by=c("node_from","layer_from","Plot"))
  
  cz_edge_module_info_i$sum_aux_c[is.na(cz_edge_module_info_i$sum_aux_c)] <- 0
  cz_edge_module_info_i$c[is.na(cz_edge_module_info_i$c)] <- 1
  return(cz_edge_module_info_i)
}
