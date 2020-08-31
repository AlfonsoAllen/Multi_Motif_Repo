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

cz_values_function <- function(edge_module_info_i){

  # k.is <- edge_module_info_i %>% filter(module_to==module_from) %>% 
  #   group_by(node_from,module_to) %>% count(wt=weight) %>%
  #   rename(k.is = n)
  
  k.is_out <- edge_module_info_i %>% filter(module_to==module_from) %>% 
    group_by(node_from,module_to) %>% count(wt=weight) %>%
    rename(k.is_out = n, node=node_from, module=module_to)
  
  k.is_in <- edge_module_info_i %>% filter(module_to==module_from) %>% 
    group_by(node_to,module_from) %>% count(wt=weight) %>% 
    rename(k.is_in = n, node=node_to, module=module_from)
  
  k.is <- k.is_in %>% left_join(k.is_out,by=c("node","module")) %>%
    mutate(k.is=k.is_out+k.is_in)
  
  # ks.bar<- k.is %>% group_by(module_to) %>% summarise(ks.bar=mean(k.is,na.rm = T))
  # SD.ks <- k.is %>% group_by(module_to) %>% summarise(SD.ks=sd(k.is,na.rm = T))
  
  ks.bar<- k.is %>% group_by(module) %>% summarise(ks.bar=mean(k.is,na.rm = T))
  SD.ks <- k.is %>% group_by(module) %>% summarise(SD.ks=sd(k.is,na.rm = T))
  
  # k.it <- edge_module_info_i %>% filter(module_to!=module_from) %>% 
  #   group_by(node_from,module_to) %>% count(wt=weight) %>% rename(k.it = n)
  
  k.it_out <- edge_module_info_i %>% filter(module_to!=module_from) %>% 
    group_by(node_from,module_to) %>% count(wt=weight) %>% 
    rename(k.it_out = n, node=node_from, module=module_to)
  
  k.it_in <- edge_module_info_i %>% filter(module_to!=module_from) %>% 
    group_by(node_to,module_from) %>% count(wt=weight) %>% 
    rename(k.it_in = n, node=node_to, module=module_from)
  
  k.it <- k.it_in %>% left_join(k.it_out,by=c("node","module")) %>%
    mutate(k.it=k.it_out+k.it_in)
  
  # k.i <- edge_module_info_i %>% group_by(node_from) %>% count(wt=weight) %>% rename(k.i = n)
  
  k.i_out <- edge_module_info_i %>% group_by(node_from) %>% count(wt=weight) %>% 
    rename(k.i_out = n,node=node_from)
  k.i_in <- edge_module_info_i %>% group_by(node_to) %>% count(wt=weight) %>% 
    rename(k.i_in = n,node=node_to)
  k.i <- k.i_in %>% left_join(k.i_out,by=c("node")) %>%
    mutate(k.i=k.i_out+k.i_in)
  
  cz_edge_module_info_i <- k.is %>% 
    left_join(ks.bar,by=c("module")) %>% 
    left_join(SD.ks,by=c("module")) %>% 
    left_join(k.it,by=c("node","module"))%>% 
    left_join(k.i,by=c("node")) 
  
  cz_edge_module_info_i$k.is[is.na(cz_edge_module_info_i$k.is)] <- 0
  cz_edge_module_info_i$k.it[is.na(cz_edge_module_info_i$k.it)] <- 0
  
  cz_edge_module_info_i <- cz_edge_module_info_i %>%
    mutate(z = (k.is - ks.bar) / SD.ks, aux_c =  (k.it/k.i)^2)
  
  c <- cz_edge_module_info_i %>% group_by(node) %>% 
    count(wt=aux_c) %>% rename(sum_aux_c = n) %>% mutate(c=1-sum_aux_c)
  
  cz_edge_module_info_i <- cz_edge_module_info_i %>% left_join(c,by="node")

  return(cz_edge_module_info_i)
}
