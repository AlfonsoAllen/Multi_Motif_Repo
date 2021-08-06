library(tidyverse)
library(igraph)
library(RColorBrewer)

# Load multilayer----------
mult_i <- readRDS("Processed_data/NN_networks/Plot_8_NN_intra_inter.rds")

# get weighted adj_matrix
as_adjacency_matrix(mult_i,attr = "weight",sparse = T)


plot.igraph(mult_i, vertex.label = F)

nodes <- V(mult_i)$name
nodes_prop <- tibble(name = nodes) %>% separate(name, c("node", "layer")," ")
nodes_prop$type <- NA
nodes_prop$type[nchar(nodes_prop$node) > 2] <- "visitor"
nodes_prop$type[!nchar(nodes_prop$node) > 2] <- "plant"

nodes_prop$shape <- NA
nodes_prop$shape[nodes_prop$type == "visitor"] <- "circle"
nodes_prop$shape[!nodes_prop$type == "visitor"] <- "rectangle"

nodes_prop$name <- paste(nodes_prop$node,nodes_prop$layer, sep = " ") 

nodes_prop <- nodes_prop %>% arrange(layer,type)

# Add color to nodes-------

fitness_data2 <- read_csv2("Raw_Data/final_Pollinators_2020.csv")
plant_list <- fitness_data2$Plant %>% unique()
plant_list <- plant_list[!plant_list %in% c("OUT","0")] 
plant_list <- plant_list[!is.na(plant_list)] 

colors_plants <- tibble(layer = plant_list,
                        color = RColorBrewer::brewer.pal(length(plant_list),"Paired"))

nodes_prop <- nodes_prop %>% left_join(colors_plants, by = "layer")


# Add coordinates to nodes-------

nodes_prop_plant_LEMA <- nodes_prop %>% filter(type == "plant",
                                               layer == "LEMA") %>% select(name) 
nodes_prop_pollinator_LEMA <- nodes_prop %>% filter(type != "plant",
                                                    layer == "LEMA") %>% select(name)
nodes_prop_plant_OTHER <- nodes_prop %>% filter(type == "plant",
                                                layer != "LEMA") %>% select(name) 
nodes_prop_pollinator_OTHER <- nodes_prop %>% filter(type != "plant",
                                                     layer != "LEMA") %>% select(name)

col_plant_LEMA <- nrow(nodes_prop_plant_LEMA)
col_pollinator_LEMA <- nrow(nodes_prop_pollinator_LEMA)
col_plant_OTHER <- nrow(nodes_prop_plant_OTHER)
col_pollinator_OTHER <- nrow(nodes_prop_pollinator_OTHER)

col_max <- max(col_plant_LEMA,col_pollinator_LEMA,
               col_plant_OTHER,col_pollinator_OTHER)

wide_min <- 2/ col_max

equidistr_nodes_x <- function(wide_min,number_nodes){
  
  if(number_nodes %% 2 != 0){
    
    aux <- wide_min*(0:(number_nodes-1))
    result <- rep(-wide_min*0.5*(number_nodes-1),number_nodes) + aux
    return(result)
  }else{
    
    aux <- wide_min*(0:(number_nodes-1))
    result <- rep(-wide_min*0.5*(number_nodes),number_nodes) + aux + wide_min*0.5
    return(result)
    
  }
  
}

nodes_prop_plant_LEMA <-nodes_prop_plant_LEMA %>%
  mutate(x = equidistr_nodes_x(wide_min,col_plant_LEMA),
         y = -1) %>% select(name, x, y)

nodes_prop_pollinator_LEMA <-nodes_prop_pollinator_LEMA %>%
  mutate(x = equidistr_nodes_x(wide_min,col_pollinator_LEMA),
         y = -1/6) %>% select(name, x, y)

nodes_prop_plant_OTHER <-nodes_prop_plant_OTHER %>%
  mutate(x = equidistr_nodes_x(wide_min,col_plant_OTHER),
         y = 1) %>% select(name, x, y)

nodes_prop_pollinator_OTHER <-nodes_prop_pollinator_OTHER %>%
  mutate(x = equidistr_nodes_x(wide_min,col_pollinator_OTHER),
         y = 1/6) %>% select(name, x, y)

coordinates_nodes <- bind_rows(nodes_prop_plant_LEMA,nodes_prop_pollinator_LEMA,
                               nodes_prop_plant_OTHER,nodes_prop_pollinator_OTHER)

data_net <- tibble(name = nodes) %>% left_join(nodes_prop, by = "name") %>%
  left_join(coordinates_nodes, by = "name")

l <- as.matrix(data_net[,c(7,8)])


# Add color to edges---------
# get weighted edge_list
edge_prop <- cbind( get.edgelist(mult_i) , round( E(mult_i)$weight, 3 ))  %>% 
  as_tibble() %>% rename(from = V1, to = V2, weight = V3) %>% 
  mutate(weight = as.numeric(weight)) %>%
  separate(from,c("node_from","layer_from")," ") %>%
  separate(to,c("node_to","layer_to")," ")

edge_prop$layer <- NA
edge_prop$type <- NA

edge_prop$type[edge_prop$layer_from==edge_prop$layer_to] <- "intra"
edge_prop$type[edge_prop$layer_from!=edge_prop$layer_to] <- "inter"
edge_prop$layer[edge_prop$type=="intra"] <- edge_prop$layer_from[edge_prop$type=="intra"]
edge_prop$layer[edge_prop$type!="intra"] <- "inter"

edge_prop_color <- edge_prop %>% left_join(colors_plants, by = "layer")
edge_prop_color$color[is.na(edge_prop_color$color)] <- "gray20"


V(mult_i)$type[data_net$type=="plant"] <- TRUE
V(mult_i)$type[data_net$type!="plant"] <- FALSE

plot.igraph(mult_i,
            vertex.label = NA, vertex.label.color = data_net$color,
            vertex.size = 5,vertex.size2 = 5,
            vertex.color = data_net$color, vertex.frame.color = "gray20",
            vertex.shape = data_net$shape,
            edge.arrow.size=0.1, edge.color = edge_prop_color$color, 
            edge.width = E(mult_i)$weight,
            edge.curved = F, 
            layout = l)+
  legend("right",legend=nodes_prop$layer %>% unique,
         fill=nodes_prop$color %>% unique,
         bty = "n")
