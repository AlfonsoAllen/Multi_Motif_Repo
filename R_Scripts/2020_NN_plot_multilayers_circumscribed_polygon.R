library(tidyverse)
library(igraph)
library(RColorBrewer)
library(installr)

min_between_poll_space_finder <- function(plot_id){
  
  # Load multilayer----------
  file_i <- paste0("Processed_data/NN_networks/Plot_",plot_id,"_NN_intra_inter.rds")
  mult_i <- readRDS(file_i)
  
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
  
  plant_per_species <- nodes_prop %>% filter(type == "plant") %>% group_by(layer) %>% 
    count() %>% rename(number_nodes = n)
  
  # is it possible to circumscribe a polygon?
  between <- 0
  while(candidates_for_between(plant_per_species,between)==F){
    between <- between + 1
  }
  
  return(between)
}

candidates_for_between <- function(plant_per_species,between){
  
  list_lengths <- NULL
  
  for(i in 1:nrow(plant_per_species)){
    
    list_lengths <- c(list_lengths,plant_per_species$number_nodes[i],between)
    
  }
  
  max_edge_aux <- max(list_lengths)
  opt_radius <- function (R) sum(2*asin(0.5*list_lengths/R))-2*pi
  
  first_limit <- opt_radius(max_edge_aux)
  
  test_denominator <- cumsum(c(0.01,rep(0.001,5000)))
  
  second_limit_candidates <- NULL
  for(j in test_denominator){
    second_limit_candidates <- c(second_limit_candidates,opt_radius(max_edge/j))
  }
  
  if(first_limit < 0){
    candidates <- second_limit_candidates[second_limit_candidates>0]
  }else{
    candidates <- second_limit_candidates[second_limit_candidates<0]
  }
  
  if((sum(is.na(candidates))+sum(is.nan(candidates))==length(candidates)) == T){
    result = F
  }else{result = T}

  
  return(result)
}


circumscribe_polygon <- function(list_nodes){
  
  # the sum of all the angles of the edge (chords) must be equal to 2pi
  max_edge <- max(list_nodes)
  
  # Select proper limits----------
  
  opt_radius <- function (R) sum(2*asin(0.5*list_nodes/R))-2*pi
  
  first_limit <- opt_radius(max_edge)
  
  test_denominator <- cumsum(c(0.01,rep(0.001,5000)))
  
  second_limit_candidates <- NULL
  for(j in test_denominator){
    second_limit_candidates <- c(second_limit_candidates,opt_radius(max_edge/j))
  }
  
  if(first_limit < 0){
    candidates <- second_limit_candidates[second_limit_candidates>0]
  }else{
    candidates <- second_limit_candidates[second_limit_candidates<0]
  }
  
  candidates <- candidates[!is.na(candidates)]
  candidates <- candidates[!is.nan(candidates)]
  
  index_candidate <- which(second_limit_candidates==candidates[1])
  
  second_limit <- max_edge/test_denominator[index_candidate]
  
  if(max_edge<second_limit){
    
    r_cycle <- uniroot(function (R) sum(2*asin(0.5*list_nodes/R))-2*pi,
                       lower = max_edge,
                       upper = second_limit,
                       tol = 0.0001)$root
    
  }else{
    
    r_cycle <- uniroot(function (R) sum(2*asin(0.5*list_nodes/R))-2*pi,
                       lower = second_limit,
                       upper = max_edge,
                       tol = 0.0001)$root
    
  }

  
  # distribute the vertices of the polygon over the circle with r_cycle
  
  angles_cycle <- 2*asin(0.5*list_nodes/r_cycle)
  angles_cycle_accumulated <- cumsum(angles_cycle)
  x_vertex <- r_cycle*cos(angles_cycle_accumulated)
  y_vertex <- r_cycle*sin(angles_cycle_accumulated)
  
  # Calculating centroid
  x_vertex_aux <- c(x_vertex,x_vertex[1])
  y_vertex_aux <- c(y_vertex,y_vertex[1])
  
  x_vertex_aux2 <- c(x_vertex2,x_vertex2[1])
  y_vertex_aux2 <- c(y_vertex2,y_vertex2[1])
  
  Cx <- 0
  Cy <- 0
  A <- 0
  
  for(i in 1:length(x_vertex)){
    
    Cx <- Cx + (1/6)*(x_vertex_aux[i] + x_vertex_aux[i+1])*(x_vertex_aux[i]*y_vertex_aux[i+1]-x_vertex_aux[i+1]*y_vertex_aux[i])
    Cy <- Cy + (1/6)*(y_vertex_aux[i] + y_vertex_aux[i+1])*(x_vertex_aux[i]*y_vertex_aux[i+1]-x_vertex_aux[i+1]*y_vertex_aux[i])
    A <- A + 0.5*(x_vertex_aux[i]*y_vertex_aux[i+1]-x_vertex_aux[i+1]*y_vertex_aux[i])
    
  }
  
  Cx <- Cx/A
  Cy <- Cy/A
  
  return(list(r_cycle,angles_cycle,Cx,Cy))
  
}

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

plot_multilayer <- function(plot_id,between_poll_space, r_ext){
  # Load multilayer----------
  file_i <- paste0("Processed_data/NN_networks/Plot_",plot_id,"_NN_intra_inter.rds")
  mult_i <- readRDS(file_i)
  
  # get weighted adj_matrix
  as_adjacency_matrix(mult_i,attr = "weight",sparse = T)
  
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
  
  ############################################################
  # Add coordinates to nodes-------
  
  # Circunscribe polygon--------
  # The length of its segments will be proportional to the number of plant nodes
  # of each species
  
  plant_per_species <- nodes_prop %>% filter(type == "plant") %>% group_by(layer) %>% 
    count() %>% rename(number_nodes = n)
  
  # is it possible to circumscribe a polygon?
  
  list_lengths <- NULL
  
  for(i in 1:nrow(plant_per_species)){
    
    list_lengths <- c(list_lengths,plant_per_species$number_nodes[i],between_poll_space)
    
  }
  
  # If the sides of the polygon are not long enough to circumscribe it, we add extra
  # distance to the "between pollinator spaces"
  
  if(sum(list_lengths)-max(list_lengths) < 0){
    
    max_edge <- max(list_lengths)
    
    extra_nodes <- (max_edge-sum(list_lengths))/nrow(plant_per_species)
    
    list_lengths[c(F,T)] <- list_lengths[c(F,T)] + extra_nodes
    
    
  }
  
  data_polygon <- circumscribe_polygon(list_lengths)
  r_cycle <- data_polygon[[1]]
  angles_cycle <- data_polygon[[2]]
  Cx <- data_polygon[[3]]
  Cy <- data_polygon[[4]]
  
  # distribute the vertices of the polygon over the circle with r_cycle
  
  angles_cycle_accumulated <- cumsum(angles_cycle)
  x_vertex <- r_cycle*cos(angles_cycle_accumulated)
  y_vertex <- r_cycle*sin(angles_cycle_accumulated)
  
  # Segments
  
  x_vertex_aux <- c(x_vertex[length(x_vertex)],x_vertex)
  y_vertex_aux <- c(y_vertex[length(y_vertex)],y_vertex)
  x_vertex_end <- x_vertex_aux[1:(length(x_vertex_aux)-1)]
  y_vertex_end <- y_vertex_aux[1:(length(y_vertex_aux)-1)]
  
  #Sanity check
  sqrt((x_vertex-x_vertex_end)^2+(y_vertex-y_vertex_end)^2)
  
  # Now we distribute the pollinator nodes over the corresponding segments---
  
  pollinator_per_species <- nodes_prop %>% filter(type != "plant") %>% group_by(layer) %>% 
    count() %>% rename(number_nodes = n)
  
  pollinator_per_species$length = list_lengths[c(T,F)]
  pollinator_per_species$xini = x_vertex[c(T,F)]
  pollinator_per_species$xend = x_vertex_end[c(T,F)]
  pollinator_per_species$yini = y_vertex[c(T,F)]
  pollinator_per_species$yend = y_vertex_end[c(T,F)]
  
  # Sanity_check
  pollinator_per_species %>% mutate(dis =
                                      sqrt((xini-xend)^2+(yini-yend)^2))
  
  min_separation <- min(pollinator_per_species$length/(pollinator_per_species$number_nodes-1))
  
  coordinates_nodes_pollinator <- NULL
  
  for(i in 1:nrow(pollinator_per_species)){
    
    coordinates_aux <- nodes_prop %>% 
      filter(layer == pollinator_per_species$layer[i],type == "visitor") %>%
      select(name)
    
    hor_points <- equidistr_nodes_x(min_separation,pollinator_per_species$number_nodes[i])
    rot_angle <- atan((pollinator_per_species$yend[i]-pollinator_per_species$yini[i])/
                        (pollinator_per_species$xend[i]-pollinator_per_species$xini[i]))
    real_x <- hor_points*cos(rot_angle)+
      0.5*(pollinator_per_species$xend[i]+pollinator_per_species$xini[i])
    real_y <- hor_points*sin(rot_angle)+
      0.5*(pollinator_per_species$yend[i]+pollinator_per_species$yini[i])
    
    coordinates_aux$x <- real_x
    coordinates_aux$y <- real_y
    
    coordinates_nodes_pollinator <- bind_rows(coordinates_nodes_pollinator,
                                              coordinates_aux)
  }
  
  # # Sanity check 
  # circleFun <- function(center = c(0,0),diameter = 2*r_cycle, npoints = 100){
  #   r = diameter / 2
  #   tt <- seq(0,2*pi,length.out = npoints)
  #   xx <- center[1] + r * cos(tt)
  #   yy <- center[2] + r * sin(tt)
  #   return(data.frame(x = xx, y = yy))
  # }
  # 
  # circle_plot <- circleFun(c(0,0),2*r_cycle,npoints = 100)
  # 
  # all_segments <- data.frame(xini = x_vertex,
  #                            xend = x_vertex_end,
  #                            yini = y_vertex,
  #                            yend = y_vertex_end)
  # 
  # ggplot(circle_plot,aes(x,y)) + geom_path()+
  #   geom_point(data = all_segments,aes(x=xini,y=yini))+
  #   geom_segment(data = all_segments,
  #                aes(x=xini,y=yini,xend=xend,yend=yend))+
  #   geom_point(data = pollinator_per_species,aes(x=xini,y=yini,color=layer))+
  #   geom_point(data = pollinator_per_species,aes(x=xend,y=yend,color=layer))+
  #   geom_point(data = coordinates_nodes_pollinator,aes(x,y))
    
  # Now we distribute the plant nodes over the arc of each cord (polygon segment)
  # We project those positions to a circunference of radius r_ext
  
  # distribute the vertices of the polygon over the circle with r_cycle
  
  angles_cycle_accumulated <- cumsum(angles_cycle)
  x_vertex_ext <- r_ext*cos(angles_cycle_accumulated)
  y_vertex_ext <-  r_ext*sin(angles_cycle_accumulated)
  
  # Segments
  
  x_vertex_ext_aux <- c(x_vertex_ext[length(x_vertex_ext)],x_vertex_ext)
  y_vertex_ext_aux <- c(y_vertex_ext[length(y_vertex_ext)],y_vertex_ext)
  x_vertex_ext_end <- x_vertex_ext_aux[1:(length(x_vertex_ext_aux)-1)]
  y_vertex_ext_end <- y_vertex_ext_aux[1:(length(y_vertex_ext_aux)-1)]
  
  
  # Now we distribute the pollinator nodes over the corresponding segments---
  
  plant_per_species <- nodes_prop %>% filter(type == "plant") %>% group_by(layer) %>% 
    count() %>% rename(number_nodes = n)
  plant_per_species$angle_arc = angles_cycle[c(T,F)]
  plant_per_species$xini = x_vertex_ext[c(T,F)]
  plant_per_species$xend = x_vertex_ext_end[c(T,F)]
  plant_per_species$yini = y_vertex_ext[c(T,F)]
  plant_per_species$yend = y_vertex_ext_end[c(T,F)]
  
  coordinates_nodes_plant <- NULL
  
  for(i in 1:nrow(plant_per_species)){
    
    coordinates_aux <- nodes_prop %>% 
      filter(layer == plant_per_species$layer[i],type == "plant") %>%
      select(name)
    
    arc_incr <- plant_per_species$angle_arc[i] /(plant_per_species$number_nodes[i]-1)
    arc_x_ini <- atan(plant_per_species$yend[i]/plant_per_species$xend[i])
    
    if(plant_per_species$yend[i]*plant_per_species$xend[i] < 0 & 
       plant_per_species$yend[i] > plant_per_species$xend[i] &
       arc_x_ini < 0){
      arc_x_ini <- arc_x_ini + pi
    }else if(plant_per_species$yend[i]*plant_per_species$xend[i] > 0 & 
        plant_per_species$yend[i] < 0 &
        arc_x_ini > 0){
      arc_x_ini <- arc_x_ini + pi
    }
    
    arc_aux <- c(arc_x_ini,rep(arc_incr,(plant_per_species$number_nodes[i]-1)))
    arc_aux_acc <- cumsum(arc_aux)
    
    coordinates_aux$x <- (r_ext*cos(arc_aux_acc))
    coordinates_aux$y <-  (r_ext*sin(arc_aux_acc))
    
    coordinates_nodes_plant <- bind_rows(coordinates_nodes_plant,
                                              coordinates_aux)
  }
  
  
  # # Sanity check
  # circleFun <- function(center = c(0,0),diameter = 2*r_ext, npoints = 100){
  #   r = diameter / 2
  #   tt <- seq(0,2*pi,length.out = npoints)
  #   xx <- center[1] + r * cos(tt)
  #   yy <- center[2] + r * sin(tt)
  #   return(data.frame(x = xx, y = yy))
  # }
  # 
  # circle_plot <- circleFun(c(0,0),2*r_ext,npoints = 100)
  # 
  # all_segments <- data.frame(xini = x_vertex_ext,
  #                            xend = x_vertex_ext_end,
  #                            yini = y_vertex_ext,
  #                            yend = y_vertex_ext_end)
  # 
  # ggplot(circle_plot,aes(x,y)) + geom_path()+
  #   geom_point(data = all_segments,aes(x=xini,y=yini))+
  #   geom_segment(data = all_segments,
  #                aes(x=xini,y=yini,xend=xend,yend=yend))+
  #   geom_point(data = plant_per_species,aes(x=xini,y=yini,color=layer))+
  #   geom_point(data = plant_per_species,aes(x=xend,y=yend,color=layer))+
  #   geom_point(data = coordinates_nodes_plant,aes(x,y))

  
  coordinates_nodes <- bind_rows(coordinates_nodes_plant,coordinates_nodes_pollinator)
  
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
  
 result <- plot.igraph(mult_i,
              vertex.label = NA, vertex.label.color = data_net$color,
              vertex.size = 5,vertex.size2 = 5,
              vertex.color = data_net$color, vertex.frame.color = "gray20",
              vertex.shape = data_net$shape,
              edge.arrow.size=0.2, edge.color = edge_prop_color$color, 
              edge.width = E(mult_i)$weight,
              edge.curved = F, 
              layout = l)+
   legend("right",legend=nodes_prop$layer %>% unique,
         fill=nodes_prop$color %>% unique,
         bty = "n")
 
 return(result)
 #
}


#dev.off()
#par(mfrow=c(3,3))
par(mar=c(0,0,0,0)+.1)
min_between_poll_space_finder(plot_id=1)
multi_1 <- plot_multilayer(1,between_poll_space=11, r_ext= 40)

min_between_poll_space_finder(plot_id=2)
multi_2 <- plot_multilayer(2,between_poll_space=5, r_ext= 40)

min_between_poll_space_finder(plot_id=3)
multi_3 <- plot_multilayer(3,between_poll_space=4, r_ext= 40)

min_between_poll_space_finder(plot_id=4)
multi_4 <- plot_multilayer(4,between_poll_space=18, r_ext= 40)

min_between_poll_space_finder(plot_id=5)
multi_5 <- plot_multilayer(5,between_poll_space=13, r_ext= 40)

min_between_poll_space_finder(plot_id=6)
multi_6 <- plot_multilayer(6,between_poll_space=9, r_ext= 40)

min_between_poll_space_finder(plot_id=7)
multi_7 <- plot_multilayer(7,between_poll_space=7, r_ext= 40)

min_between_poll_space_finder(plot_id=8)
multi_8 <- plot_multilayer(8,between_poll_space=5, r_ext= 40)

min_between_poll_space_finder(plot_id=9)
multi_9 <- plot_multilayer(9,between_poll_space=5, r_ext= 40)
