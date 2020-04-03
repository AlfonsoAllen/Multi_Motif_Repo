# Load libraries
library(attempt)
library(igraph)
library(bipartite)
library(tidyverse)
library(betalink)


#Access layers files
dir_ini <- getwd()

folder_base <- paste(dir_ini,"/Processed_data/Multilayer_WEEK",sep="")

files_base <- list.files(folder_base)

setwd(folder_base)

# Select the plot number to build up the multilayer
Plot_i <- 8

# Extract layer files for Plot_i

list_files_field_level <- files_base[grepl(paste("Plot_",Plot_i,sep = ""), files_base)]

# Generate layer metadata (layer_id,name[week])
layer_metadata_i <- tibble(layer_id=1:length(list_files_field_level), layer_name=NA)

# Extract edge_list for each layer
for (i in 1:length(list_files_field_level)){
  # Get the week (layer) information and add it to layer_metadata_i
  week <- strsplit(list_files_field_level[i],".csv")
  week <- strsplit(week[[1]][1],".WEEK_")
  week <- as.numeric(week[[1]][2])
  layer_metadata_i$layer_name[i] <- week
}


Plot_matrices_i <- NULL
# Extract edge_list for each layer
for (i in 1:length(list_files_field_level)){

  print(list_files_field_level[i])
  # Extract the incidence matrix
  inc_matrix <- read.csv(list_files_field_level[i], header=T, row.names=1)
  
  # Create a graph for each layer
  g_i <- graph_from_incidence_matrix(inc_matrix, directed = FALSE, weighted = T)
  
  # Extract the adjacency matrix for each layer
  g_i_matrix <- as_adjacency_matrix(g_i, names = T, sparse = F, attr = 'weight')
  
  #g_i_matrix[g_i_matrix > 0] <- 1
  
  # Add the adjacency matrix to our list
  print(paste("WEEK_",layer_metadata_i$layer_name[i],sep=""))
  Plot_matrices_i[[paste("WEEK_",layer_metadata_i$layer_name[i],sep="")]] <- g_i_matrix
}

# Sort matrices by week number
layer_metadata_i <- layer_metadata_i %>% 
  mutate(layer_full_name=paste("WEEK_",layer_name,sep="")) %>%
  arrange(layer_name)

Plot_matrices_i_sorted <- Plot_matrices_i[layer_metadata_i$layer_full_name]

# Set the base dir again
setwd(dir_ini)


# 3. Variation in connectivity across layers and beta-diversity
# In ecology, turnover is many times quantified using beta-diversity. 
# We will calculate beta diversity between layers following the
# Poisot et al 2012 framework. This is particularly useful for spatial data,
# We are working without interlinlks

# Differences in interactions between networks (bWN) originate from
# differences in species composition (bST, dissimilarity in 
# interaction structure introduced by dissimilarity in species composition),
# and because shared species between the two realisations may interact 
# differently (bOS, dissimilarity of interactions in co-occuring species).
# This leads to an additive view of network dissimilarity, wherein:
#  bWN = bST + bOS


# bST takes values between 0 (dissimilarity between two networks is
# entirely explained by shared species interacting differently), and
# bWN (the shared species interact in the same way, and all the difference
# between the two networks is explained by species turnover).

#https://cran.r-project.org/web/packages/betalink/betalink.pdf

# Prepare the matrices we have for betalink

# "prepare_networks": Note that a single non-directed edge between a pair of nodes is
# decodified as two directed edges
igraph_layers <- prepare_networks(Plot_matrices_i_sorted, directed = F)

# BETALINK
#While interpreting the output, it is important to consider
#that ST is strongly constrained by the values of S
#(the species composition dissimilarity). ST is
#only really meaningful when the values of S are "intermediate";
#a good example is when the networks have been sampled along a gradient,
#and a more or less equal proportion of the species show
#turnover from one step to the next. In the situations 
#where S is either really high or really low, the
#values of ST are constrained and should no be given importance. 
#The values of OS and WN, and
#how they relate to S, have more informative value.

##############################################################
##############################################################
# Meassure of dissimilarity: bw (Whittaker, 1960) [See Eq. 2 in Poisot et al. 2012]
# If the elements we compare in both networks are present in both systems bw = 0,
# whereas b_w = 1  there is no overlap.

# An example with two layers
betalink(igraph_layers[["WEEK_11"]], igraph_layers[["WEEK_17"]]) #S=0.9393939

#We compare the species
row.names(Plot_matrices_i[["WEEK_17"]]) %in% row.names(Plot_matrices_i[["WEEK_11"]])

plot(igraph_layers[["WEEK_11"]])
plot(igraph_layers[["WEEK_17"]])

#species in both netwoks
a <-  sum(row.names(Plot_matrices_i[["WEEK_11"]])%in% row.names(Plot_matrices_i[["WEEK_17"]]))
#species only in network WEEK_11
b <-  length(row.names(Plot_matrices_i[["WEEK_11"]])) - a
#species only in network WEEK_17
c <-  length(row.names(Plot_matrices_i[["WEEK_17"]])) - a

b_W <- 2*(a+b+c)/(2*a+b+c)-1 #S=0.9393939 (Check!)

#---------------------
# We compare all the connections

betalink(igraph_layers[["WEEK_11"]], igraph_layers[["WEEK_17"]]) #WN = 1

# Calculate the intersection graph edges:
E(intersection(igraph_layers[["WEEK_11"]],igraph_layers[["WEEK_17"]])) #0/0 edges

#links only in both networks
a <-  length(E(intersection(igraph_layers[["WEEK_11"]],igraph_layers[["WEEK_17"]])))
#links only in network WEEK_11
b <-   sum(E(igraph_layers[["WEEK_11"]]) %in% E(igraph_layers[["WEEK_11"]])) - a
#links only in network WEEK_17
c <-   sum(E(igraph_layers[["WEEK_17"]]) %in% E(igraph_layers[["WEEK_17"]])) - a

b_W <- 2*(a+b+c)/(2*a+b+c)-1 #WS=1 (Check!)
###########################################################
###########################################################
# Calculate all pairwise layers 

net_dissim <- as_tibble(network_betadiversity(igraph_layers))

# reorder factors
net_dissim$i <- factor(net_dissim$i, levels = layer_metadata_i$layer_full_name)
net_dissim$j <- factor(net_dissim$j, levels = layer_metadata_i$layer_full_name)


# It is easier to examine results using heatmaps
##################
#species composition dissimilarity between networks

net_dissim %>% ggplot()+
  geom_tile(aes(x=i, y=j, fill=S))+
  theme_bw()+labs(x ="",y="",title=paste("Plot ",Plot_i,': Dissimilarty of species composition',sep=""))+
  theme(axis.text.x = element_text(angle=-90))

# Dissimilarty of INTERACTIONS

net_dissim %>% ggplot()+
  geom_tile(aes(x=i, y=j, fill=WN))+theme_bw()+labs(x ="",y="",title=paste("Plot ",Plot_i,': Dissimilarty of interactions',sep=""))+
  theme(axis.text.x = element_text(angle=-90))



ggplot(net_dissim)+
  geom_point(aes(x=S, y=WN),alpha=0.2)+theme_bw()+
  labs(x ="",y="",
       title=paste("Plot ",Plot_i,': Pairwise dissimilarty of species composition (beta_S) VS Dissimilarty of interactions (beta_WN)',sep=""))



#########################
# OS <> component of the beta-diversity WHEN BOTH SPECIES OCCURS IN BOTH PATCHES BUT WEIGHT INTERACTIONS ARE DIFFERENT


net_dissim %>% ggplot()+
  geom_tile(aes(x=i, y=j, fill=OS))+theme_bw()+
  labs(x ="",y="",title=paste("Plot ",Plot_i,': Dissimilarty of interactions between shared species',sep=""))+
  theme(axis.text.x = element_text(angle=-90))



# NaN appears because terms a, b and c in equation 2 are equal to zero and 0/0=NaN appears
# (see below the new columns added to net_dissim)

for (i in 1:nrow(net_dissim)){
  m1 <- net_dissim$i[i]
  m2 <- net_dissim$j[i]
  net_dissim$s_i[i] <- length(row.names(Plot_matrices_i[[m1]]))
  net_dissim$s_j[i] <- length(row.names(Plot_matrices_i[[m2]]))
  net_dissim$s_comunes[i] <-  sum(row.names(Plot_matrices_i[[m1]]) %in% row.names(Plot_matrices_i[[m2]]))
  
  #Subgraph of layer m1 that only contains shared nodes with m2
  int_graph_1 <- induced_subgraph(igraph_layers[[m1]],
                                  names(V(igraph_layers[[m1]]))[names(V(igraph_layers[[m1]]))%in% names(V(igraph_layers[[m2]]))])
  
  #Subgraph of layer m2 that only contains shared nodes with m1
  int_graph_2 <- induced_subgraph(igraph_layers[[m2]],
                                  names(V(igraph_layers[[m2]]))[names(V(igraph_layers[[m2]]))%in% names(V(igraph_layers[[m1]]))])
  
  net_dissim$links_s_com_i[i] <- length(E(int_graph_1))
  
  net_dissim$links_s_com_j[i] <- length(E(int_graph_2)) 
  
}

