
# Hypothesis: Ratio do not change significantly when isolated nodes are added to the plot
# multilayer


library(tidyverse)

# Load data


fitness_orig <- load_data_models_2020()

fitness.data <- subset(fitness_orig,Seeds_GF > 0)

##################



# Multilayer----
Plot_i=1
graph_Plot_i_all_links <- readRDS(file = paste0("Processed_data/NN_networks/Plot_",Plot_i,"_NN_intra_inter.rds"))
page_rank_i <- igraph::page_rank(graph_Plot_i_all_links,
                                 directed = TRUE, damping = 0.85,
                                 personalized = NULL, weights = NULL, options = NULL)

page_rank_i_nodes_all_links <-
  data.frame(species=names(page_rank_i[[1]]),
             Real_PR_Multi=page_rank_i[[1]], row.names=NULL)

# Add extra nodes----

# Number of extra nodes

isolated_nodes <- fitness.data %>% filter(Plot==as.character(Plot_i),DegreeIn==0) %>%
  mutate(isolated_nodes = paste0(Plant_Simple," ",Subplot)) %>% ungroup() %>%
  select(isolated_nodes) %>% nrow()

graph_Plot_i_all_links_extra_nodes <- graph_Plot_i_all_links %>%
  igraph::add_vertices(isolated_nodes, color = "red")

page_rank_i_extra_nodes <- igraph::page_rank(graph_Plot_i_all_links_extra_nodes,
                                 directed = TRUE, damping = 0.85,
                                 personalized = NULL, weights = NULL, options = NULL)

page_rank_i_nodes_all_links_extra_nodes <-
  data.frame(species=names(page_rank_i_extra_nodes[[1]]),
             Real_PR_Multi_extra=page_rank_i_extra_nodes[[1]], row.names=NULL)

pagerank_isolated <- page_rank_i_nodes_all_links_extra_nodes$Real_PR_Multi_extra[
  is.na(page_rank_i_nodes_all_links_extra_nodes$species)] %>% unique()

# Compare
page_rank_i_nodes_all_links <- page_rank_i_nodes_all_links %>%
  left_join(page_rank_i_nodes_all_links_extra_nodes,by="species")

# Uncoupled layers----

graph_Plot_i_intra_links <- readRDS(file =
                                      paste0("Processed_data/NN_networks/Plot_",Plot_i,"_NN_intra_only.rds"))
page_rank_i <- igraph::page_rank(graph_Plot_i_intra_links,
                                 directed = TRUE, damping = 0.85,
                                 personalized = NULL, weights = NULL, options = NULL)

page_rank_i_nodes_intra_links <-
  data.frame(species=names(page_rank_i[[1]]),
             Real_PR_Layer=page_rank_i[[1]], row.names=NULL)

# Add extra nodes----

# Number of extra nodes

graph_Plot_i_intra_links_extra_nodes <- graph_Plot_i_intra_links %>%
  igraph::add_vertices(isolated_nodes, color = "red")

page_rank_i_extra_nodes <- igraph::page_rank(graph_Plot_i_intra_links_extra_nodes,
                                             directed = TRUE, damping = 0.85,
                                             personalized = NULL, weights = NULL, options = NULL)

page_rank_i_nodes_intra_links_extra_nodes <-
  data.frame(species=names(page_rank_i_extra_nodes[[1]]),
             Real_PR_Layer_extra=page_rank_i_extra_nodes[[1]], row.names=NULL)

# Compare
page_rank_i_nodes_intra_links <- page_rank_i_nodes_intra_links %>%
  left_join(page_rank_i_nodes_intra_links_extra_nodes,by="species")

# Ratio for extra nodes is 1!!!

pagee_rank_plot <- page_rank_i_nodes_all_links %>% 
  left_join(page_rank_i_nodes_intra_links,by="species") %>%
  mutate(Ratio=Real_PR_Multi/Real_PR_Layer,
         Ratio_extra=Real_PR_Multi_extra/Real_PR_Layer_extra,
         Diff=Ratio-Ratio_extra,
         Multi_ratio=Real_PR_Multi/Real_PR_Multi_extra,
         Layer_ratio=Real_PR_Layer/Real_PR_Layer_extra)

# Ratio is preserved


