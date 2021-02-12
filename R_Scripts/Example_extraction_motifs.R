# load libraries
library(tidyverse)
library(bipartite)
library(igraph)

sessionInfo()

source("R_Scripts/functions.R")

####################################################################
# Loadind example dataset

example_data <- read_csv("Example/example.csv")

# 
aggregate_total <- NULL

for (week_i in unique(example_data$Week)){

  print("WEEK")
  print(week_i)
  
  
 example_week_i <- example_data %>% filter(Week==week_i)
 
 # Aggregate visits by week 
 
 example_week_i <- example_week_i %>% group_by(Plot,Subplot,Plant,ID) %>%
    count(wt=Visits) %>% rename(Visits_tot = n)
  
 example_week_i$Subplot_Plant_Label <- paste(example_week_i$Subplot,example_week_i$Plant,sep = " ")
  
  aggregate_week_i <- example_week_i %>% ungroup() %>% 
    select(Plot,ID,Subplot_Plant_Label,Visits_tot)
  
  
  
  aggregate_week_i <- homo_hete_motifs(aggregate_week_i)
  aggregate_week_i <- aggregate_week_i %>% mutate(Week=week_i)
  
  aggregate_total <-  bind_rows(aggregate_total, aggregate_week_i) 
  
}

aggregate_total %>% separate(Subplot_Plant_Label,c("Subplot","Plant"), " ")
