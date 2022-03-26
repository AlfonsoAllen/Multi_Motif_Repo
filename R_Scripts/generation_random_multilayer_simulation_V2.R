
library(tidyverse)

source("R_Scripts/aux_functions_for_simulations_V2.R")
source("R_Scripts/data2generate_plant_phenologies_simulations.R")


##################
# Create a matrix

number_plant_sp <- 5 # Do not change
number_intralinks <- 100
number_interlinks <- 30 # Max: number_pollinator_sp x 10 # Min: 0
number_plant_individuals <- 50
number_pollinator_sp <- 10
average_intralink_strenght <- 5

#------------------------------------------------

interaction_list_with_phen <- generate_interaction_list_with_phen(number_plant_individuals,
                                                                  number_pollinator_sp,
                                                                  number_intralinks,
                                                                  average_intralink_strenght,
                                                                  number_interlinks,
                                                                  mean_phenology_plant_sp_1,
                                                                  sd_phenology_plant_sp_1,
                                                                  mean_phenology_plant_sp_2,
                                                                  sd_phenology_plant_sp_2,
                                                                  mean_phenology_plant_sp_3,
                                                                  sd_phenology_plant_sp_3,
                                                                  mean_phenology_plant_sp_4,
                                                                  sd_phenology_plant_sp_4,
                                                                  mean_phenology_plant_sp_5,
                                                                  sd_phenology_plant_sp_5)


repetitions <- 100

# SIMULATIONS CHANGING INTRALINKS----------------

# Variation in the number of intralinks while keeping constant the number_interlinks and
# total number of plant individuals

list_number_intralinks <- c(100,200,300)
number_interlinks <- 20

#The total number of plant individuals should be between (intralinks - interlinks) and intralinks
# number_pollinator_sp <- number_intralinks - number_interlinks

number_plant_individuals <- 50
number_pollinator_sp <- 25

set.seed(121)

data_changing_intralinks <- NULL

for(number_intralinks in list_number_intralinks){
  
  for(rep in 1:repetitions){
    
    source("R_Scripts/data2generate_plant_phenologies_simulations.R")
    interaction_list_with_phen <- generate_interaction_list_with_phen(number_plant_individuals,
                                                                      number_pollinator_sp,
                                                                      number_intralinks,
                                                                      average_intralink_strenght,
                                                                      number_interlinks,
                                                                      mean_phenology_plant_sp_1,
                                                                      sd_phenology_plant_sp_1,
                                                                      mean_phenology_plant_sp_2,
                                                                      sd_phenology_plant_sp_2,
                                                                      mean_phenology_plant_sp_3,
                                                                      sd_phenology_plant_sp_3,
                                                                      mean_phenology_plant_sp_4,
                                                                      sd_phenology_plant_sp_4,
                                                                      mean_phenology_plant_sp_5,
                                                                      sd_phenology_plant_sp_5)
    
    interaction_list_with_phen$Plot <- paste0(rep,"_",number_intralinks,"_intralinks",
                                              number_interlinks,"_interlinks",
                                              number_plant_individuals,"_individuals")
    
    interaction_list_with_phen <- interaction_list_with_phen[,c("Plot","Subplot","Plant","ID",
                                                                "Visits","Week")]
    
    data_changing_intralinks <- bind_rows(data_changing_intralinks,interaction_list_with_phen)
  }
  
}

# Commented for security reasons
# write_csv(data_changing_intralinks,"Processed_data/Data_simulation/data_changing_intralinks_V2.csv")
                 


                          
# SIMULATIONS CHANGING INTERLINKS----------------

# Variation in the number of interlinks while keeping constant the number_intralinks and
# total number of plant individuals

list_number_interlinks <- c(5,20,60)
number_intralinks <- 200
number_plant_individuals <- 50
number_pollinator_sp <- 25

set.seed(121)

data_changing_interlinks <- NULL

for(number_interlinks in list_number_interlinks){
  
  for(rep in 1:repetitions){
    
    source("R_Scripts/data2generate_plant_phenologies_simulations.R")
    interaction_list_with_phen <- generate_interaction_list_with_phen(number_plant_individuals,
                                                                      number_pollinator_sp,
                                                                      number_intralinks,
                                                                      average_intralink_strenght,
                                                                      number_interlinks,
                                                                      mean_phenology_plant_sp_1,
                                                                      sd_phenology_plant_sp_1,
                                                                      mean_phenology_plant_sp_2,
                                                                      sd_phenology_plant_sp_2,
                                                                      mean_phenology_plant_sp_3,
                                                                      sd_phenology_plant_sp_3,
                                                                      mean_phenology_plant_sp_4,
                                                                      sd_phenology_plant_sp_4,
                                                                      mean_phenology_plant_sp_5,
                                                                      sd_phenology_plant_sp_5)
    
    interaction_list_with_phen$Plot <- paste0(rep,"_",number_intralinks,"_intralinks",
                                              number_interlinks,"_interlinks",
                                              number_plant_individuals,"_individuals")
    
    interaction_list_with_phen <- interaction_list_with_phen[,c("Plot","Subplot","Plant","ID",
                                                                "Visits","Week")]
    
    data_changing_interlinks <- bind_rows(data_changing_interlinks,interaction_list_with_phen)
  }
  
}

# Commented for security reasons
# write_csv(data_changing_interlinks,"Processed_data/Data_simulation/data_changing_interlinks_V2.csv")

# SIMULATIONS CHANGING PLANT INDIVIDUALS--------------

number_interlinks <- 20
number_intralinks <- 200
number_pollinator_sp <- 25

list_number_plant_individuals <- c(15,50,100)

set.seed(121)

data_changing_individuals <- NULL

for(number_plant_individuals in list_number_plant_individuals){
  
  for(rep in 1:repetitions){
    
    source("R_Scripts/data2generate_plant_phenologies_simulations.R")
    interaction_list_with_phen <- generate_interaction_list_with_phen(number_plant_individuals,
                                                                      number_pollinator_sp,
                                                                      number_intralinks,
                                                                      average_intralink_strenght,
                                                                      number_interlinks,
                                                                      mean_phenology_plant_sp_1,
                                                                      sd_phenology_plant_sp_1,
                                                                      mean_phenology_plant_sp_2,
                                                                      sd_phenology_plant_sp_2,
                                                                      mean_phenology_plant_sp_3,
                                                                      sd_phenology_plant_sp_3,
                                                                      mean_phenology_plant_sp_4,
                                                                      sd_phenology_plant_sp_4,
                                                                      mean_phenology_plant_sp_5,
                                                                      sd_phenology_plant_sp_5)
    
    interaction_list_with_phen$Plot <- paste0(rep,"_",number_intralinks,"_intralinks",
                                              number_interlinks,"_interlinks",
                                              number_plant_individuals,"_individuals")
    
    interaction_list_with_phen <- interaction_list_with_phen[,c("Plot","Subplot","Plant","ID",
                                                                "Visits","Week")]
    
    data_changing_individuals <- bind_rows(data_changing_individuals,interaction_list_with_phen)
  }
  
}

# Commented for security reasons
# write_csv(data_changing_individuals,"Processed_data/Data_simulation/data_changing_individuals_V2.csv")

# SIMULATIONS CHANGING POLLINATOR SPECIES--------------

number_interlinks <- 20
number_intralinks <- 200
number_plant_individuals <- 50

list_number_pollinator_sp <- c(5,25,75)

set.seed(112)

data_changing_pollinators <- NULL

for(number_pollinator_sp in list_number_pollinator_sp){
  
  for(rep in 1:repetitions){
    
    source("R_Scripts/data2generate_plant_phenologies_simulations.R")
    interaction_list_with_phen <- generate_interaction_list_with_phen(number_plant_individuals,
                                                                      number_pollinator_sp,
                                                                      number_intralinks,
                                                                      average_intralink_strenght,
                                                                      number_interlinks,
                                                                      mean_phenology_plant_sp_1,
                                                                      sd_phenology_plant_sp_1,
                                                                      mean_phenology_plant_sp_2,
                                                                      sd_phenology_plant_sp_2,
                                                                      mean_phenology_plant_sp_3,
                                                                      sd_phenology_plant_sp_3,
                                                                      mean_phenology_plant_sp_4,
                                                                      sd_phenology_plant_sp_4,
                                                                      mean_phenology_plant_sp_5,
                                                                      sd_phenology_plant_sp_5)
    
    interaction_list_with_phen$Plot <- paste0(rep,"_",number_intralinks,"_intralinks",
                                              number_interlinks,"_interlinks",
                                              number_plant_individuals,"_individuals",
                                              number_pollinator_sp,"_pollinators")
    
    interaction_list_with_phen <- interaction_list_with_phen[,c("Plot","Subplot","Plant","ID",
                                                                "Visits","Week")]
    
    data_changing_pollinators <- bind_rows(data_changing_pollinators,interaction_list_with_phen)
  }
  
}

# Commented for security reasons
# write_csv(data_changing_pollinators,"Processed_data/Data_simulation/data_changing_pollinators_V2.csv")

