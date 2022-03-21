#####################################
# Data to generate plant phenologies 

Mean_of_means <- 16
sd_of_means <- 5
mean_phenology_plant <- round(rnorm(5, Mean_of_means, sd_of_means),0) %>% sort()


mean_phenology_plant_sp_1 <- mean_phenology_plant[1]
sd_phenology_plant_sp_1 <- 2

mean_phenology_plant_sp_2 <- mean_phenology_plant[2]
sd_phenology_plant_sp_2 <- 2

mean_phenology_plant_sp_3 <- mean_phenology_plant[3]
sd_phenology_plant_sp_3 <- 2

mean_phenology_plant_sp_4 <- mean_phenology_plant[4]
sd_phenology_plant_sp_4 <- 2

mean_phenology_plant_sp_5 <- mean_phenology_plant[5]
sd_phenology_plant_sp_5 <- 2