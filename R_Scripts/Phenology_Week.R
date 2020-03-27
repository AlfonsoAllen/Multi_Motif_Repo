# load libraries
library(tidyverse)


####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################

pollination <- read_csv("Raw_data/Metadata_Pollinators_Abundances_Seeds_2019.csv")

pollination_19 <- pollination %>% filter(Year==2019,Subplot!="OUT")


####################################################################
# Estimating phenology from pollinators'visits
####################################################################

pollination_dates_2019 <- pollination_19 %>% select(Plot,Plant_Simple,Month,Day) %>%
  group_by(Plot,Plant_Simple,Month,Day) %>%
  count() %>% mutate(Week=NA)

for (i in 1:nrow(pollination_dates_2019)){
  
  date_raw <- paste(pollination_dates_2019$Day[i],pollination_dates_2019$Month[i],"2019",sep="/")
  fecha <- as.Date(date_raw, "%d/%m/%Y")
  pollination_dates_2019$Week[i] <- format(fecha, "%V")
}

write_csv(pollination_dates_2019, "Processed_data/Phenology/phenology_2019.csv")

#Phenology by plot

pdf("Processed_data/Phenology/estimated_phenology_plot.pdf",
    width = 11.69, # The width of the plot in inches
    height = 8.27)

ggplot(pollination_dates_2019)+
  geom_point(aes(x=Week,y=Plant_Simple,size=n,color=n))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  labs(color = "#visits",size = "#visits")

dev.off()

#Phenology: CARACOLES
pdf("Processed_data/Phenology/estimated_phenology_Caracoles.pdf",
    width = 11.69, # The width of the plot in inches
    height = 8.27)

ggplot(pollination_dates_2019)+
  geom_point(aes(x=Week,y=Plant_Simple,size=n,color=n))+
  labs(color = "#visits",size = "#visits")
dev.off()
