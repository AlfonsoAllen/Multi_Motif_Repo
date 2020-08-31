# load libraries
library(tidyverse)


####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################

pollination <- read_csv2("Raw_data/raw_Pollinators_2020_1.csv")

pollination_20 <- pollination %>%  filter(!is.na(Plant),Plant!="0",Subplot!="OUT",Plant!="Ground")


####################################################################
# Estimating phenology from pollinators'visits
####################################################################

pollination_dates_2020 <- pollination_20 %>% select(Plot,Plant,Month,Day) %>%
  group_by(Plot,Plant,Month,Day) %>%
  count() %>% mutate(Week=NA,Line=NA)

for (i in 1:nrow(pollination_dates_2020)){
  
  date_raw <- paste(pollination_dates_2020$Day[i],pollination_dates_2020$Month[i],"2020",sep="/")
  fecha <- as.Date(date_raw, "%d/%m/%Y")
  pollination_dates_2020$Week[i] <- format(fecha, "%V")
  
  if(pollination_dates_2020$Plot[i] %in% c(1,2,3)){pollination_dates_2020$Line[i] <- 1}
  else if(pollination_dates_2020$Plot[i] %in% c(4,5,6)){pollination_dates_2020$Line[i] <- 2}
  else{pollination_dates_2020$Line[i] <- 3}
  
}

write_csv(pollination_dates_2020, "Processed_data/Phenology/phenology_2020.csv")

#Phenology by plot

plot_labs <-c(
  "Plot 1",
  "Plot 2",
  "Plot 3",
  "Plot 4",
  "Plot 5",
  "Plot 6",
  "Plot 7",
  "Plot 8",
  "Plot 9"
)
names(plot_labs) <- c(
  '1',
  '2',
  '3',
  '4',
  '5',
  "6",
  '7',
  '8',
  "9"
)


pdf("Processed_data/Phenology/estimated_phenology_plot_2020.pdf",
    width = 11.69, # The width of the plot in inches
    height = 8.27)

ggplot(pollination_dates_2020)+
  geom_point(aes(x=Week,y=Plant,size=n,color=n))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  ylab("Plant Species") +
  labs(color = "#visits",size = "#visits")+ theme(legend.position="bottom")+
  theme_bw()

dev.off()

pdf("Processed_data/Phenology/estimated_phenology_line.pdf",
    width = 11.69, # The width of the plot in inches
    height = 8.27)

ggplot(pollination_dates_2020)+
  geom_point(aes(x=Week,y=Plant,size=n,color=n))+
  facet_wrap(vars(Line),nrow = 3,ncol = 3)+
  labs(color = "#visits",size = "#visits")

dev.off()

#Phenology: CARACOLES
pdf("Processed_data/Phenology/estimated_phenology_Caracoles.pdf",
    width = 11.69, # The width of the plot in inches
    height = 8.27)

ggplot(pollination_dates_2020)+
  geom_point(aes(x=Week,y=Plant,size=n,color=n))+
  labs(color = "#visits",size = "#visits")
dev.off()



