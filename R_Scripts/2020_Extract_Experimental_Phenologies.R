# load libraries
library(tidyverse)
library(RColorBrewer)

####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################


pollination <- read_csv2("Raw_Data/final_Pollinators_2020.csv")
# Remove points from ID names
pollination$ID <- sub("\\.", "", pollination$ID)
pollination$ID_Simple <- sub("\\.", "", pollination$ID_Simple)

# filter tabanidae
pollination <- pollination %>% filter(ID != "Tabanidae")

pollination_20 <- pollination %>%  filter(!is.na(Plant),Plant!="0",Subplot!="OUT",Plant!="Ground")

####################################################################
# Estimating phenology from pollinators'visits
####################################################################

pollination_dates_2020 <- pollination_20 %>% dplyr::select(Plot,Plant,Month,Day) %>%
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

pollination_dates_2020 %>% ungroup() %>% select(Month,Day) %>% unique()

write_csv(pollination_dates_2020, "Processed_data/Plant_experimental_phenologies/phenology_2020.csv")

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

pollination_dates_2020$Plant[pollination_dates_2020$Plant == "BEMA"] <- "B. macrocarpa"
pollination_dates_2020$Plant[pollination_dates_2020$Plant == "CETE"] <- "C. tenuiflorum"
pollination_dates_2020$Plant[pollination_dates_2020$Plant == "CHFU"] <- "C. fuscatum"
pollination_dates_2020$Plant[pollination_dates_2020$Plant == "CHMI"] <- "C. mixtum"
pollination_dates_2020$Plant[pollination_dates_2020$Plant == "LEMA"] <- "L. maroccanus"
pollination_dates_2020$Plant[pollination_dates_2020$Plant == "MESU"] <- "M. sulcatus"
pollination_dates_2020$Plant[pollination_dates_2020$Plant == "PUPA"] <- "P. paludosa"
pollination_dates_2020$Plant[pollination_dates_2020$Plant == "SCLA"] <- "S. laciniata"
pollination_dates_2020$Plant[pollination_dates_2020$Plant == "SOAS"] <- "S. asper"
pollination_dates_2020$Plant[pollination_dates_2020$Plant == "SPRU"] <- "S. rubra"

pdf("Processed_data/Plant_experimental_phenologies/estimated_phenology_plot_2020.pdf",
    width = 11.69, # The width of the plot in inches
    height = 8.27)

ggplot(pollination_dates_2020)+
  geom_point(aes(x=Week,y=Plant,size=n,color=n))+
  scale_color_distiller(palette = "Spectral")+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3,labeller=labeller(Plot= plot_labs))+
  ylab("Plant Species") +
  theme_bw()+
  labs(color = "# visits",size = "# visits")+ theme(legend.position="bottom")+
  theme(axis.text.y = element_text(face = "italic"))
  

dev.off()

pdf("Processed_data/Plant_experimental_phenologies/estimated_phenology_line.pdf",
    width = 11.69, # The width of the plot in inches
    height = 8.27)

ggplot(pollination_dates_2020)+
  geom_point(aes(x=Week,y=Plant,size=n,color=n))+
  facet_wrap(vars(Line),nrow = 3,ncol = 3)+
  labs(color = "#visits",size = "#visits")

dev.off()

#Phenology: CARACOLES
pdf("Processed_data/Plant_experimental_phenologies/estimated_phenology_Caracoles.pdf",
    width = 11.69, # The width of the plot in inches
    height = 8.27)

pollination_dates_2020 %>% group_by(Plant,Week) %>% count(wt=n) %>%
ggplot()+
  geom_point(aes(x=Week,y=Plant,size=n,color=n))+
  scale_color_distiller(palette = "Spectral")+
  labs(color = "#visits",size = "#visits")+
  theme_bw()+
  labs(color = "# visits",size = "# visits")+ 
  theme(axis.text.y = element_text(face = "italic"))
dev.off()


#############################
# Phenology by plant and G_F

####################################################################
# Estimating phenology from pollinators'visits
####################################################################

pollination_dates_2020_GF <- pollination_20 %>% dplyr::select(Plot,Plant,Month,Day,G_F) %>%
  group_by(Plot,Plant,Month,Day,G_F) %>%
  count() %>% mutate(Week=NA,Line=NA)

for (i in 1:nrow(pollination_dates_2020_GF)){
  
  date_raw <- paste(pollination_dates_2020_GF$Day[i],pollination_dates_2020_GF$Month[i],"2020",sep="/")
  fecha <- as.Date(date_raw, "%d/%m/%Y")
  pollination_dates_2020_GF$Week[i] <- format(fecha, "%V")
  
  if(pollination_dates_2020_GF$Plot[i] %in% c(1,2,3)){pollination_dates_2020_GF$Line[i] <- 1}
  else if(pollination_dates_2020_GF$Plot[i] %in% c(4,5,6)){pollination_dates_2020_GF$Line[i] <- 2}
  else{pollination_dates_2020_GF$Line[i] <- 3}
  
}

pollination_dates_2020_GF$Plant[pollination_dates_2020_GF$Plant == "BEMA"] <- "B. macrocarpa"
pollination_dates_2020_GF$Plant[pollination_dates_2020_GF$Plant == "CETE"] <- "C. tenuiflorum"
pollination_dates_2020_GF$Plant[pollination_dates_2020_GF$Plant == "CHFU"] <- "C. fuscatum"
pollination_dates_2020_GF$Plant[pollination_dates_2020_GF$Plant == "CHMI"] <- "C. mixtum"
pollination_dates_2020_GF$Plant[pollination_dates_2020_GF$Plant == "LEMA"] <- "L. maroccanus"
pollination_dates_2020_GF$Plant[pollination_dates_2020_GF$Plant == "MESU"] <- "M. sulcatus"
pollination_dates_2020_GF$Plant[pollination_dates_2020_GF$Plant == "PUPA"] <- "P. paludosa"
pollination_dates_2020_GF$Plant[pollination_dates_2020_GF$Plant == "SCLA"] <- "S. laciniata"
pollination_dates_2020_GF$Plant[pollination_dates_2020_GF$Plant == "SOAS"] <- "S. asper"
pollination_dates_2020_GF$Plant[pollination_dates_2020_GF$Plant == "SPRU"] <- "S. rubra"


pdf("Processed_data/Plant_experimental_phenologies/estimated_phenology_Caracoles_by_GF.pdf",
    width = 11.69, # The width of the plot in inches
    height = 8.27)
pollination_dates_2020_GF %>% group_by(Plant,Week,G_F) %>% count(wt=n) %>%
  ggplot()+
  geom_point(aes(x=as.factor(as.numeric(Week)),y=Plant,size=n,color=n))+
  scale_color_distiller(palette = "Spectral")+
  labs(color = "#visits",size = "#visits", x = "Week")+
  facet_wrap(~G_F)+
  theme_bw()+
  labs(color = "# visits",size = "# visits")+ 
  theme(axis.text.y = element_text(face = "italic"))

dev.off()

# ####################################################################
# Estimating phenology from pollinators'visits by species (solitary bees, flower beetles and house flies)
####################################################################

pollination_dates_2020_sp <- pollination_20 %>% filter(G_F %in% c("Flower_beetles",
                                                                  "Small_flies",
                                                                  "Solitary_bees")) %>%
  dplyr::select(Plot,Plant,Month,Day,ID_Simple) %>%
  group_by(Plot,Plant,Month,Day,ID_Simple) %>%
  count() %>% mutate(Week=NA,Line=NA)

for (i in 1:nrow(pollination_dates_2020_sp)){
  
  date_raw <- paste(pollination_dates_2020_sp$Day[i],pollination_dates_2020_sp$Month[i],"2020",sep="/")
  fecha <- as.Date(date_raw, "%d/%m/%Y")
  pollination_dates_2020_sp$Week[i] <- format(fecha, "%V")
  
  if(pollination_dates_2020_sp$Plot[i] %in% c(1,2,3)){pollination_dates_2020_sp$Line[i] <- 1}
  else if(pollination_dates_2020_sp$Plot[i] %in% c(4,5,6)){pollination_dates_2020_sp$Line[i] <- 2}
  else{pollination_dates_2020_sp$Line[i] <- 3}
  
}

pdf("Processed_data/Plant_experimental_phenologies/estimated_phenology_Caracoles_by_sp.pdf",
    width = 11.69, # The width of the plot in inches
    height = 8.27)
pollination_dates_2020_sp %>% group_by(Plant,Week,ID_Simple) %>% count(wt=n) %>%
  ggplot()+
  geom_point(aes(x=Week,y=Plant,size=n,color=n))+
  scale_color_distiller(palette = "Spectral")+
  labs(color = "#visits",size = "#visits")+
  facet_wrap(~ID_Simple)+
  theme_bw()
dev.off()
