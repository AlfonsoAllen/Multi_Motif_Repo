# load libraries
library(tidyverse)


####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################

metadata <- read_csv2("Raw_Data/Metadata_Pollinators_Plants_fruits_abundancesmodified.csv")

metadata_19 <- metadata %>% filter(Year==2019,Subplot!="OUT",!is.na(ID)) %>% #QUITAR OUTS
  select(Day,Month,Year,Plot,Subplot,G_F,Order,Superfamily, Family, Species, ID,Plant_Simple,Visits,num.plantas,Fruit,Seed)



metadata_19 <- metadata_19 %>% filter(!ID%in%c("Formicidae","Mirini_sp.","Chrysididae",
                                            "Diplazon_sp.","Culicidae","Melanopangonius_sp.",
                                            "Nephrotoma_sp.","Larva"))

abundant_species <- c("Lasiocampa_trifolii", "Psilothrix_viridicoerulea",
                      "Eupeodes_corollae","Sphaerophoria_scripta","Nemotelus_sp.",
                      "Dilophus_sp.","Meligethes_sp.","Odontomyia_sp.","Scaeva_sp.")

for (i in 1:nrow(metadata_19)){
  
  if (metadata_19$Species[i] %in% abundant_species){
    metadata_19$ID[i] <- metadata_19$Species[i]
  }else if (!is.na(metadata_19$Family[i])){
    metadata_19$ID[i] <- metadata_19$Family[i]
  }else{
    metadata_19$ID[i] <- NA
  }
}


metadata_19 %>% filter(is.na(ID)) %>% count()
metadata_19 %>% filter(ID=="Na") %>% count()
metadata_19$ID[metadata_19$ID=="Na"] <- NA

write_csv(metadata_19,"Raw_Data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")
