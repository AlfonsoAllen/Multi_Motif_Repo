# load libraries
library(tidyverse)


####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################

metadata <- read_csv2("Raw_Data/Metadata_Pollinators_Plants_fruits_abundancesmodified.csv")
head(as.data.frame(metadata))

metadata_19 <- metadata %>% filter(Year==2019,Subplot!="OUT",!is.na(ID)) %>% #QUITAR OUTS
  select(Day,Month,Year,Plot,Subplot,G_F,Order,Superfamily, Family, Species, ID,Plant_Simple,Visits,num.plantas,Fruit,Seed)



metadata_19 <- metadata_19 %>% filter(!ID%in%c("Formicidae","Mirini_sp.","Chrysididae",
                                            "Diplazon_sp.","Culicidae","Melanopangonius_sp.",
                                            "Nephrotoma_sp.","Larva"))

sort(table(metadata_19$ID), decreasing = TRUE)
temp <- unique(metadata_19$ID)
temp[order(temp)]
#la idea es no perder ninguna visita de a) polinizadores o b) comedores de pollen (que no hay tantas)
#cosas que yo haria:
#"Apidae", "Apoidea" -> viendo que solo hay cuatro abejas, y una es muy conspicua yo asumaria que son Andrena sp.       
#"Apis mellifera" -> Habla con Maria, pero vimos en la colección que no lo era. 
#"Coccinella_septempunctata" -> Yo la excluiria, ya que no es pollinator, ni pollen feeder.
#"Chrysomelidae" -> check... comen pollen? son pollinators?
#"Coleoptera" -> dificil, pero a lo mejor el guild nos da pistas. a malas hay solo 5.
#"Diptera" y "Lepidoptera" -> idem. Diptera hay muchos! Ver comment en syrphidae.
#"Episyrphus_balteatus" Y  "Episyrphus_sp." Unificar. Unlikely there are two episyrphus species
#"Eupeodes_corollae"         "Eupeodes_sp." Unificar
#"Malachius_bipustulatus"   "Malachius_sp." Unificar
#"Hymenoptera" -> esto son No abejas? Other hymenoptera?
#"Musca_sp."  "Muscidae" Unificar
# "nemotelus_sp."             "Nemotelus_sp."   Unificar
#"Sphaerophoria_scripta"   "Sphaerophoria_sp." Unificar.
#"Syrphidae"-> de los más dificiles por que hay varias especies, pero no quiero perder 40 visitas. 
  #Seria una locura asignarlos a la especie que se haya pillado ese dia en esa planta?

metadata_19$Species == metadata_19$ID 

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
