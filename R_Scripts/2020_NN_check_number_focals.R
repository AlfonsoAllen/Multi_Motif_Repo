
# Check that the number of focals in LEMA is correct

pol_data <- read_csv2("Raw_Data/final_Pollinators_2020.csv")
# Remove points from ID names
pol_data$ID <- sub("\\.", "", pol_data$ID)
pol_data$ID_Simple <- sub("\\.", "", pol_data$ID_Simple)

# filter tabanidae
pol_data <- pol_data %>% filter(ID != "Tabanidae")


pol_data_final <- pol_data %>% filter(!is.na(Plant),Plant!="0",Subplot!="OUT",Plant!="Ground")

pol_data_final <- pol_data_final %>% dplyr::select(-Order,-Family,-Genus,-ID) %>% rename(ID=ID_Simple) %>%
  mutate(date_raw=as.Date(paste(Day,Month,Year,sep="/"), "%d/%m/%Y"),
         Week=as.numeric(format(date_raw, "%V")))

pol_data_final %>% ungroup() %>% dplyr::select(Plot,Subplot,Plant) %>% unique() %>%
  group_by(Plant) %>% count()

competition <- read_csv2("Raw_Data/competition_caracoles2020.csv") %>%
  filter(!is.na(Seeds)) %>% 
  filter(Sp.Focal %in% c("BEMA","CETE","CHFU","CHMI","LEMA","MESU","PUPA",
                         "SCLA","SOAS","SPRU")) %>% rename(Plant = Sp.Focal)

competition%>% dplyr::select(Plot,Subplot,Plant) %>% unique() %>%
  group_by(Plant) %>% count()

focal_pollinator <- pol_data_final %>% ungroup() %>% dplyr::select(Plot,Subplot,Plant) %>% unique()

focal_competition <- competition%>% dplyr::select(Plot,Subplot,Plant) %>% unique()

dplyr::intersect(focal_pollinator, focal_competition) %>% filter(Plant == "LEMA") 
#225: number of LEMAS that are common in focal_pollinator and focal_competition

focal_pollinator %>% filter(Plant == "LEMA") #249
focal_competition %>% filter(Plant == "LEMA") #269

# 293 LEMAS 
249-225
24+269
