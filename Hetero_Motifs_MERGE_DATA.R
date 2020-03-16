# load libraries
library(tidyverse)


####################################################################
# Loadind datasets (Caracoles) for 2019: Plot	G_F	Subplot_Plant_Label,
# Visits_tot,	hetero_motif,	hete_motif

####################################################################

data_real <- read_csv("Motifs/Caracoles.csv")

data_real <- data_real %>% mutate(Subplot=NA,Plant_Label=NA)

for(i in 1:nrow(data_real)){
  
  x <- strsplit(data_real$Subplot_Plant_Label[i], " ")
  data_real$Subplot[i] <- x[[1]][1]
  data_real$Plant_Label[i] <- x[[1]][2]

  }

GF_list <- data_real %>% group_by(G_F) %>% count() %>% select(G_F) %>% ungroup()
GF_list <- GF_list[[1]]
Subplot_list <- data_real %>% group_by(Subplot) %>% count() %>% select(Subplot)%>% ungroup()
Subplot_list <- Subplot_list[[1]]
Plant_list <- data_real %>% group_by(Plant_Label) %>% count() %>% select(Plant_Label)%>% ungroup()
Plant_list <- Plant_list[[1]]
tot_rows <- 9*nrow(Subplot_list)*nrow(GF_list)*nrow(Plant_list)


hetero_table <- as_tibble(expand.grid(Plot=1:9,Subplot=Subplot_list,
                                    G_F=GF_list,Plant_Label=Plant_list,
                                    stringsAsFactors = FALSE)) %>%
  arrange(Plot,Subplot,G_F)%>%
  mutate(motifs_Caracoles=NA)

for(i in 1:nrow(hetero_table)){
  
  h_m_i <- data_real %>% filter(Plot==hetero_table$Plot[i],
                                G_F==hetero_table$G_F[i],
                                Subplot_Plant_Label==paste(hetero_table$Subplot[i],
                                                           hetero_table$Plant_Label[i],sep=" "))
  
  if (length(h_m_i$hetero_motif)!=0){
    hetero_table$motifs_Caracoles[i] <- as.numeric(h_m_i$hetero_motif)}else{
      hetero_table$motifs_Caracoles[i] <- 0
    }
  
}




#########################################################################################
# Function: autofill_hetero_motifs
# INPUTS: a) list_rd_files: vector with the names of the files with random data
# b) data_table: Table/tibble that will be completed with the information in each rd_file
# OUTPUT: updated data_table
#########################################################################################

autofill_hetero_motifs <- function(list_rd_files,data_table){
  
  #data_table_or <- hetero_table #TEST
  
  #We removed from our list of file names: "Caracoles.csv" (it contains the real data)
  files_simpl <- list_rd_files[!list_rd_files %in% "Caracoles.csv"] 
  
  #files_simpl <- files_simpl[1:2]  #TEST
  
  initial_number_columns <- ncol(data_table)
  
  # New columns are added to data_table (one for each rd_file)
  
  columnsToAdd <- paste("rd_", 1:length(files_simpl),sep="") 
  data_table[,columnsToAdd] <- NA
  
  for(i in 1:length(files_simpl)){

    #i= 1  #TEST
    print(files_simpl[i])

    data_read <- read_csv(paste("Motifs/",files_simpl[i],sep=""))
    
    #data_read <- unique(data_read)  #TEST
    cont =0
    for (j in 1:nrow(data_table)){
      
      #j=19  #TEST
      #print(j)
    
      h_m_i <- data_read %>% filter(Plot==data_table$Plot[j],
                                    G_F==data_table$G_F[j],
                                    Subplot_Plant_Label==paste(data_table$Subplot[j],
                                                               data_table$Plant_Label[j],
                                                               sep=" "))
      
      
      if (length(h_m_i$hetero_motif)!=0){
        data_table[j,initial_number_columns+i] <- as.numeric(h_m_i$hetero_motif)
        cont=cont+1}else{
          data_table[j,initial_number_columns+i] <- 0
          
        }
      
      }
    
  }
  return(data_table)
}

################################################

list_rd_files <- list.files(paste(getwd(),"/Motifs",sep="")) 
hetero_data_subplot <- autofill_hetero_motifs(list_rd_files,hetero_table)

write_csv(hetero_data_subplot, paste("Total_hetero_data_subplot.csv",sep = ""))
