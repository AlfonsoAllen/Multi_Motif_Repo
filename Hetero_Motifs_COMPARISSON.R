# load libraries
library(tidyverse)

#####################################
# CONFIDENCE AND PREDICTION INTERVALS
#####################################

# Loading data for Caracoles and random samples

Homo_data_subplot <- read_csv("Total_Homo_data_subplot.csv")

Prediction_homo_subplot <- Homo_data_subplot %>%
  select(Plot,Subplot,G_F,Plant_Label,motifs_Caracoles) %>%
  mutate(Mean_PI=NA,lwr_PI=NA,upr_PI=NA,Caracoles_In=NA)


for (i in 1:nrow(Homo_data_subplot)){

  sample_i <- gather(Homo_data_subplot[i,6:ncol(data_table)], key = "sample", value = "homo")

  PI<-predict(lm(sample_i$homo~ 1), interval="predict")
  #PI[1,]
  
  Prediction_homo_subplot$Mean_PI[i] <- PI[1,1]
  Prediction_homo_subplot$lwr_PI[i] <- PI[1,2]
  Prediction_homo_subplot$upr_PI[i] <- PI[1,3]
  Prediction_homo_subplot$Caracoles_In[i] <- (PI[1,2] <= Homo_data_subplot[i,5])&
    (Homo_data_subplot[i,5] <= PI[1,3])
  

}

Prediction_homo_subplot %>% filter(Caracoles_In==FALSE,motifs_Caracoles==0)
Prediction_homo_subplot %>% filter(Caracoles_In==FALSE,motifs_Caracoles>0)

###############################################3
# It seems that homotifs per (plot,subplot) do not follow a normal distribution
################################################
i=13
sample_i <- gather(Homo_data_subplot[i,6:ncol(data_table)], key = "sample", value = "homo")

PI<-predict(lm(sample_i$homo~ 1), interval="predict")
x.test <- shapiro.test(sample_i$homo)
print(x.test)

plotn <- function(x,main="Histograma de frecuencias \ny distribución normal",
                  xlab="X",ylab="Densidad") {
  min <- min(sample_i$homo)
  max <- max(sample_i$homo)
  media <- mean(sample_i$homo)
  dt <- sd(sample_i$homo)
  hist(sample_i$homo,freq=F,main=main,xlab=xlab,ylab=ylab)
  curve(dnorm(x,media,dt), min, max,add = T,col="blue")
}

plotn(sample_i$homo,main="Distribución normal")#Grafico de x
