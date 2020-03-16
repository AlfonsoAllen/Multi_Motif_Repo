# load libraries
library(tidyverse)

#####################################
# CONFIDENCE AND PREDICTION INTERVALS
#####################################

# Loading data for Caracoles and random samples

hetero_data_subplot <- read_csv("Total_hetero_data_subplot.csv")

Prediction_hetero_subplot <- hetero_data_subplot %>%
  select(Plot,Subplot,G_F,Plant_Label,motifs_Caracoles) %>%
  mutate(Mean_PI=NA,lwr_PI=NA,upr_PI=NA,Caracoles_In=NA)


for (i in 1:nrow(hetero_data_subplot)){

  sample_i <- gather(hetero_data_subplot[i,6:ncol(hetero_data_subplot)], key = "sample", value = "hetero")

  PI<-predict(lm(sample_i$hetero~ 1), interval="predict")
  #PI[1,]
  
  Prediction_hetero_subplot$Mean_PI[i] <- PI[1,1]
  Prediction_hetero_subplot$lwr_PI[i] <- PI[1,2]
  Prediction_hetero_subplot$upr_PI[i] <- PI[1,3]
  Prediction_hetero_subplot$Caracoles_In[i] <- (PI[1,2] <= hetero_data_subplot[i,5])&
    (hetero_data_subplot[i,5] <= PI[1,3])
  

}

Prediction_hetero_subplot %>% filter(Caracoles_In==FALSE,motifs_Caracoles==0)
Prediction_hetero_subplot %>% filter(Caracoles_In==FALSE,motifs_Caracoles>0)

###############################################3
# It seems that heterotifs per (plot,subplot) do not follow a normal distribution 
################################################
i=11
sample_i <- gather(hetero_data_subplot[i,6:ncol(hetero_data_subplot)], key = "sample", value = "hetero")

PI<-predict(lm(sample_i$hetero~ 1), interval="predict")
x.test <- shapiro.test(sample_i$hetero)
print(x.test)

plotn <- function(x,main="Histograma de frecuencias \ny distribución normal",
                  xlab="X",ylab="Densidad") {
  min <- min(sample_i$hetero)
  max <- max(sample_i$hetero)
  media <- mean(sample_i$hetero)
  dt <- sd(sample_i$hetero)
  hist(sample_i$hetero,freq=F,main=main,xlab=xlab,ylab=ylab)
  curve(dnorm(x,media,dt), min, max,add = T,col="blue")
}

plotn(sample_i$hetero,main="Distribución normal")#Grafico de x
