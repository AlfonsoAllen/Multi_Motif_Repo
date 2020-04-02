#This scripts introduce multivariate analysis
library(tidyverse)
library(vegan)

#NMDS
#The data

abundances <- read_csv2("Raw_Data/abundances_2019.csv")
abundances_19 <- abundances %>% filter(year==2019) %>% 
  select(plot,subplot,species,individuals) %>% arrange(plot,subplot)

abundances_19 %>% group_by(plot) %>% count(wt = individuals)
# Plot 7, 8, 9 have almost one third of plant individuals

Caracoles <- abundances_19 %>% spread(species,individuals)
head(Caracoles)
str(Caracoles)
summary(Caracoles)
hist(Caracoles$plot,breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5))


#simplify objects to use
Plot <- Caracoles$plot
Subplot <- Caracoles$subplot
#select the community
Car_community <- Caracoles[,3:ncol(Caracoles)]

#cargar librerias----
#install.packages("vegan") #si no tienes la libreria, instalala primero
library("vegan")

#The basic is the distance measure you use:
#distance selection!----
?dist
?vegdist
?betadiver

v <- vegdist(Car_community, "horn")
head(v, 20)
Car_community[c(1,13),]
vegdist(Car_community[c(1,13),], "horn")

#A few words:
#binary:
#Jackard
#Sorensen (This coefficient weights matches in species composition
#between the two samples more heavily than mismatches)
#quantitative
#Euclidian: simple distance, good for e.g. distance between sites
#bray: 0-1 The Bray-Curtis measure ignores cases in which the species
#is absent in both community samples, and it is dominated
#by the abundant species so that rare species add very little to the
#value of the coefficient
#morisita: 0-1 independent of sample size. Only for counts. Recomended.
#kulczynski: Weigth more rare species.
#gower (allows factors!)

#Best is Legendre book numerical ecology
#(only found Krebs online): http://www.zoology.ubc.ca/~krebs/downloads/krebs_chapter_12_2014.pdf

# NMDS
# Addtional Function Documentation: https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/metaMDS

set.seed(123)
Car_community.mds <- metaMDS(comm = Car_community,
                             distance = "horn", trymax=100,maxit=1000,noshare=T,k=2,
                             trace = TRUE, autotransform = FALSE)

#Rerun the following command till convergence is achieved
Car_community.mds <- metaMDS(comm = Car_community,previous.best=Car_community.mds,
                              distance = "horn", trymax=100,maxit=5000,noshare=T,k=2,
                              trace = TRUE, autotransform = FALSE)


#Herb_community.mds <- metaMDS(v) #idem
plot(Car_community.mds$points, col = Plot, pch = 16)


MDS_test <- tibble(MDS1=Car_community.mds$points[,1],
       MDS2=Car_community.mds$points[,2],
       Plot=Caracoles$plot,
       Subplot=Caracoles$subplot)

ggplot(MDS_test)+
  geom_point(aes( MDS1,MDS2,color=as.factor(Subplot)))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  labs(title="NMDS Caracoles' abundances for 2019",
       x ="MDS1", y = "MDS2",color="Subplot")

###########################################
###########################################

set.seed(123)
Car_community.mds_3D <- metaMDS(comm = Car_community,
                             distance = "horn", trymax=100,maxit=1000,noshare=T,k=3,
                             trace = TRUE, autotransform = FALSE)


#devtools::install_github("AckerDWM/gg3D")

MDS_test_3D <- tibble(MDS1=Car_community.mds_3D$points[,1],
                   MDS2=Car_community.mds_3D$points[,2],
                   MDS3=Car_community.mds_3D$points[,3],
                   Plot=Caracoles$plot,
                   Subplot=Caracoles$subplot)

library("gg3D")
# http://htmlpreview.github.io/?https://github.com/AckerDWM/gg3D/blob/master/gg3D-vignette.html
theta=0 
phi=20

ggplot(MDS_test_3D,aes(x=MDS1,y=MDS2,z=MDS3,color=as.factor(Subplot))) + 
  geom_point() +
  axes_3D() +
  stat_3D()+
  labs_3D(labs=c("MDS1", "MDS2", "MDS3"))+
facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  labs(title="NMDS Caracoles' abundances for 2019",color="Subplot")+
  theme_void()

#######################################
#assumptions:
Car_community.mds$stress

#If the stress value is greater than 0.2, it is advisable to include an
#additional dimension, but remember that human brains are not very well
#equipped to visualise objects in more than 2-dimensions.

#Transformation and standardisation. Transforming data sets prior to
#creating a MDS plot is often desirable, not to meet assumptions of normality,
#but to reduce the influence of extreme values. For example,


Car_community.sq <- sqrt(Car_community)
set.seed(123)
Car_community.sq.mds <- metaMDS(comm = Car_community.sq,
                                distance = "horn", trymax=100,maxit=1000,k=3,
                                noshare=T,trace = 1,
                                autotransform = FALSE)
Car_community.sq.mds <- metaMDS(comm = Car_community.sq,previous.best=Car_community.sq.mds,
                                distance = "horn", trymax=100,maxit=1000,k=3,
                                noshare=T,trace = 1,
                                autotransform = FALSE)

MDS_test_3D_sq <- tibble(MDS1=Car_community.sq.mds$points[,1],
                      MDS2=Car_community.sq.mds$points[,2],
                      MDS3=Car_community.sq.mds$points[,3],
                      Plot=Caracoles$plot,
                      Subplot=Caracoles$subplot)

ggplot(MDS_test_3D_sq,aes(x=MDS1,y=MDS2,z=MDS3,color=as.factor(Subplot))) + 
  geom_point() +
  axes_3D() +
  stat_3D()+
  labs_3D(labs=c("MDS1", "MDS2", "MDS3"))+
  facet_wrap(vars(Plot),nrow = 3,ncol = 3)+
  labs(title="NMDS Caracoles' abundances for 2019",color="Subplot")+
  theme_void()

#alternative plot
ordiplot(Car_community.mds)
points(Car_community.mds$points, cex = 0.5,
       pch = 16,
       col = "darkgrey")
ordihull(Car_community.mds,groups=Plot,draw="polygon",col="grey90",
         label=FALSE)
orditorp(Car_community.mds,display="species",col="red",air=0.01)
ordiellipse(ord = Car_community.mds,
            groups = Plot,
            col = "red")


#Testing hypothesis PERMANOVA----
# PERMANOVA is used to compare groups of objects and test the
# null hypothesis that the centroids and dispersion of the groups
# as defined by measure space are equivalent for all groups. A rejection
# of the null hypothesis means that either the centroid and/or the spread
# of the objects is different between the groups.


#First test: Are centroids different?
adonis(Car_community ~ Plot, method = "horn")
#Second test: Is the spread different?
b <- betadisper(vegdist(Car_community, method = "horn"), group = Plot)
anova(b)
boxplot(b)

# ANOVA test, a significant p-value indicates that some of the group means
# are different, but we don’t know which pairs of groups are different.
# We can compute Tukey HSD (Tukey Honest Significant Differences,
# R function: TukeyHSD()) for performing multiple pairwise-comparison
# between the means of groups.

TukeyHSD(b)


