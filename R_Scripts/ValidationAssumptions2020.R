d <- read.csv("Raw_Data/raw_Pollinators_2020_1.csv", sep = ";")
head(d)
f <- read.csv("Raw_Data/Fitness_2020.csv", sep = ";")
head(f)

#mean visits of each guild per plant

library(reshape2)

d2 <- dcast(d, Plot + Subplot + Plant ~ Group, fun.aggregate = sum, value.var = "Visits")
head(d2)
f2 <- dcast(f, Plot + Subplot + Plant ~ ., fun.aggregate = mean, value.var = "Seeds.Fruit")
colnames(f2)[4] <- "Seeds"
head(f2)
d2$link <- paste0(d2$Plot,d2$Subplot,d2$Plant)
f2$link <- paste0(f2$Plot,f2$Subplot,f2$Plant)
unique(d2$link)
unique(f2$link)
dim(d2);dim(f2)
dat <-  merge(f2,d2, by = "link")
dim(dat) #only half match!!!

head(dat)
#all
scatter.smooth(dat$Seeds ~ dat$Bee)
scatter.smooth(dat$Seeds ~ dat$Beetle)
scatter.smooth(dat$Seeds ~ dat$Butterfly)
scatter.smooth(dat$Seeds ~ dat$Fly)

#CHUFU
chfu <- subset(dat, Plant.x == "CHFU")
scatter.smooth(chfu$Seeds ~ chfu$Bee) #sparse, but good
scatter.smooth(chfu$Seeds ~ chfu$Beetle) #no pattern
#scatter.smooth(chfu$Seeds ~ chfu$Butterfly)
scatter.smooth(chfu$Seeds ~ chfu$Fly) #no pattern

#LEMA
lema <- subset(dat, Plant.x == "LEMA")
scatter.smooth(lema$Seeds ~ lema$Bee) # no pattert (outlyer
scatter.smooth(lema$Seeds ~ lema$Beetle) #largelly positive!
#scatter.smooth(lema$Seeds ~ lema$Butterfly)
scatter.smooth(lema$Seeds ~ lema$Fly) #slightly positive?

#PUPA
pupa <- subset(dat, Plant.x == "PUPA")
scatter.smooth(pupa$Seeds ~ pupa$Bee) # positive, despite outlyer
#scatter.smooth(pupa$Seeds ~ pupa$Beetle) 
scatter.smooth(pupa$Seeds ~ pupa$Butterfly) #sparse
scatter.smooth(pupa$Seeds ~ pupa$Fly) #sparse

#CETE
cete <- subset(dat, Plant.x == "CETE")
scatter.smooth(cete$Seeds ~ cete$Bee) # zero seeds?
#scatter.smooth(cete$Seeds ~ cete$Beetle) 
#scatter.smooth(cete$Seeds ~ cete$Butterfly) #sparse
#scatter.smooth(cete$Seeds ~ cete$Fly) #sparse

#CHMI
chmi <- subset(dat, Plant.x == "CHMI")
#scatter.smooth(chmi$Seeds ~ chmi$Bee) # positive, despite outlyer
#scatter.smooth(chmi$Seeds ~ chmi$Beetle) 
#scatter.smooth(chmi$Seeds ~ chmi$Butterfly) #sparse
scatter.smooth(chmi$Seeds ~ chmi$Fly) #sparse

#SPRU
spru <- subset(dat, Plant.x == "SPRU")
#scatter.smooth(spru$Seeds ~ spru$Bee) # positive, despite outlyer
#scatter.smooth(spru$Seeds ~ spru$Beetle) 
#scatter.smooth(spru$Seeds ~ spru$Butterfly) #sparse
#scatter.smooth(spru$Seeds ~ spru$Fly) #sparse

#all
scatter.smooth(dat$Seeds ~ dat$Bee)
scatter.smooth(dat$Seeds ~ dat$Beetle)
scatter.smooth(dat$Seeds ~ dat$Butterfly)
scatter.smooth(dat$Seeds ~ dat$Fly)
