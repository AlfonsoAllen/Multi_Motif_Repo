# load libraries
library(tidyverse)

####################################################################
# Loadind Plant-pollinator dataset (Caracoles) for 2019: visits, abundances, seeds
####################################################################

fitness_data2 <- read_csv("Raw_data/Metadata_Pollinators_Abundances_Seeds_2019_ID.csv")

fitness2 <- fitness_data2 %>% filter(Year==2019)

fitness2 <- fitness2 %>% select(-Order,-Family,-Superfamily,-ID) %>% rename(ID=ID_Simple) %>%
  mutate(date_raw=as.Date(paste(Day,Month,Year,sep="/"), "%d/%m/%Y"),
         Week=as.numeric(format(date_raw, "%V")))

fitness2$Line <- NA

for (i in 1:nrow(fitness2 )){
  if(fitness2$Plot[i] %in% c(1,2,3)){fitness2$Line[i] <- 1}
  else if(fitness2$Plot[i] %in% c(4,5,6)){fitness2$Line[i] <- 2}
  else{fitness2$Line[i] <- 3}
} 

fitness3 <- fitness2 %>% group_by(Line,Plot,Subplot,Plant_Simple,Week,ID) %>%
  summarize(Seeds_GF = mean(Seed),
            Fruit_GF = mean(Fruit),
            visits_GF = sum(Visits))

#####################################
# Uploading motifs data
#####################################

caracoles_motif <- read_csv("Processed_data/Motifs_WEEK/Caracoles_WEEK_SPECIES.csv")

#Adding Subplot and Plant_Simple
for (i in 1:nrow(caracoles_motif)){
  x <- strsplit(caracoles_motif$Subplot_Plant_Label[i], " ")
  caracoles_motif$Subplot[i] <- x[[1]][1]
  caracoles_motif$Plant_Simple[i] <- x[[1]][2]
}

caracoles_motif <- caracoles_motif %>% filter(!Plant_Simple%in%c("HOMA","Lysimachia_arvensis")) %>%
  select(-Visits_tot,-Subplot_Plant_Label)

#####################################
# Merging motifs data and fitness
#####################################

fitness <- caracoles_motif %>% left_join(fitness3, by=c("Plot","Subplot","Plant_Simple","ID","Week"))

########################################################################
# SOME MOTIF VS SEEDS GRAPHS
########################################################################

fitness_SUM <- fitness %>% group_by(Line,Plot,Subplot,Plant_Simple,ID) %>%
  summarize(Seeds_GF = mean(Seeds_GF),
            Fruit_GF = mean(Fruit_GF),
            visits_GF = sum(visits_GF),
            Homo_Sum=sum(homo_motif),
            Hete_Sum=sum(hete_motif))

fitness_Spread <- fitness_SUM %>% mutate(ID_homo=paste0(ID,"_homo"),ID_hete=paste0(ID,"_hete")) %>%
  spread(ID_homo,Homo_Sum) %>%spread(ID_hete,Hete_Sum)

fitness_Spread[,9:ncol(fitness_Spread)][is.na(fitness_Spread[,9:ncol(fitness_Spread)])] <- 0
fitness_Spread$Homo_total <- rowSums(fitness_Spread[,9:45],na.rm=TRUE)
fitness_Spread$Hete_total <- rowSums(fitness_Spread[,46:82],na.rm=TRUE)
#############################################################################
#################################################################################
# PRINCIPAL COMPONENT ANALYSIS
# TUTORIAL: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
#################################################################################
#################################################################################

fitness_PCA <- fitness_Spread
row.names(fitness_PCA) <- paste(fitness_PCA$Plot,fitness_PCA$Subplot,
                                fitness_PCA$Plant_Simple,fitness_PCA$ID,sep="_")

#plot_i <- 8

fitness_PCA_i <- fitness_PCA %>% select(Line,Plot,Subplot,Plant_Simple,
                                        Seeds_GF,visits_GF,Fruit_GF,Psilothrix_viridicoerulea_homo,Coleoptera_homo,Homo_total)

library("FactoMineR")
#res.pca <- PCA(fitness_PCA[,c(4:8,10:ncol(fitness_PCA))], ncp = 5, graph = FALSE)
res.pca <- PCA(fitness_PCA_i[,c(6:ncol(fitness_PCA_i))], ncp = 5, graph = FALSE)

print(res.pca)

library("factoextra")
eig.val <- get_eigenvalue(res.pca)
eig.val
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

var <- get_pca_var(res.pca)
var

# Correlation circle
# Coordinates of variables
head(var$coord, 4)

fviz_pca_var(res.pca, col.var = "black")

# Quality of representation
head(var$cos2, 4)

library("corrplot")
corrplot(var$cos2, is.corr=FALSE)

# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(res.pca, choice = "var", axes = 1:2)


# Color by cos2 values: quality on the factor map
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

# Contributions of variables to PCs
# The contributions of variables in accounting for the variability in a given principal 
# component are expressed in percentage.
head(var$contrib, 5)
corrplot(var$contrib, is.corr=FALSE)  

# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)

#The most important (or, contributing) variables can be highlighted on the correlation plot
#as follow:
  
fviz_pca_var(res.pca, col.var = "contrib",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
  )

###################################################
# Graph of individuals
# Results

ind <- get_pca_ind(res.pca)
ind

# Coordinates of individuals
head(ind$coord)
# Quality of individuals
head(ind$cos2)
# Contributions of individuals
head(ind$contrib)

#Plots: quality and contribution

fviz_pca_ind(res.pca, col.ind = "cos2", 
             geom.ind = "point", # show points only (nbut not "text")
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE # Avoid text overlapping (slow if many points)
  )

#Color by groups
library(randomcoloR)

rcolors <- randomColor(36)

fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = as.factor(fitness_PCA_i$Line), # color by groups
             palette = rcolors,
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Line"
)

fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = as.factor(fitness_PCA_i$Plant_Simple), # color by groups
             palette = rcolors,
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Plant"
)

fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = as.factor(fitness_PCA_i$Plot), # color by groups
             palette = rcolors,
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Plots",
  
)

