
#________________________________________#
#________________________________________#
##########  PACKAGE UPLOADING   ##########
#________________________________________#
#________________________________________#

# Data treatment and analysis
library(nlme)
library(vegan)
library(car)
library(mgcv)
library(Rmisc)
library(gamm4)
library(MuMIn)
library(usdm) # Collinearity
library(randomForestSRC) # RF
library(ggRandomForests) # RF
library(gbm) # BRT
library(dismo) # BRT
library(caret)

# Ploting and image 
library(ggplot2)
library(ggforce)
library(gtable)    
library(grid)
library(gridExtra) 
library(png)
library(viridis)
library(cowplot)
library(magick)

# Geographicla treatment
library(geosphere) # Charging package to caluclate distances
library(sp) # Other packages to deal with geospatial data 
library(stars) # Other packages to deal with geospatial data 
library(ggspatial) # Other packages to deal with geospatial data 
library(foreign) # Work with dbf files
library(rgdal) # upload and treat spatial data

# Plots to obtain a gg network plot
library(GGally)
library(network)
library(ggnetwork)

# Network treatment and analysis
library(sna)
library(igraph)
library(dplyr)

#________________________________________#
#________________________________________#
##########    DATA UPLOADING    ##########
#________________________________________#
#________________________________________#

# NOTES regarding data correction____________________####

# env: "Chl.a"     "Cond"      "lake_area" "Altitude"  - all log transformed
# rare taxa are excluded for phyto (min occurence 2 samples), 
# 16S & 18S (min both 10 samples)
# for phyto 16S & 18S abundances are relative and cubic root transformed
# for zoo cub root tranformed

##### 16S bacterial plankton
load("Database/com.16S.Rdata")
load("Database/env.16S.Rdata")
com.16S <- com
env.16S <- env 
# select OTUS with min occurence of 10 sites (=exclude rare OTUs)
seltax<-which(apply(com.16S,2,sum)>10)
com.16S<-com.16S[,seltax]

##### 18S bacterial plankton
load("Database/S18.prot.Rdata")
load("Database/env.18S.Rdata")
com.18S <- com
env.18S <- env 
# select OTUS with min occurence of 10 sites (=exclude rare OTUs)
com.18S <-com.18S[-c(9,13,25,28,39),] 
env.18S <- env.18S[-c(1,39),] 
seltax<-which(apply(com.18S,2,sum)>10)
com.18S<-com.18S[,seltax]

##### Phytoplankton bacterial plankton
load("Database/com.phyt.Rdata")
load("Database/env.phyt.Rdata")
com.phy <- com
env.phy <- env 
# select spp with min occurence of 3 sites (=exclude rare OTUs)
seltax<--which(apply(com.phy,2,sum)>3)
com.phy<-com.phy[,seltax]

##### Zooplankton 
load("Database/com.zoo.Rdata")
load("Database/env.zoo.Rdata")
com.zoo <- com
env.zoo <- env 
# select spp with min occurence of 3 sites (=exclude rare OTUs)
seltax<-which(apply(com.zoo,2,sum)>3)
com.zoo<-com.zoo[,seltax]

##### 18S bacterial zooplankton
load("Database/S18.zoo.Rdata")
load("Database/env.18S.Rdata")
com.18S.zoo <- S18.zoo
env.18S.zoo <- env 
# select OTUS with min occurence of 10 sites (=exclude rare OTUs)
seltax<-which(apply(ifelse(com.18S.zoo>0,1,0),2,sum)>10)
com.18S.zoo<-com.18S.zoo[,seltax]

#_______________________________________________________________###
# For the CCA with abundance values we keep the abundances of the available groups
#### WARNING: ONLY FOR THE CCA!
list.names <- c("S16", "S18", "phyto", "zoo", "18S.zoo")
comm_data_QUANTIT<- list(com.16S, com.18S, com.phy, com.zoo, com.18S.zoo)
names(comm_data_QUANTIT) <-list.names 






library(ggfortify)
PCA_network_results <- list()
PCA_network_plot<- list()
for (r in 1:length(network_data)) {
  PCA_result <- prcomp(network_data[[r]], center = T, scale. = T)
  PCA_network_plot[[r]] <- PCA_result
  PCA_network_results[[r]] <-network_data[[r]][,1] #PCA_result$x[,1] 
  hist(PCA_network_results[[r]])
}


# PCA 
names(PCA_network_results) <- c("600 km","300 km"," 100 km","60 km","6 km")

#grid.arrange(autoplot(PCA_network_plot[[1]],loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
#                      loadings.colour = 'red', loadings.label.colour="black")+
#               geom_point(size=2, alpha=0.1,color=CUNILLERA_cols("black"))+labs(title="~600km")+
#               theme_bw(), 
#             autoplot(PCA_network_plot[[2]],loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
#                      loadings.colour = 'red', loadings.label.colour="black")+
#               geom_point(size=2, alpha=0.1,color=CUNILLERA_cols("black"))+labs(title="~300km")+
#               theme_bw(),
#             autoplot(PCA_network_plot[[3]],loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
#                      loadings.colour = 'red', loadings.label.colour="black")+
#               geom_point(size=2, alpha=0.1,color=CUNILLERA_cols("black"))+labs(title="~100km")+
#               theme_bw(),
#             autoplot(PCA_network_plot[[4]],loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
#                      loadings.colour = 'red', loadings.label.colour="black")+
#               geom_point(size=2, alpha=0.1,color=CUNILLERA_cols("black"))+labs(title="~60km")+
#               theme_bw(),
#             autoplot(PCA_network_plot[[5]],loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
#                      loadings.colour = 'red', loadings.label.colour="black")+
#               geom_point(size=2, alpha=0.1,color=CUNILLERA_cols("black"))+labs(title="~6km")+
#               theme_bw(),  
#             ncol = 2, nrow=3)
#
#PCA_network_results[[1]] <- PCA_network_results[[1]]*-1
#PCA_network_results[[2]] <- PCA_network_results[[2]]*-1
#PCA_network_results[[5]] <- PCA_network_results[[5]]*-1

#autoplot(PCA_fluvial_network_plot[[1]],loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
#         loadings.colour = 'red', loadings.label.colour="black")+
#  geom_point(size=2, alpha=0.1,color=CUNILLERA_cols("black"))+labs(title="~600km")+
#  theme_bw()



# Calculate the fitted values according to "The method" 
### FITTINGS OF OBSERVED VS FITTEDS WITH CCA AND DBRDA

# CCA
tst_coeficient_out <- list()
for (r in 1:5) {
  require(pscl)
  cca1<-cca(comm_data[[r]]~as.matrix(env_data[[r]][,4:8]))
  fitteds <- fitted(cca1, model="CCA", type= "response")
  coeffici <- c()
  for (e in 1:nrow(comm_data[[r]])) {
    model <- glm(as.numeric(comm_data[[r]][e,])~fitteds[e,],family = "binomial")
    coeffici[e] <- pR2(model)["McFadden"]
  }
  #coeffici <- ifelse(coeffici<0,0,coeffici)
  tst_coeficient_out[[r]] <- coeffici
}

community_indices[[1]] <- tst_coeficient_out
names(community_indices[[1]]) <- c("EnvTrack_CCA_S16","EnvTrack_CCA_S18","EnvTrack_CCA_PHY","EnvTrack_CCA_ZOO", "EnvTrack_CCA_18S.ZOO") 



##CCA ___________________________________________________________________________________________
#p.val[1] <-summary.gam(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
#output_results[[1]] <- summary.gam(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
#output_model_results[[1]] <- gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
#preds_1 <- predict(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)


