
##########################################
##########                      ##########
##########  PACKAGE UPLOADING   ##########
##########                      ##########
##########################################

library(vegan)
library(betareg)
library(corrplot)
library(mgcv)
library(MuMIn)

source("C:/Users/Cunilleramontcusi/Dropbox/DAVID DOC/LLAM al DIA/CUNILLERA_palette.R")

##########################################
# LOAD DATASET from ALPINE LAKES _______________________________####
##### 16S bacterial plankton
load("S16-values/com.16S.Rdata")
load("S16-values/env.16S.Rdata")
com.16S <- com
env.16S <- env 
# select OTUS with min occurence of 10 sites (=exclude rare OTUs)
seltax<-which(colSums(com.16S>0)>9)
com_S16<-com.16S[,seltax]

##### 18S bacterial plankton
load("S16-values/com.18S.Rdata")
load("S16-values/env.18S.Rdata")
com.18S <- com
env.18S <- env 
# select OTUS with min occurence of 10 sites (=exclude rare OTUs)
seltax<-which(colSums(com.18S>0)>9)
com.18S<-com.18S[,seltax]

##### Phytoplankton bacterial plankton
load("S16-values/com.phyt.Rdata")
load("S16-values/env.phyt.Rdata")
com.phy <- com
env.phy <- env 
# select OTUS with min occurence of 2 sites (=exclude rare OTUs)
seltax<-which(colSums(com.phy>0)>1)
com.phy<-com.phy[,seltax]

##### Zooplankton bacterial plankton
load("S16-values/com.zoo.Rdata")
load("S16-values/env.zoo.Rdata")
com.zoo <- com
env.zoo <- env 

#_______________________________________________________________###
# IF you want a check with PA DATA run the following: 
#### WARNING: Change dissimilarity to jaccard if you are doing it!
com.16S <- ifelse(com.16S>0,1,0)
com.18S <- ifelse(com.18S>0,1,0)
com.phy <- ifelse(com.phy>0,1,0)
com.zoo <- ifelse(com.zoo>0,1,0)
#_______________________________________________________________###

# Data preparation and list building 
list.names <- c("S16", "S18", "phyto", "zoo")
comm_data<- list(com.16S, com.18S, com.phy, com.zoo)
names(comm_data) <-list.names 

list.names.env <- c("env.S16", "env.S18", "env.phyto", "env.zoo")
env_data<- list(env.16S, env.18S, env.phy, env.zoo)
names(env_data) <-list.names.env 

# Give individual names to each community to better further treatmetn  
ALP_S16 <- comm_data$S16
apply(ALP_S16,1,sum)

ALP_S18 <- comm_data$S18
apply(ALP_S18,1,sum)

ALP_phyto <- comm_data$phyto
apply(ALP_phyto,1,sum)

ALP_Zoo <- comm_data$zoo
apply(ALP_Zoo,1,sum)

# The same for environmental variables (the ones selected)
ALP_S16_env <- as.matrix(env_data$env.S16[,4:6])
ALP_S18_env <- as.matrix(env_data$env.S18[,4:6])
ALP_phyto_env <- as.matrix(env_data$env.phyto[,4:6])
ALP_Zoo_env <- as.matrix(env_data$env.zoo[,4:6])
#_______________________________________________________________###

# LOAD DATASET from NORWAY _______________________________####

load("C:/Users/Cunilleramontcusi/Dropbox/DAVID DOC/LLAM al DIA/1. Lunz al DIA/com.fish.Rdata")
load("C:/Users/Cunilleramontcusi/Dropbox/DAVID DOC/LLAM al DIA/1. Lunz al DIA/com.zoo.Rdata")
load("C:/Users/Cunilleramontcusi/Dropbox/DAVID DOC/LLAM al DIA/1. Lunz al DIA/com.phyto.Rdata")

Norw.fish <- fish
apply(Norw.fish,1,sum)

Norw.zoo <- ifelse(zoo>0,1,0)
apply(Norw.zoo,1,sum)
Norw.zoo <- Norw.zoo[-5,] 

Norw.phyto <- ifelse(com>0,1,0)
apply(Norw.phyto,1,sum)

load("C:/Users/Cunilleramontcusi/Dropbox/DAVID DOC/LLAM al DIA/1. Lunz al DIA/env.phyto.Rdata")
Norw.env <- env[,c(7,8,14,16:23)]
Norw.env[is.na(Norw.env)] <- 0
Norw.env <- as.matrix(log(Norw.env+1))

Norw.env.zoo <- Norw.env[-5,] 

# LOAD DATASET from SODA PANS _______________________________####
setwd("C:/Users/Cunilleramontcusi/Dropbox/DAVID DOC/LLAM al DIA/1. Lunz al DIA/S16-values/SodaPANS_Zsófia")
load("Seewinkel_predictors.RData")
load("Seewinkel_zoopl.RData")
setwd("C:/Users/Cunilleramontcusi/Dropbox/DAVID DOC/LLAM al DIA/1. Lunz al DIA/ALPINE_Lakes/Alpine_Lakes")

Sods.env <- cbind(pred$lZs,pred$lCond,pred$lTSS)
Sods.env <- as.matrix(Sods.env)

Sods.com <- ifelse(zoopl0>0,1,0)
Sods.com <- Sods.com[,which(apply(Sods.com,2,sum)>2)]
ncol(Sods.com)

##########################################
# Building lists with all the datasets and charging function  ####

# Community data
com_list <- list(ALP_S16, ALP_S18, ALP_phyto, ALP_Zoo, Norw.fish, Norw.zoo, Norw.phyto, Sods.com)
# Environmental data
env_list <- list(ALP_S16_env, ALP_S18_env, ALP_phyto_env, ALP_Zoo_env, Norw.env, Norw.env.zoo, Norw.env, Sods.env)

# Functions uploading
### Warning: distance in the dbRDA_Env_Track need to be written as character " "
dbRDA_Env_Track <- function(communnity_dataset,environmental_dataset, distance){
  # Calculate the dbRDA with the community data and the environmental and the corresponding distance
  # P/A - Jaccard distance
  dbRDA1<-capscale(communnity_dataset~environmental_dataset, distance = distance,add = "lingoes")
  
  # Fitteds - we predict according to the model where our samples should be
  fitteds <- fitted(dbRDA1, model="CCA", type= "response")
  
  tst_coeficient <- c()# Output to drop the coefficients
  # Correlation between predicted and observed values
  for (e in 1:nrow(communnity_dataset)) {
    tst_coeficient[e] <- cor(as.matrix(fitteds)[e,-e], # As matrix the fitteds to select the corresponding row, minus the column corresponding to themselves
                             as.matrix(vegdist(communnity_dataset, method = distance))[e,-e], method = "spearman") # Same for the fitteds
  }
tst_coeficient
}
CCA_Env_Track <- function(communnity_dataset,environmental_dataset){
  # Calculate the CCA with the community data and the environmental and the corresponding distance
  # P/A - Jaccard distance
  cca1<-cca(communnity_dataset~environmental_dataset)
  
  # Fitteds - we predict according to the model where our samples should be
  fitteds <- fitted(cca1, model="CCA", type= "response")
  
  tst_coeficient <- c()# Output to drop the coefficients
  # Correlation between predicted and observed values
  for (e in 1:nrow(communnity_dataset)) {
    tst_coeficient[e] <- cor(fitteds[e,], # As matrix the fitteds to select the corresponding row, minus the column corresponding to themselves
                             as.numeric(communnity_dataset[e,]), method = "spearman")
  }
tst_coeficient
}

# Test the dataset ####
# dbRDA test for data
output_dbRDA <- list()
for (e in 1:length(com_list)) {
output_dbRDA[[e]] <- dbRDA_Env_Track(communnity_dataset = com_list[[e]] ,
                               environmental_dataset = env_list[[e]] ,
                               distance = "jaccard")
}
# CCA test for data
output_CCA <- list()
for (e in 1:length(com_list)) {
  output_CCA[[e]] <- CCA_Env_Track(communnity_dataset = com_list[[e]] ,
                                       environmental_dataset = env_list[[e]])
}


# Results treatment
# As in Wind dispersal results in a gradient of dispersal limitation and environmental match among discrete aquatic habitats
# DOI: 10.1111/ecog.01685
Results_treatment <- function(Env_Track, lon, lat,number_colors){
Best_Models <- list()
Plots <- list()

value.results <- as.data.frame(cbind(Env_Track,lon, lat))
colnames(value.results) <- c("EnvTrack","lon","lat")

gm1<-gam(EnvTrack~s(lon, lat), data =value.results) 
gm2<-gam(EnvTrack~s(lon) + lat, data =value.results) 
gm3<-gam(EnvTrack~s(lat) + lon, data =value.results) 
lm3<-lm(EnvTrack~lon+lat, data =value.results) 

gamModels <- list(lm3, gm1, gm2, gm3)
Best_Models[[1]] <- model.sel(gamModels)

# Add some more place on the map around peripheral points 
xy<-expand.grid(seq(min(value.results$lon)-0.03, max(value.results$lon)+0.03, l=30), 
                seq(min(value.results$lat)-0.03, max(value.results$lat)+0.03, l=30))
xy_layout <- data.frame(lon=xy[,1], lat=xy[,2])
pred2<-predict(gamModels[[as.numeric(rownames(model.sel(gamModels)[1]))]], newdata=xy_layout)
  
# Choose colours
my.cols<-CUNILLERA_pal("LGTBI")(number_colors)

# Project environmental mismatch on the map 
filled.contour(seq(min(value.results$lon, na.rm=TRUE)-0.03, max(value.results$lon,na.rm=TRUE)+0.03, l=30),
               seq(min(value.results$lat, na.rm=TRUE)-0.03, max(value.results$lat, na.rm=TRUE)+0.03, l=30),
               matrix(pred2, nrow=30), 
               col=my.cols, 
               nlevels=30, 
               plot.axes = { points(value.results$lat~value.results$lon,
                                    bg=my.cols[cut(value.results$EnvTrack, length(my.cols), lab=F)], pch=21, cex=1.2, col="black");
                 axis(1); axis(2); title(main="Env_Track", xlab="Longitude", ylab="Latitude")
               })
Best_Models[[1]]
}

# ALPINE RESULTS ####
# S16 ________________________________________________________________________________________#
#dbRDA
Results_treatment(Env_Track =output_dbRDA[[1]] ,lon = env.16S$lon,lat =env.16S$lat, number_colors=35)
#CCA
Results_treatment(Env_Track =output_CCA[[1]] ,lon = env.16S$lon,lat =env.16S$lat, number_colors=35)

# S18 ________________________________________________________________________________________#
#dbRDA
Results_treatment(Env_Track =output_dbRDA[[2]] ,lon = env.18S$lon,lat =env.18S$lat, number_colors=35)
#CCA
Results_treatment(Env_Track =output_CCA[[2]] ,lon = env.18S$lon,lat =env.18S$lat, number_colors=35)

# PHYT ________________________________________________________________________________________#
#dbRDA
Results_treatment(Env_Track =output_dbRDA[[3]] ,lon = env.phy$lon,lat =env.phy$lat, number_colors=35)
#CCA
Results_treatment(Env_Track =output_CCA[[3]] ,lon = env.phy$lon,lat =env.phy$lat, number_colors=20)

# ZOO ________________________________________________________________________________________#
#dbRDA
Results_treatment(Env_Track =output_dbRDA[[4]] ,lon = env.zoo$lon,lat =env.zoo$lat, number_colors=35)
#CCA
Results_treatment(Env_Track =output_CCA[[4]] ,lon = env.zoo$lon,lat =env.zoo$lat, number_colors=25)

# NORW RESULTS ####
# Fish ________________________________________________________________________________________#
#dbRDA
Results_treatment(Env_Track =output_dbRDA[[5]] ,lon = env$Longitude,lat =env$Latitude, number_colors=40)
#CCA
Results_treatment(Env_Track =output_CCA[[5]] ,lon = env$Longitude,lat =env$Latitude, number_colors=40)

# Zoo ________________________________________________________________________________________#
#dbRDA
Results_treatment(Env_Track =output_dbRDA[[6]] ,lon = env$Longitude[-5],lat =env$Latitude[-5], number_colors=45)
#CCA
Results_treatment(Env_Track =output_CCA[[6]] ,lon = env$Longitude[-5],lat =env$Latitude[-5], number_colors=35)

# Phy ________________________________________________________________________________________#
#dbRDA
Results_treatment(Env_Track =output_dbRDA[[7]],lon = env$Longitude,lat =env$Latitude, number_colors=45)
#CCA
Results_treatment(Env_Track =output_CCA[[7]],lon = env$Longitude,lat =env$Latitude, number_colors=30)

# SODA PANS RESULTS ####
# dbRDA __________________________________________________________________________________#
#dbRDA
Results_treatment(Env_Track =output_dbRDA[[8]],lon = pred$lon,lat =pred$lat, number_colors=40)

#CCA
Results_treatment(Env_Track =output_CCA[[8]],lon = pred$lon,lat =pred$lat, number_colors=45)












