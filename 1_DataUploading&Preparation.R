
#________________________________________#
#________________________________________#
##########  PACKAGE UPLOADING   ##########
#________________________________________#
#________________________________________#

# Data treatment and analysis
library(nlme); library(vegan);library(car);library(mgcv);library(Rmisc);library(gamm4);
library(MuMIn);library(usdm);library(randomForestSRC);library(ggRandomForests);library(gbm)
library(dismo);library(caret)

# Ploting and image 
library(ggplot2);library(ggforce);library(gtable);library(grid);
library(gridExtra);library(png);library(cowplot);library(magick)

# Geographicla treatment
library(geosphere) # Charging package to caluclate distances
library(sp) # Other packages to deal with geospatial data 
library(stars) # Other packages to deal with geospatial data 
library(ggspatial) # Other packages to deal with geospatial data 
library(foreign) # Work with dbf files
library(rgdal) # upload and treat spatial data

# Plots to obtain a gg network plot
library(GGally);library(network);library(ggnetwork)

# Network treatment and analysis
library(sna);library(igraph);library(dplyr);library(tidyr);library(stringr)


# Nice colours? CUNILLERA_palette is what you need
source("C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/CUNILLERA_palette.R")
library(viridis)
source("Alpine_Lakes_NetMetrics_functions.R")

#________________________________________#
#________________________________________#
##########    DATA UPLOADING    ##########
#________________________________________#
#________________________________________#

Data_directory <- c("C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/1. Lunz al DIA/ALPINE_Lakes/Alpine_Lakes")

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
com.18S <- as.data.frame(S18.prot) %>% tibble::rownames_to_column()
env.18S <- env 

correct_tot <- left_join(com.18S,env.18S, by=c("rowname"="lake")) %>% drop_na()

# select OTUS with min occurence of 3 sites (=exclude rare OTUs)
com.18S <-correct_tot %>% select(-colnames(env.18S)[2:8]) %>% tibble::column_to_rownames("rowname")
env.18S <- correct_tot %>% select(colnames(env.18S)[2:8]) %>% mutate(lake=correct_tot$rowname,.before=lon)%>% `rownames<-`(.[,1])

seltax<-which(apply(com.18S,2,sum)>3)
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
# IF you want a check with PA DATA run the following: 
#### WARNING: Change dissimilarity to jaccard if you are doing it!
com.16S <- ifelse(com.16S>0,1,0)
com.18S <- ifelse(com.18S>0,1,0)
com.phy <- ifelse(com.phy>0,1,0)
com.zoo <- ifelse(com.zoo>0,1,0)
com.18S.zoo <- ifelse(com.18S.zoo>0,1,0)
#_______________________________________________________________###

# Data preparation and list building 
list.names <- c("S16", "S18", "phyto", "zoo", "18S.zoo")
comm_data<- list(com.16S, com.18S, com.phy, com.zoo, com.18S.zoo)
names(comm_data) <-list.names 

list.names.env <- c("env.S16", "env.S18", "env.phyto", "env.zoo", "env.18S.zoo")
env_data<- list(env.16S, env.18S, env.phy, env.zoo, env.18S.zoo)
names(env_data) <-list.names.env 
#_______________________________________________________________###


#_______________________________________________________________###
# Geographic position loading and treatment

corr_LongLat <- read.table("Database/Lat_Long_correction.txt",header = F)
colnames(corr_LongLat) <- c("lake", "lat", "lon")
#_______________________________________________________________###
#_______________________________________________________________###
# IMPORTANT! RUN THE FOLLOWING LINES !!!!!!!!!!!!!!!!!!!!!!!!!!!
# Do it individually for each value... For 1 for 2 for 3 and for 4.... 
# It has to be done manually because for some "s" there is an error and the "for" stops running. This error is acutally 
#due to the fact that not all env_data have the third lake.
s <- 1
out <- c()
out[1] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[1])
out[2] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[2])
out[3] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[3])
out[4] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[4])
out[5] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[5])
out[6] <- NA

par(mfrow=c(1,2))
plot(env_data[[s]][out[-which(is.na(out))],2],env_data[[s]][out[-which(is.na(out))],3])
env_data[[s]][out[-which(is.na(out))],2] <- corr_LongLat$lon[which(is.na(out)==FALSE)]
env_data[[s]][out[-which(is.na(out))],3] <- corr_LongLat$lat[which(is.na(out)==FALSE)]
plot(env_data[[s]][out[-which(is.na(out))],2],env_data[[s]][out[-which(is.na(out))],3], col="red")

s <- 2
out <- c()
out[1] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[1])
out[2] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[2])
out[3] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[3])
out[4] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[4])
out[5] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[5])
out[6] <- NA

par(mfrow=c(1,2))
plot(env_data[[s]][out[-which(is.na(out))],2],env_data[[s]][out[-which(is.na(out))],3])
env_data[[s]][out[-which(is.na(out))],2] <- corr_LongLat$lon[which(is.na(out)==FALSE)]
env_data[[s]][out[-which(is.na(out))],3] <- corr_LongLat$lat[which(is.na(out)==FALSE)]
plot(env_data[[s]][out[-which(is.na(out))],2],env_data[[s]][out[-which(is.na(out))],3], col="red")

s <- 3
out <- c()
out[1] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[1])
out[2] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[2])
out[3] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[3])
out[4] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[4])
out[5] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[5])
out[6] <- NA

par(mfrow=c(1,2))
plot(env_data[[s]][out[-which(is.na(out))],2],env_data[[s]][out[-which(is.na(out))],3])
env_data[[s]][out[-which(is.na(out))],2] <- corr_LongLat$lon[which(is.na(out)==FALSE)]
env_data[[s]][out[-which(is.na(out))],3] <- corr_LongLat$lat[which(is.na(out)==FALSE)]
plot(env_data[[s]][out[-which(is.na(out))],2],env_data[[s]][out[-which(is.na(out))],3], col="red")

s <- 4
out <- c()
out[1] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[1])
out[2] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[2])
out[3] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[3])
out[4] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[4])
out[5] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[5])
out[6] <- NA

par(mfrow=c(1,2))
plot(env_data[[s]][out[-which(is.na(out))],2],env_data[[s]][out[-which(is.na(out))],3])
env_data[[s]][out[-which(is.na(out))],2] <- corr_LongLat$lon[which(is.na(out)==FALSE)]
env_data[[s]][out[-which(is.na(out))],3] <- corr_LongLat$lat[which(is.na(out)==FALSE)]
plot(env_data[[s]][out[-which(is.na(out))],2],env_data[[s]][out[-which(is.na(out))],3], col="red")

s <- 5
out <- c()
out[1] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[1])
out[2] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[2])
out[3] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[3])
out[4] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[4])
out[5] <- which(as.character(env_data[[s]][,1])==as.character(corr_LongLat$lake)[5])
out[6] <- NA

par(mfrow=c(1,2))
plot(env_data[[s]][out[-which(is.na(out))],2],env_data[[s]][out[-which(is.na(out))],3])
env_data[[s]][out[-which(is.na(out))],2] <- corr_LongLat$lon[which(is.na(out)==FALSE)]
env_data[[s]][out[-which(is.na(out))],3] <- corr_LongLat$lat[which(is.na(out)==FALSE)]
plot(env_data[[s]][out[-which(is.na(out))],2],env_data[[s]][out[-which(is.na(out))],3], col="red")

#_________________________________________________________________________________________________####
#_____________________________________________________________________________________________________
#_____________________________________________________________________________________________________

save.image("Database.RData")

