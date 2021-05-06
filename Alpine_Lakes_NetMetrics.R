
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


# Ploting and image 
library(ggplot2)
library(ggforce)
library(gtable)    
library(grid)
library(gridExtra) 
library(png)
library(viridis)

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
load("S16-values/com.16S.Rdata")
load("S16-values/env.16S.Rdata")
com.16S <- com
env.16S <- env 
# select OTUS with min occurence of 10 sites (=exclude rare OTUs)
seltax<-which(apply(com.16S,2,sum)>10)
com.16S<-com.16S[,seltax]

##### 18S bacterial plankton
load("S16-values/S18.prot.Rdata")
load("S16-values/env.18S.Rdata")
com.18S <- com
env.18S <- env 
# select OTUS with min occurence of 10 sites (=exclude rare OTUs)
com.18S <-com.18S[-c(9,13,25,28,39),] 
env.18S <- env.18S[-c(1,39),] 
seltax<-which(apply(com.18S,2,sum)>10)
com.18S<-com.18S[,seltax]

##### Phytoplankton bacterial plankton
load("S16-values/com.phyt.Rdata")
load("S16-values/env.phyt.Rdata")
com.phy <- com
env.phy <- env 
# select spp with min occurence of 3 sites (=exclude rare OTUs)
seltax<--which(apply(com.phy,2,sum)>3)
com.phy<-com.phy[,seltax]

##### Zooplankton 
load("S16-values/com.zoo.Rdata")
load("S16-values/env.zoo.Rdata")
com.zoo <- com
env.zoo <- env 
# select spp with min occurence of 3 sites (=exclude rare OTUs)
seltax<-which(apply(com.zoo,2,sum)>3)
com.zoo<-com.zoo[,seltax]

##### 18S bacterial zooplankton
load("S16-values/S18.zoo.Rdata")
load("S16-values/env.18S.Rdata")
com.18S.zoo <- com
env.18S.zoo <- env 
# select OTUS with min occurence of 10 sites (=exclude rare OTUs)
com.18S.zoo <-com.18S.zoo[-c(9,13,25,28,39),]
env.18S.zoo <- env.18S.zoo[-c(1,39),] 
seltax<-which(apply(com.18S.zoo,2,sum)>3)
com.18S.zoo<-com.18S.zoo[,seltax]

#_______________________________________________________________###
# For the CCA with abundance values we keep the abundances of the available groups
#### WARNING: ONLY FOR THE CCA!
list.names <- c("S16", "S18", "phyto", "zoo", "18S.zoo")
comm_data_QUANTIT<- list(com.16S, com.18S, com.phy, com.zoo, com.18S.zoo)
names(comm_data_QUANTIT) <-list.names 

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

corr_LongLat <- read.table("S16-values/Lat_Long_correction.txt",header = F)
colnames(corr_LongLat) <- c("lake", "lat", "lon")
#_______________________________________________________________###
#_______________________________________________________________###
# IMPORTANT! READ THE FOLLOWING LINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Do it individually for each value... For 1 for 2 for 3 and for 4.... 
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

# Nice colours? CUNILLERA_palette is what you need
source("C:/Users/Cunilleramontcusi/Dropbox/DAVID DOC/LLAM al DIA/CUNILLERA_palette.R")

#________________________________________#
#________________________________________#
##########  NETWORK BUILDING    ##########
#________________________________________#
#________________________________________#


# Uploading 
source("C:/Users/Cunilleramontcusi/Dropbox/DAVID DOC/LLAM al DIA/CUNILLERA_maxcomp_gradiente.R")

# Checking only the network of sampled lakes (not accounting with the lakes in between)
dist.percol <- list()
for (r in 1:4) {
  situacio <- env_data[[r]][,2:3]
  library(ggplot2)
  library(geosphere)
  library(sp)
  library(stars)
  library(ggspatial)
  # Transforming to sf format
  lakes_sw_small <- st_as_sf(situacio, coords = c("lon","lat"))
  # Convert to "old" format to caluclate to distances
  lakes_sw_sp_small <- as(lakes_sw_small, "Spatial")
  # Calculate the distance matrix (in METRES!!!)
  dist_lakes_small <- distm(lakes_sw_sp_small, fun = distGeo)
  
  max.comp_gradiente(dist_lakes_small,10)->aa # funci? "max.comp_gradiente" calcula de totes les dist?ncies entre tots els nodes la dist?ncia 
  # i el nombre de nodes connectats entre ells, per tots els nodes 
  which(aa[,2]==ncol(dist_lakes_small))[1]->k #De la llista generada selecciona la primera dist?ncia en que tots els nodes s'ajunten (la xarx de percolaci?)
  aa[k,1]->d.percol #li donem el nom a i a?llem aquesta dist?ncia de percolaci?
  dist.percol[[r]] <- d.percol
}

# MAPS_ALPINE LAKES DATA
maps_list <- list()
for (r in 1:4) {
  comunitats<- comm_data[[r]]
  situacio <- env_data[[r]][,2:3]
  library(ggplot2)
  library(geosphere)
  library(sp)
  library(stars)
  library(ggspatial)
  # Transforming to sf format
  lakes_sw_small <- st_as_sf(situacio, coords = c("lon","lat"))
  # Convert to "old" format to caluclate to distances
  lakes_sw_sp_small <- as(lakes_sw_small, "Spatial")
  # Calculate the distance matrix (in METRES!!!)
  dist_lakes_small <- distm(lakes_sw_sp_small, fun = distGeo)
  ALP_xarxa <- ifelse(dist_lakes_small>dist.percol[[r]],0,1)
  library(igraph)
  # In order to plot the network in a ggplot way do the following
  n<- network(ALP_xarxa, directed=F, diag=F) 
  #n %v% "family" <- xy.guils[,3] # Family is an standard name for the categortical variable that we are creating 
  detach("package:igraph", unload = TRUE)
  # ALL THESE LINES MUST BE RUNNED!!!
  #from here_____________________________________________________________________
  cols <- CUNILLERA_pal(palette = "sea")(4)[as.numeric(cut(env_data[[r]][,7],breaks = 4))]
  #valors_factors <- which(factors==1)
  
  maps_list[[r]] <- ggplot(n, layout=as.matrix(situacio) ,
                           aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(color = "grey40", size=0.5)+
    geom_nodes(size=5 ,color="black", shape=21, fill=cols, alpha=.85)+
    labs(x="",y="")+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),legend.position = "none")
}  
grid.newpage()
grid.arrange(maps_list[[1]],maps_list[[2]],maps_list[[3]],maps_list[[4]], ncol = 1, nrow=4)


#________________________________________#
# Building the whole network based on lks.global _____________________________________________________________####
#________________________________________#

load("lks.global.Rdata")

# Obtain the longitude and latitude of half europe (enough to built a buffer)
sel<-which(lks.global$lon>-6 & lks.global$lon<26 & lks.global$lat>37 & lks.global$lat<55)
dat<-data.frame(lks.global[sel, c("lon", "lat")])

# Extract the coordinates of all the sampled lakes (55 different lakes)
all_lakes_coord <- unique(rbind(env_data[[1]][,2:3],env_data[[2]][,2:3],env_data[[3]][,2:3],env_data[[4]][,2:3]))
dat1<-data.frame(lon=all_lakes_coord[,1], lat=all_lakes_coord[,2], sampled=1)

# Lakes buffer extracts a buffer of lakes surrounding the sampled lakes 
### EU_lakes: The lks.global dataset
### Samp_lakes: the lakes that have been sampled in the context of our study. 
### Buffer distance: This value will define the distance of the buffer. 
###  by default the bufffer is set at the maximum distance between the sampled lakes (buffer_distance = 1). In our case 600km
###  each value that we add there is divides this maximum value (buffer_distance = 2 means a distance of 600/2=300km).

# The output of this function is a list with: 
## [[1]]= A matrix with the coordinates of all the lakes being the sampled ones the las 55 lakes.
## [[2]]= A distance matrix of all the lakes beteen each other in meters 

Lakes_buffer <- function(EU_lakes, Samp_lakes, buffer_distance, check_EU_Sampl_plots=F, check_BUFFER_plots=F,check_ELIMINATION_goodLAKES_plots=F){
  
  # Calculating the corresponding sites of the global lakes database 
  dist.between.lakes <- dist(rbind(EU_lakes,Samp_lakes[,1:2]),method = "euclidean" )
  # Extraction of the closest point to each lake
  corr.lakes <- c()
  checking.value<- c()
  for (i in 1:nrow(Samp_lakes)) {
    corresponding.lake <- which(as.matrix(dist.between.lakes)[nrow(EU_lakes)+i,]==
                                  min(as.matrix(dist.between.lakes)[nrow(EU_lakes)+i,][1:nrow(EU_lakes)]))
    checking.value[i] <- min(as.matrix(dist.between.lakes)[nrow(EU_lakes)+i,][1:nrow(EU_lakes)])
    corr.lakes[i]<- rownames(EU_lakes[corresponding.lake,])
  }
  corr.lakes
  checking.value
  the.good.lakes <- c()
  for (e in 1:length(corr.lakes)) {
    the.good.lakes[e] <- which(rownames(EU_lakes)==corr.lakes[e])
  }
  #The good lakes are the rows of the "EU_lakes" dataset that correspond to the lakes
  
  ###________________________
  # Plots checking if the coordinates of the selected lakes are right or not. 
  if (check_EU_Sampl_plots==T) {
    par(mfrow=c(1,2))
    plot(EU_lakes[the.good.lakes,1],Samp_lakes[,1])
    abline(a=0,b=1, col="red")
    plot(EU_lakes[the.good.lakes,2],Samp_lakes[,2])
    abline(a=0,b=1, col="red")
    par(mfrow=c(1,1))  
    # Plot the whole dataset of points
    #plot(EU_lakes$lon, EU_lakes$lat)    
  }
  
  # Transforming to sf format
  lakes_sw <- st_as_sf(EU_lakes, coords = c("lon","lat"))
  # Convert to "old" format to caluclate to distances
  lakes_sw_sp <- as(lakes_sw, "Spatial")
  # Calculate the distance matrix (in METRES!!!)
  dist_lakes <- distm(lakes_sw_sp, fun = distGeo)
  
  #Calculate all the distances between the smapled lakes
  dist.between.the.good.lakes <- dist_lakes[the.good.lakes,the.good.lakes]
  #dist.between.the.good.lakes <- dist.between.the.good.lakes/2
  # Transform EU_lakes into a new object (just security)
  dat_0 <- EU_lakes
  # Add a third column informing the current order for being able to detect known nodes
  dat_0[,3] <- seq(from=1, to = nrow(EU_lakes), by = 1)
  
  # Let's extract the nodes identities that are wihtin a range of the sampled lakes
  # We build a buffer based on the maximum distance between the sampled lakes: max(dist.between.the.good.lakes)
  datum <- list()
  dat_total <- NULL
  for (f in 1:length(the.good.lakes)) {
    datum[[f]]<- dat_0[-c(which(dist_lakes[the.good.lakes[f],]>max(dist.between.the.good.lakes)/buffer_distance)),]
    dat_total <- rbind(dat_total,datum[[f]])
  }
  # There are duplicates so we collapse these values with "unique"
  dat_final <- unique(dat_total)
  
  if (check_BUFFER_plots==T) {
    #Check versus other unique plots (built with the "buffer" around an specific lake)
    par(mfrow=c(2,3))
    plot(dat_0[-c(which(dist_lakes[the.good.lakes[1],]>max(dist.between.the.good.lakes)/buffer_distance)),1:2])
    plot(dat_0[-c(which(dist_lakes[the.good.lakes[52],]>max(dist.between.the.good.lakes)/buffer_distance)),1:2])
    plot(dat_final[,1:2])
    
    #Check versus unique/total. The same plot but with less nodes should be observed
    plot(dat_final[,1:2])
    plot(dat_total[,1:2])
    
    #Check the distances, shoulg be shorter in "final" 
    barplot(c(nrow(dat_final), nrow(dat_total)))
    par(mfrow=c(1,1))
  }
  
  # Out is the "new" rows of "the.good.lakes"
  out <- c()
  for (ff in 1:length(the.good.lakes)) {
    out[ff] <- which(dat_final[,3]==dat_0[the.good.lakes,3][ff])
  }
  out
  
  if (check_ELIMINATION_goodLAKES_plots==T) {
    # We check that we are eliminating some of the lakes (the sampled ones)
    par(mfrow=c(1,2))
    plot(dat_final[,1:2])
    plot(dat_final[-out,1:2])
    par(mfrow=c(1,1))
  }
  
  # We deffinetly extract the good lakes values (lack of precision for them)
  dat_buffer <- dat_final[-c(out),]
  
  #Adding the real coordinates of the lakes (this is done because some of them were too close to one another)
  dat_corr <- rbind(dat_buffer[,1:2], all_lakes_coord[,1:2])
  #plot(dat_corr)
  # Transforming to sf format
  lakes_sw <- st_as_sf(dat_corr, coords = c("lon","lat"))
  # Convert to "old" format to caluclate to distances
  lakes_sw_sp <- as(lakes_sw, "Spatial")
  # Calculate the distance matrix (in METRES!!!)
  dist_lakes <- distm(lakes_sw_sp, fun = distGeo)
  
  function_result <- list(dat_corr,dist_lakes)
}

#
#
# TOPOGRÀFIC DISTANCES TRIES
#
#library(raster)
#library(topoDistance)
#
#a <- c('C:/Users/Cunilleramontcusi/Desktop/GIS - Alpine/dem_suiss.tif')
#
#template <- raster(a)
#
#plot(template)
#
#pitit <- matrix(ncol = 2, byrow = TRUE,
#                c(as.numeric(dat1[42,1:2]),
#                  as.numeric(dat1[6,1:2]),
#                  as.numeric(dat1[54,1:2])
#                  ))
#
#tdist <- topoDist(template$dem_suiss, pitit, paths = TRUE)
#
#topoPathMap(template$test_reproj, pitit, 
#            topoPaths = tdist, type = "hillshade",
#            pathWidth = 4, cex = 2, bg = "blue")
#
#
#

#________________________________________#
# Maxmimum distance between lakes -- 600 km approximately
## This has been uploaded to the supercomputer to calculate d.percol__________________
max_distance <- Lakes_buffer(EU_lakes = dat, Samp_lakes = dat1, buffer_distance = 1, 
                             check_EU_Sampl_plots = F, check_BUFFER_plots = T,check_ELIMINATION_goodLAKES_plots = F )
# We extract the distance matrix and save it to run it externally (accelerate computing times)
max_distance_down <- max_distance[[2]]
#save(max_distance_down, file = "S16-values/AlPS_lakes_percol_MAX_DISTANCE.RData")

# Run the following chunk in some "supercomputer" that can last some days running and keep the result!
# save.image with the "dist_lakes" object 
# Also remember that you have to upload the max.comp.gradiente_NORW
#
# Something like this but well written... if not do it manually deleting everything except "dist_lakes" it is not that terrible...  
#save(list =dist_lakes, file = "S16-values/AlPS_lakes_percol.RData")
#
#library(doParallel)
#registerDoParallel(cores = detectCores())
#out<-foreach(i =T , .combine=rbind) %dopar% {
#  library("sna")
#  max.comp_gradiente_NORW(min_distancias = 140000, max_distancias = 150000,dist_lakes,1000)->aa
#}

#aa <- out#taula amb aquests valors
#which(aa[,2]==ncol(dist_lakes))[1]->k #De la llista generada selecciona la primera dist?ncia en que tots els nodes s'ajunten (la xarx de percolaci?)
#aa[k,1]->d.percol #li donem el nom a i a?llem aquesta dist?ncia de percolaci?
#d.percol #ladist?ncia en concret

#________________________________________#
# Middle distance between lakes -- 300 km approximately
## This has been uploaded to the supercomputer to calculate d.percol__________________
mid_distance <- Lakes_buffer(EU_lakes = dat, Samp_lakes = dat1, buffer_distance = 2, 
                             check_EU_Sampl_plots = F, check_BUFFER_plots = T,check_ELIMINATION_goodLAKES_plots = F )
# We extract the distance matrix and save it to run it externally (accelerate computing times)
mid_distance_down <- mid_distance[[2]]
#save(mid_distance_down, file = "S16-values/AlPS_lakes_percol_MID_DISTANCE.RData")

# Run the following chunk in some "supercomputer" that can last some days running and keep the result!
# save.image with the "dist_lakes" object 
# Also remember that you have to upload the max.comp.gradiente_NORW
#
# Something like this but well written... if not do it manually deleting everything except "dist_lakes" it is not that terrible...  
#save(list =dist_lakes, file = "S16-values/AlPS_lakes_percol.RData")
#
#library(doParallel)
#registerDoParallel(cores = detectCores())
#out<-foreach(i =T , .combine=rbind) %dopar% {
#  library("sna")
#  max.comp_gradiente_NORW(min_distancias = 50000, max_distancias = 150000,dist_lakes,100000)->aa
#}

#aa <- out#taula amb aquests valors
#which(aa[,2]==ncol(dist_lakes))[1]->k #De la llista generada selecciona la primera dist?ncia en que tots els nodes s'ajunten (la xarx de percolaci?)
#aa[k,1]->d.percol #li donem el nom a i a?llem aquesta dist?ncia de percolaci?
#d.percol #ladist?ncia en concret

#________________________________________#
# Medium middle (MID_MID) distance between lakes -- 100 km approximately
## This has been uploaded to the supercomputer to calculate d.percol__________________
mid_mid_distance <- Lakes_buffer(EU_lakes = dat, Samp_lakes = dat1, buffer_distance = 6, 
                                 check_EU_Sampl_plots = F, check_BUFFER_plots = T,check_ELIMINATION_goodLAKES_plots = F )
# We extract the distance matrix and save it to run it externally (accelerate computing times)
mid_mid_distance_down <- mid_mid_distance[[2]]
#save(mid_mid_distance_down, file = "S16-values/AlPS_lakes_percol_MID_MID_DISTANCE.RData")

# Run the following chunk in some "supercomputer" that can last some days running and keep the result!
# save.image with the "dist_lakes" object 
# Also remember that you have to upload the max.comp.gradiente_NORW
#
# Something like this but well written... if not do it manually deleting everything except "dist_lakes" it is not that terrible...  
#save(list =dist_lakes, file = "S16-values/AlPS_lakes_percol.RData")
#
#library(doParallel)
#registerDoParallel(cores = detectCores())
#out<-foreach(i =T , .combine=rbind) %dopar% {
#  library("sna")
#  max.comp_gradiente_NORW(min_distancias = 20000, max_distancias = 150000,dist_lakes,10000)->aa
#}

#aa <- out#taula amb aquests valors
#which(aa[,2]==ncol(dist_lakes))[1]->k #De la llista generada selecciona la primera dist?ncia en que tots els nodes s'ajunten (la xarx de percolaci?)
#aa[k,1]->d.percol #li donem el nom a i a?llem aquesta dist?ncia de percolaci?
#d.percol #ladist?ncia en concret

#________________________________________#
# Small distance, A tenth of the maximum distance -- 60km approximately
small_distance <- Lakes_buffer(EU_lakes = dat, Samp_lakes = dat1, buffer_distance = 10, 
                               check_EU_Sampl_plots = F, check_BUFFER_plots = T,check_ELIMINATION_goodLAKES_plots = F )
### WARNING! ###
# Small computing times and the fact that the two "sub-networks" are now disconnected made us to treat them differently: 

## First: We extract the two "subnetworks" 
small_distance_down <- small_distance[[2]]
position_node_above10_SMALL <- which(small_distance[[1]][,1]>10)
position_node_below10_SMALL <- which(small_distance[[1]][,1]<10)

## Second: We calculate the percolation distance of each of these subnetworks. 
aa <- max.comp_gradiente(distancias =small_distance_down[position_node_above10_SMALL,position_node_above10_SMALL], br =1000)
k <- which(aa[,2]==ncol(small_distance_down[position_node_above10_SMALL,position_node_above10_SMALL]))[1] #De la llista generada selecciona la primera dist?ncia en que tots els nodes s'ajunten (la xarx de percolaci?)
aa[k,1]->d.percol_FIRST #li donem el nom a i a?llem aquesta dist?ncia de percolaci?
d.percol_FIRST #ladist?ncia en concret

aa <- max.comp_gradiente(distancias =small_distance_down[position_node_below10_SMALL,position_node_below10_SMALL], br =1000)
k <- which(aa[,2]==ncol(small_distance_down[position_node_below10_SMALL,position_node_below10_SMALL]))[1] #De la llista generada selecciona la primera dist?ncia en que tots els nodes s'ajunten (la xarx de percolaci?)
aa[k,1]->d.percol_SECOND #li donem el nom a i a?llem aquesta dist?ncia de percolaci?
d.percol_SECOND #ladist?ncia en concret

## Third: We change the distances according to each of the two percolation distances (greater=0, smaller=1). 
## Finally, all other distances (those between the two subnetworks are changed to 0)
ALP_SMALL_xarxa <- small_distance_down
ALP_SMALL_xarxa[position_node_above10_SMALL,position_node_above10_SMALL] <- ifelse(ALP_SMALL_xarxa[position_node_above10_SMALL,position_node_above10_SMALL]>d.percol_FIRST,0,1)
ALP_SMALL_xarxa[position_node_below10_SMALL,position_node_below10_SMALL] <- ifelse(ALP_SMALL_xarxa[position_node_below10_SMALL,position_node_below10_SMALL]>d.percol_SECOND,0,1)
ALP_SMALL_xarxa <- ifelse(ALP_SMALL_xarxa>1,0,ALP_SMALL_xarxa)

#________________________________________#
# The smallest distance (minimum), basically individual lakes or nearby (almost touching) lakes -- 6km approximately
min_distance <- Lakes_buffer(EU_lakes = dat, Samp_lakes = dat1, buffer_distance = 100, 
                             check_EU_Sampl_plots = F, check_BUFFER_plots = T,check_ELIMINATION_goodLAKES_plots = F )
### WARNING! ###
# Small computing times and the fact that the two "sub-networks" are now disconnected made us to treat them differently: 
min_distance_down <- min_distance[[2]]
position_node_above10_MIN <- which(min_distance[[1]][,1]>10)
position_node_below10_MIN <- which(min_distance[[1]][,1]<10)

## Second: We calculate the percolation distance of each of these subnetworks. 
aa <- max.comp_gradiente(distancias =min_distance_down[position_node_above10_MIN,position_node_above10_MIN], br =1000)
k <- which(aa[,2]==ncol(min_distance_down[position_node_above10_MIN,position_node_above10_MIN]))[1] #De la llista generada selecciona la primera dist?ncia en que tots els nodes s'ajunten (la xarx de percolaci?)
aa[k,1]->d.percol_FIRST #li donem el nom a i a?llem aquesta dist?ncia de percolaci?
d.percol_FIRST #ladist?ncia en concret

aa <- max.comp_gradiente(distancias =min_distance_down[position_node_below10_MIN,position_node_below10_MIN], br =1000)
k <- which(aa[,2]==ncol(min_distance_down[position_node_below10_MIN,position_node_below10_MIN]))[1] #De la llista generada selecciona la primera dist?ncia en que tots els nodes s'ajunten (la xarx de percolaci?)
aa[k,1]->d.percol_SECOND #li donem el nom a i a?llem aquesta dist?ncia de percolaci?
d.percol_SECOND #ladist?ncia en concret

## Third: We change the distances according to each of the two percolation distances (greater=0, smaller=1). 
## Finally, all other distances (those between the two subnetworks are changed to 0)
ALP_MIN_xarxa <- min_distance_down
ALP_MIN_xarxa[position_node_above10_MIN,position_node_above10_MIN] <- ifelse(ALP_MIN_xarxa[position_node_above10_MIN,position_node_above10_MIN]>d.percol_FIRST,0,1)
ALP_MIN_xarxa[position_node_below10_MIN,position_node_below10_MIN] <- ifelse(ALP_MIN_xarxa[position_node_below10_MIN,position_node_below10_MIN]>d.percol_SECOND,0,1)
ALP_MIN_xarxa <- ifelse(ALP_MIN_xarxa>1,0,ALP_MIN_xarxa)

#________________________________________#
# I already calculated the percolation distance for each network and already noted it in the following lines. 
#there is no need to upload the values in the three major cases. 
#MAX= 1422222
#load("S16-values/PC_cluster_out/ALPINE_MAX_distance.RData")
#MID= 68732.87
#load("S16-values/PC_cluster_out/ALPINE_MID_distance.RData")
#MID_MID= 72036.2
#load("S16-values/PC_cluster_out/ALPINE_MID_MID_distance.RData")

ALP_MAX_xarxa <- max_distance_down
ALP_MAX_xarxa <- ifelse(max_distance_down>142222.2,0,1)

ALP_MID_MID_xarxa <- mid_mid_distance_down
ALP_MID_MID_xarxa <- ifelse(ALP_MID_MID_xarxa>72036.2,0,1)

ALP_MID_xarxa <- mid_distance_down
ALP_MID_xarxa <- ifelse(ALP_MID_xarxa>68732.87,0,1)

cordenades_xarxes <-list(max_distance[[1]],mid_distance[[1]],mid_mid_distance[[1]], small_distance[[1]], min_distance[[1]]) 
MAPS_xarxes <- list(ALP_MAX_xarxa, ALP_MID_xarxa, ALP_MID_MID_xarxa, ALP_SMALL_xarxa, ALP_MIN_xarxa)

#________________________________________#
#________________________________________#
######## RIVER NETWORK BUILDING ##########
#________________________________________#
#________________________________________#

# Print de coordinates of the two smaller scales to posteriorly trim the "river" network in Arcgis
#write.dbf(max_distance[[1]], "GIS_data/max_distance.dbf")
#write.dbf(mid_distance[[1]], "GIS_data/mid_distance.dbf")
#write.dbf(mid_mid_distance[[1]], "GIS_data/mid_mid_distance.dbf")
#write.dbf(small_distance[[1]], "GIS_data/small_distance.dbf")
#write.dbf(min_distance[[1]], "GIS_data/min_distance.dbf")

max_distance_BASINS <- read.dbf("GIS_data/max_distance_BASINS.dbf")
mid_distance_BASINS <- read.dbf("GIS_data/mid_distance_BASINS.dbf")
mid_mid_distance_BASINS <- read.dbf("GIS_data/mid_mid_distance_BASINS.dbf")
small_distance_BASINS <- read.dbf("GIS_data/small_distance_BASINS.dbf")
min_distance_BASINS <- read.dbf("GIS_data/min_distance_BASINS.dbf")

cordenades_xarxes_BASINS <- list(max_distance_BASINS,mid_distance_BASINS, mid_mid_distance_BASINS,small_distance_BASINS,min_distance_BASINS)

GRAPH_xarxes_fluvial <- list()
all_lakes_BASINS_fluvial <- list()
correspondence_BASINS_fluvial <- list()

detach("package:sna", unload = TRUE)
library(igraph)  
library(shp2graph)

correspondence_BASINS <- c()
for (a in 1:55) {
  correspondence_BASINS[a] <- which(round(cordenades_xarxes[[1]][(nrow(cordenades_xarxes[[1]])-55+a),1],4)==round(cordenades_xarxes_BASINS[[1]][,1],4))  
}
correspondence_BASINS_fluvial[[1]] <- correspondence_BASINS
# We rewrite "ptsxy" with the converted coordinates to lat/long
ptsxy <- coordinates(cordenades_xarxes_BASINS[[1]])

#   Load the river shapefile. This shapefil is the entire basin of the Ebre river. 
# Extracted from https://www.hydrosheds.org/downloads where we obtained the catchments and rivers
# shapefile and filtered the corresponding Ebre basin. 
shape_Rivers <- st_read("GIS_data/Rivers_BASINS.shp")
shape_Rivers2 <- as(shape_Rivers, "Spatial")
#plot(shape_Rivers2)

# With the points2network we find which points of the shapefile are closer to our points (ptsxy).
# Thus, we generate a vector of IDs of which nodes of the shape correspond to each point in the ptsxy
# approach=1 searches the closer "node" to the desired point. 
# approach=2 relocates the points as new network nodes (adds new nodes)-> Not recomended as then it is difficult to 
#work with the values related with the graph... there are no equivalences...
res.nv <- points2network(ntdata = shape_Rivers2, pointsxy= ptsxy, ELComputed = TRUE,approach = 1)
# Checking the correspondence between points and their corresponding nodes
#ptsinnt.view(ntdata = shape_Rivers2, nodelist = res.nv[[1]], pointsxy = ptsxy, CoorespondIDs= res.nv[[3]])

# Generate the vector with the nodes "IDs" that correspond to our sampling points
# Check this vector
all_lakes_BASINS_fluvial[[1]] <- unlist(res.nv[[3]])

# Now that we have the list of nodes that correspond to our sampling points. Let's built de graph and play a bit: 
# Shape to igraph -> every segment is converted to a central "node" and their corresponding "edges"
rtNEL1 <-readshpnw(shape_Rivers2,ELComputed = T)

# Now we generate the strict graph
GRAPH_xarxes_fluvial[[1]] <- nel2igraph(nodelist = rtNEL1[[2]],edgelist = rtNEL1[[3]],Directed = T,weight = rtNEL1[[4]])

ConComp = components(GRAPH_xarxes_fluvial[[1]])
cols <-ConComp$membership 
cols[unlist(res.nv[[3]])[correspondence_BASINS]] <- "red"
cols <- ifelse(cols==1, "green",ifelse(cols==2,"blue","red"))
sizes <- ifelse(cols=="red",1.2,0.4)

png(filename = "C:/Users/Cunilleramontcusi/Alpine_fluvial.png", width = 20000, height = 20000, res=500)
par(mar=c(0,0,0,0))
plot(GRAPH_xarxes_fluvial[[1]], vertex.label = NA, vertex.size = sizes, vertex.size2 = sizes, vertex.color=cols,
     edge.width=0.01,edge.size=0.01, edge.color="grey50")
dev.off()
detach("package:shp2graph", unload = TRUE)

#________________________________________#
#________________________________________#
########    CENTRALITY VALUES   ##########
#________________________________________#
#________________________________________#

# Calculation of network data for the 55 lakes (closeness, degree, between, evcent, subgraph centrality)
network_data <- list()
detach("package:sna", unload = TRUE)
library(igraph)
for (i in 1:length(MAPS_xarxes)) {
  Xarx_grpah <- graph.adjacency(MAPS_xarxes[[i]], mode = "undirected",diag = F)
  output <- matrix(nrow = nrow(MAPS_xarxes[[i]]), ncol = 3)
  
  network_data[[i]] <- output
  
  network_data[[i]][,1] <- closeness(Xarx_grpah, mode = "all")
  network_data[[i]][,2] <- degree(Xarx_grpah)
  network_data[[i]][,3] <- betweenness(Xarx_grpah, directed = F)
  
  if(components(Xarx_grpah)$no>1){
    # Closenness
    a <- which(components(Xarx_grpah)$membership==1)
    network_data[[i]][a,1]<- network_data[[1]][a,1]/max(network_data[[1]][a,1])
    b <- which(components(Xarx_grpah)$membership==2)
    network_data[[i]][b,1] <- network_data[[1]][b,1]/max(network_data[[1]][b,1])
    
    # Betweenness
    a <- which(components(Xarx_grpah)$membership==1)
    network_data[[i]][a,3]<- network_data[[1]][a,3]/max(network_data[[1]][a,3])
    b <- which(components(Xarx_grpah)$membership==2)
    network_data[[i]][b,3] <- network_data[[1]][b,3]/max(network_data[[1]][b,3])
  }
  colnames(network_data[[i]]) <- c("clo_ALPS","deg_ALPS","bet_ALPS")
}

# Flulvial
library(igraph)
fluvial_network_data <- list()
output <- matrix(nrow = length(V(GRAPH_xarxes_fluvial[[1]])), ncol = 3)

fluvial_network_data[[1]] <- output

fluvial_network_data[[1]][,1] <- closeness(GRAPH_xarxes_fluvial[[1]], mode = "out")
fluvial_network_data[[1]][,2] <- degree(GRAPH_xarxes_fluvial[[1]],mode = "out")
fluvial_network_data[[1]][,3] <- betweenness(GRAPH_xarxes_fluvial[[1]])

colnames(fluvial_network_data[[1]]) <- c("clo_ALPS","deg_ALPS","bet_ALPS")

#Homogeneizing dataset according to catchments
#Closeness
a <- which(ConComp$membership==1)
fluvial_network_data[[1]][a,1]<- fluvial_network_data[[1]][a,1]/max(fluvial_network_data[[1]][a,1])
b <- which(ConComp$membership==2)
fluvial_network_data[[1]][b,1] <- fluvial_network_data[[1]][b,1]/max(fluvial_network_data[[1]][b,1])

# PCA and plots

library(ggfortify)
PCA_network_results <- list()
PCA_network_plot<- list()
for (r in 1:length(network_data)) {
  PCA_result <- prcomp(network_data[[r]], center = T, scale. = T)
  PCA_network_plot[[r]] <- PCA_result
  PCA_network_results[[r]] <-network_data[[r]][,1] #PCA_result$x[,1]
}
names(PCA_network_results) <- c("max_PCA_network","mid_PCA_network","mid_mid_PCA_network",
                                "small_PCA_network","min_PCA_network")

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

#Plots closeness maps - pallete viridis
maps_Alps <- list()
detach("package:igraph", unload = TRUE)
library(sna)
for (e in 1:length(cordenades_xarxes)) {
  factors <- rep("No_Sampled",nrow(MAPS_xarxes[[e]]))
  factors[c(nrow(MAPS_xarxes[[e]])-54):nrow(MAPS_xarxes[[e]])] <- "Sampled"
  CC_values <- PCA_network_results[[e]]
  
  n<- network(MAPS_xarxes[[e]], directed=F, diag=F)
  n %v% "family" <- factors # Family is an standard name for the categortical variable that we are creating
  n %v% "CC_values" <- CC_values
  
  maps_Alps[[e]] <- ggplot(n, layout=as.matrix(cordenades_xarxes[[e]][,1:2]),
                           aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( color = "grey40", size=0.1, alpha=0.4) +
    geom_nodes(aes(fill=CC_values,alpha=family, size=family),color ="grey20" ,shape=21, alpha=.75)+
    scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
    labs(x="",y="")+
    scale_alpha_manual(values = c(0.2,5))+
    scale_size_manual(values = c(1,3))+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_rect(colour="black"))
}

png(filename ="C:/Users/Cunilleramontcusi/Network_maps_Alps.png",width =745 ,height =742 ,units ="px",res = 100)
grid.arrange(maps_Alps[[1]],maps_Alps[[2]],maps_Alps[[3]],maps_Alps[[4]],maps_Alps[[5]])
dev.off()


library(ggfortify)
PCA_fluvial_network_results <- list()
PCA_fluvial_network_plot<- list()

PCA_fluvial_result <- prcomp(fluvial_network_data[[1]], center = T, scale. = T)
PCA_fluvial_network_plot[[1]] <- PCA_fluvial_result
PCA_fluvial_network_results[[1]] <-fluvial_network_data[[1]][,1]

names(PCA_fluvial_network_results) <- c("PCA_fluvial_network")

#autoplot(PCA_fluvial_network_plot[[1]],loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
#         loadings.colour = 'red', loadings.label.colour="black")+
#  geom_point(size=2, alpha=0.1,color=CUNILLERA_cols("black"))+labs(title="~600km")+
#  theme_bw()

#Plot centrality values in river network
cols <-ConComp$membership 
cols[unlist(res.nv[[3]])[correspondence_BASINS]] <- "red"
cols <- ifelse(cols==1, "green",ifelse(cols==2,"blue","red"))
sizes <- ifelse(cols=="red",1.5,0.5)

rbPal <- viridis_pal(1,0,1)
#Plots closeness "all"
# Relevant that "png" output size will be big (nice arrows plotted)
#This adds a column of color values based on the y values
Col <- rbPal(length(fluvial_network_data[[1]][,1]))[as.numeric(cut(fluvial_network_data[[1]][,1],breaks = length(fluvial_network_data[[1]][,1])))]
Col[which(cols=="red")] <- "red"
png(filename = "C:/Users/Cunilleramontcusi/Alpine_fluvial_CLOS.png", width = 20000, height = 20000, res=500)
par(mar=c(0,0,0,0))
plot(GRAPH_xarxes_fluvial[[1]], vertex.label = NA, vertex.size = 1, vertex.size2 = sizes, vertex.color=Col,
     edge.width=0.01, edge.color="grey50")
dev.off()

#________________________________________#
#________________________________________#
########    DIVERSITY VALUES    ##########
#________________________________________#
#________________________________________#

#### Community indices list to gather everything 

community_indices<- list()

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

# ALPHA DIVERSITY
# Richness
tst_coeficient_out <- list()
for (r in 1:5) {
  ### WARNING: CHANGE TO JACCARD IF WORKING WITH PA
  rich <- apply(comm_data[[r]],1,sum)
  tst_coeficient_out[[r]] <- rich
}
community_indices[[2]] <- tst_coeficient_out
names(community_indices[[2]])<- c("Rich_S16","Rich_S18","Rich_PHY","Rich_ZOO","Rich_ZOO.18S") 

# BETA DIVERSITY
# LCBD
library(tidyr)
library(betareg)
library(vegan)
library(adespatial)
library(ade4)
tst_coeficient_out <- list()
for (r in 1:5) {
  ### WARNING: CHANGE TO JACCARD IF WORKING WITH PA
  LCBD_values <- beta.div(comm_data[[r]], method = "jaccard", nperm = 9999)
  tst_coeficient_out[[r]] <- LCBD_values$LCBD
}
community_indices[[3]] <- tst_coeficient_out
names(community_indices[[3]])<- c("LCBD_S16","LCBD_S18","LCBD_PHY","LCBD_ZOO","LCBD_18S.ZOO") 

# Betadiversity components
tst_coeficient_out <- list()
for (r in 1:5) {
  ### WARNING: CHANGE TO JACCARD IF WORKING WITH PA
  taxa.q_Ruziska_AB <- beta.div.comp(comm_data[[r]], coef = "J", quant = TRUE)
  replacemnet <- apply(as.matrix(taxa.q_Ruziska_AB$repl),1,mean)
  rich_diff <- apply(as.matrix(taxa.q_Ruziska_AB$rich),1,mean)
  beta_div_components <- cbind(replacemnet,rich_diff)
  tst_coeficient_out[[r]] <- beta_div_components
}

community_indices[[4]]<- list(tst_coeficient_out[[1]][,1],tst_coeficient_out[[2]][,1],
                              tst_coeficient_out[[3]][,1], tst_coeficient_out[[4]][,1],tst_coeficient_out[[5]][,1])
names(community_indices[[4]])<- c("Repl_S16","Repl_S18","Repl_PHY","Repl_ZOO","Repl_18S.ZOO") 
community_indices[[5]]<- list(tst_coeficient_out[[1]][,2],tst_coeficient_out[[2]][,2],
                              tst_coeficient_out[[3]][,2], tst_coeficient_out[[4]][,2],tst_coeficient_out[[5]][,2])
names(community_indices[[5]])<- c("RicDif_S16","RicDif_S18","RicDif_PHY","RicDif_ZOO","RicDif_18S.ZOO") 

biod <- list()
for (r in 1:5) {
  biod[[r]] <- cbind(community_indices[[1]][[r]],community_indices[[2]][[r]],community_indices[[3]][[r]],
                     community_indices[[4]][[r]],community_indices[[5]][[r]])
  noms <-c(names(community_indices[[1]])[r],names(community_indices[[2]])[r],names(community_indices[[3]])[r],
           names(community_indices[[4]])[r],names(community_indices[[5]])[r]) 
  colnames(biod[[r]]) <- noms
} 
biod

#________________________________________#
#________________________________________#
##########  NET vs DIV VALUES   ##########
#________________________________________#
#________________________________________#

# Extracting the values of which lakes have been sampled for each taxonomic group
# S16 = 52 lakes
# S18 = 48 lakes
# Phy = 50 lakes
# Zoo = 52 lakes
# 18S.Zoo = 48 lakes

biod_names <- c("S16","S18","phy","zoo", "zoo.18S")

coincidence <- c()
coincidence_values <- list()
for (r in 1:5) {
  for (coinc in 1:nrow(all_lakes_coord)) {
    coincidence_value <- 0
    coincidence_value<-which(rownames(all_lakes_coord)==rownames(env_data[[r]])[coinc])
    coincidence_value<-ifelse(length(coincidence_value)==0, 0 ,coincidence_value)
    coincidence[coinc] <- coincidence_value
  }  
  coincidence_values[[r]] <- coincidence
}


# GAM models______________________________________ ####
Names_Networks <- c("600 km", "300 km","100 km","60 km","6 km")

GAMmodel_resutls_total <- list()
GAM_model_resutls <- list()
GAM_direct_model <- list()
GAM_direct_model_total <- list()

# For Euclidean network
for (groups in 1:5) {
  color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")
  for (net in 1:5) {
    coin <- PCA_network_results[[net]][c(length(PCA_network_results[[net]])-54):length(PCA_network_results[[net]])]
    dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                     biod[[groups]][,1:5])
    colnames(dataset)[1] <-c("Network")
    plots_grups <- list()
    output_results <- list()
    output_model_results <- list()
    
    p.val <- c()
    #CCA
    p.val[1] <-summary.gam(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
    output_results[[1]] <- summary.gam(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
    output_model_results[[1]] <- gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
    preds_1 <- predict(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
    #Richness
    p.val[2] <- summary.gam(gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
    output_results[[2]] <- summary.gam(gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
    output_model_results[[2]] <- gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
    preds_2<- predict(gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
    #LCBD
    p.val[3] <-summary.gam(gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
    output_results[[3]] <- summary.gam(gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
    output_model_results[[3]] <- gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
    preds_3<- predict(gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
    #Turn
    p.val[4] <-summary.gam(gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
    output_results[[4]] <- summary.gam(gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
    output_model_results[[4]] <- gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
    preds_4<- predict(gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
    #RichDiff
    p.val[5] <-summary.gam(gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
    output_results[[5]] <- summary.gam(gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
    output_model_results[[5]] <- gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
    preds_5<- predict(gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
    
    GAM.pred <- list(preds_1, preds_2, preds_3, preds_4, preds_5)
    
    for(var in 1:5){
      if(p.val[var]>0.05){
        my_data <- data.frame(cbind(dataset[,var+1],dataset[,1]),
                              mu   = GAM.pred[[var]]$fit,
                              low  = GAM.pred[[var]]$fit - 1.96 * GAM.pred[[var]]$se.fit,
                              high = GAM.pred[[var]]$fit + 1.96 * GAM.pred[[var]]$se.fit)
        plots_grups[[var]] <-
          ggplot(my_data, aes(x = X2, y = X1)) +
          geom_jitter(alpha=0.2, shape=21, size=3, colour="black", aes(fill=X2))+
          scale_fill_continuous(type = "viridis")+
          geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="black",linetype=2, size=2)+
          labs(title=colnames(dataset)[var+1])+ylab(colnames(dataset)[var+1])+xlab("Centrality-Isolation")+
          theme_classic()+
          theme(panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)),
                legend.position = "none") 
      }else{
        my_data <- data.frame(cbind(dataset[,var+1],dataset[,1]),
                              mu   = GAM.pred[[var]]$fit,
                              low  = GAM.pred[[var]]$fit - 1.96 * GAM.pred[[var]]$se.fit,
                              high = GAM.pred[[var]]$fit + 1.96 * GAM.pred[[var]]$se.fit)
        plots_grups[[var]] <-ggplot(my_data, aes(x = X2, y = X1)) +
          geom_jitter(alpha=0.9, shape=21, size=3, colour="black", aes(fill=X2))+
          scale_fill_continuous(type = "viridis")+
          geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="black",linetype=1, size=2)+
          labs(title=colnames(dataset)[var+1])+ylab(colnames(dataset)[var+1])+xlab("Centrality-Isolation")+
          theme_classic()+
          theme(panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)),
                legend.position = "none")  
      }
    }
    
    GAM_model_resutls[[net]] <-output_results
    GAM_direct_model[[net]]<- output_model_results
    
    png(filename =paste("C:/Users/Cunilleramontcusi/","GAM_Divers",biod_names[[groups]],"_",names(PCA_network_results)[[net]],".png"),
        width =582*2 ,height =629*2 ,units ="px",res = 200)
    grid.arrange(plots_grups[[1]],plots_grups[[2]],
                 plots_grups[[3]],plots_grups[[4]],
                 plots_grups[[5]],
                 ncol=2,nrow=3, top=names(PCA_network_results)[[net]])
    dev.off()
  }
  GAMmodel_resutls_total[[groups]] <- GAM_model_resutls
  GAM_direct_model_total[[groups]]<- GAM_direct_model
}             

# For fluvial network
GAMmodel_resutls_fluvial <- list()
GAM_direct_model_fluvial<- list()

for (groups in 1:5) {
  color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")
  coin <- PCA_fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
  dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                   biod[[groups]][,1:5])
  colnames(dataset)[1] <-c("Network")
  plots_grups <- list()
  output_results <- list()
  output_model_results <- list()
  
  p.val <- c()
  #CCA
  p.val[1] <-summary.gam(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
  output_results[[1]] <- summary.gam(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
  output_model_results[[1]] <- gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
  preds_1 <- predict(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
  #Richness
  p.val[2] <- summary.gam(gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
  output_results[[2]] <- summary.gam(gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
  output_model_results[[2]] <- gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
  preds_2<- predict(gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
  #LCBD
  p.val[3] <-summary.gam(gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
  output_results[[3]] <- summary.gam(gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
  output_model_results[[3]] <- gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
  preds_3<- predict(gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
  #Turn
  p.val[4] <-summary.gam(gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
  output_results[[4]] <- summary.gam(gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
  output_model_results[[4]] <- gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
  preds_4<- predict(gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
  #RichDiff
  p.val[5] <-summary.gam(gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
  output_results[[5]] <- summary.gam(gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
  output_model_results[[5]] <- gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
  preds_5<- predict(gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
  
  GAM.pred <- list(preds_1, preds_2, preds_3, preds_4, preds_5)
  
  for(var in 1:5){
    if(p.val[var]>0.05){
      my_data <- data.frame(cbind(dataset[,var+1],dataset[,1]),
                            mu   = GAM.pred[[var]]$fit,
                            low  = GAM.pred[[var]]$fit - 1.96 * GAM.pred[[var]]$se.fit,
                            high = GAM.pred[[var]]$fit + 1.96 * GAM.pred[[var]]$se.fit)
      plots_grups[[var]] <-
        ggplot(my_data, aes(x = X2, y = X1)) +
        geom_jitter(alpha=0.2, shape=21, size=3, colour="black", aes(fill=X2))+
        scale_fill_continuous(type = "viridis")+
        geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="black",linetype=2, size=2)+
        labs(title=colnames(dataset)[var+1])+ylab(colnames(dataset)[var+1])+xlab("Centrality-Isolation")+
        theme_classic()+
        theme(panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)),
              legend.position = "none")  
    }else{
      my_data <- data.frame(cbind(dataset[,var+1],dataset[,1]),
                            mu   = GAM.pred[[var]]$fit,
                            low  = GAM.pred[[var]]$fit - 1.96 * GAM.pred[[var]]$se.fit,
                            high = GAM.pred[[var]]$fit + 1.96 * GAM.pred[[var]]$se.fit)
      plots_grups[[var]] <-ggplot(my_data, aes(x = X2, y = X1)) +
        geom_jitter(alpha=0.9, shape=21, size=3, colour="black", aes(fill=X2))+
        scale_fill_continuous(type = "viridis")+
        geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="black",linetype=1, size=2)+
        labs(title=colnames(dataset)[var+1])+ylab(colnames(dataset)[var+1])+xlab("Centrality-Isolation")+
        theme_classic()+
        theme(panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)),
              legend.position = "none")  
    }
  }
  
  GAMmodel_resutls_fluvial[[groups]] <-output_results
  GAM_direct_model_fluvial[[groups]] <-output_model_results
  
  png(filename =paste("C:/Users/Cunilleramontcusi/","GAM_Divers",biod_names[[groups]],"_Fluvial",".png"),
      width =582*2 ,height =629*2 ,units ="px",res = 200)
  grid.arrange(plots_grups[[1]],plots_grups[[2]],
               plots_grups[[3]],plots_grups[[4]],
               plots_grups[[5]],
               ncol=2,nrow=3, top="Fluvial network")
  dev.off()
}



# GAM models result in table format - Supplementary like_####
biod_names <- c("S16","S18","phy","zoo", "zoo.18S")
Names_Networks <- c("600 km", "300 km","100 km","60 km","6 km", "Fluvial")

# Add the fluvial as a sixth network
for (t in 1:5) {
GAMmodel_resutls_total[[t]][[6]] <- GAMmodel_resutls_fluvial[[t]]
}

table_groups <- list()
for (groups in 1:5) {
  a <- matrix(nrow =length(GAMmodel_resutls_total)*6,ncol = 10 )
  row_reference <- 0
  colnames(a) <- c("Network","Variable",
                   "Intercept","Std.Err.","t-value","p-value","Smooth F-value","Smooth p-value","Adj R-sqr","Expl.Deviance")
    for (net in 1:6) {
      for (var in 1:5) {
        row_reference <- row_reference+1
        
        Netw_value <- Names_Networks[net]
        Variable_value<- colnames(biod[[groups]])[var]
        
        a[row_reference,] <- c(Netw_value,Variable_value,
                               round(GAMmodel_resutls_total[[groups]][[net]][[var]]$p.coeff,2),
                               round(GAMmodel_resutls_total[[groups]][[net]][[var]]$se[1],2),
                               round(GAMmodel_resutls_total[[groups]][[net]][[var]]$p.t,2),
                               round(GAMmodel_resutls_total[[groups]][[net]][[var]]$p.pv,2),
                               round(GAMmodel_resutls_total[[groups]][[net]][[var]]$chi.sq,2),
                               round(GAMmodel_resutls_total[[groups]][[net]][[var]]$s.pv,2),
                               round(GAMmodel_resutls_total[[groups]][[net]][[var]]$r.sq,2),
                               round(GAMmodel_resutls_total[[groups]][[net]][[var]]$dev.expl,2)) 
  }
 }
  table_groups[[groups]] <- as.data.frame(a)
  write.table(table_groups[[groups]], file = paste(biod_names[groups],"GAM_Results",".txt",sep = ""), sep = ",", quote = FALSE, row.names = F)
}



# Summary plot GAM models______________________________________ ####
sign_netw <- list()
sign_groups <- list()
for (group in 1:5) {
  for (netw in 1:5) {
    sign <- c()
    sign[1] <-GAMmodel_resutls_total[[group]][[netw]][[1]][[8]]
    sign[2] <-GAMmodel_resutls_total[[group]][[netw]][[2]][[8]]
    sign[3] <-GAMmodel_resutls_total[[group]][[netw]][[3]][[8]]
    sign[4] <-GAMmodel_resutls_total[[group]][[netw]][[4]][[8]]
    sign[5] <-GAMmodel_resutls_total[[group]][[netw]][[5]][[8]] 
    sign_netw[[netw]] <- sign
  }
  sign_groups[[group]] <- sign_netw
}

flu_sign_groups <- list()
for (group in 1:5) {
  sign <- c()
  sign[1] <-GAMmodel_resutls_fluvial[[group]][[1]][[8]]
  sign[2] <-GAMmodel_resutls_fluvial[[group]][[2]][[8]]
  sign[3] <-GAMmodel_resutls_fluvial[[group]][[3]][[8]]
  sign[4] <-GAMmodel_resutls_fluvial[[group]][[4]][[8]]
  sign[5] <-GAMmodel_resutls_fluvial[[group]][[5]][[8]]
  flu_sign_groups[[group]] <- sign
}


Names_Variab <- c("Environmental tracking", "Species richness", "LCBD", "Replacement", "Richness difference")
plots_significance <- list()
for (variable in 1:5) {
  
  max_netw <- cbind(c(sign_groups[[1]][[1]][[variable]], sign_groups[[2]][[1]][[variable]], sign_groups[[3]][[1]][[variable]], 
                      sign_groups[[4]][[1]][[variable]], sign_groups[[5]][[1]][[variable]]),
                    rep("600 km", 5), 
                    c("S16","S18","Phy","Zoo", "S18zoo"))
  
  mid_netw <- cbind(c(sign_groups[[1]][[2]][[variable]], sign_groups[[2]][[2]][[variable]], sign_groups[[3]][[2]][[variable]], 
                      sign_groups[[4]][[2]][[variable]], sign_groups[[5]][[2]][[variable]]),
                    rep("300 km", 5), 
                    c("S16","S18","Phy","Zoo", "S18zoo"))
  
  mid_mid_netw <- cbind(c(sign_groups[[1]][[3]][[variable]], sign_groups[[2]][[3]][[variable]], sign_groups[[3]][[3]][[variable]], 
                          sign_groups[[4]][[3]][[variable]], sign_groups[[5]][[3]][[variable]]),
                        rep("100 km", 5), 
                        c("S16","S18","Phy","Zoo", "S18zoo"))
  
  small_netw <- cbind(c(sign_groups[[1]][[4]][[variable]], sign_groups[[2]][[4]][[variable]], sign_groups[[3]][[4]][[variable]], 
                        sign_groups[[4]][[4]][[variable]], sign_groups[[5]][[4]][[variable]]),
                      rep("60 km", 5), 
                      c("S16","S18","Phy","Zoo", "S18zoo"))
  
  min_netw <- cbind(c(sign_groups[[1]][[5]][[variable]], sign_groups[[2]][[5]][[variable]], sign_groups[[3]][[5]][[variable]], 
                      sign_groups[[4]][[5]][[variable]], sign_groups[[5]][[5]][[variable]]),
                    rep("6 km", 5), 
                    c("S16","S18","Phy","Zoo", "S18zoo"))
  
  fluv_netw <- cbind(c(flu_sign_groups[[1]][[variable]], flu_sign_groups[[2]][[variable]], flu_sign_groups[[3]][[variable]], 
                       flu_sign_groups[[4]][[variable]], flu_sign_groups[[5]][[variable]]),
                     rep("Fluvial", 5), 
                     c("S16","S18","Phy","Zoo", "S18zoo"))
  
  
  dataset_pval <- as.data.frame(rbind(max_netw,mid_netw,mid_mid_netw,small_netw,min_netw,fluv_netw))
  colnames(dataset_pval) <- c("val","Network","Group")
  dataset_pval$val <-as.numeric(dataset_pval$val)
  dataset_pval$Network <- factor(dataset_pval$Network,
                                 levels = c("600 km", "300 km","100 km","60 km","6 km", "Fluvial"))
  dataset_pval$Group <- factor(dataset_pval$Group,
                               levels = c("S16","S18","Phy","Zoo", "S18zoo"))
  significants <- rep("NoSign",5*6)
  significants[which(dataset_pval$val<0.05)]<- "Sign"
  dataset_pval$Sign <- factor(significants) 
  
  color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")
  
  plots_significance[[variable]] <-  ggplot(dataset_pval, aes(x=Network, y=as.numeric(val)))+
    geom_abline(slope = 0,intercept = 0.05, colour="black", linetype=2,size=1)+
    geom_jitter(aes(fill=Group, alpha=Sign, size=Sign),shape=21,width = 0.5)+
    scale_alpha_manual(values = c(0.2,0.9))+
    scale_size_manual(values = c(2,6))+
    scale_fill_manual(values=c(color_groups[1],color_groups[2],color_groups[3],
                               color_groups[4],color_groups[5],color_groups[6]))+
    scale_y_continuous(expand = c(0.2,0.01),
                       breaks =c(0.2,0.4,0.6,0.8,1) )+
    geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5), size=1, colour="grey70")+
    labs(title=Names_Variab[variable])+ylab("P-values")+xlab("Network")+
    theme_classic()
}

png(filename =paste("C:/Users/Cunilleramontcusi/All_Significance_Groups.png"),
    width =900*5 ,height =700*5 ,units ="px",res = 400)
grid.arrange(plots_significance[[1]],plots_significance[[2]],
             plots_significance[[3]],plots_significance[[4]],
             plots_significance[[5]],
             ncol=2,nrow=3, top="Significance values")
dev.off()      




# GAM plot significant______________________________________ ####
# For Euclidean network
GAM_Sign_plots_total <- list()
ref_value <- 0
for (groups in 1:5) {
  color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")
  for (net in 1:5) {
    coin <- PCA_network_results[[net]][c(length(PCA_network_results[[net]])-54):length(PCA_network_results[[net]])]
    dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                     biod[[groups]][,1:5])
    colnames(dataset)[1] <-c("Network")
    plots_grups <- list()
    
    p.val <- c()
    #CCA
    p.val[1] <-summary.gam(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
    preds_1 <- predict(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
    #Richness
    p.val[2] <- summary.gam(gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
    preds_2<- predict(gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
    #LCBD
    p.val[3] <-summary.gam(gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
    preds_3<- predict(gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
    #Turn
    p.val[4] <-summary.gam(gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
    preds_4<- predict(gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
    #RichDiff
    p.val[5] <-summary.gam(gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
    preds_5<- predict(gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
    
    GAM.pred <- list(preds_1, preds_2, preds_3, preds_4, preds_5)
    
    select_p.val <- which(p.val<0.05)
    
    if(length(select_p.val)>0){
      for(var in 1:length(select_p.val)){
        ref_value <- ref_value+1
        
        my_data <- data.frame(cbind(dataset[,select_p.val[var]+1],dataset[,1]),
                              mu   = GAM.pred[[select_p.val[var]]]$fit,
                              low  = GAM.pred[[select_p.val[var]]]$fit - 1.96 * GAM.pred[[select_p.val[var]]]$se.fit,
                              high = GAM.pred[[select_p.val[var]]]$fit + 1.96 * GAM.pred[[select_p.val[var]]]$se.fit)
        GAM_Sign_plots_total[[ref_value]] <-ggplot(my_data, aes(x = X2, y = X1)) +
          geom_jitter(alpha=0.9, shape=21, size=3, colour="black", aes(fill=X2))+
          scale_fill_continuous(type = "viridis")+
          geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="black",linetype=1, size=2)+
          labs(title=paste(Names_Networks[[net]],colnames(dataset)[select_p.val[var]+1]))+
          ylab(colnames(dataset)[var+1])+xlab("Centrality-Isolation")+
          theme_classic()+
          theme(panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)),
                legend.position = "none")  
      }
    }
  }
}             
# For Fluvial network
GAM_Sign_plots_total_Fluvial <- list()
ref_value <- 0
for (groups in 1:5) {
  color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")
  coin <- PCA_fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
  dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                   biod[[groups]][,1:5])
  colnames(dataset)[1] <-c("Network")
  plots_grups <- list()
  output_results <- list()
  p.val <- c()
  
  #CCA
  p.val[1] <-summary.gam(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
  preds_1 <- predict(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
  #Richness
  p.val[2] <- summary.gam(gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
  preds_2<- predict(gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
  #LCBD
  p.val[3] <-summary.gam(gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
  preds_3<- predict(gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
  #Turn
  p.val[4] <-summary.gam(gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
  preds_4<- predict(gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
  #RichDiff
  p.val[5] <-summary.gam(gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
  preds_5<- predict(gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
  
  GAM.pred <- list(preds_1, preds_2, preds_3, preds_4, preds_5)
  
  select_p.val <- which(p.val<0.05)
  
  if(length(select_p.val)>0){
    for(var in 1:length(select_p.val)){
      ref_value <- ref_value+1
      
      my_data <- data.frame(cbind(dataset[,select_p.val[var]+1],dataset[,1]),
                            mu   = GAM.pred[[select_p.val[var]]]$fit,
                            low  = GAM.pred[[select_p.val[var]]]$fit - 1.96 * GAM.pred[[select_p.val[var]]]$se.fit,
                            high = GAM.pred[[select_p.val[var]]]$fit + 1.96 * GAM.pred[[select_p.val[var]]]$se.fit)
      
      GAM_Sign_plots_total_Fluvial[[ref_value]] <-ggplot(my_data, aes(x = X2, y = X1)) +
        geom_jitter(alpha=0.9, shape=21, size=3, colour="black", aes(fill=X2))+
        scale_fill_continuous(type = "viridis")+
        geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="black",linetype=1, size=2)+
        labs(title=paste("Fluvial",colnames(dataset)[select_p.val[var]+1]))+
        ylab(colnames(dataset)[select_p.val[var]+1])+xlab("Centrality-Isolation")+
        theme_classic()+
        theme(panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)),
              legend.position = "none")  
      
    }
  }
}

# Plot extraction
plot_plot_sign_plot <- list()
for(t in 1:length(GAM_Sign_plots_total)){
  plot_plot_sign_plot[[t]] <- GAM_Sign_plots_total[[t]]
}
for(d in c(length(GAM_Sign_plots_total)+1):c(length(GAM_Sign_plots_total)+length(GAM_Sign_plots_total_Fluvial))){
  plot_plot_sign_plot[[d]] <- GAM_Sign_plots_total_Fluvial[[d-length(GAM_Sign_plots_total)]]
}

# Printing for Env_Tracking
png(filename ="C:/Users/Cunilleramontcusi/GAM_Sign_EnvTrack.png",
    width =582*4 ,height =629*2 ,units ="px",res = 300)
grid.arrange(plot_plot_sign_plot[[1]],plot_plot_sign_plot[[2]],
             ncol=2,nrow=1, top="Assembly (Environmental Tracking)")
dev.off()

# Printing for Diverse
png(filename ="C:/Users/Cunilleramontcusi/GAM_Sign_Diverse.png",
    width =629*6 ,height =629*6 ,units ="px",res = 300)
grid.arrange(plot_plot_sign_plot[[3]],plot_plot_sign_plot[[5]], plot_plot_sign_plot[[6]],
             plot_plot_sign_plot[[7]],plot_plot_sign_plot[[8]],plot_plot_sign_plot[[9]],plot_plot_sign_plot[[19]],
             plot_plot_sign_plot[[10]],plot_plot_sign_plot[[12]],
             plot_plot_sign_plot[[13]],plot_plot_sign_plot[[14]],plot_plot_sign_plot[[15]],
             plot_plot_sign_plot[[16]],plot_plot_sign_plot[[17]],plot_plot_sign_plot[[18]],
             ncol=4,nrow=4, top="Diversity")
dev.off()



# NMDS plots______________________________________ ####
biod_names <- c("S16","S18","phy","zoo", "zoo.18S")
net_names <- c("600 km", "300 km","100 km","60 km","6 km", "Fluvial")
color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")

# For Euclidean network
plots_NMDS <- list()
plots_NMDS_total <- list()

plots_NMDS_model_result <- list()
plots_NMDS_total_model_result <- list()

for (groups in 1:5) {
  for (net in 1:5) {
coin <- PCA_network_results[[net]][c(length(PCA_network_results[[net]])-54):length(PCA_network_results[[net]])]
centr_iso <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]])

spe.abu.jac<-vegdist(comm_data[[groups]],method = "jaccard") #Càlcul de la matriu de similitud amb Bray-Curtis
spe.abu.MDS<-metaMDS(spe.abu.jac, distance="jaccard", k=2, try = 999) #MDS amb similitud bray-curtis

x <- spe.abu.MDS$points[,1]
y <- spe.abu.MDS$points[,2]

dataset <- data.frame(x,y,centr_iso)
NMDS_model <- ordisurf(spe.abu.MDS ~ dataset[,3], plot = F,)
NMDS_model_results <- summary(NMDS_model)

plots_NMDS_model_result[[net]] <- NMDS_model_results

extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("x", "y", "z")
  return(xyz)
}
contour.vals <- extract.xyz(obj = NMDS_model)

    plots_NMDS[[net]] <- ggplot(dataset, aes(x=x,y=y))+
                  geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+
                  geom_jitter(shape=21, size=5, alpha=0.8, aes(fill=centr_iso))+
                  geom_contour(data=contour.vals, aes(x, y, z = z, colour = ..level..))+
                  #geom_path(data = df_ellipse, aes(x=x, y=y, colour=Group), size=2, show.legend = FALSE)+
                  scale_colour_viridis(alpha = 1,begin = 1,end = 0)+
                  #manual(values = viridis_pal(0.9,1,0,direction = -1)(length(unique(df_ellipse$Group))))+
                  scale_fill_viridis(alpha = 1,begin = 1,end = 0)+
                  labs(title = paste(biod_names[groups], net_names[net]),
                       subtitle = paste("Stress=",round(spe.abu.MDS$stress,2), 
                                        #"Rsqr=", round(NMDS_model_results$r.sq,2),
                                        "Expl. Dev.=",round(NMDS_model_results$dev.expl,2)))+
                  xlab("NMDS1")+ylab("NMDS2")+
                  theme_classic()+
                  theme(legend.position = "none",
                        panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)))
  }
  plots_NMDS_total[[groups]] <- plots_NMDS
  plots_NMDS_total_model_result[[groups]] <- plots_NMDS_model_result
}

# For Fluvial network
plots_NMDS_fluvial<- list()
plots_NMDS_total_fluvial <- list()

plots_NMDS_fluvial_model_result <- list()
plots_NMDS_fluvial_total_model_result <- list()

for (groups in 1:5) {
    coin <- PCA_fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
    centr_iso <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]])
    
    spe.abu.jac<-vegdist(comm_data[[groups]],method = "jaccard") #Càlcul de la matriu de similitud amb Bray-Curtis
    spe.abu.MDS<-metaMDS(spe.abu.jac, distance="jaccard", k=2, try = 999) #MDS amb similitud bray-curtis
    
    x <- spe.abu.MDS$points[,1]
    y <- spe.abu.MDS$points[,2]
    
    dataset <- data.frame(x,y,centr_iso)
    NMDS_model <- ordisurf(spe.abu.MDS ~ dataset[,3], plot = F)
    NMDS_model_results <- summary(NMDS_model)
    
    plots_NMDS_fluvial_model_result[[groups]] <- NMDS_model_results
    
    extract.xyz <- function(obj) {
      xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
      xyz <- cbind(xy, c(obj$grid$z))
      names(xyz) <- c("x", "y", "z")
      return(xyz)
    }
    contour.vals <- extract.xyz(obj = NMDS_model)
    
    plots_NMDS_fluvial[[groups]] <- ggplot(dataset, aes(x=x,y=y))+
      geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+
      geom_jitter(shape=21, size=5, alpha=0.8, aes(fill=centr_iso))+
      geom_contour(data=contour.vals, aes(x, y, z = z, colour = ..level..))+
      scale_colour_viridis(alpha = 1,begin = 1,end = 0)+
      scale_fill_viridis(alpha = 1,begin = 1,end = 0)+
      labs(title = paste(biod_names[groups], net_names[net]),
           subtitle = paste("Stress=",round(spe.abu.MDS$stress,2), 
                            "Rsqr=", round(NMDS_model_results$r.sq,2),
                            "Expl. Dev.=",round(NMDS_model_results$dev.expl,2)))+
      xlab("NMDS1")+ylab("NMDS2")+
      theme_classic()+
      theme(legend.position = "none",
            panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)))
    
  plots_NMDS_total_fluvial[[groups]] <- plots_NMDS_fluvial
  plots_NMDS_fluvial_total_model_result[[groups]] <- plots_NMDS_fluvial_model_result
}


# Print NMDS
png(filename ="C:/Users/Cunilleramontcusi/NMDS_Diverse.png",
    width =729*10, height =629*10 ,units ="px",res = 300)
grid.arrange(plots_NMDS_total[[1]][[1]],plots_NMDS_total[[1]][[2]],plots_NMDS_total[[1]][[3]],
             plots_NMDS_total[[1]][[4]],plots_NMDS_total[[1]][[5]],plots_NMDS_total_fluvial[[1]][[1]],
             
             plots_NMDS_total[[2]][[1]],plots_NMDS_total[[2]][[2]],plots_NMDS_total[[2]][[3]],
             plots_NMDS_total[[2]][[4]],plots_NMDS_total[[2]][[5]],plots_NMDS_total_fluvial[[2]][[2]],
             
             plots_NMDS_total[[3]][[1]],plots_NMDS_total[[3]][[2]],plots_NMDS_total[[3]][[3]],
             plots_NMDS_total[[3]][[4]],plots_NMDS_total[[3]][[5]],plots_NMDS_total_fluvial[[3]][[3]],
             
             plots_NMDS_total[[4]][[1]],plots_NMDS_total[[4]][[2]],plots_NMDS_total[[4]][[3]],
             plots_NMDS_total[[4]][[4]],plots_NMDS_total[[4]][[5]],plots_NMDS_total_fluvial[[4]][[4]],
             
             plots_NMDS_total[[5]][[1]],plots_NMDS_total[[5]][[2]],plots_NMDS_total[[5]][[3]],
             plots_NMDS_total[[5]][[4]],plots_NMDS_total[[5]][[5]],plots_NMDS_total_fluvial[[5]][[5]],
             
             ncol=6,nrow=5, top="NMDS groups")
dev.off()


# NMDS plots significant Ordisurfs______________________________________ ####

# For Euclidean network
plots_NMDS_sign <- list()

ref_value <- 0

for (groups in 1:5) {
  for (net in 1:5) {
    coin <- PCA_network_results[[net]][c(length(PCA_network_results[[net]])-54):length(PCA_network_results[[net]])]
    centr_iso <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]])
    
    spe.abu.jac<-vegdist(comm_data[[groups]],method = "jaccard") #Càlcul de la matriu de similitud amb Bray-Curtis
    spe.abu.MDS<-metaMDS(spe.abu.jac, distance="jaccard", k=2, try = 999) #MDS amb similitud bray-curtis
    
    x <- spe.abu.MDS$points[,1]
    y <- spe.abu.MDS$points[,2]
    
    dataset <- data.frame(x,y,centr_iso)
    NMDS_model <- ordisurf(spe.abu.MDS ~ dataset[,3], plot = F)
    NMDS_model_results <- summary(NMDS_model)
    
    plots_NMDS_model_pval <- NMDS_model_results$s.pv
    
    if (plots_NMDS_model_pval<0.05) {
      
      ref_value <- ref_value+1
      
      extract.xyz <- function(obj) {
      xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
      xyz <- cbind(xy, c(obj$grid$z))
      names(xyz) <- c("x", "y", "z")
      return(xyz)
    }
    contour.vals <- extract.xyz(obj = NMDS_model)
    
    plots_NMDS_sign[[ref_value]] <- ggplot(dataset, aes(x=x,y=y))+
      geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+
      geom_jitter(shape=21, size=5, alpha=0.8, aes(fill=centr_iso))+
      geom_contour(data=contour.vals, aes(x, y, z = z, colour = ..level..))+
      scale_colour_viridis(alpha = 1,begin = 1,end = 0)+
      scale_fill_viridis(alpha = 1,begin = 1,end = 0)+
      labs(title = paste(biod_names[groups], net_names[net]),
           subtitle = paste("Stress=",round(spe.abu.MDS$stress,2), 
                            "Rsqr=", round(NMDS_model_results$r.sq,2),
                            "Expl. Dev.=",round(NMDS_model_results$dev.expl,2)))+
      xlab("NMDS1")+ylab("NMDS2")+
      theme_classic()+
      theme(legend.position = "none",
            panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1))) 
    }
  }
  plots_NMDS_sign
}

# For Fluvial network
plots_NMDS_fluvial_sign<- list()

ref_value <- 0

for (groups in 1:5) {
  coin <- PCA_fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
  centr_iso <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]])
  
  spe.abu.jac<-vegdist(comm_data[[groups]],method = "jaccard") #Càlcul de la matriu de similitud amb Bray-Curtis
  spe.abu.MDS<-metaMDS(spe.abu.jac, distance="jaccard", k=2, try = 999) #MDS amb similitud bray-curtis
  
  x <- spe.abu.MDS$points[,1]
  y <- spe.abu.MDS$points[,2]
  
  dataset <- data.frame(x,y,centr_iso)
  NMDS_model <- ordisurf(spe.abu.MDS ~ dataset[,3], plot = F)
  NMDS_model_results <- summary(NMDS_model)
  
  plots_NMDS_model_pval <- NMDS_model_results$s.pv
  
  if (plots_NMDS_model_pval<0.05) {
    
    ref_value <- ref_value+1
    
    extract.xyz <- function(obj) {
      xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
      xyz <- cbind(xy, c(obj$grid$z))
      names(xyz) <- c("x", "y", "z")
      return(xyz)
    }
    contour.vals <- extract.xyz(obj = NMDS_model)
  
    plots_NMDS_fluvial_sign[[ref_value]] <- ggplot(dataset, aes(x=x,y=y))+
    geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+
    geom_jitter(shape=21, size=5, alpha=0.8, aes(fill=centr_iso))+
    geom_contour(data=contour.vals, aes(x, y, z = z, colour = ..level..))+
    scale_colour_viridis(alpha = 1,begin = 1,end = 0)+
    #manual(values = viridis_pal(0.9,1,0,direction = -1)(length(unique(df_ellipse$Group))))+
    scale_fill_viridis(alpha = 1,begin = 1,end = 0)+
    labs(title = paste(biod_names[groups], "Fluvial"),
         subtitle = paste("Stress=",round(spe.abu.MDS$stress,2), 
                          "Rsqr=", round(NMDS_model_results$r.sq,2),
                          "Expl. Dev.=",round(NMDS_model_results$dev.expl,2)))+
    xlab("NMDS1")+ylab("NMDS2")+
    theme_classic()+
    theme(legend.position = "none",
          panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)))
  }
plots_NMDS_fluvial_sign
}

png(filename ="C:/Users/Cunilleramontcusi/NMDS_Diverse_Sign.png",
    width =429*10, height =629*10 ,units ="px",res = 300)
grid.arrange(plots_NMDS_sign[[1]],plots_NMDS_sign[[2]],plots_NMDS_sign[[3]],plots_NMDS_sign[[4]],
             plots_NMDS_sign[[5]],plots_NMDS_sign[[6]],plots_NMDS_sign[[7]],plots_NMDS_sign[[8]],
             plots_NMDS_sign[[9]],
             plots_NMDS_fluvial_sign[[1]],plots_NMDS_fluvial_sign[[2]],
             plots_NMDS_fluvial_sign[[3]],plots_NMDS_fluvial_sign[[4]],
             ncol=3,nrow=5, top="Composition (NMDS)")
dev.off()

# GAM NMDS models result in table format - Supplementary like_####
biod_names <- c("S16","S18","phy","zoo", "zoo.18S")
Names_Networks <- c("600 km", "300 km","100 km","60 km","6 km", "Fluvial")

# Add the fluvial as a sixth network
for (t in 1:5) {
  plots_NMDS_total_model_result[[t]][[6]] <- plots_NMDS_fluvial_total_model_result[[5]][[t]]
}

# 
table_groups <- list()
for (groups in 1:5) {
  a <- matrix(nrow =length(plots_NMDS_total_model_result[[groups]]),ncol = 10 )
  row_reference <- 0
  colnames(a) <- c("Network","Variable",
                   "Intercept","Std.Err.","t-value","p-value","Smooth F-value","Smooth p-value","Adj R-sqr","Expl.Deviance")
  for (net in 1:6) {
      row_reference <- row_reference+1
      
      Netw_value <- Names_Networks[net]
      Variable_value<- "NMDS X"
      
      a[row_reference,] <- c(Netw_value,Variable_value,
                             round(plots_NMDS_total_model_result[[groups]][[net]]$p.coeff,2),
                             round(plots_NMDS_total_model_result[[groups]][[net]]$se[1],2),
                             round(plots_NMDS_total_model_result[[groups]][[net]]$p.t,2),
                             round(plots_NMDS_total_model_result[[groups]][[net]]$p.pv,2),
                             round(plots_NMDS_total_model_result[[groups]][[net]]$chi.sq,2),
                             round(plots_NMDS_total_model_result[[groups]][[net]]$s.pv,2),
                             round(plots_NMDS_total_model_result[[groups]][[net]]$r.sq,2),
                             round(plots_NMDS_total_model_result[[groups]][[net]]$dev.expl,2)) 
  }
  table_groups[[groups]] <- as.data.frame(a)
  write.table(table_groups[[groups]], file = paste(biod_names[groups],"GAM_NMDS_Results",".txt",sep = ""), sep = ",", quote = FALSE, row.names = F)
}


##########        END           ##########
##########                      ##########
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################

# Old versions

## Print NMDS x
#png(filename ="C:/Users/Cunilleramontcusi/X_NMDS_Diverse.png",
#    width =729*10, height =629*10 ,units ="px",res = 300)
#grid.arrange(plots_NMDS_total_x[[1]][[1]],plots_NMDS_total_x[[1]][[2]],plots_NMDS_total_x[[1]][[3]],
#             plots_NMDS_total_x[[1]][[4]],plots_NMDS_total_x[[1]][[5]],plots_NMDS_total_x_fluvial[[1]][[1]],
#             
#             plots_NMDS_total_x[[2]][[1]],plots_NMDS_total_x[[2]][[2]],plots_NMDS_total_x[[2]][[3]],
#             plots_NMDS_total_x[[2]][[4]],plots_NMDS_total_x[[2]][[5]],plots_NMDS_total_x_fluvial[[2]][[2]],
#             
#             plots_NMDS_total_x[[3]][[1]],plots_NMDS_total_x[[3]][[2]],plots_NMDS_total_x[[3]][[3]],
#             plots_NMDS_total_x[[3]][[4]],plots_NMDS_total_x[[3]][[5]],plots_NMDS_total_x_fluvial[[3]][[3]],
#             
#             plots_NMDS_total_x[[4]][[1]],plots_NMDS_total_x[[4]][[2]],plots_NMDS_total_x[[4]][[3]],
#             plots_NMDS_total_x[[4]][[4]],plots_NMDS_total_x[[4]][[5]],plots_NMDS_total_x_fluvial[[4]][[4]],
#             
#             plots_NMDS_total_x[[5]][[1]],plots_NMDS_total_x[[5]][[2]],plots_NMDS_total_x[[5]][[3]],
#             plots_NMDS_total_x[[5]][[4]],plots_NMDS_total_x[[5]][[5]],plots_NMDS_total_x_fluvial[[5]][[5]],
#             
#             ncol=6,nrow=5, top="NMDS groups")
#dev.off()

## Print NMDS y
#png(filename ="C:/Users/Cunilleramontcusi/Y_NMDS_Diverse.png",
#    width =729*10, height =629*10 ,units ="px",res = 300)
#grid.arrange(plots_NMDS_total_y[[1]][[1]],plots_NMDS_total_y[[1]][[2]],plots_NMDS_total_y[[1]][[3]],
#             plots_NMDS_total_y[[1]][[4]],plots_NMDS_total_y[[1]][[5]],plots_NMDS_total_y_fluvial[[1]][[1]],
#             
#             plots_NMDS_total_y[[2]][[1]],plots_NMDS_total_y[[2]][[2]],plots_NMDS_total_y[[2]][[3]],
#             plots_NMDS_total_y[[2]][[4]],plots_NMDS_total_y[[2]][[5]],plots_NMDS_total_y_fluvial[[2]][[2]],
#             
#             plots_NMDS_total_y[[3]][[1]],plots_NMDS_total_y[[3]][[2]],plots_NMDS_total_y[[3]][[3]],
#             plots_NMDS_total_y[[3]][[4]],plots_NMDS_total_y[[3]][[5]],plots_NMDS_total_y_fluvial[[3]][[3]],
#             
#             plots_NMDS_total_y[[4]][[1]],plots_NMDS_total_y[[4]][[2]],plots_NMDS_total_y[[4]][[3]],
#             plots_NMDS_total_y[[4]][[4]],plots_NMDS_total_y[[4]][[5]],plots_NMDS_total_y_fluvial[[4]][[4]],
#             
#             plots_NMDS_total_y[[5]][[1]],plots_NMDS_total_y[[5]][[2]],plots_NMDS_total_y[[5]][[3]],
#             plots_NMDS_total_y[[5]][[4]],plots_NMDS_total_y[[5]][[5]],plots_NMDS_total_y_fluvial[[5]][[5]],
#             
#             ncol=6,nrow=5, top="NMDS groups")
#dev.off()

# NMDS significant plots______________________________________ ####
biod_names <- c("S16","S18","phy","zoo", "zoo.18S")
net_names <- c("600 km", "300 km","100 km","60 km","6 km", "Fluvial")
color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")

# For Euclidean network
plots_NMDS_x_sign <- list()
plots_NMDS_total_x_sign <- list()

plots_NMDS_y_sign <- list()
plots_NMDS_total_y_sign <- list()

ref_value_X <- 0
ref_value_Y <- 0

for (groups in 1:5) {
  for (net in 1:5) {
    coin <- PCA_network_results[[net]][c(length(PCA_network_results[[net]])-54):length(PCA_network_results[[net]])]
    centr_iso <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]])
    
    spe.abu.jac<-vegdist(comm_data[[groups]],method = "jaccard") #Càlcul de la matriu de similitud amb Bray-Curtis
    spe.abu.MDS<-metaMDS(spe.abu.jac, distance="jaccard", k=2, try = 999) #MDS amb similitud bray-curtis
    
    x <- spe.abu.MDS$points[,1]
    y <- spe.abu.MDS$points[,2]
    
    dataset <- data.frame(x,y,centr_iso)
    
    sign_x <- summary.gam(gam(dataset[,1]~ s(dataset[,3], k=2, bs="cr"), method = "REML"))[[8]]
    pred_x <- predict(gam(dataset[,1]~ s(dataset[,3], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
    my_data <- data.frame(cbind(dataset[,1],dataset[,3]),
                          mu   = pred_x$fit,
                          low  = pred_x$fit - 1.96 * pred_x$se.fit,
                          high = pred_x$fit + 1.96 * pred_x$se.fit)
    if(sign_x<0.05){
      
      ref_value_X <- ref_value_X+1
      
      plots_NMDS_x_sign[[ref_value_X]]<-ggplot(my_data, aes(x = X2, y = X1)) +
        geom_jitter(alpha=0.8,shape=21, size=5, aes(fill=X2))+
        geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity",  colour="black",size=2,se = TRUE)+
        scale_fill_continuous(type = "viridis")+
        labs(title = biod_names[groups], 
             subtitle = paste(net_names[net]," ", "Stress=",round(spe.abu.MDS$stress,2)))+
        xlab("Centrality-Isolation")+ylab("NMDS1")+
        theme_classic()+
        theme(legend.position = "none",
              panel.background=element_rect(colour="black", fill=alpha(color_groups[groups],0.1))) 
    }
    
    sign_y <- summary.gam(gam(dataset[,2]~ s(dataset[,3], k=2, bs="cr"), method = "REML"))[[8]]
    pred_y <- predict(gam(dataset[,2]~ s(dataset[,3], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
    my_data <- data.frame(cbind(dataset[,2],dataset[,3]),
                          mu   = pred_y$fit,
                          low  = pred_y$fit - 1.96 * pred_y$se.fit,
                          high = pred_y$fit + 1.96 * pred_y$se.fit)
    if(sign_y<0.05){
      
      ref_value_Y <- ref_value_Y+1
      
      plots_NMDS_y_sign[[ref_value_Y]] <- ggplot(my_data, aes(x = X2, y = X1)) +
        geom_jitter(alpha=0.8,shape=21, size=5, aes(fill=X2))+
        geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="black",size=2,se = T)+
        scale_fill_continuous(type = "viridis")+
        labs(title = biod_names[groups], 
             subtitle = paste(net_names[net]," ", "Stress=",round(spe.abu.MDS$stress,2)))+
        xlab("Centrality-Isolation")+ylab("NMDS2")+
        theme_classic()+
        theme(legend.position = "none",
              panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1))) 
    }
  }
  plots_NMDS_total_x_sign[[groups]] <- plots_NMDS_x_sign
  plots_NMDS_total_y_sign[[groups]] <- plots_NMDS_y_sign
}

# For Fluvial network
plots_NMDS_x_fluvial_sign <- list()
plots_NMDS_total_x_fluvial_sign <- list()

plots_NMDS_y_fluvial_sign <- list()
plots_NMDS_total_y_fluvial_sign <- list()

ref_value_X <- 0
ref_value_Y <- 0

for (groups in 1:5) {
  coin <- PCA_fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
  centr_iso <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]])
  
  spe.abu.jac<-vegdist(comm_data[[groups]],method = "jaccard") #Càlcul de la matriu de similitud amb Bray-Curtis
  spe.abu.MDS<-metaMDS(spe.abu.jac, distance="jaccard", k=2, try = 999) #MDS amb similitud bray-curtis
  
  x <- spe.abu.MDS$points[,1]
  y <- spe.abu.MDS$points[,2]
  
  dataset <- data.frame(x,y,centr_iso)
  
  sign_x <- summary.gam(gam(dataset[,1]~ s(dataset[,3], k=2, bs="cr"), method = "REML"))[[8]]
  pred_x <- predict(gam(dataset[,1]~ s(dataset[,3], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
  my_data <- data.frame(cbind(dataset[,1],dataset[,3]),
                        mu   = pred_x$fit,
                        low  = pred_x$fit - 1.96 * pred_x$se.fit,
                        high = pred_x$fit + 1.96 * pred_x$se.fit)
  if(sign_x<0.05){
    
    ref_value_X <- ref_value_X+1
    
    plots_NMDS_x_fluvial_sign[[ref_value_X]] <-ggplot(my_data, aes(x = X2, y = X1)) +
      geom_jitter(alpha=0.8,shape=21, size=5, aes(fill=X2))+
      geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity",  colour="black",size=2,se = TRUE)+
      scale_fill_continuous(type = "viridis")+
      labs(title = biod_names[groups], 
           subtitle = paste(net_names[6]," ", "Stress=",round(spe.abu.MDS$stress,2)))+
      xlab("Centrality-Isolation")+ylab("NMDS1")+
      theme_classic()+
      theme(legend.position = "none",
            panel.background=element_rect(colour="black", fill=alpha(color_groups[groups],0.1))) 
  }
  
  sign_y <- summary.gam(gam(dataset[,2]~ s(dataset[,3], k=2, bs="cr"), method = "REML"))[[8]]
  pred_y <- predict(gam(dataset[,2]~ s(dataset[,3], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
  my_data <- data.frame(cbind(dataset[,2],dataset[,3]),
                        mu   = pred_y$fit,
                        low  = pred_y$fit - 1.96 * pred_y$se.fit,
                        high = pred_y$fit + 1.96 * pred_y$se.fit)
  if(sign_y<0.05){
    
    ref_value_Y <- ref_value_Y+1
    
    plots_NMDS_y_fluvial_sign[[ref_value_Y]] <- ggplot(my_data, aes(x = X2, y = X1)) +
      geom_jitter(alpha=0.8,shape=21, size=5, aes(fill=X2))+
      geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="black",size=2,se = T)+
      scale_fill_continuous(type = "viridis")+
      labs(title = biod_names[groups], 
           subtitle = paste(net_names[6]," ", "Stress=",round(spe.abu.MDS$stress,2)))+
      xlab("Centrality-Isolation")+ylab("NMDS2")+
      theme_classic()+
      theme(legend.position = "none",
            panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1))) 
  }
  plots_NMDS_total_x_fluvial_sign[[groups]] <- plots_NMDS_x_fluvial_sign 
  plots_NMDS_total_y_fluvial_sign[[groups]] <- plots_NMDS_y_fluvial_sign
}

# Plot extraction
plot_plot_sign_plot_NMDS_x <- list()
for(t in 1:length(plots_NMDS_total_x_sign[[5]])){
  plot_plot_sign_plot_NMDS_x[[t]] <- plots_NMDS_total_x_sign[[5]][[t]]
}
for(d in c(length(plots_NMDS_total_x_sign[[5]])+1):c(length(plots_NMDS_total_x_sign[[5]])+length(plots_NMDS_total_x_fluvial_sign[[5]]))){
  plot_plot_sign_plot_NMDS_x[[d]] <- plots_NMDS_total_x_fluvial_sign[[5]][[d-length(plots_NMDS_total_x_sign[[5]])]]
}

plot_plot_sign_plot_NMDS_y <- list()
for(t in 1:length(plots_NMDS_total_y_sign[[5]])){
  plot_plot_sign_plot_NMDS_y[[t]] <- plots_NMDS_total_y_sign[[5]][[t]]
}
for(d in c(length(plots_NMDS_total_y_sign[[5]])+1):c(length(plots_NMDS_total_y_sign[[5]])+length(plots_NMDS_total_y_fluvial_sign[[5]]))){
  plot_plot_sign_plot_NMDS_y[[d]] <- plots_NMDS_total_y_fluvial_sign[[5]][[d-length(plots_NMDS_total_y_sign[[5]])]]
}

# Print together
png(filename ="C:/Users/Cunilleramontcusi/GAM_Sign_NMDS.png",
    width =629*6 ,height =629*6 ,units ="px",res = 300)
grid.arrange(plot_plot_sign_plot_NMDS_x[[1]],plot_plot_sign_plot_NMDS_x[[2]],plot_plot_sign_plot_NMDS_x[[3]],
             plot_plot_sign_plot_NMDS_x[[4]],plot_plot_sign_plot_NMDS_x[[5]],plot_plot_sign_plot_NMDS_x[[6]],
             plot_plot_sign_plot_NMDS_x[[7]],plot_plot_sign_plot_NMDS_x[[8]],
             plot_plot_sign_plot_NMDS_y[[1]],plot_plot_sign_plot_NMDS_y[[2]],plot_plot_sign_plot_NMDS_y[[3]],
             plot_plot_sign_plot_NMDS_y[[4]],plot_plot_sign_plot_NMDS_y[[5]],
             ncol=3,nrow=5, top="Sign_Ttla")
dev.off()


centr_iso_fact <- cut(dataset$centr_iso, breaks = 4,labels=c("A","B","C","D"))
dataset[,4] <- centr_iso_fact

# Ordiplots in ggplot from: https://rdrr.io/github/jfq3/ggordiplots/src/R/gg_ordiplot.R
#' Make Ordination Axis Labels
#' Makes ordination axis labels that include, if apprpriate, the \% total variance explained by each axis.
ord_labels <-function(ord){
  ev <- vegan::eigenvals(ord)
  if (!is.na(ev)[1]) {
    tol <- -(1e-07)*ev[1]
    ord.labels <- rep("", length(ev))
    if ((any(is.na(ev))) | (any(ev < tol))) {
      for ( i in 1:length(ev)) {
        ord.labels[i] <- paste("DIM", i, sep = "")
      }
    }
    else {
      ev.pc <- round(100*(ev/sum(ev)), 2)
      axis.names <- names(ev)
      if (is.null(axis.names)) {
        for ( i in 1:length(ev.pc)) {
          ord.labels[i] <- paste("DIM", i, " ", sprintf(ev.pc[i], fmt = '%#.1f'), "%", sep="")
        }
      } else {
        for (i in 1:length(ev.pc)){
          ord.labels[i] <- paste(axis.names[i], " ", ev.pc[i],"%", sep="")
        }
      }
    }
  } else {
    ord.labels <- colnames(vegan::scores(ord))
  }
  
  return(ord.labels)
}

# Get site coordinates to plot.
df_ord <- vegan::scores(spe.abu.MDS, display = "sites", scaling=1, choices=c(1,2))
axis.labels <- ord_labels(spe.abu.MDS)[c(1,2)]
df_ord <- data.frame(x=df_ord[,1], y=df_ord[,2], Group= dataset[,4])
# Get parameters from the ordiellipse function.
rslt <- vegan::ordiellipse(spe.abu.MDS, groups=dataset[,4], display = "sites", 
                           scaling=1, choices=c(1,2), kind = "sd", show.groups = as.vector(levels(dataset[,4])), draw = "none", label = F)
# Get points to plot for the ellipses.
df_ellipse <- data.frame()
for(g in as.vector(levels(dataset[,4]))) {
  df_ellipse <- rbind(df_ellipse, 
                      cbind(as.data.frame(with(df_ord[df_ord$Group==g,],
                                               vegan:::veganCovEllipse(rslt[[g]]$cov,rslt[[g]]$center, rslt[[g]]$scale))),Group=g))
}





# Linear models __________________________________________________________________ ####
    
    model_resutls_total <- list()
    linear_model_resutls <- list()
    linear_direct_model <- list()
    linear_direct_model_total <- list()
    
    for (groups in 1:5) {
      color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")
      for (net in 1:5) {
        coin <- PCA_network_results[[net]][c(length(PCA_network_results[[net]])-54):length(PCA_network_results[[net]])]
        dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                         biod[[groups]][,1:5])
        colnames(dataset)[1] <-c("Network")
        plots_grups <- list()
        output_results <- list()
        output_model_results <- list()
        
        require(betareg)
        p.val <- c()
        #CCA
        p.val[1] <-summary(betareg(abs(dataset[,2])~dataset[,1]))[[1]]$mean[2,4]
        output_results[[1]] <- summary(betareg(abs(dataset[,2])~dataset[,1]))
        output_model_results[[1]] <- betareg(abs(dataset[,2])~dataset[,1])
        #Richness
        p.val[2] <- summary(lm(log(dataset[,3]+1)~dataset[,1]))[[4]][2,4]
        output_results[[2]] <- summary(lm(log(dataset[,3]+1)~dataset[,1]))
        output_model_results[[2]] <- lm(log(dataset[,3]+1)~dataset[,1])
        #LCBD
        p.val[3] <-summary(betareg(dataset[,4]~dataset[,1]))[[1]]$mean[2,4]
        output_results[[3]] <- summary(betareg(dataset[,4]~dataset[,1]))
        output_model_results[[3]] <- betareg(dataset[,4]~dataset[,1])
        #Turn
        p.val[4] <-summary(betareg(dataset[,5]~dataset[,1]))[[1]]$mean[2,4]
        output_results[[4]] <- summary(betareg(dataset[,5]~dataset[,1]))
        output_model_results[[4]] <- betareg(dataset[,5]~dataset[,1])
        #RichDiff
        p.val[5] <-summary(betareg(dataset[,6]~dataset[,1]))[[1]]$mean[2,4]
        output_results[[5]] <- summary(betareg(dataset[,6]~dataset[,1]))
        output_model_results[[5]] <- betareg(dataset[,6]~dataset[,1])
        
        for(var in 1:5){
          if(p.val[var]>0.05){
            plots_grups[[var]] <-ggplot(dataset, aes_string(x=as.data.frame(dataset)[,1],
                                                            y=as.data.frame(dataset)[,var+1]))+
              geom_jitter(alpha=0.2, shape=21, size=3, fill=color_groups[groups], colour="black")+
              geom_smooth(method = "lm", se=F, colour="black",linetype=2, size=2)+
              labs(title=colnames(dataset)[var+1])+ylab(colnames(dataset)[var+1])+
              theme_classic()
          }else{
            plots_grups[[var]] <-ggplot(dataset, aes_string(x=as.data.frame(dataset)[,1],
                                                            y=as.data.frame(dataset)[,var+1]))+
              geom_jitter(alpha=0.9, shape=21, size=3, fill=color_groups[groups], colour="black")+
              geom_smooth(method = "lm", se=F, colour="black",linetype=1, size=2)+
              labs(title=colnames(dataset)[var+1])+ylab(colnames(dataset)[var+1])+
              theme_classic()
          }
        }
        
        linear_model_resutls[[net]] <-output_results
        linear_direct_model[[net]] <- output_model_results
        
        png(filename =paste("C:/Users/Cunilleramontcusi/","Divers",biod_names[[groups]],"_",names(PCA_network_results)[[net]],".png"),
            width =582*2 ,height =629*2 ,units ="px",res = 200)
        grid.arrange(plots_grups[[1]],plots_grups[[2]],
                     plots_grups[[3]],plots_grups[[4]],
                     plots_grups[[5]],
                     ncol=2,nrow=3, top=names(PCA_network_results)[[net]])
        dev.off()
      }
      model_resutls_total[[groups]] <- linear_model_resutls
      linear_direct_model_total[[groups]] <- linear_direct_model
    }   
    
    #_____________________________________________________________________________________________________________________________________________#
    model_resutls_fluvial <- list()
    direct_model_resutls_fluvial <- list()
    
    for (groups in 1:5) {
      color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")
      coin <- PCA_fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
      dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                       biod[[groups]][,1:5])
      colnames(dataset)[1] <-c("Network")
      plots_grups <- list()
      output_results <- list()
      output_model_results <- list()
      
      require(betareg)
      p.val <- c()
      #CCA
      p.val[1] <-summary(betareg(abs(dataset[,2])~dataset[,1]))[[1]]$mean[2,4]
      output_results[[1]] <- summary(betareg(abs(dataset[,2])~dataset[,1]))
      output_model_results[[1]] <- betareg(abs(dataset[,2])~dataset[,1])
      #Richness
      p.val[2] <- summary(lm(log(dataset[,3]+1)~dataset[,1]))[[4]][2,4]
      output_results[[2]] <- summary(lm(log(dataset[,3]+1)~dataset[,1]))
      output_model_results[[2]] <- lm(log(dataset[,3]+1)~dataset[,1])
      #LCBD
      p.val[3] <-summary(betareg(dataset[,4]~dataset[,1]))[[1]]$mean[2,4]
      output_results[[3]] <- summary(betareg(dataset[,4]~dataset[,1]))
      output_model_results[[3]] <- betareg(dataset[,4]~dataset[,1])
      #Turn
      p.val[4] <-summary(betareg(dataset[,5]~dataset[,1]))[[1]]$mean[2,4]
      output_results[[4]] <- summary(betareg(dataset[,5]~dataset[,1]))
      output_model_results[[4]] <- betareg(dataset[,5]~dataset[,1])
      #RichDiff
      p.val[5] <-summary(betareg(dataset[,6]~dataset[,1]))[[1]]$mean[2,4]
      output_results[[5]] <- summary(betareg(dataset[,6]~dataset[,1]))
      output_model_results[[5]] <- betareg(dataset[,6]~dataset[,1])
      
      for(var in 1:5){
        if(p.val[var]>0.05){
          plots_grups[[var]] <-ggplot(dataset, aes_string(x=as.data.frame(dataset)[,1],
                                                          y=as.data.frame(dataset)[,var+1]))+
            geom_jitter(alpha=0.2, shape=21, size=3, fill=color_groups[groups], colour="black")+
            geom_smooth(method = "lm", se=F, colour="black",linetype=2, size=2)+
            labs(title=colnames(dataset)[var+1])+ylab(colnames(dataset)[var+1])+
            theme_classic()
        }else{
          plots_grups[[var]] <-ggplot(dataset, aes_string(x=as.data.frame(dataset)[,1],
                                                          y=as.data.frame(dataset)[,var+1]))+
            geom_jitter(alpha=0.9, shape=21, size=3, fill=color_groups[groups], colour="black")+
            geom_smooth(method = "lm", se=F, colour="black",linetype=1, size=2)+
            labs(title=colnames(dataset)[var+1])+ylab(colnames(dataset)[var+1])+
            theme_classic()
        }
      }
      
      model_resutls_fluvial[[groups]] <-output_results
      direct_model_resutls_fluvial[[groups]] <- 
        
        png(filename =paste("C:/Users/Cunilleramontcusi/","Divers",biod_names[[groups]],"_Fluvial",".png"),
            width =582*2 ,height =629*2 ,units ="px",res = 200)
      grid.arrange(plots_grups[[1]],plots_grups[[2]],
                   plots_grups[[3]],plots_grups[[4]],
                   plots_grups[[5]],
                   ncol=2,nrow=3, top="Fluvial network")
      dev.off()
    }
    
    
    sign_netw <- list()
    sign_groups <- list()
    for (group in 1:5) {
      for (netw in 1:5) {
        sign <- c()
        sign[1] <-model_resutls_total[[group]][[netw]][[1]][[1]]$mean[2,4]
        sign[2] <-model_resutls_total[[group]][[netw]][[2]][[4]][2,4]
        sign[3] <-model_resutls_total[[group]][[netw]][[3]][[1]]$mean[2,4]
        sign[4] <-model_resutls_total[[group]][[netw]][[4]][[1]]$mean[2,4]
        sign[5] <-model_resutls_total[[group]][[netw]][[5]][[1]]$mean[2,4]  
        sign_netw[[netw]] <- sign
      }
      sign_groups[[group]] <- sign_netw
    }
    
    flu_sign_groups <- list()
    for (group in 1:5) {
      sign <- c()
      sign[1] <-model_resutls_fluvial[[group]][[1]][[1]]$mean[2,4]
      sign[2] <-model_resutls_fluvial[[group]][[2]][[4]][2,4]
      sign[3] <-model_resutls_fluvial[[group]][[3]][[1]]$mean[2,4]
      sign[4] <-model_resutls_fluvial[[group]][[4]][[1]]$mean[2,4]
      sign[5] <-model_resutls_fluvial[[group]][[5]][[1]]$mean[2,4]  
      flu_sign_groups[[group]] <- sign
    }
    
    Names_Variab <- c("Environmental tracking", "Species richness", "LCBD", "Replacement", "Richness difference")
    plots_significance <- list()
    for (variable in 1:5) {
      
      max_netw <- cbind(c(sign_groups[[1]][[1]][[variable]], sign_groups[[2]][[1]][[variable]], sign_groups[[3]][[1]][[variable]], 
                          sign_groups[[4]][[1]][[variable]], sign_groups[[5]][[1]][[variable]]),
                        rep("600 km", 5), 
                        c("S16","S18","Phy","Zoo", "S18zoo"))
      
      mid_netw <- cbind(c(sign_groups[[1]][[2]][[variable]], sign_groups[[2]][[2]][[variable]], sign_groups[[3]][[2]][[variable]], 
                          sign_groups[[4]][[2]][[variable]], sign_groups[[5]][[2]][[variable]]),
                        rep("300 km", 5), 
                        c("S16","S18","Phy","Zoo", "S18zoo"))
      
      mid_mid_netw <- cbind(c(sign_groups[[1]][[3]][[variable]], sign_groups[[2]][[3]][[variable]], sign_groups[[3]][[3]][[variable]], 
                              sign_groups[[4]][[3]][[variable]], sign_groups[[5]][[3]][[variable]]),
                            rep("100 km", 5), 
                            c("S16","S18","Phy","Zoo", "S18zoo"))
      
      small_netw <- cbind(c(sign_groups[[1]][[4]][[variable]], sign_groups[[2]][[4]][[variable]], sign_groups[[3]][[4]][[variable]], 
                            sign_groups[[4]][[4]][[variable]], sign_groups[[5]][[4]][[variable]]),
                          rep("60 km", 5), 
                          c("S16","S18","Phy","Zoo", "S18zoo"))
      
      min_netw <- cbind(c(sign_groups[[1]][[5]][[variable]], sign_groups[[2]][[5]][[variable]], sign_groups[[3]][[5]][[variable]], 
                          sign_groups[[4]][[5]][[variable]], sign_groups[[5]][[5]][[variable]]),
                        rep("6 km", 5), 
                        c("S16","S18","Phy","Zoo", "S18zoo"))
      
      fluv_netw <- cbind(c(flu_sign_groups[[1]][[variable]], flu_sign_groups[[2]][[variable]], flu_sign_groups[[3]][[variable]], 
                           flu_sign_groups[[4]][[variable]], flu_sign_groups[[5]][[variable]]),
                         rep("Fluvial", 5), 
                         c("S16","S18","Phy","Zoo", "S18zoo"))
      
      
      dataset_pval <- as.data.frame(rbind(max_netw,mid_netw,mid_mid_netw,small_netw,min_netw,fluv_netw))
      colnames(dataset_pval) <- c("val","Network","Group")
      dataset_pval$val <-as.numeric(dataset_pval$val)
      dataset_pval$Network <- factor(dataset_pval$Network,
                                     levels = c("600 km", "300 km","100 km","60 km","6 km", "Fluvial"))
      dataset_pval$Group <- factor(dataset_pval$Group,
                                   levels = c("S16","S18","Phy","Zoo", "S18zoo"))
      significants <- rep("NoSign",5*6)
      significants[which(dataset_pval$val<0.05)]<- "Sign"
      dataset_pval$Sign <- factor(significants) 
      
      color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")
      
      plots_significance[[variable]] <-  ggplot(dataset_pval, aes(x=Network, y=as.numeric(val)))+
        geom_abline(slope = 0,intercept = 0.05, colour="black", linetype=2,size=1)+
        geom_jitter(aes(fill=Group, alpha=Sign, size=Sign), shape=21,width = 0.5)+
        scale_size_manual(values = c(2,6))+
        scale_alpha_manual(values = c(0.2,0.9))+
        scale_fill_manual(values=c(color_groups[1],color_groups[2],color_groups[3],
                                   color_groups[4],color_groups[5],color_groups[6]))+
        scale_y_continuous(expand = c(0.2,0.01),
                           breaks =c(0.2,0.4,0.6,0.8,1) )+
        
        geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5), size=1, colour="grey70")+
        labs(title=Names_Variab[variable])+ylab("P-values")+xlab("Network")+
        theme_classic()
    }
    
    png(filename =paste("C:/Users/Cunilleramontcusi/Significance_Groups.png"),
        width =900*5 ,height =700*5 ,units ="px",res = 400)
    grid.arrange(plots_significance[[1]],plots_significance[[2]],
                 plots_significance[[3]],plots_significance[[4]],
                 plots_significance[[5]],
                 ncol=2,nrow=3, top="Significance values")
    dev.off()
    
    #_____________________________________________________________________________________________________________________________________________#


# Linear sign models ____________________####
    Sign_plots_total <- list()
    ref_value <- 0
    Names_Networks <- c("600 km", "300 km","100 km","60 km","6 km")
    for (groups in 1:5) {
      color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")
      for (net in 1:5) {
        coin <- PCA_network_results[[net]][c(length(PCA_network_results[[net]])-54):length(PCA_network_results[[net]])]
        dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                         biod[[groups]][,1:5])
        colnames(dataset)[1] <-c("Network")
        plots_grups <- list()
        output_results <- list()
        
        require(betareg)
        p.val <- c()
        #CCA
        p.val[1] <-summary(betareg(abs(dataset[,2])~dataset[,1]))[[1]]$mean[2,4]
        #Richness
        p.val[2] <- summary(lm(log(dataset[,3]+1)~dataset[,1]))[[4]][2,4]
        #LCBD
        p.val[3] <-summary(betareg(dataset[,4]~dataset[,1]))[[1]]$mean[2,4]
        #Turn
        p.val[4] <-summary(betareg(dataset[,5]~dataset[,1]))[[1]]$mean[2,4]
        #RichDiff
        p.val[5] <-summary(betareg(dataset[,6]~dataset[,1]))[[1]]$mean[2,4]
        
        select_p.val <- which(p.val<0.05)
        
        if(length(select_p.val)>0){
          for(var in 1:length(select_p.val)){
            ref_value <- ref_value+1
            Sign_plots_total[[ref_value]] <- ggplot(dataset, aes_string(x=as.data.frame(dataset)[,1],
                                                                        y=as.data.frame(dataset)[,select_p.val[var]+1]))+
              geom_jitter(alpha=0.9, shape=21, size=3, fill=color_groups[groups], colour="black")+
              geom_smooth(method = "lm", se=F, colour="black",linetype=1, size=2)+
              labs(title=paste(Names_Networks[net],colnames(dataset)[select_p.val[var]+1]))+
              ylab(colnames(dataset)[select_p.val[var]+1])+
              theme_classic()
          }
        }
      }
    }   
    
    png(filename ="C:/Users/Cunilleramontcusi/Sign_Lineal_Diverse.png",
        width =582*4 ,height =629*4 ,units ="px",res = 300)
    grid.arrange(Sign_plots_total[[1]],Sign_plots_total[[2]],
                 Sign_plots_total[[3]],Sign_plots_total[[4]],
                 Sign_plots_total[[5]],Sign_plots_total[[6]],
                 Sign_plots_total[[7]],Sign_plots_total[[8]],
                 Sign_plots_total[[9]],Sign_plots_total[[10]],
                 Sign_plots_total[[11]],Sign_plots_total[[12]],
                 ncol=3,nrow=4, top="Sign_Ttla")
    dev.off()
    
    
## MODEL SELECTION AND PLOTTING OF SIGNIFICANT BEST MODELS___####
#_____________________________________________________________________________________________________________________________________________#
Result_output <- list()
Result_output_complete <- list()
Result_output_total <- list()

Model_result_output <- list()
Model_result_output_complete <- list()
Model_result_output_total <- list()

Result_pval <- list()
Result_pval_complete <- list()
Result_pval_total <- list()

for (groups in 1:5) {
  color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")
  for (net in 1:5) {
    coin <- PCA_network_results[[net]][c(length(PCA_network_results[[net]])-54):length(PCA_network_results[[net]])]
    dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                     biod[[groups]][,1:5])
    colnames(dataset)[1] <-c("Network")
    plots_grups <- list()
    
    # Linear model
    LINEAR_output_results <- list()
    LINEAR_output_model_results <- list()
    
    require(betareg)
    LINEAR_p.val <- c()
    #CCA
    LINEAR_p.val[1] <-summary(betareg(abs(dataset[,2])~dataset[,1]))[[1]]$mean[2,4]
    LINEAR_output_results[[1]] <- summary(betareg(abs(dataset[,2])~dataset[,1]))
    LINEAR_output_model_results[[1]] <- betareg(abs(dataset[,2])~dataset[,1])
    #Richness
    LINEAR_p.val[2] <- summary(lm(log(dataset[,3]+1)~dataset[,1]))[[4]][2,4]
    LINEAR_output_results[[2]] <- summary(lm(log(dataset[,3]+1)~dataset[,1]))
    LINEAR_output_model_results[[2]] <- lm(log(dataset[,3]+1)~dataset[,1])
    #LCBD
    LINEAR_p.val[3] <-summary(betareg(dataset[,4]~dataset[,1]))[[1]]$mean[2,4]
    LINEAR_output_results[[3]] <- summary(betareg(dataset[,4]~dataset[,1]))
    LINEAR_output_model_results[[3]] <- betareg(dataset[,4]~dataset[,1])
    #Turn
    LINEAR_p.val[4] <-summary(betareg(dataset[,5]~dataset[,1]))[[1]]$mean[2,4]
    LINEAR_output_results[[4]] <- summary(betareg(dataset[,5]~dataset[,1]))
    LINEAR_output_model_results[[4]] <- betareg(dataset[,5]~dataset[,1])
    #RichDiff
    LINEAR_p.val[5] <-summary(betareg(dataset[,6]~dataset[,1]))[[1]]$mean[2,4]
    LINEAR_output_results[[5]] <- summary(betareg(dataset[,6]~dataset[,1]))
    LINEAR_output_model_results[[5]] <- betareg(dataset[,6]~dataset[,1])
    
    # GAM model
    GAM_output_results <- list()
    GAM_output_model_results <- list()
    
    GAM_p.val <- c()
    #CCA
    GAM_p.val[1] <-summary.gam(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
    GAM_output_results[[1]] <- summary.gam(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
    GAM_output_model_results[[1]] <- gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
    GAM_preds_1 <- predict(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
    #Richness
    GAM_p.val[2] <- summary.gam(gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
    GAM_output_results[[2]] <- summary.gam(gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
    GAM_output_model_results[[2]] <- gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
    GAM_preds_2<- predict(gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
    #LCBD
    GAM_p.val[3] <-summary.gam(gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
    GAM_output_results[[3]] <- summary.gam(gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
    GAM_output_model_results[[3]] <- gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
    GAM_preds_3<- predict(gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
    #Turn
    GAM_p.val[4] <-summary.gam(gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
    GAM_output_results[[4]] <- summary.gam(gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
    GAM_output_model_results[[4]] <- gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
    GAM_preds_4<- predict(gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
    #RichDiff
    GAM_p.val[5] <-summary.gam(gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
    GAM_output_results[[5]] <- summary.gam(gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
    GAM_output_model_results[[5]] <- gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
    GAM_preds_5<- predict(gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
    
    GAM.pred <- list(GAM_preds_1, GAM_preds_2, GAM_preds_3, GAM_preds_4, GAM_preds_5)
    
    for (best_bar in 1:5) {
      gam_mode <- LINEAR_output_model_results[[best_bar]]
      lin_mod<- GAM_output_model_results [[best_bar]]
      
      gamModels <- list(lin_mod, gam_mode)
      a <- model.sel(gamModels)[1]
      
      # Plots for GAM relationships  
      if(a$class=="gam"){
        
        Result_output[[best_bar]] <- GAM_output_results[[best_bar]]
        Model_result_output[[best_bar]] <- GAM_output_model_results[[best_bar]]
        Result_pval[[best_bar]] <-GAM_p.val[best_bar] 
        
        if(GAM_p.val[best_bar]>0.05){
          my_data <- data.frame(cbind(dataset[,best_bar+1],dataset[,1]),
                                mu   = GAM.pred[[best_bar]]$fit,
                                low  = GAM.pred[[best_bar]]$fit - 1.96 * GAM.pred[[best_bar]]$se.fit,
                                high = GAM.pred[[best_bar]]$fit + 1.96 * GAM.pred[[best_bar]]$se.fit)
          plots_grups[[best_bar]] <-
            ggplot(my_data, aes(x = X2, y = X1)) +
            geom_jitter(alpha=0.2, shape=21, size=3, fill=color_groups[groups], colour="black")+
            geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="#4d004d",linetype=2, size=2)+
            labs(title=colnames(dataset)[best_bar+1])+ylab(colnames(dataset)[best_bar+1])+
            theme_classic()
        }else{
          my_data <- data.frame(cbind(dataset[,best_bar+1],dataset[,1]),
                                mu   = GAM.pred[[best_bar]]$fit,
                                low  = GAM.pred[[best_bar]]$fit - 1.96 * GAM.pred[[best_bar]]$se.fit,
                                high = GAM.pred[[best_bar]]$fit + 1.96 * GAM.pred[[best_bar]]$se.fit)
          plots_grups[[best_bar]] <-ggplot(my_data, aes(x = X2, y = X1)) +
            geom_jitter(alpha=0.9, shape=21, size=3, fill=color_groups[groups], colour="black")+
            geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="#4d004d",linetype=1, size=2)+
            labs(title=colnames(dataset)[best_bar+1])+ylab(colnames(dataset)[best_bar+1])+
            theme_classic()
        }
        
        # Plots for linear relationships  
      }else{
        
        Result_output[[best_bar]] <- LINEAR_output_results[[best_bar]]
        Model_result_output[[best_bar]] <- LINEAR_output_model_results[[best_bar]]
        Result_pval[[best_bar]] <-LINEAR_p.val[best_bar] 
        
        if(LINEAR_p.val[best_bar]>0.05){
          plots_grups[[best_bar]] <-ggplot(dataset, aes_string(x=as.data.frame(dataset)[,1],
                                                               y=as.data.frame(dataset)[,best_bar+1]))+
            geom_jitter(alpha=0.2, shape=21, size=3, fill=color_groups[groups], colour="black")+
            geom_smooth(method = "lm", se=F, colour="black",linetype=2, size=2)+
            labs(title=colnames(dataset)[best_bar+1])+ylab(colnames(dataset)[best_bar+1])+
            theme_classic()
        }else{
          plots_grups[[best_bar]] <-ggplot(dataset, aes_string(x=as.data.frame(dataset)[,1],
                                                               y=as.data.frame(dataset)[,best_bar+1]))+
            geom_jitter(alpha=0.9, shape=21, size=3, fill=color_groups[groups], colour="black")+
            geom_smooth(method = "lm", se=F, colour="black",linetype=1, size=2)+
            labs(title=colnames(dataset)[best_bar+1])+ylab(colnames(dataset)[best_bar+1])+
            theme_classic()
          
        }
      }
    }
    
    Result_output_complete[[net]] <- Result_output
    Model_result_output_complete[[net]]<- Model_result_output
    Result_pval_complete[[net]] <- Result_pval
    
    png(filename =paste("C:/Users/Cunilleramontcusi/","All_Divers",biod_names[[groups]],"_",names(PCA_network_results)[[net]],".png"),
        width =582*2 ,height =629*2 ,units ="px",res = 200)
    grid.arrange(plots_grups[[1]],plots_grups[[2]],
                 plots_grups[[3]],plots_grups[[4]],
                 plots_grups[[5]],
                 ncol=2,nrow=3, top=names(PCA_network_results)[[net]])
    dev.off()
  }
  Result_output_total[[groups]] <- Result_output_complete
  Model_result_output_total[[groups]] <- Model_result_output_complete
  Result_pval_total[[groups]] <- Result_pval_complete
}             

# Fluvial
Result_output_fluvial <- list()
Result_output_total_fluvial <- list()

Model_result_output_fluvial <- list()
Model_result_output_total_fluvial <- list()

Result_pval_fluvial <- list()
Result_pval_total_fluvial <- list()

for (groups in 1:5) {
  color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")
  
  coin <- PCA_fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
  dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                   biod[[groups]][,1:5])
  colnames(dataset)[1] <-c("Network")
  plots_grups <- list()
  
  # Linear model
  LINEAR_output_results <- list()
  LINEAR_output_model_results <- list()
  
  require(betareg)
  LINEAR_p.val <- c()
  #CCA
  LINEAR_p.val[1] <-summary(betareg(abs(dataset[,2])~dataset[,1]))[[1]]$mean[2,4]
  LINEAR_output_results[[1]] <- summary(betareg(abs(dataset[,2])~dataset[,1]))
  LINEAR_output_model_results[[1]] <- betareg(abs(dataset[,2])~dataset[,1])
  #Richness
  LINEAR_p.val[2] <- summary(lm(log(dataset[,3]+1)~dataset[,1]))[[4]][2,4]
  LINEAR_output_results[[2]] <- summary(lm(log(dataset[,3]+1)~dataset[,1]))
  LINEAR_output_model_results[[2]] <- lm(log(dataset[,3]+1)~dataset[,1])
  #LCBD
  LINEAR_p.val[3] <-summary(betareg(dataset[,4]~dataset[,1]))[[1]]$mean[2,4]
  LINEAR_output_results[[3]] <- summary(betareg(dataset[,4]~dataset[,1]))
  LINEAR_output_model_results[[3]] <- betareg(dataset[,4]~dataset[,1])
  #Turn
  LINEAR_p.val[4] <-summary(betareg(dataset[,5]~dataset[,1]))[[1]]$mean[2,4]
  LINEAR_output_results[[4]] <- summary(betareg(dataset[,5]~dataset[,1]))
  LINEAR_output_model_results[[4]] <- betareg(dataset[,5]~dataset[,1])
  #RichDiff
  LINEAR_p.val[5] <-summary(betareg(dataset[,6]~dataset[,1]))[[1]]$mean[2,4]
  LINEAR_output_results[[5]] <- summary(betareg(dataset[,6]~dataset[,1]))
  LINEAR_output_model_results[[5]] <- betareg(dataset[,6]~dataset[,1])
  
  # GAM model
  GAM_output_results <- list()
  GAM_output_model_results <- list()
  
  GAM_p.val <- c()
  #CCA
  GAM_p.val[1] <-summary.gam(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
  GAM_output_results[[1]] <- summary.gam(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
  GAM_output_model_results[[1]] <- gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
  GAM_preds_1 <- predict(gam(dataset[,2]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
  #Richness
  GAM_p.val[2] <- summary.gam(gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
  GAM_output_results[[2]] <- summary.gam(gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
  GAM_output_model_results[[2]] <- gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
  GAM_preds_2<- predict(gam(dataset[,3]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
  #LCBD
  GAM_p.val[3] <-summary.gam(gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
  GAM_output_results[[3]] <- summary.gam(gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
  GAM_output_model_results[[3]] <- gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
  GAM_preds_3<- predict(gam(dataset[,4]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
  #Turn
  GAM_p.val[4] <-summary.gam(gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
  GAM_output_results[[4]] <- summary.gam(gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
  GAM_output_model_results[[4]] <- gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
  GAM_preds_4<- predict(gam(dataset[,5]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
  #RichDiff
  GAM_p.val[5] <-summary.gam(gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))[[8]]
  GAM_output_results[[5]] <- summary.gam(gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML"))
  GAM_output_model_results[[5]] <- gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML")
  GAM_preds_5<- predict(gam(dataset[,6]~ s(dataset[,1], k=2, bs="cr"), method = "REML"), se.fit = TRUE)
  
  GAM.pred <- list(GAM_preds_1, GAM_preds_2, GAM_preds_3, GAM_preds_4, GAM_preds_5)
  
  for (best_bar in 1:5) {
    gam_mode <- LINEAR_output_model_results[[best_bar]]
    lin_mod<- GAM_output_model_results [[best_bar]]
    
    gamModels <- list(lin_mod, gam_mode)
    a <- model.sel(gamModels)[1]
    
    # Plots for GAM relationships  
    if(a$class=="gam"){
      
      Result_output_fluvial[[best_bar]] <- GAM_output_results[[best_bar]]
      Model_result_output_fluvial[[best_bar]] <- GAM_output_model_results[[best_bar]]
      Result_pval_fluvial[[best_bar]] <- GAM_p.val[best_bar]
      
      if(GAM_p.val[best_bar]>0.05){
        my_data <- data.frame(cbind(dataset[,best_bar+1],dataset[,1]),
                              mu   = GAM.pred[[best_bar]]$fit,
                              low  = GAM.pred[[best_bar]]$fit - 1.96 * GAM.pred[[best_bar]]$se.fit,
                              high = GAM.pred[[best_bar]]$fit + 1.96 * GAM.pred[[best_bar]]$se.fit)
        plots_grups[[best_bar]] <-
          ggplot(my_data, aes(x = X2, y = X1)) +
          geom_jitter(alpha=0.2, shape=21, size=3, fill=color_groups[groups], colour="black")+
          geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="#4d004d",linetype=2, size=2)+
          labs(title=colnames(dataset)[best_bar+1])+ylab(colnames(dataset)[best_bar+1])+
          theme_classic()
      }else{
        my_data <- data.frame(cbind(dataset[,best_bar+1],dataset[,1]),
                              mu   = GAM.pred[[best_bar]]$fit,
                              low  = GAM.pred[[best_bar]]$fit - 1.96 * GAM.pred[[best_bar]]$se.fit,
                              high = GAM.pred[[best_bar]]$fit + 1.96 * GAM.pred[[best_bar]]$se.fit)
        plots_grups[[best_bar]] <-ggplot(my_data, aes(x = X2, y = X1)) +
          geom_jitter(alpha=0.9, shape=21, size=3, fill=color_groups[groups], colour="black")+
          geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="#4d004d",linetype=1, size=2)+
          labs(title=colnames(dataset)[best_bar+1])+ylab(colnames(dataset)[best_bar+1])+
          theme_classic()
      }
      
      # Plots for linear relationships  
    }else{
      
      Result_output_fluvial[[best_bar]] <- LINEAR_output_results[[best_bar]]
      Model_result_output_fluvial[[best_bar]] <- LINEAR_output_model_results[[best_bar]]
      Result_pval_fluvial[[best_bar]] <- LINEAR_p.val[best_bar]
      
      if(LINEAR_p.val[best_bar]>0.05){
        plots_grups[[best_bar]] <-ggplot(dataset, aes_string(x=as.data.frame(dataset)[,1],
                                                             y=as.data.frame(dataset)[,best_bar+1]))+
          geom_jitter(alpha=0.2, shape=21, size=3, fill=color_groups[groups], colour="black")+
          geom_smooth(method = "lm", se=F, colour="black",linetype=2, size=2)+
          labs(title=colnames(dataset)[best_bar+1])+ylab(colnames(dataset)[best_bar+1])+
          theme_classic()
      }else{
        plots_grups[[best_bar]] <-ggplot(dataset, aes_string(x=as.data.frame(dataset)[,1],
                                                             y=as.data.frame(dataset)[,best_bar+1]))+
          geom_jitter(alpha=0.9, shape=21, size=3, fill=color_groups[groups], colour="black")+
          geom_smooth(method = "lm", se=F, colour="black",linetype=1, size=2)+
          labs(title=colnames(dataset)[best_bar+1])+ylab(colnames(dataset)[best_bar+1])+
          theme_classic()
        
      }
    }
  }
  
  Result_output_total_fluvial[[groups]] <- Result_output_fluvial
  Model_result_output_total_fluvial[[groups]] <- Model_result_output_fluvial
  Result_pval_total_fluvial[[groups]] <- Result_pval_fluvial 
  
  png(filename =paste("C:/Users/Cunilleramontcusi/","All_Fluvial_Divers",biod_names[[groups]],"_",names(PCA_network_results)[[net]],".png"),
      width =582*2 ,height =629*2 ,units ="px",res = 200)
  grid.arrange(plots_grups[[1]],plots_grups[[2]],
               plots_grups[[3]],plots_grups[[4]],
               plots_grups[[5]],
               ncol=2,nrow=3, top=names(PCA_network_results)[[net]])
  dev.off()
}             



# PCA_approach __________________________________________________________________________________________####


PCA_biod_results <- list()
PCA_biod_plots <- list()
for (r in 1:length(biod)) {
  PCA_result <- prcomp(biod[[r]][,c(3,4,6,7)], center = T, scale. = T)
  PCA_biod_plots[[r]] <- PCA_result
  PCA_biod_results[[r]] <-PCA_result$x[,1:2]
}
names(PCA_biod_results) <- c("S16","S18","phy","zoo")

grid.arrange(autoplot(PCA_biod_plots[[1]],loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
                      loadings.colour = 'red', loadings.label.colour="black")+
               geom_point(size=2, alpha=0.4,color=CUNILLERA_cols("yellow"))+labs(title="Bacteria S16")+
               theme_bw(),
             autoplot(PCA_biod_plots[[2]],loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
                      loadings.colour = 'red', loadings.label.colour="black")+
               geom_point(size=2, alpha=0.4,color=CUNILLERA_cols("blue"))+labs(title="Protista S18")+
               theme_bw(),
             autoplot(PCA_biod_plots[[3]],loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
                      loadings.colour = 'red', loadings.label.colour="black")+
               geom_point(size=2, alpha=0.4,color=CUNILLERA_cols("green"))+labs(title="Phytoplankton")+
               theme_bw(),
             autoplot(PCA_biod_plots[[4]],loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
                      loadings.colour = 'red', loadings.label.colour="black")+
               geom_point(size=2, alpha=0.4,color=CUNILLERA_cols("red"))+labs(title="Zooplankton")+
               theme_bw(), 
             ncol = 2, nrow=2)
# Extracting the values of which lakes have been sampled for each taxonomic group
# S16 = 52 lakes
# S18 = 48 lakes
# Phy = 50 lakes
# Zoo = 52 lakes

coincidence <- c()
coincidence_values <- list()
for (r in 1:4) {
  for (coinc in 1:nrow(all_lakes_coord)) {
    coincidence_value <- 0
    coincidence_value<-which(rownames(all_lakes_coord)==rownames(env_data[[r]])[coinc])
    coincidence_value<-ifelse(length(coincidence_value)==0, 0 ,coincidence_value)
    coincidence[coinc] <- coincidence_value
  }  
  coincidence_values[[r]] <- coincidence
}


for (groups in 1:4) {
  plots_grups <- list()
  color_groups <- CUNILLERA_cols("yellow","blue","green","red")
  names_networks <- c("~600km","~300km","~100km","~60km","~6km")
  noms_groups <- c("Bacteria S16","Protista S18","Phytoplankton","Zooplankton")
  for (net in 1:5) {
    coin <- PCA_network_results[[net]][c(length(PCA_network_results[[net]])-54):length(PCA_network_results[[net]])]
    dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                     PCA_biod_results[[groups]][,1:2])
    colnames(dataset) <-c("Network","Div_PC1","Div_PC2")
    
    p.val_1 <- cor.test(dataset[,1],dataset[,2])[3]
    if (p.val_1>0.05) {
      plots_grups[[net]] <-ggplot(dataset, aes(x=Network, y=Div_PC1))+
        geom_jitter(alpha=0.2, shape=21, size=3, fill=color_groups[groups], colour="black")+
        geom_smooth(method = "lm", se=F, colour="black",linetype=2, size=2)+labs(title=names_networks[net])+
        theme_bw()
    }else{
      plots_grups[[net]] <-ggplot(dataset, aes(x=Network, y=Div_PC1))+
        geom_jitter(alpha=0.9, shape=21, size=3, fill=color_groups[groups], colour="black")+
        geom_smooth(method = "lm", se=F, colour="black",linetype=1, size=2)+labs(title=names_networks[net])+
        theme_bw()
    }
    
    p.val_2 <- cor.test(dataset[,1],dataset[,3])[3]
    if (p.val_2>0.05) {
      plots_grups[[net+5]] <-ggplot(dataset, aes(x=Network, y=Div_PC2))+
        geom_jitter(alpha=0.2, shape=21, size=3, fill=color_groups[groups], colour="black")+
        geom_smooth(method = "lm", se=F, colour="black",linetype=2, size=2)+labs(title=names_networks[net])+
        theme_bw()
    }else{
      plots_grups[[net+5]] <-ggplot(dataset, aes(x=Network, y=Div_PC2))+
        geom_jitter(alpha=0.9, shape=21, size=3, fill=color_groups[groups], colour="black")+
        geom_smooth(method = "lm", se=F, colour="black",linetype=1, size=2)+labs(title=names_networks[net])+
        theme_bw()
    }
    
  }
  png(filename =paste("C:/Users/Cunilleramontcusi/","PCA",names(PCA_biod_results)[[groups]],".png"),
      width =582*2 ,height =729*2 ,units ="px",res = 200)
  grid.arrange(plots_grups[[1]],plots_grups[[6]],
               plots_grups[[2]],plots_grups[[7]],
               plots_grups[[3]],plots_grups[[8]],
               plots_grups[[4]],plots_grups[[9]],
               plots_grups[[5]],plots_grups[[10]],
               ncol=2,nrow=5,top=noms_groups[groups])
  dev.off()
}             




















