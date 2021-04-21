##########################################
##########################################
##########                      ##########
##########  PACKAGE UPLOADING   ##########
##########                      ##########
##########################################
##########################################

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
library(gtable)    
library(grid)
library(gridExtra) 
library(png)

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

##########################################
##########################################
##########                      ##########
##########  DATA UPLOADING      ##########
##########                      ##########
##########################################
##########################################

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

##########################################
##########################################
##########                      ##########
##########  NETWORK BUILDING    ##########
##########                      ##########
##########################################
##########################################

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

####__________________________________________________________________________________________________________###
####__________________________________________________________________________________________________________###
# Building the whole network based on lks.global _____________________________________________________________####
####__________________________________________________________________________________________________________###

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
#library(raster)
#
#a <- c('C:/Users/Cunilleramontcusi/Desktop/495044851ffb5dc430127353e17ee99eec4e04e4/eu_dem_v11_E40N20/eu_dem_v11_E40N20.TIF.tif')
#e <- extent(min(dat1[,1]),max(dat1[,1]), min(dat1[,2]), max(dat1[,2]))
#template <- raster(e)
#proj4string(template) <- CRS(taiwan_elevation)
#writeRaster(template, file="MiamiWatershed.tif", format="GTiff")
#mosaic_rasters(gdalfile=a,dst_dataset="MiamiWatershed.tif",of="GTiff")
#gdalinfo("MiamiWatershed.tif")
#
#
#library(topoDistance)
#
#taiwan_elevation <- raster("C:/Users/Cunilleramontcusi/Desktop/495044851ffb5dc430127353e17ee99eec4e04e4/eu_dem_v11_E40N20/eu_dem_v11_E40N20.tif")
#e <- extent(min(dat1[,1]),max(dat1[,1]), min(dat1[,2]), max(dat1[,2]))
#plot(dat1[,1],dat1[,2])
#
#writeRaster(template, file="MiamiWatershed.tif", format="GTiff")
#
#
#plot(taiwan_elevation)
#plot(dat1)
#
#td <- topoDist(taiwan_elevation, dat1[,1:2], paths = TRUE)
#
#tlcp <- topoLCP(Yosemite$DEM, Yosemite$SDM, xy, 
#                paths = TRUE) 
#topoPathMap(Yosemite$DEM, xy, tlcp, costSurface = 
#              Yosemite$SDM, 
#            type = "hillshade", pathColor = "purple")
#

#_______________________________________________________________________________________________________________###
# Maxmimum distance between lakes -- 600 km approximately
max_distance <- Lakes_buffer(EU_lakes = dat, Samp_lakes = dat1, buffer_distance = 1, 
                             check_EU_Sampl_plots = F, check_BUFFER_plots = T,check_ELIMINATION_goodLAKES_plots = F )
# We extract the distance matrix and save it to run it externally (accelerate computing times)
max_distance_down <- max_distance[[2]]
#save(max_distance_down, file = "S16-values/AlPS_lakes_percol_MAX_DISTANCE.RData")

##_____
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
##____________________________________________________________ This has been uploaded to the supercomputer to calculate d.percol

#_______________________________________________________________________________________________________________###
# Middle distance between lakes -- 300 km approximately
mid_distance <- Lakes_buffer(EU_lakes = dat, Samp_lakes = dat1, buffer_distance = 2, 
                             check_EU_Sampl_plots = F, check_BUFFER_plots = T,check_ELIMINATION_goodLAKES_plots = F )
# We extract the distance matrix and save it to run it externally (accelerate computing times)
mid_distance_down <- mid_distance[[2]]
#save(mid_distance_down, file = "S16-values/AlPS_lakes_percol_MID_DISTANCE.RData")

##_____
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
##____________________________________________________________ This has been uploaded to the supercomputer to calculate d.percol

#_______________________________________________________________________________________________________________###
# Medium middle (MID_MID) distance between lakes -- 100 km approximately
mid_mid_distance <- Lakes_buffer(EU_lakes = dat, Samp_lakes = dat1, buffer_distance = 6, 
                                 check_EU_Sampl_plots = F, check_BUFFER_plots = T,check_ELIMINATION_goodLAKES_plots = F )
# We extract the distance matrix and save it to run it externally (accelerate computing times)
mid_mid_distance_down <- mid_mid_distance[[2]]
#save(mid_mid_distance_down, file = "S16-values/AlPS_lakes_percol_MID_MID_DISTANCE.RData")

##_____
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
##____________________________________________________________ This has been uploaded to the supercomputer to calculate d.percol

#_______________________________________________________________________________________________________________###
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
##____________________________________________________________

#_______________________________________________________________________________________________________________###
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
##____________________________________________________________ 


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

maps_Alps <- list()
detach("package:igraph", unload = TRUE)
for (e in 1:length(cordenades_xarxes)) {
  factors <- rep("No_Sampled",nrow(MAPS_xarxes[[e]]))
  factors[c(nrow(MAPS_xarxes[[e]])-54):nrow(MAPS_xarxes[[e]])] <- "Sampled"
  
  n<- network(MAPS_xarxes[[e]], directed=F, diag=F)
  n %v% "family" <- factors # Family is an standard name for the categortical variable that we are creating 
  maps_Alps[[e]] <- ggplot(n, layout=as.matrix(cordenades_xarxes[[e]][,1:2]),
                           aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( color = "grey45", alpha= 0.4) +
    geom_nodes(aes(fill=family, alpha=family), size=2 ,color="black", shape=21)+
    labs(x="",y="")+
    scale_fill_manual(values = c("grey35","red"))+
    scale_alpha_manual(values = c(0.5,1))+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none")
}

png(filename ="C:/Users/Cunilleramontcusi/Network_maps_Alps.png",width =745 ,height =742 ,units ="px",res = 100)
grid.arrange(maps_Alps[[1]],maps_Alps[[2]],maps_Alps[[3]],maps_Alps[[4]],maps_Alps[[5]])
dev.off()

############################################
############################################
##########                        ##########
########## RIVER NETWORK BUILDING ##########
##########                        ##########
############################################
############################################
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
sizes <- ifelse(cols=="red",10,1)

png(filename = "C:/Users/Cunilleramontcusi/Alpine_fluvial.png", width = 20000, height = 20000, res=1000)
par(mar=c(0,0,0,0))
plot(GRAPH_xarxes_fluvial[[1]], vertex.label = NA, vertex.size = sizes, vertex.size2 = sizes, vertex.color=cols,
     edge.width=0.01, edge.color="grey50")
dev.off()
detach("package:shp2graph", unload = TRUE)


##########################################
##########################################
##########                      ##########
##########  CENTRALITY VALUES   ##########
##########                      ##########
##########################################
##########################################

# Calculation of network data for the 55 lakes (closeness, degree, between, evcent, subgraph centrality)
network_data <- list()
for (i in 1:length(MAPS_xarxes)) {
  Xarx_grpah <- graph.adjacency(MAPS_xarxes[[i]], mode = "undirected",diag = F)
  output <- matrix(nrow = nrow(MAPS_xarxes[[i]]), ncol = 5)
  
  network_data[[i]] <- output
  
  network_data[[i]][,1] <- closeness(Xarx_grpah, mode = "all")
  network_data[[i]][,2] <- subgraph.centrality(Xarx_grpah)
  network_data[[i]][,3] <- degree(Xarx_grpah)
  network_data[[i]][,4] <- evcent(Xarx_grpah)[[1]]
  network_data[[i]][,5] <- betweenness(Xarx_grpah, directed = F)
  
  
  colnames(network_data[[i]]) <- c("clo_ALPS","Sub_centr_ALPS","deg_ALPS","Eigvec_cent_ALPS","bet_ALPS")
}


fluvial_network_data <- list()
output <- matrix(nrow = length(V(GRAPH_xarxes_fluvial[[1]])), ncol = 3)

fluvial_network_data[[1]] <- output

fluvial_network_data[[1]][,1] <- closeness(GRAPH_xarxes_fluvial[[1]], mode = "out")
fluvial_network_data[[1]][,2] <- degree(GRAPH_xarxes_fluvial[[1]],mode = "out")
fluvial_network_data[[1]][,3] <- betweenness(GRAPH_xarxes_fluvial[[1]])

colnames(fluvial_network_data[[1]]) <- c("clo_ALPS","deg_ALPS","bet_ALPS")

# Saving the list with the FIVE matrices with the network values 
save(network_data, file = "S16-values/network_dataset.RData")
save(fluvial_network_data, file = "S16-values/fluvial_network_dataset.RData")

library(ggfortify)
PCA_network_results <- list()
PCA_network_plot<- list()
for (r in 1:length(MAPS_xarxes)) {
  PCA_result <- prcomp(network_data[[r]], center = T, scale. = T)
  PCA_network_plot[[r]] <- PCA_result
  PCA_network_results[[r]] <-PCA_result$x[,1]
}
names(PCA_network_results) <- c("max_PCA_network","mid_PCA_network","mid_mid_PCA_network",
                                "small_PCA_network","min_PCA_network")

PCA_network_results[[1]] <- PCA_network_results[[1]]*-1

grid.arrange(autoplot(PCA_network_plot[[1]],loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
                      loadings.colour = 'red', loadings.label.colour="black")+
               geom_point(size=2, alpha=0.1,color=CUNILLERA_cols("black"))+labs(title="~600km")+
               theme_bw(), 
             autoplot(PCA_network_plot[[2]],loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
                      loadings.colour = 'red', loadings.label.colour="black")+
               geom_point(size=2, alpha=0.1,color=CUNILLERA_cols("black"))+labs(title="~300km")+
               theme_bw(),
             autoplot(PCA_network_plot[[3]],loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
                      loadings.colour = 'red', loadings.label.colour="black")+
               geom_point(size=2, alpha=0.1,color=CUNILLERA_cols("black"))+labs(title="~100km")+
               theme_bw(),
             autoplot(PCA_network_plot[[4]],loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
                      loadings.colour = 'red', loadings.label.colour="black")+
               geom_point(size=2, alpha=0.1,color=CUNILLERA_cols("black"))+labs(title="~60km")+
               theme_bw(),
             autoplot(PCA_network_plot[[5]],loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
                      loadings.colour = 'red', loadings.label.colour="black")+
               geom_point(size=2, alpha=0.1,color=CUNILLERA_cols("black"))+labs(title="~6km")+
               theme_bw(),  
             ncol = 2, nrow=3)


library(ggfortify)
PCA_fluvial_network_results <- list()
PCA_fluvial_network_plot<- list()

PCA_fluvial_result <- prcomp(fluvial_network_data[[1]], center = T, scale. = T)
PCA_fluvial_network_plot[[1]] <- PCA_fluvial_result
PCA_fluvial_network_results[[1]] <-PCA_fluvial_result$x[,1]

names(PCA_fluvial_network_results) <- c("PCA_fluvial_network")

autoplot(PCA_fluvial_network_plot[[1]],loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
         loadings.colour = 'red', loadings.label.colour="black")+
  geom_point(size=2, alpha=0.1,color=CUNILLERA_cols("black"))+labs(title="~600km")+
  theme_bw()


##########################################
##########################################
##########                      ##########
##########  DIVERSITY VALUES    ##########
##########                      ##########
##########################################
##########################################

##### Community indices list to gather everything 

community_indices<- list()

# Calculate the fitted values according to "The method" 
### FITTINGS OF OBSERVED VS FITTEDS WITH CCA AND DBRDA
# CCA
tst_coeficient_out <- list()
for (r in 1:5) {
  cca1<-cca(comm_data[[r]]~as.matrix(env_data[[r]][,4:6]))
  fitteds <- fitted(cca1, model="CCA", type= "response")
  coeffici <- c()
  for (e in 1:nrow(comm_data[[r]])) {
    coeffici[e] <- cor(fitteds[e,],as.numeric(comm_data[[r]][e,]),method = "pearson")
  }
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

##########################################
##########################################
##########                      ##########
##########  NET vs DIV VALUES   ##########
##########                      ##########
##########################################
##########################################

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

model_resutls_total <- list()
linear_model_resutls <- list()
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
     output_results[[1]] <- summary(betareg(abs(dataset[,2])~dataset[,1]))
    #Richness
     p.val[2] <- summary(lm(log(dataset[,3]+1)~dataset[,1]))[[4]][2,4]
     output_results[[2]] <- summary(lm(log(dataset[,3]+1)~dataset[,1]))  
    #LCBD
     p.val[3] <-summary(betareg(dataset[,4]~dataset[,1]))[[1]]$mean[2,4]
     output_results[[3]] <- summary(betareg(dataset[,4]~dataset[,1]))
    #Turn
     p.val[4] <-summary(betareg(dataset[,5]~dataset[,1]))[[1]]$mean[2,4]
     output_results[[4]] <- summary(betareg(dataset[,5]~dataset[,1]))
    #RichDiff
     p.val[5] <-summary(betareg(dataset[,6]~dataset[,1]))[[1]]$mean[2,4]
     output_results[[5]] <- summary(betareg(dataset[,6]~dataset[,1]))
      
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

    png(filename =paste("C:/Users/Cunilleramontcusi/","Divers",biod_names[[groups]],"_",names(PCA_network_results)[[net]],".png"),
        width =582*2 ,height =629*2 ,units ="px",res = 200)
    grid.arrange(plots_grups[[1]],plots_grups[[2]],
                 plots_grups[[3]],plots_grups[[4]],
                 plots_grups[[5]],
                 ncol=2,nrow=3, top=names(PCA_network_results)[[net]])
    dev.off()
  }
  model_resutls_total[[groups]] <- linear_model_resutls
}             

model_resutls_fluvial <- list()
for (groups in 1:5) {
  color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")
  coin <- PCA_fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
  dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                   biod[[groups]][,1:5])
  colnames(dataset)[1] <-c("Network")
  plots_grups <- list()
  output_results <- list()
  require(betareg)
  p.val <- c()
  #CCA
  p.val[1] <-summary(betareg(abs(dataset[,2])~dataset[,1]))[[1]]$mean[2,4]
  output_results[[1]] <- summary(betareg(abs(dataset[,2])~dataset[,1]))
  #Richness
  p.val[2] <- summary(lm(log(dataset[,3]+1)~dataset[,1]))[[4]][2,4]
  output_results[[2]] <- summary(lm(log(dataset[,3]+1)~dataset[,1]))  
  #LCBD
  p.val[3] <-summary(betareg(dataset[,4]~dataset[,1]))[[1]]$mean[2,4]
  output_results[[3]] <- summary(betareg(dataset[,4]~dataset[,1]))
  #Turn
  p.val[4] <-summary(betareg(dataset[,5]~dataset[,1]))[[1]]$mean[2,4]
  output_results[[4]] <- summary(betareg(dataset[,5]~dataset[,1]))
  #RichDiff
  p.val[5] <-summary(betareg(dataset[,6]~dataset[,1]))[[1]]$mean[2,4]
  output_results[[5]] <- summary(betareg(dataset[,6]~dataset[,1]))
  
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
  
  png(filename =paste("C:/Users/Cunilleramontcusi/","Divers",biod_names[[groups]],"_Fluvial",".png"),
      width =582*2 ,height =629*2 ,units ="px",res = 200)
  grid.arrange(plots_grups[[1]],plots_grups[[2]],
               plots_grups[[3]],plots_grups[[4]],
               plots_grups[[5]],
               ncol=2,nrow=3, top="Fluvial network")
  dev.off()
}




##########################################
##########################################
##########                      ##########
##########        END           ##########
##########                      ##########
##########################################
##########################################

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




















