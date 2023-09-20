
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
##########  NETWORK BUILDING    ##########
#________________________________________#
#________________________________________#

#________________________________________#
# Building the whole network based on lks.global and using the Lakes_buffer function found at "Alpine_Lakes_NetMetrics_functions.R"

load("lks.global.Rdata")

# Obtain the longitude and latitude of half europe (enough to built a buffer)
sel<-which(lks.global$lon>-6 & lks.global$lon<26 & lks.global$lat>37 & lks.global$lat<55)
dat<-data.frame(lks.global[sel, c("lon", "lat")])

# Extract the coordinates of all the sampled lakes (55 different lakes)
all_lakes_coord <- unique(rbind(env_data[[1]][,2:3],
                                env_data[[2]][,2:3],
                                env_data[[3]][,2:3],
                                env_data[[4]][,2:3]))
dat1<-data.frame(lon=all_lakes_coord[,1], lat=all_lakes_coord[,2], sampled=1)


# Transforming to sf format
lakes_sw <- st_as_sf(dat1, coords = c("lon","lat"))
# Convert to "old" format to caluclate to distances
lakes_sw_sp <- as(lakes_sw, "Spatial")
# Calculate the distance matrix (in METRES!!!)
dist_lakes <- distm(lakes_sw_sp, fun = distGeo)


#________________________________________#
# Maxmimum distance between lakes -- 600 km approximately
## This has been uploaded to the supercomputer to calculate d.percol__________________
max_distance <- Lakes_buffer(EU_lakes = dat, Samp_lakes = dat1, buffer_distance = 1, 
                             check_EU_Sampl_plots = F, 
                             check_BUFFER_plots = T,
                             check_ELIMINATION_goodLAKES_plots = F )
# We extract the distance matrix and save it to run it externally (accelerate computing times)
max_distance_down <- max_distance[[2]]
#save(max_distance_down, file = "Database/AlPS_lakes_percol_MAX_DISTANCE.RData")

### IMPORTANT MESSAGE: The following lines have been run on a computer cluster. 
# Run the following chunk in some "supercomputer" that can last some days running and keep the result!
# save.image with the "dist_lakes" object 
# Also remember that you have to upload the max.comp.gradiente_NORW

# Something like this but well written... if not do it manually deleting everything except "dist_lakes" it is not that terrible...  
#save(list =dist_lakes, file = "Database/AlPS_lakes_percol.RData")


##library(doParallel)
##registerDoParallel(cores = detectCores())
##out<-foreach(i =T , .combine=rbind) %dopar% {
##  library("sna")
##  max.comp_gradiente_NORW(min_distancias = 140000, max_distancias = 150000,dist_lakes,1000)->aa
##}
#
##aa <- out#taula amb aquests valors
##which(aa[,2]==ncol(dist_lakes))[1]->k #De la llista generada selecciona la primera dist?ncia en que tots els nodes s'ajunten (la xarx de percolaci?)
##aa[k,1]->d.percol #li donem el nom a i a?llem aquesta dist?ncia de percolaci?
##d.percol #ladist?ncia en concret

#________________________________________#
# Middle distance between lakes -- 300 km approximately
## This has been uploaded to the supercomputer to calculate d.percol__________________
mid_distance <- Lakes_buffer(EU_lakes = dat, Samp_lakes = dat1, buffer_distance = 2, 
                             check_EU_Sampl_plots = F, 
                             check_BUFFER_plots = T,
                             check_ELIMINATION_goodLAKES_plots = F )
# We extract the distance matrix and save it to run it externally (accelerate computing times)
mid_distance_down <- mid_distance[[2]]
#save(mid_distance_down, file = "Database/AlPS_lakes_percol_MID_DISTANCE.RData")

### IMPORTANT MESSAGE: The following lines have been run on a computer cluster.
# Run the following chunk in some "supercomputer" that can last some days running and keep the result!
# save.image with the "dist_lakes" object 
# Also remember that you have to upload the max.comp.gradiente_NORW
#
# Something like this but well written... if not do it manually deleting everything except "dist_lakes" it is not that terrible...  
#save(list =dist_lakes, file = "Database/AlPS_lakes_percol.RData")


##library(doParallel)
##registerDoParallel(cores = detectCores())
##out<-foreach(i =T , .combine=rbind) %dopar% {
##  library("sna")
##  max.comp_gradiente_NORW(min_distancias = 50000, max_distancias = 150000,dist_lakes,100000)->aa
##}
#
##aa <- out#taula amb aquests valors
##which(aa[,2]==ncol(dist_lakes))[1]->k #De la llista generada selecciona la primera dist?ncia en que tots els nodes s'ajunten (la xarx de percolaci?)
##aa[k,1]->d.percol #li donem el nom a i a?llem aquesta dist?ncia de percolaci?
##d.percol #ladist?ncia en concret

#________________________________________#
# Medium middle (MID_MID) distance between lakes -- 100 km approximately
## This has been uploaded to the supercomputer to calculate d.percol__________________
mid_mid_distance <- Lakes_buffer(EU_lakes = dat, Samp_lakes = dat1, buffer_distance = 6, 
                                 check_EU_Sampl_plots = F, 
                                 check_BUFFER_plots = T,
                                 check_ELIMINATION_goodLAKES_plots = F )
# We extract the distance matrix and save it to run it externally (accelerate computing times)
mid_mid_distance_down <- mid_mid_distance[[2]]
#save(mid_mid_distance_down, file = "Database/AlPS_lakes_percol_MID_MID_DISTANCE.RData")

### IMPORTANT MESSAGE: The following lines have been run on a computer cluster.
# Run the following chunk in some "supercomputer" that can last some days running and keep the result!
# save.image with the "dist_lakes" object 
# Also remember that you have to upload the max.comp.gradiente_NORW
#
# Something like this but well written... if not do it manually deleting everything except "dist_lakes" it is not that terrible...  
#save(list =dist_lakes, file = "Database/AlPS_lakes_percol.RData")


##library(doParallel)
##registerDoParallel(cores = detectCores())
##out<-foreach(i =T , .combine=rbind) %dopar% {
##  library("sna")
##  max.comp_gradiente_NORW(min_distancias = 20000, max_distancias = 150000,dist_lakes,10000)->aa
##}
#
##aa <- out#taula amb aquests valors
##which(aa[,2]==ncol(dist_lakes))[1]->k #De la llista generada selecciona la primera dist?ncia en que tots els nodes s'ajunten (la xarx de percolaci?)
##aa[k,1]->d.percol #li donem el nom a i a?llem aquesta dist?ncia de percolaci?
##d.percol #ladist?ncia en concret

#________________________________________#
# Small distance, A tenth of the maximum distance -- 60km approximately
small_distance <- Lakes_buffer(EU_lakes = dat, Samp_lakes = dat1, buffer_distance = 10, 
                               check_EU_Sampl_plots = F, 
                               check_BUFFER_plots = T,
                               check_ELIMINATION_goodLAKES_plots = F )
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
#MAX= 142222.2
#load("Database/PC_cluster_out/ALPINE_MAX_distance.RData")
#MID= 68732.87
#load("Database/PC_cluster_out/ALPINE_MID_distance.RData")
#MID_MID= 72036.2
#load("Database/PC_cluster_out/ALPINE_MID_MID_distance.RData")

ALP_MAX_xarxa <- max_distance_down
ALP_MAX_xarxa <- ifelse(max_distance_down>142222.2,0,1)

ALP_MID_MID_xarxa <- mid_mid_distance_down
ALP_MID_MID_xarxa <- ifelse(ALP_MID_MID_xarxa>72036.2,0,1)

ALP_MID_xarxa <- mid_distance_down
ALP_MID_xarxa <- ifelse(ALP_MID_xarxa>68732.87,0,1)

cordenades_xarxes <-list(max_distance[[1]],mid_distance[[1]],mid_mid_distance[[1]], small_distance[[1]], min_distance[[1]]) 
MAPS_xarxes <- list(ALP_MAX_xarxa, ALP_MID_xarxa, ALP_MID_MID_xarxa, ALP_SMALL_xarxa, ALP_MIN_xarxa)

save.image("Database.RData") 
