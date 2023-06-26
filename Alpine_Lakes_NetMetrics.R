
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

#___________________________________________________#
#___________________________________________________#
##########  ENVIRONMENTAL DATA CHECKING    ##########
#___________________________________________________#
#___________________________________________________#

# Printing all data as a table and retransforming area and altitude to normal data 
library("rnaturalearth")
library("rnaturalearthdata")
world <- ne_countries(scale = "medium", returnclass = "sf")

plot_loc <- world %>% 
  mutate(Coloret = case_when(
    str_detect(admin,"Austria") ~ "Yes", 
    str_detect(admin,"Switzerland") ~ "Yes",
    str_detect(admin,"Germany") ~ "Yes",
    TRUE ~ "No")) %>%               
  ggplot() +geom_sf(aes(fill=Coloret))+
  scale_fill_manual(values = c("grey90","grey20"))+
  coord_sf(xlim = c(min(env_data[[1]]$lon)-15, 
                    max(env_data[[1]]$lon)+15), 
           ylim = c(min(env_data[[1]]$lat)-15,
                    max(env_data[[1]]$lat)+15), expand = T)+
  theme_classic()+theme(legend.position = "none",axis.text = element_blank(),axis.ticks = element_blank())

TOTAL_env_data <- data.frame() 
TOTAL_env_names <- data.frame()
Plot_maps <- list()
Subplot_title <- c("16S_Bacterioplankton", "18S_Protista", "Phytoplankton", "Zooplankton","18S_Zooplankton")
color_groups <- as.character(CUNILLERA_cols("yellow","blue","green","red","cyan"))
for (group in 1:length(env_data)) {
TOTAL_env_data <-env_data[[group]] %>% mutate(Group=Subplot_title[[group]]) %>% bind_rows(TOTAL_env_data)}

TOTAL_env_data <- TOTAL_env_data %>% pivot_wider(names_from =Group,values_from = Group) %>% 
                   mutate(lake_area=exp(1)^lake_area,
                          Altitude=exp(1)^Altitude,
                          Cond=exp(1)^Cond) %>% 
                  rename(Lake_name=lake, Longitude=lon,Latitude=lat,Chlorophyl.a=Chl.a,Conductivity=Cond,Lake_area=lake_area) %>% 
                  distinct(Lake_name,.keep_all = T)

the_plot_db <-TOTAL_env_data %>% 
  mutate(A=ifelse(is.na(`18S_Zooplankton`)==T,0,1))%>%
  mutate(B=ifelse(is.na(Zooplankton)==T,0,1))%>% 
  mutate(C=ifelse(is.na(Phytoplankton)==T,0,1))%>% 
  mutate(D=ifelse(is.na(`18S_Protista`)==T,0,1))%>% 
  mutate(E=ifelse(is.na(`16S_Bacterioplankton`)==T,0,1)) %>% 
  rowwise() %>% 
  mutate(Count=sum(A,B,C,D,E)) %>% ungroup()
  
the_plot <- ggplot(data = world) +geom_sf()+
  coord_sf(xlim = c(min(TOTAL_env_data$Longitude)-1, 
                    max(TOTAL_env_data$Longitude)+1), 
           ylim = c(min(TOTAL_env_data$Latitude)-1,
                    max(TOTAL_env_data$Latitude)+1), expand = T)+
  theme_classic()+
  geom_point(data = the_plot_db, 
             aes(x = Longitude, y = Latitude,shape=as.factor(Count),fill =as.factor(Count)),size = 3, alpha=0.5)+
    scale_shape_manual(values = c(21,22,23,24,25))+
    scale_fill_brewer(type = "qual",direction = -1)+
  labs(title ="Sampled lakes",
       shape="Sampled groups in each lake",
       fill="Sampled groups in each lake")+
  theme(legend.position = "bottom")
  
the_plot <- ggdraw() +  draw_plot(the_plot)+draw_plot(plot_loc,scale = 0.3,x = 0.4,y = -0.14) 

write.csv2(TOTAL_env_data,"Results_tables/TOTAL_env_data.csv")

png(filename = "Figures/Sampl_Maps.png",
    width = 1030*3, height = 587*3, res=300)
grid.arrange(the_plot,ncol=1)
dev.off()

# Checking only the network of sampled lakes (not accounting with the lakes in between)
dist.percol <- list()
for (r in 1:4) {
  situacio <- env_data[[r]][,2:3]
  # Transforming to sf format
  lakes_sw_small <- st_as_sf(situacio, coords = c("lon","lat"))
  # Convert to "old" format to caluclate to distances
  lakes_sw_sp_small <- as(lakes_sw_small, "Spatial")
  # Calculate the distance matrix (in METRES!!!)
  dist_lakes_small <- distm(lakes_sw_sp_small, fun = distGeo)
  
  max.comp_gradiente(dist_lakes_small,100)->aa # the function "max.comp_gradiente" calculates for all distances (form the min to the max distance)
  # in "dist_lakes_small" a network and assess how many nodes are connected within it.
  which(aa[,2]==ncol(dist_lakes_small))[1]->k # Then, from the list generated we can obtain which is the first distance where only one component
  # is obtained (a network where all nodes are connected). This distance and the corresponding network are called percolation distance and network. 
  aa[k,1]->d.percol
  dist.percol[[r]] <- d.percol
}

# MAPS_ALPINE LAKES ENVIRONMENTAL DATA
Plot_title <- c("Temperature", "Chl.a", "Conductivity", "Lake area", "Altitude")
Plot_colors <- viridis(length(Plot_title),option = "A",alpha = 0.1)
Subplot_title <- c("S16- Bacterioplankton", "S18- Protista", "Phytoplankton", "Zooplankton")

Env_maps_list <- list()
for (varib in 4:8) {
maps_list <- list()
for (r in 1:4) {
  comunitats<- comm_data[[r]]
  situacio <- env_data[[r]][,2:3]
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
  n %v% "family" <- env_data[[r]][,varib]# Family is an standard name for the categortical variable that we are creating 
  detach("package:igraph", unload = TRUE)
  maps_list[[r]] <- ggplot(n, layout=as.matrix(situacio) ,
                           aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(color = "grey40", size=0.5)+
    geom_nodes( aes(fill=family),size=5 ,color="black", shape=21, alpha=.85)+
    labs(x="",y="",title = Subplot_title[r])+
    scale_fill_viridis()+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right",
          legend.title = element_blank(),
          panel.background = element_rect(fill=Plot_colors[varib-3]))
}  
Env_maps_list[[varib]] <- grid.arrange(maps_list[[1]],maps_list[[2]],maps_list[[3]],
                                       maps_list[[4]], ncol = 2, nrow=2,
                                       top=Plot_title[varib-3])

}

# Printing plot
png(filename = "Figures/Env_Maps.png",
    width = 5000, height = 5000, res=300)
grid.arrange(Env_maps_list[[4]],Env_maps_list[[5]],Env_maps_list[[6]],
             Env_maps_list[[7]],Env_maps_list[[8]])
dev.off()

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

#________________________________________#
#________________________________________#
######## RIVER NETWORK BUILDING ##########
#________________________________________#
#________________________________________#

# We first will print the xy coordinates of all the lakes that we could detect with the lakes_buffer function and will later find all these
#lakes within the shapefile of the two rivers. Therefore, we will be able to find which points of the rivers correspond to the lakes 
#from which we have information. 

# Print de coordinates of the two smaller scales to posteriorly trim the "river" network in Arcgis
#write.dbf(max_distance[[1]], "GIS_data/max_distance.dbf")
max_distance_BASINS <- read.dbf("GIS_data/max_distance_BASINS.dbf")

cordenades_xarxes_BASINS <- list(max_distance_BASINS)

GRAPH_xarxes_fluvial <- list()
all_lakes_BASINS_fluvial <- list()
correspondence_BASINS_fluvial <- list()

detach("package:sna", unload = TRUE)
library(igraph)  
library(shp2graph)

correspondence_BASINS <- c()
for (a in 1:55) {
  correspondence_BASINS[a] <- which(round(cordenades_xarxes[[1]][(nrow(cordenades_xarxes[[1]])-55+a),1],4)==
                                    round(cordenades_xarxes_BASINS[[1]][,1],4))  
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

png(filename = "Figures/Alpine_fluvial.png",
    width = 20000, height = 20000, res=500)
par(mar=c(0,0,0,0))
plot(GRAPH_xarxes_fluvial[[1]], vertex.label = NA, vertex.size = sizes, vertex.size2 = sizes, vertex.color=cols,
     edge.width=0.01,edge.size=0.01, edge.color="grey50")
dev.off()
detach("package:shp2graph", unload = TRUE)

site_cord <- matrix(unlist(rtNEL1[[2]][,2]), ncol=2, byrow=TRUE)
colnames(site_cord) <- c("lon","lat")
cordenades_xarxes[[6]] <- as.data.frame(site_cord)
sizes_points <- c(2,2,2,2,2,0.1)
scale_names <- c("650 km","325 km"," 100 km","65 km","6.5 km","Fluvial")
lakes_loc_plot <- list()
for (lake in 1:length(cordenades_xarxes)) {
  lakes_loc_plot[[lake]] <-   ggplot(data = world) +geom_sf()+
    coord_sf(xlim = c(min(max_distance[[1]][,1])-1, 
                      max(site_cord[,1])+1), 
             ylim = c(min(site_cord[,2])-1,
                      max(max_distance[[1]][,2])+1), expand = T)+
    theme_classic()+
    geom_point(data = cordenades_xarxes[[lake]],aes(x = lon, y = lat), size = sizes_points[lake],shape = 21, fill ="black", alpha=0.5)+
    labs(title =scale_names[[lake]])+xlab("Longitude")+ylab("Latitude")
}

png(filename = "Figures/AlpineLakes_scales.png",
    width = 5000, height = 5000, res=300)
grid.arrange(lakes_loc_plot[[1]],lakes_loc_plot[[2]],lakes_loc_plot[[3]],lakes_loc_plot[[4]],lakes_loc_plot[[5]],lakes_loc_plot[[6]],ncol=2)
dev.off()

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
  
  network_data[[i]][,1] <- igraph::closeness(Xarx_grpah, mode = "all")
  network_data[[i]][,2] <- igraph::degree(Xarx_grpah)
  network_data[[i]][,3] <- igraph::betweenness(Xarx_grpah, directed = F)
  
  if(components(Xarx_grpah)$no>1){
    # Closenness
    a <- which(components(Xarx_grpah)$membership==1)
    b <- which(components(Xarx_grpah)$membership==2)
    
    Xarx_grpah <- graph.adjacency(MAPS_xarxes[[i]][a,a], mode = "undirected",diag = F)
    network_data[[i]][a,1] <- igraph::closeness(Xarx_grpah, mode = "all")/max(igraph::closeness(Xarx_grpah, mode = "all"))
   
    Xarx_grpah <- graph.adjacency(MAPS_xarxes[[i]][b,b], mode = "undirected",diag = F)
    network_data[[i]][b,1] <- igraph::closeness(Xarx_grpah, mode = "all")/max(igraph::closeness(Xarx_grpah, mode = "all"))
  }

  colnames(network_data[[i]]) <- c("clo_ALPS","deg_ALPS","bet_ALPS")
}

# Flulvial
detach("package:sna", unload = TRUE)
library(igraph)
fluvial_network_data <- list()
output <- matrix(nrow = length(igraph::V(GRAPH_xarxes_fluvial[[1]])), ncol = 4)
fluvial_network_data[[1]] <- output

fluvial_network_data[[1]][,1] <- 1-igraph::harmonic_centrality(GRAPH_xarxes_fluvial[[1]], mode = "out")/max(igraph::harmonic_centrality(GRAPH_xarxes_fluvial[[1]], mode = "out"))
fluvial_network_data[[1]][,2] <- igraph::degree(GRAPH_xarxes_fluvial[[1]],mode = "out")
fluvial_network_data[[1]][,3] <- igraph::betweenness(GRAPH_xarxes_fluvial[[1]])
fluvial_network_data[[1]][,4] <- igraph::closeness(GRAPH_xarxes_fluvial[[1]],mode = "out",normalized = T)

a <- which(components(GRAPH_xarxes_fluvial[[1]])$membership==1)
b <- which(components(GRAPH_xarxes_fluvial[[1]])$membership==2)

Basins_Euro <- igraph::decompose.graph(GRAPH_xarxes_fluvial[[1]])

a_clos <-igraph::closeness(Basins_Euro[[1]],mode = "out", normalized = T)
b_clos <-igraph::closeness(Basins_Euro[[2]],mode = "out", normalized = T)

fluvial_network_data[[1]][a,4]<- a_clos/max(a_clos[which(is.na(a_clos)!=T)])
fluvial_network_data[[1]][b,4]<- b_clos/max(b_clos[which(is.na(b_clos)!=T)])

colnames(fluvial_network_data[[1]]) <- c("HarC_ALPS","deg_ALPS","bet_ALPS", "clo_ALPS")

library(ggfortify)
Network_results <- list()
for (r in 1:length(network_data)) {
  Network_results[[r]] <-network_data[[r]][,1] #PCA_result$x[,1] 
}
names(Network_results) <- c("650 km","325 km"," 100 km","65 km","6.5 km")

#Plots closeness maps - pallete viridis
maps_Alps <- list()
detach("package:igraph", unload = TRUE)
library(sna)
for (e in 1:(length(cordenades_xarxes))) {
  factors <- rep("No_Sampled",nrow(MAPS_xarxes[[e]]))
  factors[c(nrow(MAPS_xarxes[[e]])-54):nrow(MAPS_xarxes[[e]])] <- "Sampled"
  CC_values <- Network_results[[e]]
  
  n<- network(MAPS_xarxes[[e]], directed=F, diag=F)
  n %v% "family" <- factors # Family is an standard name for the categortical variable that we are creating
  n %v% "CC_values" <- CC_values
  
  maps_Alps[[e]] <- ggplot(n, layout=as.matrix(cordenades_xarxes[[e]][,1:2]),
                           aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges( color = "grey60", size=0.1, alpha=0.4) +
    geom_nodes(aes(fill=CC_values,alpha=family, size=family, color=family) ,shape=21, alpha=.75)+
    scale_fill_viridis(discrete = F,alpha = 1,begin = 1,end = 0)+
    labs(x="",y="")+
    scale_alpha_manual(values = c(0.2,5))+
    scale_size_manual(values = c(1,3))+
    scale_color_manual(values = c("grey30","red"))+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
}

png(filename ="Figures/Network_maps_Alps.png",
    width =745 ,height =742 ,units ="px",res = 100)
grid.arrange(maps_Alps[[1]],maps_Alps[[2]],maps_Alps[[3]],maps_Alps[[4]],maps_Alps[[5]])
dev.off()

png(filename ="Figures/Network_maps_Alps_max.png",
    width =600 ,height =500 ,units ="px",res = 200)
grid.arrange(maps_Alps[[1]])
dev.off()

png(filename ="Figures/Network_maps_Alps_med.png",
    width =600 ,height =500 ,units ="px",res = 200)
grid.arrange(maps_Alps[[2]])
dev.off()

png(filename ="Figures/Network_maps_Alps_med_med.png",
    width =600 ,height =500 ,units ="px",res = 200)
grid.arrange(maps_Alps[[3]])
dev.off()

png(filename ="Figures/Network_maps_Alps_small.png",
    width =600 ,height =500 ,units ="px",res = 200)
grid.arrange(maps_Alps[[4]])
dev.off()

png(filename ="Figures/Network_maps_Alps_min.png",
    width =600 ,height =500 ,units ="px",res = 200)
grid.arrange(maps_Alps[[5]])
dev.off()


library(ggfortify)
Fluvial_network_results <- list()
Fluvial_network_results[[1]] <-fluvial_network_data[[1]][,1]
names(Fluvial_network_results) <- c("Fluvial")

#Plot centrality values in river network
cols <-ConComp$membership 
cols[unlist(res.nv[[3]])[correspondence_BASINS]] <- "red"
cols <- ifelse(cols==1, "green",ifelse(cols==2,"blue","red"))
sizes <- ifelse(cols=="red",3,0.6)

rbPal <- viridis_pal(alpha = 1,direction = -1)
#Plots closeness "all"
# Relevant that "png" output size will be big (nice arrows plotted)
#This adds a column of color values based on the y values
Col <- rbPal(length(fluvial_network_data[[1]][,1]))[
     as.numeric(cut(fluvial_network_data[[1]][,1],
    breaks = length(fluvial_network_data[[1]][,1])))]
Col[which(cols=="red")] <- "red"
detach("package:sna", unload = TRUE)
library(igraph)
png(filename = "Figures/Alpine_fluvial_CLOS.png", width = 20000, height = 20000, res=500)
par(mar=c(0,0,0,0))
plot(GRAPH_xarxes_fluvial[[1]], vertex.label = NA, vertex.size = sizes, vertex.size2 = sizes, vertex.color=Col,
     edge.width=0.01, edge.color="grey50")
dev.off()


# Print network structure results_______________________________________________
detach("package:sna", unload = TRUE)
library(igraph)

out_netw_struc <- matrix(nrow = 6, ncol = 8)

for (i in 1:length(MAPS_xarxes)) {
Xarx_grpah <- graph.adjacency(MAPS_xarxes[[i]], mode = "undirected",diag = F)

out_netw_struc[i,1] <- mean(degree(Xarx_grpah))
out_netw_struc[i,2] <- edge_density(Xarx_grpah)
out_netw_struc[i,3] <- diameter(Xarx_grpah)
out_netw_struc[i,4] <- length(V(Xarx_grpah))
out_netw_struc[i,5] <- mean(Network_results[[i]])
out_netw_struc[i,6] <- sd(Network_results[[i]])
out_netw_struc[i,7] <- max(Network_results[[i]])
out_netw_struc[i,8] <- min(Network_results[[i]])
}

out_netw_struc[6,1] <- mean(degree(GRAPH_xarxes_fluvial[[1]]))
out_netw_struc[6,2] <- edge_density(GRAPH_xarxes_fluvial[[1]])
out_netw_struc[6,3] <- diameter(GRAPH_xarxes_fluvial[[1]])
out_netw_struc[6,4] <- length(V(GRAPH_xarxes_fluvial[[1]]))
out_netw_struc[6,5] <- mean(fluvial_network_data[[1]][,1])
out_netw_struc[6,6] <- sd(fluvial_network_data[[1]][,1])
out_netw_struc[6,7] <- max(fluvial_network_data[[1]][,1])
out_netw_struc[6,8] <- min(fluvial_network_data[[1]][,1])

colnames(out_netw_struc) <- c("Linkage density", "Connectance", "Diameter","Number of nodes" ,"Mean CC", "sd", "Max CC", "Min CC")
rownames(out_netw_struc) <- c(names(Network_results), "Fluvial")
out_netw_struc

write.csv2(out_netw_struc, file = "Results_Tables/out_netw_struc.csv")

#________________________________________#
#________________________________________#
########    DIVERSITY VALUES    ##########
#________________________________________#
#________________________________________#

#### Community indices list to gather everything 
community_indices<- list()

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
  taxa.q_Ruziska_AB <- beta.div.comp(comm_data[[r]], coef = "J", quant = FALSE)
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
biod[[r]] <- cbind(community_indices[[2]][[r]],community_indices[[3]][[r]],community_indices[[4]][[r]],community_indices[[5]][[r]])
noms <-c(names(community_indices[[2]])[r],names(community_indices[[3]])[r],names(community_indices[[4]])[r],names(community_indices[[5]])[r]) 
colnames(biod[[r]]) <- noms
} 
biod

image_list <- list("s16_image.png", "s18_image.png","phy_image.png","zoo_image.png","zooS18_image.png")

TOTAL_biod_plot <- list()
Subplot_title <- c("16S_Bacterioplankton", "18S_Protista", "Phytoplankton", "Zooplankton","18S_Zooplankton")
color_groups <- as.character(CUNILLERA_cols("yellow","blue","green","red","cyan"))
for (group in 1:length(biod)) {
ComStrMetr <- colnames(biod[[group]])
Plot_maps_biod <- list()
for (variab in 1:length(ComStrMetr)) {
the_plot <- as.data.frame(biod[[group]]) %>%tibble::rownames_to_column() %>% 
            left_join(env_data[[group]], by=c("rowname"="lake")) %>% 
            #Plot
            ggplot()+
            geom_point(aes(x = lon, y = lat, size=biod[[group]][,variab]),shape = 21, fill =color_groups[[group]], alpha=0.5)+
              labs(title =Subplot_title[[group]],size=ComStrMetr[[variab]])+
              theme_classic()

Plot_maps_biod[[variab]] <-ggdraw() +  draw_plot(the_plot)+ 
                                       draw_image(magick::image_read(image_list[[group]]),scale = 0.2,x = 0.4,y = 0.38) 
}
TOTAL_biod_plot[[group]] <- Plot_maps_biod
}

Out_Biod <- data.frame()
for (r in 1:length(biod)) {
cols_biod <- colnames(biod[[r]])
temp <- as.data.frame(biod[[r]]) %>%
  rename("Rich"=cols_biod[[1]],"LCBD"=cols_biod[[2]],"Repl"=cols_biod[[3]],"RicDif"=cols_biod[[4]])%>% 
  #pivot_longer(cols =c("Rich","LCBD","Repl","RicDif")) %>% 
  #group_by(name) %>%
  mutate(Taxon_Group=Subplot_title[[r]], Taxon_Lake=paste(rownames(biod[[r]]),Subplot_title[[r]]))
Out_Biod <- bind_rows(Out_Biod,temp) 
}
write.csv2(Out_Biod,"Results_Tables/Out_Biod.csv")


Out_Mean_Biod <- data.frame()
for (r in 1:length(biod)) {
cols_biod <- colnames(biod[[r]])
temp <- as.data.frame(biod[[r]]) %>% rename("Rich"=cols_biod[[1]],
                                            "LCBD"=cols_biod[[2]],
                                            "Repl"=cols_biod[[3]],
                                            "RicDif"=cols_biod[[4]])%>% 
  pivot_longer(cols =c("Rich","LCBD","Repl","RicDif")) %>% 
  group_by(name) %>% summarise(Mean=mean(value),sd=sd(value)) %>% 
  mutate(Taxon_Group=Subplot_title[[r]])
Out_Mean_Biod <- bind_rows(Out_Mean_Biod,temp) 
}
Out_Mean_Biod <- pivot_wider(Out_Mean_Biod,names_from = name,values_from = c(Mean,sd))
write.csv2(Out_Mean_Biod,"Results_Tables/Out_Mean_Biod.csv")

      
png(filename = "Figures/Res_Biod_Maps.png", width = 20000, height = 12000, res=550)
grid.arrange(
grid.arrange(grobs = TOTAL_biod_plot[[1]], ncol = 4),
grid.arrange(grobs = TOTAL_biod_plot[[2]], ncol = 4),
grid.arrange(grobs = TOTAL_biod_plot[[3]], ncol = 4),
grid.arrange(grobs = TOTAL_biod_plot[[4]], ncol = 4),
grid.arrange(grobs = TOTAL_biod_plot[[5]], ncol = 4),
nrow=5)
dev.off()

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

histo_list <- list()
for (histo in 1:5) {
coin <- Network_results[[histo]][c(length(Network_results[[histo]])-54):length(Network_results[[histo]])]
histo_list[[histo]] <- ggplot(data.frame(y=Network_results[[histo]]))+
  geom_histogram(aes(x=y), bins=10,colour="black",fill="grey60")+
  geom_vline(xintercept = coin, colour="black")+
  scale_y_continuous(expand = c(0,0))+theme_classic()+labs(x="Closeness centrality", title = names(Network_results)[histo])
}

coin <- Fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
histo_list[[6]] <- ggplot(data.frame(y=fluvial_network_data[[1]][,1]))+
  geom_histogram(aes(x=y), bins=10,colour="black",fill="grey60")+
  geom_vline(xintercept = coin, colour="black")+
  scale_y_continuous(expand = c(0,0))+theme_classic()+labs(x="Out-closeness", title = "Fluvial")


png(filename ="Figures/Clos_Histo.png",
    width =582*4 ,height =629*4 ,units ="px",res = 300)
grid.arrange(histo_list[[1]],histo_list[[2]],histo_list[[3]],histo_list[[4]],histo_list[[5]],histo_list[[6]])
dev.off()


# GAM models_______________________________________________________________________________ ####
Names_Networks <- c("650 km", "325 km","100 km","65 km","6.5 km")

GAMmodel_resutls_total <- list()
GAM_model_resutls <- list()
GAM_direct_model <- list()
GAM_direct_model_total <- list()

# For Euclidean networks with 600, 300 100, 60 and 6 km
for (groups in 1:5) {
  color_groups <- as.character(CUNILLERA_cols("yellow","blue","green","red","cyan"))
  for (net in 1:5) {
    coin <- Network_results[[net]][c(length(Network_results[[net]])-54):length(Network_results[[net]])]
    dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                     biod[[groups]][,1:4])
    colnames(dataset)[1] <-c("Network")
    
    plots_grups <- list()
    output_results <- list()
    output_model_results <- list()
    p.val <- c()
   
    #Richness ___________________________________________________________________________________________
    # Variance Inflation Factor
    vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
    Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
    # Random Forest (RF) - selection of most informative variables for each DIV metric + group
    Model_RF = cbind(dataset[,2],Sel_vari) 
    colnames(Model_RF)[1] <- "Rich"
    RF_output <- rfsrc(Rich ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
    # Checking plots (not necessary to run every time)
    #S.rf.vimp <- gg_vimp (S.rf)
    #S.rf.vimp
    #plot (S.rf.vimp) #plot variable importance
    out_selected_vari <- var.select(RF_output,conservative = "low")
    env_selected_var <- c()
    for (u in 1:length(out_selected_vari$topvars)) {
      env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
    }
    env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
    model <- lm(dataset[,2]~env_data_model[,1:length(ncol(env_data_model))])
    summary(model)
    resid_values <- model$residuals
    dataset[,2] <- resid_values
    
    #GAM model
    p.val[1] <- summary.gam(gam(dataset[,2] ~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[8]]
    output_results[[1]] <- summary.gam(gam(dataset[,2] ~ s(dataset[,1], k=3, bs="cr"), method = "REML"))
    output_model_results[[1]] <- gam(dataset[,2] ~ s(dataset[,1], k=3, bs="cr"), method = "REML")
    preds_1<- predict(gam(dataset[,2] ~ s(dataset[,1], k=3, bs="cr"), method = "REML"), se.fit = TRUE)
    
    #LCBD ___________________________________________________________________________________________
    # Variance Inflation Factor
    vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
    Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
    # Random Forest (RF) - selection of most informative variables for each DIV metric + group
    Model_RF = cbind(dataset[,3],Sel_vari) 
    colnames(Model_RF)[1] <- "LCBD"
    RF_output <- rfsrc(LCBD ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
    # Checking plots (not necessary to run every time)
    #S.rf.vimp <- gg_vimp (S.rf)
    #S.rf.vimp
    #plot (S.rf.vimp) #plot variable importance
    out_selected_vari <- var.select(RF_output,conservative = "low")
    env_selected_var <- c()
    for (u in 1:length(out_selected_vari$topvars)) {
      env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
    }
    env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
    model <- lm(dataset[,3]~env_data_model[,1:length(ncol(env_data_model))])
    summary(model)
    resid_values <- model$residuals
    dataset[,3] <- resid_values
    #GAM model
    p.val[2] <-summary.gam(gam(dataset[,3]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[8]]
    output_results[[2]] <- summary.gam(gam(dataset[,3]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))
    output_model_results[[2]] <- gam(dataset[,3]~ s(dataset[,1], k=3, bs="cr"), method = "REML")
    preds_2<- predict(gam(dataset[,3]~ s(dataset[,1], k=3, bs="cr"), method = "REML"), se.fit = TRUE)
    
    #Turn ___________________________________________________________________________________________
    # Variance Inflation Factor
    vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
    Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
    # Random Forest (RF) - selection of most informative variables for each DIV metric + group
    Model_RF = cbind(dataset[,4],Sel_vari) 
    colnames(Model_RF)[1] <- "Turn"
    RF_output <- rfsrc(Turn ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
    # Checking plots (not necessary to run every time)
    #S.rf.vimp <- gg_vimp (S.rf)
    #S.rf.vimp
    #plot (S.rf.vimp) #plot variable importance
    out_selected_vari <- var.select(RF_output,conservative = "low")
    env_selected_var <- c()
    for (u in 1:length(out_selected_vari$topvars)) {
      env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
    }
    env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
    model <- lm(dataset[,4]~env_data_model[,1:length(ncol(env_data_model))])
    summary(model)
    resid_values <- model$residuals
    dataset[,4] <- resid_values
    # GAM model
    p.val[3] <-summary.gam(gam(dataset[,4]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[8]]
    output_results[[3]] <- summary.gam(gam(dataset[,4]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))
    output_model_results[[3]] <- gam(dataset[,4]~ s(dataset[,1], k=3, bs="cr"), method = "REML")
    preds_3<- predict(gam(dataset[,4]~ s(dataset[,1], k=3, bs="cr"), method = "REML"), se.fit = TRUE)
    
    #RichDiff ___________________________________________________________________________________________
    # Variance Inflation Factor
    vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
    Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
    # Random Forest (RF) - selection of most informative variables for each DIV metric + group
    Model_RF = cbind(dataset[,5],Sel_vari) 
    colnames(Model_RF)[1] <- "RichDiff"
    RF_output <- rfsrc(RichDiff ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
    # Checking plots (not necessary to run every time)
    #S.rf.vimp <- gg_vimp (S.rf)
    #S.rf.vimp
    #plot (S.rf.vimp) #plot variable importance
    out_selected_vari <- var.select(RF_output,conservative = "low")
    env_selected_var <- c()
    for (u in 1:length(out_selected_vari$topvars)) {
      env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
    }
    env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
    model <- lm(dataset[,5]~env_data_model[,1:length(ncol(env_data_model))])
    summary(model)
    resid_values <- model$residuals
    dataset[,5] <- resid_values
    # GAM model 
    p.val[4] <-summary.gam(gam(dataset[,5]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[8]]
    output_results[[4]] <- summary.gam(gam(dataset[,5]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))
    output_model_results[[4]] <- gam(dataset[,5]~ s(dataset[,1], k=3, bs="cr"), method = "REML")
    preds_4<- predict(gam(dataset[,5]~ s(dataset[,1], k=3, bs="cr"), method = "REML"), se.fit = TRUE)
    
    GAM.pred <- list(preds_1, preds_2, preds_3, preds_4)
    
    for(var in 1:4){
      if(p.val[var]>0.05){
        my_data <- data.frame(cbind(dataset[,var+1],dataset[,1]),
                              mu   = GAM.pred[[var]]$fit,
                              low  = GAM.pred[[var]]$fit - 1.96 * GAM.pred[[var]]$se.fit,
                              high = GAM.pred[[var]]$fit + 1.96 * GAM.pred[[var]]$se.fit)
          the_plot<-
          ggplot(my_data, aes(x = X2, y = X1)) +
          geom_jitter(alpha=0.2, shape=21, size=3, colour="black", aes(fill=X2))+
          scale_fill_continuous(type = "viridis",  alpha = 1,begin = 1,end = 0)+
          geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="black",linetype=2, size=2)+
          labs(title=colnames(dataset)[var+1], subtitle = paste("R2=", round(output_results[[var]][[10]],2),
                                                                "ED=",round(output_results[[var]][[14]],2)))+
          ylab(colnames(dataset)[var+1])+xlab("Centrality-Isolation")+
          theme_classic()+
          theme(panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)),
                legend.position = "none")

          plots_grups[[var]]<- ggdraw() +  draw_plot(the_plot)+
                                            draw_image(magick::image_read(image_list[[groups]]),
                                                       scale = 0.2,x = 0.4,y = 0.38) 

      }else{
        my_data <- data.frame(cbind(dataset[,var+1],dataset[,1]),
                              mu   = GAM.pred[[var]]$fit,
                              low  = GAM.pred[[var]]$fit - 1.96 * GAM.pred[[var]]$se.fit,
                              high = GAM.pred[[var]]$fit + 1.96 * GAM.pred[[var]]$se.fit)
        the_plot<-
          ggplot(my_data, aes(x = X2, y = X1)) +
          geom_jitter(alpha=0.9, shape=21, size=3, colour="black", aes(fill=X2))+
          scale_fill_continuous(type = "viridis",  alpha = 1,begin = 1,end = 0)+
          geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="black",linetype=1, size=2)+
          labs(title=colnames(dataset)[var+1],subtitle = paste("R2=", round(output_results[[var]][[10]],2),
                                                               "ED=",round(output_results[[var]][[14]],2)))+
          ylab(colnames(dataset)[var+1])+xlab("Centrality-Isolation")+
          theme_classic()+
          theme(panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)),
                legend.position = "none")  
        
        plots_grups[[var]]<- ggdraw() +  draw_plot(the_plot)+
                            draw_image(magick::image_read(image_list[[groups]]),
                                                  scale = 0.2,x = 0.4,y = 0.38) 
        
      }
    }
    
    GAM_model_resutls[[net]] <-output_results
    GAM_direct_model[[net]]<- output_model_results
    
    png(filename =paste("Figures/","GAM_Divers",biod_names[[groups]],"_",Names_Networks[[net]],".png"),
        width =582*4 ,height =629*4 ,units ="px",res = 300)
    grid.arrange(plots_grups[[1]],plots_grups[[2]],plots_grups[[3]],plots_grups[[4]],
                 ncol=2,nrow=2, top=Names_Networks[[net]])
    dev.off()
  }
  GAMmodel_resutls_total[[groups]] <- GAM_model_resutls
  GAM_direct_model_total[[groups]]<- GAM_direct_model
}             

# For fluvial network
GAMmodel_resutls_fluvial <- list()
GAM_direct_model_fluvial<- list()

for (groups in 1:5) {
  color_groups <- as.character(CUNILLERA_cols("yellow","blue","green","red","cyan"))
  coin <- Fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
  dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                   biod[[groups]][,1:4])
  colnames(dataset)[1] <-c("Network")
  plots_grups <- list()
  output_results <- list()
  output_model_results <- list()
  
  p.val <- c()
  
  #Richness ___________________________________________________________________________________________
  # Variance Inflation Factor
  vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
  Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
  # Random Forest (RF) - selection of most informative variables for each DIV metric + group
  Model_RF = cbind(dataset[,2],Sel_vari) 
  colnames(Model_RF)[1] <- "Rich"
  RF_output <- rfsrc(Rich ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
  # Checking plots (not necessary to run every time)
  #S.rf.vimp <- gg_vimp (S.rf)
  #S.rf.vimp
  #plot (S.rf.vimp) #plot variable importance
  out_selected_vari <- var.select(RF_output,conservative = "low")
  env_selected_var <- c()
  for (u in 1:length(out_selected_vari$topvars)) {
    env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
  }
  env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
  model <- lm(dataset[,2]~env_data_model[,1:length(ncol(env_data_model))])
  summary(model)
  resid_values <- model$residuals
  dataset[,2] <- resid_values
  #GAM model
  p.val[1] <- summary.gam(gam(dataset[,2] ~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[8]]
  output_results[[1]] <- summary.gam(gam(dataset[,2] ~ s(dataset[,1], k=3, bs="cr"), method = "REML"))
  output_model_results[[1]] <- gam(dataset[,2] ~ s(dataset[,1], k=3, bs="cr"), method = "REML")
  preds_1<- predict(gam(dataset[,2] ~ s(dataset[,1], k=3, bs="cr"), method = "REML"), se.fit = TRUE)
  
  #LCBD ___________________________________________________________________________________________
  # Variance Inflation Factor
  vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
  Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
  # Random Forest (RF) - selection of most informative variables for each DIV metric + group
  Model_RF = cbind(dataset[,3],Sel_vari) 
  colnames(Model_RF)[1] <- "LCBD"
  RF_output <- rfsrc(LCBD ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
  # Checking plots (not necessary to run every time)
  #S.rf.vimp <- gg_vimp (S.rf)
  #S.rf.vimp
  #plot (S.rf.vimp) #plot variable importance
  out_selected_vari <- var.select(RF_output,conservative = "low")
  env_selected_var <- c()
  for (u in 1:length(out_selected_vari$topvars)) {
    env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
  }
  env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
  model <- lm(dataset[,3]~env_data_model[,1:length(ncol(env_data_model))])
  summary(model)
  resid_values <- model$residuals
  dataset[,3] <- resid_values
  #GAM model
  p.val[2] <-summary.gam(gam(dataset[,3]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[8]]
  output_results[[2]] <- summary.gam(gam(dataset[,3]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))
  output_model_results[[2]] <- gam(dataset[,3]~ s(dataset[,1], k=3, bs="cr"), method = "REML")
  preds_2<- predict(gam(dataset[,3]~ s(dataset[,1], k=3, bs="cr"), method = "REML"), se.fit = TRUE)
  
  #Turn ___________________________________________________________________________________________
  # Variance Inflation Factor
  vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
  Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
  # Random Forest (RF) - selection of most informative variables for each DIV metric + group
  Model_RF = cbind(dataset[,4],Sel_vari) 
  colnames(Model_RF)[1] <- "Turn"
  RF_output <- rfsrc(Turn ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
  # Checking plots (not necessary to run every time)
  #S.rf.vimp <- gg_vimp (S.rf)
  #S.rf.vimp
  #plot (S.rf.vimp) #plot variable importance
  out_selected_vari <- var.select(RF_output,conservative = "low")
  env_selected_var <- c()
  for (u in 1:length(out_selected_vari$topvars)) {
    env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
  }
  env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
  model <- lm(dataset[,4]~env_data_model[,1:length(ncol(env_data_model))])
  summary(model)
  resid_values <- model$residuals
  dataset[,4] <- resid_values
  # GAM model
  p.val[3] <-summary.gam(gam(dataset[,4]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[8]]
  output_results[[3]] <- summary.gam(gam(dataset[,4]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))
  output_model_results[[3]] <- gam(dataset[,4]~ s(dataset[,1], k=3, bs="cr"), method = "REML")
  preds_3<- predict(gam(dataset[,4]~ s(dataset[,1], k=3, bs="cr"), method = "REML"), se.fit = TRUE)
  
  #RichDiff ___________________________________________________________________________________________
  # Variance Inflation Factor
  vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
  Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
  # Random Forest (RF) - selection of most informative variables for each DIV metric + group
  Model_RF = cbind(dataset[,5],Sel_vari) 
  colnames(Model_RF)[1] <- "RichDiff"
  RF_output <- rfsrc(RichDiff ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
  # Checking plots (not necessary to run every time)
  #S.rf.vimp <- gg_vimp (S.rf)
  #S.rf.vimp
  #plot (S.rf.vimp) #plot variable importance
  out_selected_vari <- var.select(RF_output,conservative = "low")
  env_selected_var <- c()
  for (u in 1:length(out_selected_vari$topvars)) {
    env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
  }
  env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
  model <- lm(dataset[,5]~env_data_model[,1:length(ncol(env_data_model))])
  summary(model)
  resid_values <- model$residuals
  dataset[,5] <- resid_values
  # GAM model 
  p.val[4] <-summary.gam(gam(dataset[,5]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[8]]
  output_results[[4]] <- summary.gam(gam(dataset[,5]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))
  output_model_results[[4]] <- gam(dataset[,5]~ s(dataset[,1], k=3, bs="cr"), method = "REML")
  preds_4<- predict(gam(dataset[,5]~ s(dataset[,1], k=3, bs="cr"), method = "REML"), se.fit = TRUE)
  
  GAM.pred <- list(preds_1, preds_2, preds_3, preds_4)
  
  for(var in 1:4){
    if(p.val[var]>0.05){
      my_data <- data.frame(cbind(dataset[,var+1],dataset[,1]),
                            mu   = GAM.pred[[var]]$fit,
                            low  = GAM.pred[[var]]$fit - 1.96 * GAM.pred[[var]]$se.fit,
                            high = GAM.pred[[var]]$fit + 1.96 * GAM.pred[[var]]$se.fit)
     the_plot <- 
        ggplot(my_data, aes(x = X2, y = X1)) +
        geom_jitter(alpha=0.2, shape=21, size=3, colour="black", aes(fill=X2))+
        scale_fill_continuous(type = "viridis",  alpha = 1,begin = 1,end = 0)+
        geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="black",linetype=2, size=2)+
        labs(title=colnames(dataset)[var+1],subtitle = paste("R2=", round(output_results[[var]][[10]],2),
                                                             "ED=",round(output_results[[var]][[14]],2)))+
        ylab(colnames(dataset)[var+1])+xlab("Centrality-Isolation")+
        theme_classic()+
        theme(panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)),
              legend.position = "none")  
     
        plots_grups[[var]] <- ggdraw() +  draw_plot(the_plot)+
          draw_image(magick::image_read(image_list[[groups]]),
                     scale = 0.2,x = 0.4,y = 0.38) 
    }else{
      my_data <- data.frame(cbind(dataset[,var+1],dataset[,1]),
                            mu   = GAM.pred[[var]]$fit,
                            low  = GAM.pred[[var]]$fit - 1.96 * GAM.pred[[var]]$se.fit,
                            high = GAM.pred[[var]]$fit + 1.96 * GAM.pred[[var]]$se.fit)
      the_plot <- 
        ggplot(my_data, aes(x = X2, y = X1)) +
        geom_jitter(alpha=0.9, shape=21, size=3, colour="black", aes(fill=X2))+
        scale_fill_continuous(type = "viridis",  alpha = 1,begin = 1,end = 0)+
        geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="black",linetype=1, size=2)+
        labs(title=colnames(dataset)[var+1],subtitle =paste("R2=", round(output_results[[var]][[10]],2),
                                                            "ED=",round(output_results[[var]][[14]],2)))+
        ylab(colnames(dataset)[var+1])+xlab("Centrality-Isolation")+
        theme_classic()+
        theme(panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)),
              legend.position = "none") 
      
      plots_grups[[var]] <- ggdraw() +  draw_plot(the_plot)+
        draw_image(magick::image_read(image_list[[groups]]),
                   scale = 0.2,x = 0.4,y = 0.38) 
    }
  }
  
  GAMmodel_resutls_fluvial[[groups]] <-output_results
  GAM_direct_model_fluvial[[groups]] <-output_model_results
  
  png(filename =paste("Figures/","GAM_Divers",biod_names[[groups]],"_Fluvial",".png"),
      width =582*4 ,height =629*4 ,units ="px",res = 300)
  grid.arrange(plots_grups[[1]],plots_grups[[2]],
               plots_grups[[3]],plots_grups[[4]],
               ncol=2,nrow=2, top="Fluvial network")
  dev.off()
}


# GAM models result in table format - Supplementary like___________________________________ ####
biod_names <- c("S16","S18","phy","zoo", "zoo.18S")
Names_Networks <- c("650 km", "325 km","100 km","65 km","6.5 km", "Fluvial")

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
      for (var in 1:4) {
        row_reference <- row_reference+1
        
        Netw_value <- Names_Networks[net]
        Variable_value<- colnames(biod[[groups]])[var]
        
        a[row_reference,] <- c(Netw_value,Variable_value,
                               GAMmodel_resutls_total[[groups]][[net]][[var]]$p.coeff,
                               GAMmodel_resutls_total[[groups]][[net]][[var]]$se[1],
                               GAMmodel_resutls_total[[groups]][[net]][[var]]$p.t,
                               GAMmodel_resutls_total[[groups]][[net]][[var]]$p.pv,
                               GAMmodel_resutls_total[[groups]][[net]][[var]]$chi.sq,
                               GAMmodel_resutls_total[[groups]][[net]][[var]]$s.pv,
                               GAMmodel_resutls_total[[groups]][[net]][[var]]$r.sq,
                               GAMmodel_resutls_total[[groups]][[net]][[var]]$dev.expl) 
  }
 }
  table_groups[[groups]] <- as.data.frame(a)
  write.csv2(table_groups[[groups]], file = paste("Results_Tables/",biod_names[groups],"GAM_Results",".csv",sep = ""),dec = ",")
  #write.table(table_groups[[groups]], file = paste("Results_Tables/",biod_names[groups],"GAM_Results",".txt",sep = ""),
  #            sep = ",", quote = FALSE, row.names = F,dec = ".")
}


# GAM plot significant_____________________________________________________________________ ####
# For Euclidean network
GAM_Sign_plots_total <- list()
ref_value <- 0
for (groups in 1:5) {
  color_groups <- as.character(CUNILLERA_cols("yellow","blue","green","red","cyan"))
  for (net in 1:5) {
    coin <- Network_results[[net]][c(length(Network_results[[net]])-54):length(Network_results[[net]])]
    dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                     biod[[groups]][,1:4])
    colnames(dataset)[1] <-c("Network")
    plots_grups <- list()
    
    p.val <- c()
    r_sqr <- c()
    dev_expl <- c()
    
    #Richness ___________________________________________________________________________________________
    # Variance Inflation Factor
    vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
    Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
    # Random Forest (RF) - selection of most informative variables for each DIV metric + group
    Model_RF = cbind(dataset[,2],Sel_vari) 
    colnames(Model_RF)[1] <- "Rich"
    RF_output <- rfsrc(Rich ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
    # Checking plots (not necessary to run every time)
    #S.rf.vimp <- gg_vimp (S.rf)
    #S.rf.vimp
    #plot (S.rf.vimp) #plot variable importance
    out_selected_vari <- var.select(RF_output,conservative = "low")
    env_selected_var <- c()
    for (u in 1:length(out_selected_vari$topvars)) {
      env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
    }
    env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
    model <- lm(dataset[,2]~env_data_model[,1:length(ncol(env_data_model))])
    summary(model)
    resid_values <- model$residuals
    dataset[,2] <- resid_values
    #GAM model
    p.val[1] <- summary.gam(gam(dataset[,2]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[8]]
    r_sqr[1] <- summary.gam(gam(dataset[,2]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[10]]
    dev_expl[1] <- summary.gam(gam(dataset[,2]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[14]]
    preds_1<- predict(gam(dataset[,2]~ s(dataset[,1], k=3, bs="cr"), method = "REML"), se.fit = TRUE)
    
    #LCBD ___________________________________________________________________________________________
    # Variance Inflation Factor
    vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
    Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
    # Random Forest (RF) - selection of most informative variables for each DIV metric + group
    Model_RF = cbind(dataset[,3],Sel_vari) 
    colnames(Model_RF)[1] <- "LCBD"
    RF_output <- rfsrc(LCBD ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
    # Checking plots (not necessary to run every time)
    #S.rf.vimp <- gg_vimp (S.rf)
    #S.rf.vimp
    #plot (S.rf.vimp) #plot variable importance
    out_selected_vari <- var.select(RF_output,conservative = "low")
    env_selected_var <- c()
    for (u in 1:length(out_selected_vari$topvars)) {
      env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
    }
    env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
    model <- lm(dataset[,3]~env_data_model[,1:length(ncol(env_data_model))])
    summary(model)
    resid_values <- model$residuals
    dataset[,3] <- resid_values
    #GAM model
    p.val[2] <-summary.gam(gam(dataset[,3]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[8]]
    r_sqr[2] <- summary.gam(gam(dataset[,3]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[10]]
    dev_expl[2] <- summary.gam(gam(dataset[,3]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[14]]
    preds_2<- predict(gam(dataset[,3]~ s(dataset[,1], k=3, bs="cr"), method = "REML"), se.fit = TRUE)
    
    #Turn ___________________________________________________________________________________________
    # Variance Inflation Factor
    vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
    Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
    # Random Forest (RF) - selection of most informative variables for each DIV metric + group
    Model_RF = cbind(dataset[,4],Sel_vari) 
    colnames(Model_RF)[1] <- "Turn"
    RF_output <- rfsrc(Turn ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
    # Checking plots (not necessary to run every time)
    #S.rf.vimp <- gg_vimp (S.rf)
    #S.rf.vimp
    #plot (S.rf.vimp) #plot variable importance
    out_selected_vari <- var.select(RF_output,conservative = "low")
    env_selected_var <- c()
    for (u in 1:length(out_selected_vari$topvars)) {
      env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
    }
    env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
    model <- lm(dataset[,4]~env_data_model[,1:length(ncol(env_data_model))])
    summary(model)
    resid_values <- model$residuals
    dataset[,4] <- resid_values
    # GAM model
    p.val[3] <-summary.gam(gam(dataset[,4]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[8]]
    r_sqr[3] <- summary.gam(gam(dataset[,4]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[10]]
    dev_expl[3] <- summary.gam(gam(dataset[,4]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[14]]
    preds_3<- predict(gam(dataset[,4]~ s(dataset[,1], k=3, bs="cr"), method = "REML"), se.fit = TRUE)
    
    #RichDiff ___________________________________________________________________________________________
    # Variance Inflation Factor
    vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
    Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
    # Random Forest (RF) - selection of most informative variables for each DIV metric + group
    Model_RF = cbind(dataset[,5],Sel_vari) 
    colnames(Model_RF)[1] <- "RichDiff"
    RF_output <- rfsrc(RichDiff ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
    # Checking plots (not necessary to run every time)
    #S.rf.vimp <- gg_vimp (S.rf)
    #S.rf.vimp
    #plot (S.rf.vimp) #plot variable importance
    out_selected_vari <- var.select(RF_output,conservative = "low")
    env_selected_var <- c()
    for (u in 1:length(out_selected_vari$topvars)) {
      env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
    }
    env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
    model <- lm(dataset[,5]~env_data_model[,1:length(ncol(env_data_model))])
    summary(model)
    resid_values <- model$residuals
    dataset[,5] <- resid_values
    # GAM model 
    p.val[4] <-summary.gam(gam(dataset[,5]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[8]]
    r_sqr[4] <- summary.gam(gam(dataset[,5]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[10]]
    dev_expl[4] <- summary.gam(gam(dataset[,5]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[14]]
    preds_4<- predict(gam(dataset[,5]~ s(dataset[,1], k=3, bs="cr"), method = "REML"), se.fit = TRUE)
    
    GAM.pred <- list(preds_1, preds_2, preds_3, preds_4)
    
    select_p.val <- which(p.val<0.05)
    
    if(length(select_p.val)>0){
      for(var in 1:length(select_p.val)){
        ref_value <- ref_value+1
        
        my_data <- data.frame(cbind(dataset[,select_p.val[var]+1],dataset[,1]),
                              mu   = GAM.pred[[select_p.val[var]]]$fit,
                              low  = GAM.pred[[select_p.val[var]]]$fit - 1.96 * GAM.pred[[select_p.val[var]]]$se.fit,
                              high = GAM.pred[[select_p.val[var]]]$fit + 1.96 * GAM.pred[[select_p.val[var]]]$se.fit)
         the_plot <-
          ggplot(my_data, aes(x = X2, y = X1)) +
          geom_jitter(alpha=0.9, shape=21, size=3, colour="black", aes(fill=X2))+
          scale_fill_continuous(type = "viridis",  alpha = 1,begin = 1,end = 0)+
          geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="black",linetype=1, size=2)+
          labs(title=paste(Names_Networks[[net]],colnames(dataset)[select_p.val[var]+1]),
               subtitle = paste("R2=", round(r_sqr[select_p.val[var]],2),
                                "ED=",round(dev_expl[select_p.val[var]],2)))+
          ylab(colnames(dataset)[select_p.val[var]+1])+xlab("Centrality-Isolation")+
          theme_classic()+
          theme(panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)),
                legend.position = "none") 
         
        GAM_Sign_plots_total[[ref_value]] <-ggdraw() +  draw_plot(the_plot)+
          draw_image(magick::image_read(image_list[[groups]]),
                     scale = 0.2,x = 0.4,y = 0.38) 
      }
    }
  }
}             
# For Fluvial network
GAM_Sign_plots_total_Fluvial <- list()
ref_value <- 0
for (groups in 1:5) {
  color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")
  coin <- Fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
  dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                   biod[[groups]][,1:4])
  colnames(dataset)[1] <-c("Network")
  plots_grups <- list()
  output_results <- list()
  p.val <- c()

  #Richness ___________________________________________________________________________________________
  # Variance Inflation Factor
  vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
  Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
  # Random Forest (RF) - selection of most informative variables for each DIV metric + group
  Model_RF = cbind(dataset[,2],Sel_vari) 
  colnames(Model_RF)[1] <- "Rich"
  RF_output <- rfsrc(Rich ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
  # Checking plots (not necessary to run every time)
  #S.rf.vimp <- gg_vimp (S.rf)
  #S.rf.vimp
  #plot (S.rf.vimp) #plot variable importance
  out_selected_vari <- var.select(RF_output,conservative = "low")
  env_selected_var <- c()
  for (u in 1:length(out_selected_vari$topvars)) {
    env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
  }
  env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
  model <- lm(dataset[,2]~env_data_model[,1:length(ncol(env_data_model))])
  summary(model)
  resid_values <- model$residuals
  dataset[,2] <- resid_values
  #GAM model
  p.val[1] <- summary.gam(gam(dataset[,2]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[8]]
  r_sqr[1] <- summary.gam(gam(dataset[,2]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[10]]
  dev_expl[1] <- summary.gam(gam(dataset[,2]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[14]]
  preds_1<- predict(gam(dataset[,2]~ s(dataset[,1], k=3, bs="cr"), method = "REML"), se.fit = TRUE)
  
  #LCBD ___________________________________________________________________________________________
  # Variance Inflation Factor
  vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
  Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
  # Random Forest (RF) - selection of most informative variables for each DIV metric + group
  Model_RF = cbind(dataset[,3],Sel_vari) 
  colnames(Model_RF)[1] <- "LCBD"
  RF_output <- rfsrc(LCBD ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
  # Checking plots (not necessary to run every time)
  #S.rf.vimp <- gg_vimp (S.rf)
  #S.rf.vimp
  #plot (S.rf.vimp) #plot variable importance
  out_selected_vari <- var.select(RF_output,conservative = "low")
  env_selected_var <- c()
  for (u in 1:length(out_selected_vari$topvars)) {
    env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
  }
  env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
  model <- lm(dataset[,3]~env_data_model[,1:length(ncol(env_data_model))])
  summary(model)
  resid_values <- model$residuals
  dataset[,3] <- resid_values
  #GAM model
  p.val[2] <-summary.gam(gam(dataset[,3]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[8]]
  r_sqr[2] <- summary.gam(gam(dataset[,3]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[10]]
  dev_expl[2] <- summary.gam(gam(dataset[,3]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[14]]
  preds_2<- predict(gam(dataset[,3]~ s(dataset[,1], k=3, bs="cr"), method = "REML"), se.fit = TRUE)
  
  #Turn ___________________________________________________________________________________________
  # Variance Inflation Factor
  vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
  Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
  # Random Forest (RF) - selection of most informative variables for each DIV metric + group
  Model_RF = cbind(dataset[,4],Sel_vari) 
  colnames(Model_RF)[1] <- "Turn"
  RF_output <- rfsrc(Turn ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
  # Checking plots (not necessary to run every time)
  #S.rf.vimp <- gg_vimp (S.rf)
  #S.rf.vimp
  #plot (S.rf.vimp) #plot variable importance
  out_selected_vari <- var.select(RF_output,conservative = "low")
  env_selected_var <- c()
  for (u in 1:length(out_selected_vari$topvars)) {
    env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
  }
  env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
  model <- lm(dataset[,4]~env_data_model[,1:length(ncol(env_data_model))])
  summary(model)
  resid_values <- model$residuals
  dataset[,4] <- resid_values
  # GAM model
  p.val[3] <-summary.gam(gam(dataset[,4]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[8]]
  r_sqr[3] <- summary.gam(gam(dataset[,4]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[10]]
  dev_expl[3] <- summary.gam(gam(dataset[,4]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[14]]
  preds_3<- predict(gam(dataset[,4]~ s(dataset[,1], k=3, bs="cr"), method = "REML"), se.fit = TRUE)
  
  #RichDiff ___________________________________________________________________________________________
  # Variance Inflation Factor
  vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
  Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
  # Random Forest (RF) - selection of most informative variables for each DIV metric + group
  Model_RF = cbind(dataset[,5],Sel_vari) 
  colnames(Model_RF)[1] <- "RichDiff"
  RF_output <- rfsrc(RichDiff ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
  # Checking plots (not necessary to run every time)
  #S.rf.vimp <- gg_vimp (S.rf)
  #S.rf.vimp
  #plot (S.rf.vimp) #plot variable importance
  out_selected_vari <- var.select(RF_output,conservative = "low")
  env_selected_var <- c()
  for (u in 1:length(out_selected_vari$topvars)) {
    env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
  }
  env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
  model <- lm(dataset[,5]~env_data_model[,1:length(ncol(env_data_model))])
  summary(model)
  resid_values <- model$residuals
  dataset[,5] <- resid_values
  # GAM model 
  p.val[4] <-summary.gam(gam(dataset[,5]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[8]]
  r_sqr[4] <- summary.gam(gam(dataset[,5]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[10]]
  dev_expl[4] <- summary.gam(gam(dataset[,5]~ s(dataset[,1], k=3, bs="cr"), method = "REML"))[[14]]
  preds_4<- predict(gam(dataset[,5]~ s(dataset[,1], k=3, bs="cr"), method = "REML"), se.fit = TRUE)
  
  GAM.pred <- list(preds_1, preds_2, preds_3, preds_4)
  
  select_p.val <- which(p.val<0.05)
  
  if(length(select_p.val)>0){
    for(var in 1:length(select_p.val)){
      ref_value <- ref_value+1
      
      my_data <- data.frame(cbind(dataset[,select_p.val[var]+1],dataset[,1]),
                            mu   = GAM.pred[[select_p.val[var]]]$fit,
                            low  = GAM.pred[[select_p.val[var]]]$fit - 1.96 * GAM.pred[[select_p.val[var]]]$se.fit,
                            high = GAM.pred[[select_p.val[var]]]$fit + 1.96 * GAM.pred[[select_p.val[var]]]$se.fit)
      
       the_plot <-ggplot(my_data, aes(x = X2, y = X1)) +
        geom_jitter(alpha=0.9, shape=21, size=3, colour="black", aes(fill=X2))+
        scale_fill_continuous(type = "viridis",  alpha = 1,begin = 1,end = 0)+
        geom_smooth(aes(ymin = low, ymax = high, y = mu), stat = "identity", colour="black",linetype=1, size=2)+
        labs(title=paste("Fluvial",colnames(dataset)[select_p.val[var]+1]),
                         subtitle = paste("R2=", round(r_sqr[select_p.val[var]],2),
                                          "ED=",round(dev_expl[select_p.val[var]],2)))+
        ylab(colnames(dataset)[select_p.val[var]+1])+xlab("Centrality-Isolation")+
        theme_classic()+
        theme(panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)),
              legend.position = "none")

      GAM_Sign_plots_total_Fluvial[[ref_value]] <-ggdraw() +  draw_plot(the_plot)+
         draw_image(magick::image_read(image_list[[groups]]),
                    scale = 0.2,x = 0.4,y = 0.38) 
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

length(plot_plot_sign_plot)
# Printing for Diverse
png(filename ="Figures/GAM_Sign_Diverse.png",
    width =629*6 ,height =850*6 ,units ="px",res = 300)
grid.arrange(plot_plot_sign_plot[[1]],plot_plot_sign_plot[[2]],plot_plot_sign_plot[[3]],
             plot_plot_sign_plot[[4]],plot_plot_sign_plot[[5]],plot_plot_sign_plot[[6]],
             plot_plot_sign_plot[[7]],plot_plot_sign_plot[[8]],plot_plot_sign_plot[[9]],
             plot_plot_sign_plot[[10]],plot_plot_sign_plot[[11]],plot_plot_sign_plot[[12]],
             plot_plot_sign_plot[[13]],plot_plot_sign_plot[[14]],
             ncol=3,nrow=5, top="Community Diversity (alpha & beta)")
dev.off()


# NMDS plots_______________________________________________________________________________ ####
biod_names <- c("S16","S18","phy","zoo", "zoo.18S")
net_names <- c("650 km", "325 km","100 km","65 km","6.5 km", "Fluvial")
color_groups <- as.character(CUNILLERA_cols("yellow","blue","green","red","cyan"))

# For Euclidean network
plots_NMDS <- list()
plots_NMDS_total <- list()

plots_NMDS_model_result <- list()
plots_NMDS_total_model_result <- list()

for (groups in 1:5) {
  for (net in 1:5) {
coin <- Network_results[[net]][c(length(Network_results[[net]])-54):length(Network_results[[net]])]
centr_iso <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]])

dbRDA <- capscale(comm_data[[groups]]~
                  env_data[[groups]][,4]+env_data[[groups]][,5]+env_data[[groups]][,6]+env_data[[groups]][,7]+env_data[[groups]][,8],
                  distance = "jaccard",add = "lingoes")

x<- dbRDA$CA$u[,1]
y <- dbRDA$CA$u[,2]

dataset <- data.frame(x,y,centr_iso)
NMDS_model <- ordisurf(dbRDA$CA$u ~ dataset[,3], plot = F,)
NMDS_model_results <- summary(NMDS_model)

plots_NMDS_model_result[[net]] <- NMDS_model_results

contour.vals <- extract.xyz(obj = NMDS_model)

    the_plot <-ggplot(dataset, aes(x=x,y=y))+
                  geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+
                  geom_jitter(shape=21, size=5, alpha=0.8, aes(fill=centr_iso))+
                  geom_contour(data=contour.vals, aes(x, y, z = z, colour = ..level..))+
                  #geom_path(data = df_ellipse, aes(x=x, y=y, colour=Group), size=2, show.legend = FALSE)+
                  scale_colour_viridis(alpha = 1,begin = 1,end = 0)+
                  #manual(values = viridis_pal(0.9,1,0,direction = -1)(length(unique(df_ellipse$Group))))+
                  scale_fill_viridis(alpha = 1,begin = 1,end = 0)+
                  labs(title = paste(biod_names[groups], net_names[net]),
                       subtitle = paste("Unc. Inert.=",round(dbRDA$CA$tot.chi*1/dbRDA$tot.chi,2), 
                                        "R2=", round(NMDS_model_results$r.sq,2),
                                        "ED=",round(NMDS_model_results$dev.expl,2)))+
                  xlab("dbRDA1")+ylab("dbRDA2")+
                  theme_classic()+
                  theme(legend.position = "none",
                        panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)))
    
    plots_NMDS[[net]] <- ggdraw() +  draw_plot(the_plot)+
      draw_image(magick::image_read(image_list[[groups]]),
                 scale = 0.2,x = 0.4,y = 0.38) 
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
    coin <- Fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
    centr_iso <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]])

    dbRDA <- capscale(comm_data[[groups]]~
                      env_data[[groups]][,4]+env_data[[groups]][,5]+env_data[[groups]][,6]+env_data[[groups]][,7]+env_data[[groups]][,8],
                      distance = "jaccard",add = "lingoes")
    
    x<- dbRDA$CA$u[,1]
    y <- dbRDA$CA$u[,2]
    
    dataset <- data.frame(x,y,centr_iso)
    NMDS_model <- ordisurf(dbRDA$CA$u ~ dataset[,3], plot = F)
    NMDS_model_results <- summary(NMDS_model)
    
    plots_NMDS_fluvial_model_result[[groups]] <- NMDS_model_results
    
    contour.vals <- extract.xyz(obj = NMDS_model)
    
    the_plot<- ggplot(dataset, aes(x=x,y=y))+
      geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+
      geom_jitter(shape=21, size=5, alpha=0.8, aes(fill=centr_iso))+
      geom_contour(data=contour.vals, aes(x, y, z = z, colour = ..level..))+
      scale_colour_viridis(alpha = 1,begin = 1,end = 0)+
      scale_fill_viridis(alpha = 1,begin =1 ,end = 0)+
      labs(title = paste(biod_names[groups], net_names[net]),
           subtitle = paste("Unc. Inert.=",round(dbRDA$CA$tot.chi*1/dbRDA$tot.chi,2), 
                            "R2=", round(NMDS_model_results$r.sq,2),
                            "ED=",round(NMDS_model_results$dev.expl,2)))+
      xlab("dbRDA1")+ylab("dbRDA2")+
      theme_classic()+
      theme(legend.position = "none",
            panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)))
    
    plots_NMDS_fluvial[[groups]]<-ggdraw() +  draw_plot(the_plot)+
      draw_image(magick::image_read(image_list[[groups]]),
                 scale = 0.2,x = 0.4,y = 0.38) 
    
  plots_NMDS_total_fluvial[[groups]] <- plots_NMDS_fluvial
  plots_NMDS_fluvial_total_model_result[[groups]] <- plots_NMDS_fluvial_model_result
}


# Print NMDS
png(filename ="Figures/NMDS_Diverse.png",
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
             
             ncol=6,nrow=5, top="NMDS")
dev.off()


# NMDS plots significant Ordisurfs_________________________________________________________ ####

# For Euclidean network
plots_NMDS_sign <- list()
ref_value <- 0
for (groups in 1:5) {
  for (net in 1:5) {
    coin <- Network_results[[net]][c(length(Network_results[[net]])-54):length(Network_results[[net]])]
    centr_iso <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]])
    
    dbRDA <- capscale(comm_data[[groups]]~
                      env_data[[groups]][,4]+env_data[[groups]][,5]+env_data[[groups]][,6]+env_data[[groups]][,7]+env_data[[groups]][,8],
                      distance = "jaccard",add = "lingoes")
    
    x<- dbRDA$CA$u[,1]
    y <- dbRDA$CA$u[,2]
    
    dataset <- data.frame(x,y,centr_iso)
    NMDS_model <- ordisurf(dbRDA$CA$u ~ dataset[,3], plot = F)
    NMDS_model_results <- summary(NMDS_model)
    
    plots_NMDS_model_pval <- NMDS_model_results$s.pv
    
    if (plots_NMDS_model_pval<0.05) {
      
    ref_value <- ref_value+1
    contour.vals <- extract.xyz(obj = NMDS_model)
    
    the_plot<- ggplot(dataset, aes(x=x,y=y))+
      geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+
      geom_jitter(shape=21, size=5, alpha=0.8, aes(fill=centr_iso))+
      geom_contour(data=contour.vals, aes(x, y, z = z, colour = ..level..))+
      scale_colour_viridis(alpha = 1,begin = 1,end = 0)+
      scale_fill_viridis(alpha = 1,begin = 1,end = 0)+
      labs(title = paste(biod_names[groups], net_names[net]),
           subtitle = paste("Unc. Inert.=",round(dbRDA$CA$tot.chi*1/dbRDA$tot.chi,2), 
                            "R2=", round(NMDS_model_results$r.sq,2),
                            "ED=",round(NMDS_model_results$dev.expl,2)))+
      xlab("dbRDA1")+ylab("dbRDA2")+
      theme_classic()+
      theme(legend.position = "none",
            panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)))
    
    plots_NMDS_sign[[ref_value]] <- ggdraw() +  draw_plot(the_plot)+
      draw_image(magick::image_read(image_list[[groups]]),
                 scale = 0.2,x = 0.4,y = 0.38) 
    
    }
  }
  plots_NMDS_sign
}

# For Fluvial network
plots_NMDS_fluvial_sign<- list()
ref_value <- 0
for (groups in 1:5) {
  coin <- Fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
  centr_iso <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]])
  
  dbRDA <- capscale(comm_data[[groups]]~
                    env_data[[groups]][,4]+env_data[[groups]][,5]+env_data[[groups]][,6]+env_data[[groups]][,7]+env_data[[groups]][,8],
                    distance = "jaccard",add = "lingoes")
  
  x<- dbRDA$CA$u[,1]
  y <- dbRDA$CA$u[,2]
  
  dataset <- data.frame(x,y,centr_iso)
  NMDS_model <- ordisurf(dbRDA$CA$u ~ dataset[,3], plot = F)
  NMDS_model_results <- summary(NMDS_model)
  
  plots_NMDS_model_pval <- NMDS_model_results$s.pv
  
  if (plots_NMDS_model_pval<0.05) {
    
    ref_value <- ref_value+1
    contour.vals <- extract.xyz(obj = NMDS_model)
  
    the_plot <- ggplot(dataset, aes(x=x,y=y))+
    geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+
    geom_jitter(shape=21, size=5, alpha=0.8, aes(fill=centr_iso))+
    geom_contour(data=contour.vals, aes(x, y, z = z, colour = ..level..))+
    scale_colour_viridis(alpha = 1,begin = 1,end = 0)+
    #manual(values = viridis_pal(0.9,1,0,direction = -1)(length(unique(df_ellipse$Group))))+
    scale_fill_viridis(alpha = 1,begin = 1,end = 0)+
    labs(title = paste(biod_names[groups], "Fluvial"),
         subtitle = paste("Unc. Inert.=",round(dbRDA$CA$tot.chi*1/dbRDA$tot.chi,2), 
                          "R2=", round(NMDS_model_results$r.sq,2),
                          "ED=",round(NMDS_model_results$dev.expl,2)))+
    xlab("dbRDA1")+ylab("dbRDA2")+
    theme_classic()+
    theme(legend.position = "none",
          panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)))
    
    plots_NMDS_fluvial_sign[[ref_value]] <-  ggdraw() +  draw_plot(the_plot)+
      draw_image(magick::image_read(image_list[[groups]]),
                 scale = 0.2,x = 0.4,y = 0.38) 
  }
plots_NMDS_fluvial_sign
}

png(filename ="Figures/NMDS_Diverse_Sign.png",
    width =629*8 ,height =629*5 ,units ="px",res = 250)
grid.arrange(plots_NMDS_sign[[1]],plots_NMDS_sign[[2]],plots_NMDS_sign[[3]],plots_NMDS_sign[[4]],
             plots_NMDS_sign[[5]],plots_NMDS_sign[[6]],plots_NMDS_sign[[7]],
             plots_NMDS_sign[[8]],plots_NMDS_sign[[9]],plots_NMDS_sign[[10]],
             plots_NMDS_sign[[11]],
             ncol=5,nrow=3, top="Community Composition (NMDS SS)")
dev.off()

# GAM NMDS models result in table format - Supplementary like______________________________ ####
biod_names <- c("S16","S18","phy","zoo", "zoo.18S")
Names_Networks <- c("650 km", "325 km","100 km","65 km","6.5 km", "Fluvial")

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
                             plots_NMDS_total_model_result[[groups]][[net]]$p.coeff,
                             plots_NMDS_total_model_result[[groups]][[net]]$se[1],
                             plots_NMDS_total_model_result[[groups]][[net]]$p.t,
                             plots_NMDS_total_model_result[[groups]][[net]]$p.pv,
                             plots_NMDS_total_model_result[[groups]][[net]]$chi.sq,
                             plots_NMDS_total_model_result[[groups]][[net]]$s.pv,
                             plots_NMDS_total_model_result[[groups]][[net]]$r.sq,
                             plots_NMDS_total_model_result[[groups]][[net]]$dev.expl) 
  }
  table_groups[[groups]] <- as.data.frame(a)
  write.csv2(table_groups[[groups]], file = paste("Results_Tables/",biod_names[groups],"GAM_NMDS_Results",".csv",sep = ""))
  #write.table(table_groups[[groups]], file = paste("Results_Tables/",biod_names[groups],"GAM_NMDS_Results",".txt",sep = ""), 
  #            sep = ",", quote = FALSE, row.names = F, dec = ".")
}


# Summary plot GAM models__________________________________________________________________ ####
sign_netw <- list()
sign_groups <- list()

for (group in 1:5) {
  for (netw in 1:5) {
    sign <- c()
    sign[1] <-GAMmodel_resutls_total[[group]][[netw]][[1]][[8]]
    sign[2] <-GAMmodel_resutls_total[[group]][[netw]][[2]][[8]]
    sign[3] <-GAMmodel_resutls_total[[group]][[netw]][[3]][[8]]
    sign[4] <-GAMmodel_resutls_total[[group]][[netw]][[4]][[8]]
    sign[5] <-plots_NMDS_total_model_result[[group]][[netw]][[8]]
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
  sign[5] <-plots_NMDS_fluvial_total_model_result[[5]][[group]][[8]]
  flu_sign_groups[[group]] <- sign
}


Names_Variab <- c("Species richness", "LCBD", "Replacement", "Richness difference","dbRDA" )
plots_significance <- list()
for (variable in 1:5) {
  
  max_netw <- cbind(c(sign_groups[[1]][[1]][[variable]], sign_groups[[2]][[1]][[variable]], sign_groups[[3]][[1]][[variable]], 
                      sign_groups[[4]][[1]][[variable]], sign_groups[[5]][[1]][[variable]]),
                    rep("650 km", 5), 
                    c("S16","S18","Phy","Zoo", "S18zoo"))
  
  mid_netw <- cbind(c(sign_groups[[1]][[2]][[variable]], sign_groups[[2]][[2]][[variable]], sign_groups[[3]][[2]][[variable]], 
                      sign_groups[[4]][[2]][[variable]], sign_groups[[5]][[2]][[variable]]),
                    rep("325 km", 5), 
                    c("S16","S18","Phy","Zoo", "S18zoo"))
  
  mid_mid_netw <- cbind(c(sign_groups[[1]][[3]][[variable]], sign_groups[[2]][[3]][[variable]], sign_groups[[3]][[3]][[variable]], 
                          sign_groups[[4]][[3]][[variable]], sign_groups[[5]][[3]][[variable]]),
                        rep("100 km", 5), 
                        c("S16","S18","Phy","Zoo", "S18zoo"))
  
  small_netw <- cbind(c(sign_groups[[1]][[4]][[variable]], sign_groups[[2]][[4]][[variable]], sign_groups[[3]][[4]][[variable]], 
                        sign_groups[[4]][[4]][[variable]], sign_groups[[5]][[4]][[variable]]),
                      rep("65 km", 5), 
                      c("S16","S18","Phy","Zoo", "S18zoo"))
  
  min_netw <- cbind(c(sign_groups[[1]][[5]][[variable]], sign_groups[[2]][[5]][[variable]], sign_groups[[3]][[5]][[variable]], 
                      sign_groups[[4]][[5]][[variable]], sign_groups[[5]][[5]][[variable]]),
                    rep("6.5 km", 5), 
                    c("S16","S18","Phy","Zoo", "S18zoo"))
  
  fluv_netw <- cbind(c(flu_sign_groups[[1]][[variable]], flu_sign_groups[[2]][[variable]], flu_sign_groups[[3]][[variable]], 
                       flu_sign_groups[[4]][[variable]], flu_sign_groups[[5]][[variable]]),
                     rep("Fluvial", 5), 
                     c("S16","S18","Phy","Zoo", "S18zoo"))
  
  
  dataset_pval <- as.data.frame(rbind(max_netw,mid_netw,mid_mid_netw,small_netw,min_netw,fluv_netw))
  colnames(dataset_pval) <- c("pval","Network","Group")
  dataset_pval$pval <-as.numeric(dataset_pval$pval)
  dataset_pval$Network <- factor(dataset_pval$Network,
                                 levels = c("650 km", "325 km","100 km","65 km","6.5 km", "Fluvial"))
  dataset_pval$Group <- factor(dataset_pval$Group,
                               levels = c("S16","S18","Phy","Zoo", "S18zoo"))
  significants <- rep("NoSign",5*6)
  significants[which(dataset_pval$pval<0.05)]<- "Sign"
  dataset_pval$Sign <- factor(significants) 
  
  color_groups <- as.character(CUNILLERA_cols("yellow","blue","green","red","cyan"))
  
  plots_significance[[variable]] <-  ggplot(dataset_pval, aes(x=Network, y=as.numeric(pval)))+
    geom_abline(slope = 0,intercept = 0.05, colour="black", linetype=2,size=1)+
    geom_jitter(aes(fill=Group, alpha=Sign, size=Sign),shape=21,width = 0.5)+
    scale_x_discrete(limits=c("6.5 km","65 km", "100 km","325 km","650 km", "Fluvial"))+
    scale_alpha_manual(values = c(0.3,0.9))+
    scale_size_manual(values = c(2,7))+
    scale_fill_manual(values=c(color_groups[1],color_groups[2],color_groups[3],
                               color_groups[4],color_groups[5],color_groups[6]))+
    scale_y_continuous(expand = c(0.2,0.01),
                       breaks =c(0.2,0.4,0.6,0.8,1) )+
    geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5), size=1, colour="grey70")+
    labs(title="")+ylab(Names_Variab[variable])+xlab("")+
    theme_classic()+theme(legend.position = "none")
}

# Summary plots NMDS _________________________________________________________ ####

# For Euclidean network
plots_NMDS_sign <- list()
ref_value <- 0
TOTAL_pval_summary_NMDS <- list()
pval_summary_NMDS <- list()
for (groups in 1:5) {
  for (net in 1:5) {
    coin <- Network_results[[net]][c(length(Network_results[[net]])-54):length(Network_results[[net]])]
    centr_iso <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]])
    
    dbRDA <- capscale(comm_data[[groups]]~
                        env_data[[groups]][,4]+env_data[[groups]][,5]+env_data[[groups]][,6]+env_data[[groups]][,7]+env_data[[groups]][,8],
                      distance = "jaccard",add = "lingoes")
    
    x<- dbRDA$CA$u[,1]
    y <- dbRDA$CA$u[,2]
    
    dataset <- data.frame(x,y,centr_iso)
    NMDS_model <- ordisurf(dbRDA$CA$u ~ dataset[,3], plot = F)
    NMDS_model_results <- summary(NMDS_model)
    pval_summary_NMDS[[net]] <- NMDS_model_results$s.pv
  }
  TOTAL_pval_summary_NMDS[[groups]] <- pval_summary_NMDS
}


# For Fluvial network
plots_NMDS_fluvial_sign<- list()
ref_value <- 0
TOTAL_FLU_pval_summary_NMDS <- list()
for (groups in 1:5) {
  coin <- Fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
  centr_iso <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]])
  
  dbRDA <- capscale(comm_data[[groups]]~
                      env_data[[groups]][,4]+env_data[[groups]][,5]+env_data[[groups]][,6]+env_data[[groups]][,7]+env_data[[groups]][,8],
                    distance = "jaccard",add = "lingoes")
  
  x<- dbRDA$CA$u[,1]
  y <- dbRDA$CA$u[,2]
  
  dataset <- data.frame(x,y,centr_iso)
  NMDS_model <- ordisurf(dbRDA$CA$u ~ dataset[,3], plot = F)
  NMDS_model_results <- summary(NMDS_model)
  TOTAL_FLU_pval_summary_NMDS[[groups]] <- NMDS_model_results$s.pv
}

Group_names <- c("S16","S18","Phy","Zoo", "S18zoo")
output <- data.frame()
for (group in 1:length(Group_names)) {
  out <- data.frame("pval"=c(unlist(TOTAL_pval_summary_NMDS[[group]]),unlist(TOTAL_FLU_pval_summary_NMDS[[group]])),
                    "Network"=c("650 km", "325 km","100 km","65 km","6.5 km", "Fluvial"),
                    "Group"=rep(Group_names[group],6)) %>% 
    mutate(Sign=ifelse(pval<0.05,"Sign","NoSign"))
  output <- bind_rows(output,out)
}
dataset_pval <- output
dataset_pval$Network <- factor(dataset_pval$Network, levels = c("650 km", "325 km","100 km","65 km","6.5 km", "Fluvial"))
dataset_pval$Group <- factor(dataset_pval$Group,levels = c("S16","S18","Phy","Zoo", "S18zoo"))

color_groups <- as.character(CUNILLERA_cols("yellow","blue","green","red","cyan"))

plots_significance[[5]] <-  ggplot(dataset_pval, aes(x=Network, y=as.numeric(pval)))+
  geom_abline(slope = 0,intercept = 0.05, colour="black", linetype=2,size=1)+
  geom_jitter(aes(fill=Group, alpha=Sign, size=Sign),shape=21,width = 0.5)+
  scale_x_discrete(limits=c("6.5 km","65 km", "100 km","325 km","650 km", "Fluvial"))+
  scale_alpha_manual(values = c(0.3,0.9))+
  scale_size_manual(values = c(2,7))+
  scale_fill_manual(values=c(color_groups[1],color_groups[2],color_groups[3],
                             color_groups[4],color_groups[5],color_groups[6]))+
  scale_y_continuous(expand = c(0.2,0.01),
                     breaks =c(0.2,0.4,0.6,0.8,1) )+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5), size=1, colour="grey70")+
  labs(title="")+ylab("dbRDA ordisurf")+xlab("")+
  theme_classic()+theme(legend.position = "bottom")


legend_try <- cowplot::get_legend(plots_significance[[5]])


png(filename =paste("Figures/All_Significance_Groups.png"),
    width =800*3 ,height =750*5 ,units ="px",res = 400)
grid.arrange(plots_significance[[1]],
             plots_significance[[3]],plots_significance[[4]],plots_significance[[2]],
             plots_significance[[5]]+theme(legend.position = "none"),
             legend_try,
             ncol=1,nrow=6, top="Significance values")
dev.off()      


# Summary plot CORRELATION ________________________________________________________________ ####

# For Euclidean network
p.val_netw <- list()
r.value_netw <- list()

p.val_groups <- list()
r.value_groups <- list()

for (groups in 1:5) {
  color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")
  for (net in 1:5) {
    coin <- Network_results[[net]][c(length(Network_results[[net]])-54):length(Network_results[[net]])]
    dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                     biod[[groups]][,1:4])
    colnames(dataset)[1] <-c("Network")
    
    p.val <- c()
    r.value <- c()
    
    #Richness ___________________________________________________________________________________________
    # Variance Inflation Factor
    vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
    Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
    # Random Forest (RF) - selection of most informative variables for each DIV metric + group
    Model_RF = cbind(dataset[,2],Sel_vari) 
    colnames(Model_RF)[1] <- "Rich"
    RF_output <- rfsrc(Rich ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
    # Checking plots (not necessary to run every time)
    #S.rf.vimp <- gg_vimp (S.rf)
    #S.rf.vimp
    #plot (S.rf.vimp) #plot variable importance
    out_selected_vari <- var.select(RF_output,conservative = "low")
    env_selected_var <- c()
    for (u in 1:length(out_selected_vari$topvars)) {
      env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
    }
    env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
    model <- lm(dataset[,2]~env_data_model[,1:length(ncol(env_data_model))])
    summary(model)
    resid_values <- model$residuals
    dataset[,2] <- resid_values
    p.val[1] <-cor.test(dataset[,2],dataset[,1],method = "pearson")[3]
    r.value[1] <-cor.test(dataset[,2],dataset[,1],method = "pearson")[4]
    
    #LCBD ___________________________________________________________________________________________
    # Variance Inflation Factor
    vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
    Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
    # Random Forest (RF) - selection of most informative variables for each DIV metric + group
    Model_RF = cbind(dataset[,3],Sel_vari) 
    colnames(Model_RF)[1] <- "LCBD"
    RF_output <- rfsrc(LCBD ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
    # Checking plots (not necessary to run every time)
    #S.rf.vimp <- gg_vimp (S.rf)
    #S.rf.vimp
    #plot (S.rf.vimp) #plot variable importance
    out_selected_vari <- var.select(RF_output,conservative = "low")
    env_selected_var <- c()
    for (u in 1:length(out_selected_vari$topvars)) {
      env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
    }
    env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
    model <- lm(dataset[,3]~env_data_model[,1:length(ncol(env_data_model))])
    summary(model)
    resid_values <- model$residuals
    dataset[,3] <- resid_values
    p.val[2] <-cor.test(dataset[,3],dataset[,1],method = "pearson")[3]
    r.value[2] <-cor.test(dataset[,3],dataset[,1],method = "pearson")[4]
    
    #Turn ___________________________________________________________________________________________
    # Variance Inflation Factor
    vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
    Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
    # Random Forest (RF) - selection of most informative variables for each DIV metric + group
    Model_RF = cbind(dataset[,4],Sel_vari) 
    colnames(Model_RF)[1] <- "Turn"
    RF_output <- rfsrc(Turn ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
    # Checking plots (not necessary to run every time)
    #S.rf.vimp <- gg_vimp (S.rf)
    #S.rf.vimp
    #plot (S.rf.vimp) #plot variable importance
    out_selected_vari <- var.select(RF_output,conservative = "low")
    env_selected_var <- c()
    for (u in 1:length(out_selected_vari$topvars)) {
      env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
    }
    env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
    model <- lm(dataset[,4]~env_data_model[,1:length(ncol(env_data_model))])
    summary(model)
    resid_values <- model$residuals
    dataset[,4] <- resid_values
    p.val[3] <-cor.test(dataset[,4],dataset[,1],method = "pearson")[3]
    r.value[3] <-cor.test(dataset[,4],dataset[,1],method = "pearson")[4]
    
    #RichDiff ___________________________________________________________________________________________
    # Variance Inflation Factor
    vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
    Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
    # Random Forest (RF) - selection of most informative variables for each DIV metric + group
    Model_RF = cbind(dataset[,5],Sel_vari) 
    colnames(Model_RF)[1] <- "RichDiff"
    RF_output <- rfsrc(RichDiff ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
    # Checking plots (not necessary to run every time)
    #S.rf.vimp <- gg_vimp (S.rf)
    #S.rf.vimp
    #plot (S.rf.vimp) #plot variable importance
    out_selected_vari <- var.select(RF_output,conservative = "low")
    env_selected_var <- c()
    for (u in 1:length(out_selected_vari$topvars)) {
      env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
    }
    env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
    model <- lm(dataset[,5]~env_data_model[,1:length(ncol(env_data_model))])
    summary(model)
    resid_values <- model$residuals
    dataset[,5] <- resid_values
    p.val[4] <-cor.test(dataset[,5],dataset[,1],method = "pearson")[3]
    r.value[4] <-cor.test(dataset[,5],dataset[,1],method = "pearson")[4]
    
    p.val_netw[[net]] <- p.val
    r.value_netw[[net]] <- r.value
  }
  p.val_groups[[groups]] <- p.val_netw
  r.value_groups[[groups]] <- r.value_netw
}

p.val_groups_flu <- list()
r.value_groups_flu <- list()
for (groups in 1:5) {
  color_groups <- as.character(CUNILLERA_cols("yellow","blue","green","red","cyan"))
  coin <- Fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
  dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                   biod[[groups]][,1:4])
  colnames(dataset)[1] <-c("Network")
  
  p.val <- c()
  r.value <- c()

  #Richness ___________________________________________________________________________________________
  # Variance Inflation Factor
  vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
  Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
  # Random Forest (RF) - selection of most informative variables for each DIV metric + group
  Model_RF = cbind(dataset[,2],Sel_vari) 
  colnames(Model_RF)[1] <- "Rich"
  RF_output <- rfsrc(Rich ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
  # Checking plots (not necessary to run every time)
  #S.rf.vimp <- gg_vimp (S.rf)
  #S.rf.vimp
  #plot (S.rf.vimp) #plot variable importance
  out_selected_vari <- var.select(RF_output,conservative = "low")
  env_selected_var <- c()
  for (u in 1:length(out_selected_vari$topvars)) {
    env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
  }
  env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
  model <- lm(dataset[,3]~env_data_model[,1:length(ncol(env_data_model))])
  summary(model)
  resid_values <- model$residuals
  dataset[,2] <- resid_values
  p.val[1] <-cor.test(dataset[,2],dataset[,1],method = "pearson")[3]
  r.value[1] <-cor.test(dataset[,2],dataset[,1],method = "pearson")[4]
  
  #LCBD ___________________________________________________________________________________________
  # Variance Inflation Factor
  vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
  Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
  # Random Forest (RF) - selection of most informative variables for each DIV metric + group
  Model_RF = cbind(dataset[,3],Sel_vari) 
  colnames(Model_RF)[1] <- "LCBD"
  RF_output <- rfsrc(LCBD ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
  # Checking plots (not necessary to run every time)
  #S.rf.vimp <- gg_vimp (S.rf)
  #S.rf.vimp
  #plot (S.rf.vimp) #plot variable importance
  out_selected_vari <- var.select(RF_output,conservative = "low")
  env_selected_var <- c()
  for (u in 1:length(out_selected_vari$topvars)) {
    env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
  }
  env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
  model <- lm(dataset[,3]~env_data_model[,1:length(ncol(env_data_model))])
  summary(model)
  resid_values <- model$residuals
  dataset[,3] <- resid_values
  p.val[2] <-cor.test(dataset[,3],dataset[,1],method = "pearson")[3]
  r.value[2] <-cor.test(dataset[,3],dataset[,1],method = "pearson")[4]
  
  #Turn ___________________________________________________________________________________________
  # Variance Inflation Factor
  vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
  Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
  # Random Forest (RF) - selection of most informative variables for each DIV metric + group
  Model_RF = cbind(dataset[,4],Sel_vari) 
  colnames(Model_RF)[1] <- "Turn"
  RF_output <- rfsrc(Turn ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
  # Checking plots (not necessary to run every time)
  #S.rf.vimp <- gg_vimp (S.rf)
  #S.rf.vimp
  #plot (S.rf.vimp) #plot variable importance
  out_selected_vari <- var.select(RF_output,conservative = "low")
  env_selected_var <- c()
  for (u in 1:length(out_selected_vari$topvars)) {
    env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
  }
  env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
  model <- lm(dataset[,4]~env_data_model[,1:length(ncol(env_data_model))])
  summary(model)
  resid_values <- model$residuals
  dataset[,4] <- resid_values
  p.val[3] <-cor.test(dataset[,4],dataset[,1],method = "pearson")[3]
  r.value[3] <-cor.test(dataset[,4],dataset[,1],method = "pearson")[4]
  
  #RichDiff ___________________________________________________________________________________________
  # Variance Inflation Factor
  vif_value <- vifstep (env_data[[groups]][,4:8], th=7) # threshold set to VIF<7
  Sel_vari <- exclude(env_data[[groups]][,4:8], vif_value)
  # Random Forest (RF) - selection of most informative variables for each DIV metric + group
  Model_RF = cbind(dataset[,5],Sel_vari) 
  colnames(Model_RF)[1] <- "RichDiff"
  RF_output <- rfsrc(RichDiff ~ ., mtry = 6, ntree = 1000, importance="permute", 	data = Model_RF) 
  # Checking plots (not necessary to run every time)
  #S.rf.vimp <- gg_vimp (S.rf)
  #S.rf.vimp
  #plot (S.rf.vimp) #plot variable importance
  out_selected_vari <- var.select(RF_output,conservative = "low")
  env_selected_var <- c()
  for (u in 1:length(out_selected_vari$topvars)) {
    env_selected_var[u] <- which(colnames(env_data[[groups]][,4:8])==out_selected_vari$topvars[u])  
  }
  env_data_model <- as.matrix(env_data[[groups]][,c(3+env_selected_var)])
  model <- lm(dataset[,5]~env_data_model[,1:length(ncol(env_data_model))])
  summary(model)
  resid_values <- model$residuals
  dataset[,5] <- resid_values
  p.val[4] <-cor.test(dataset[,5],dataset[,1],method = "pearson")[3]
  r.value[4] <-cor.test(dataset[,5],dataset[,1],method = "pearson")[4]
  
  p.val_groups_flu[[groups]] <- p.val
  r.value_groups_flu[[groups]] <- r.value
}
  
Names_Variab <- c("Species richness", "LCBD", "Replacement", "Richness difference","NMDS SS" )
plots_significance <- list()
for (variable in 1:4) {
  #R-squared
  max_netw <- cbind(c(r.value_groups[[1]][[1]][[variable]], r.value_groups[[2]][[1]][[variable]], r.value_groups[[3]][[1]][[variable]], 
                      r.value_groups[[4]][[1]][[variable]],r.value_groups[[5]][[1]][[variable]]),
                    rep("650 km", 5), 
                    c("S16","S18","Phy","Zoo", "S18zoo"))
  
  mid_netw <- cbind(c(r.value_groups[[1]][[2]][[variable]], r.value_groups[[2]][[2]][[variable]], r.value_groups[[3]][[2]][[variable]], 
                      r.value_groups[[4]][[2]][[variable]],r.value_groups[[5]][[2]][[variable]]),
                    rep("325 km", 5), 
                    c("S16","S18","Phy","Zoo", "S18zoo"))
  
  mid_mid_netw <- cbind(c(r.value_groups[[1]][[3]][[variable]], r.value_groups[[2]][[3]][[variable]], r.value_groups[[3]][[3]][[variable]], 
                          r.value_groups[[4]][[3]][[variable]],r.value_groups[[5]][[3]][[variable]]),
                        rep("100 km", 5), 
                        c("S16","S18","Phy","Zoo", "S18zoo"))
  
  small_netw <- cbind(c(r.value_groups[[1]][[4]][[variable]], r.value_groups[[2]][[4]][[variable]], r.value_groups[[3]][[4]][[variable]], 
                        r.value_groups[[4]][[4]][[variable]],r.value_groups[[5]][[4]][[variable]]),
                      rep("65 km", 5), 
                      c("S16","S18","Phy","Zoo", "S18zoo"))
  
  min_netw <- cbind(c(r.value_groups[[1]][[5]][[variable]], r.value_groups[[2]][[5]][[variable]], r.value_groups[[3]][[5]][[variable]], 
                      r.value_groups[[4]][[5]][[variable]],r.value_groups[[5]][[5]][[variable]]),
                    rep("6.5 km", 5), 
                    c("S16","S18","Phy","Zoo", "S18zoo"))
  
  fluv_netw <- cbind(c(r.value_groups_flu[[1]][[variable]], r.value_groups_flu[[2]][[variable]], r.value_groups_flu[[3]][[variable]], 
                       r.value_groups_flu[[4]][[variable]],r.value_groups_flu[[5]][[variable]]),
                     rep("Fluvial", 5), 
                     c("S16","S18","Phy","Zoo", "S18zoo"))
  
  # P-value
  max_netw_Pval <- cbind(c(p.val_groups[[1]][[1]][[variable]], p.val_groups[[2]][[1]][[variable]], p.val_groups[[3]][[1]][[variable]], 
                      p.val_groups[[4]][[1]][[variable]],p.val_groups[[5]][[1]][[variable]]),
                    rep("650 km", 5), 
                    c("S16","S18","Phy","Zoo", "S18zoo"))
  
  mid_netw_Pval <- cbind(c(p.val_groups[[1]][[2]][[variable]], p.val_groups[[2]][[2]][[variable]], p.val_groups[[3]][[2]][[variable]], 
                      p.val_groups[[4]][[2]][[variable]],p.val_groups[[5]][[2]][[variable]]),
                    rep("325 km", 5), 
                    c("S16","S18","Phy","Zoo", "S18zoo"))
  
  mid_mid_netw_Pval <- cbind(c(p.val_groups[[1]][[3]][[variable]], p.val_groups[[2]][[3]][[variable]], p.val_groups[[3]][[3]][[variable]], 
                          p.val_groups[[4]][[3]][[variable]],p.val_groups[[5]][[3]][[variable]]),
                        rep("100 km", 5), 
                        c("S16","S18","Phy","Zoo", "S18zoo"))
  
  small_netw_Pval <- cbind(c(p.val_groups[[1]][[4]][[variable]], p.val_groups[[2]][[4]][[variable]], p.val_groups[[3]][[4]][[variable]], 
                        p.val_groups[[4]][[4]][[variable]],p.val_groups[[5]][[4]][[variable]]),
                      rep("65 km", 5), 
                      c("S16","S18","Phy","Zoo", "S18zoo"))
  
  min_netw_Pval <- cbind(c(p.val_groups[[1]][[5]][[variable]], p.val_groups[[2]][[5]][[variable]], p.val_groups[[3]][[5]][[variable]], 
                      p.val_groups[[4]][[5]][[variable]],p.val_groups[[5]][[5]][[variable]]),
                    rep("6.5 km", 5), 
                    c("S16","S18","Phy","Zoo", "S18zoo"))
  
  fluv_netw_Pval <- cbind(c(p.val_groups_flu[[1]][[variable]], p.val_groups_flu[[2]][[variable]], p.val_groups_flu[[3]][[variable]], 
                            p.val_groups_flu[[4]][[variable]],p.val_groups_flu[[5]][[variable]]),
                     rep("Fluvial", 5), 
                     c("S16","S18","Phy","Zoo", "S18zoo"))
  dataset_pval <- as.data.frame(cbind(rbind(max_netw,mid_netw,mid_mid_netw,small_netw,min_netw,fluv_netw),
                                      rbind(max_netw_Pval,mid_netw_Pval,mid_mid_netw_Pval,small_netw_Pval,min_netw_Pval,fluv_netw_Pval)))
  dataset_pval <- dataset_pval[,c(1,4,2,3)]
  colnames(dataset_pval) <- c("Rsqr","pval","Network","Group")
  dataset_pval$Rsqr <-as.numeric(dataset_pval$Rsqr)
  dataset_pval$pval <-as.numeric(dataset_pval$pval)
  dataset_pval$Network <- factor(dataset_pval$Network,
                                 levels = c("650 km", "325 km","100 km","65 km","6.5 km", "Fluvial"))
  dataset_pval$Group <- factor(dataset_pval$Group,
                               levels = c("S16","S18","Phy","Zoo", "S18zoo"))
  significants <- rep("NoSign",5*6)
  significants[which(dataset_pval$pval<0.05)]<- "Sign"
  dataset_pval$Sign <- factor(significants) 
  
  color_groups <- as.character(CUNILLERA_cols("yellow","blue","green","red","cyan"))
  
  plots_significance[[variable]] <-  ggplot(dataset_pval, aes(x=Network, y=as.numeric(Rsqr)))+
    geom_abline(slope = 0,intercept = 0.05, colour="black", linetype=2,size=1)+
    geom_jitter(aes(fill=Group, alpha=Sign, size=Sign),shape=21,width = 0.5)+
    scale_alpha_manual(values = c(0.3,0.9))+
    scale_size_manual(values = c(2,7))+
    scale_fill_manual(values=c(color_groups[1],color_groups[2],color_groups[3],
                               color_groups[4],color_groups[5],color_groups[6]))+
    scale_y_continuous(limits=c(-1,1),
                       breaks =c(-1,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1) )+
    geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5), size=1, colour="grey70")+
    labs(title=Names_Variab[variable])+ylab("Corr. Coef.")+xlab("Network")+
    theme_classic()
}

png(filename ="Figures/All_Correlation_Groups.png",
    width =900*5 ,height =700*5 ,units ="px",res = 400)
grid.arrange(plots_significance[[1]],plots_significance[[2]],
             plots_significance[[3]],plots_significance[[4]],
             ncol=2,nrow=2, top="Significance values")
dev.off()      



    
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
net_names <- c("650 km", "325 km","100 km","65 km","6.5 km", "Fluvial")
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
    coin <- Network_results[[net]][c(length(Network_results[[net]])-54):length(Network_results[[net]])]
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
        scale_fill_continuous(type = "viridis",  alpha = 1,begin = 1,end = 0)+
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
        scale_fill_continuous(type = "viridis",  alpha = 1,begin = 1,end = 0)+
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
  coin <- Fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
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
      scale_fill_continuous(type = "viridis",  alpha = 1,begin = 1,end = 0)+
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
      scale_fill_continuous(type = "viridis",  alpha = 1,begin = 1,end = 0)+
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
        coin <- Network_results[[net]][c(length(Network_results[[net]])-54):length(Network_results[[net]])]
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
        
        png(filename =paste("C:/Users/Cunilleramontcusi/","Divers",biod_names[[groups]],"_",names(Network_results)[[net]],".png"),
            width =582*2 ,height =629*2 ,units ="px",res = 200)
        grid.arrange(plots_grups[[1]],plots_grups[[2]],
                     plots_grups[[3]],plots_grups[[4]],
                     plots_grups[[5]],
                     ncol=2,nrow=3, top=names(Network_results)[[net]])
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
      coin <- Fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
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
                        rep("650 km", 5), 
                        c("S16","S18","Phy","Zoo", "S18zoo"))
      
      mid_netw <- cbind(c(sign_groups[[1]][[2]][[variable]], sign_groups[[2]][[2]][[variable]], sign_groups[[3]][[2]][[variable]], 
                          sign_groups[[4]][[2]][[variable]], sign_groups[[5]][[2]][[variable]]),
                        rep("325 km", 5), 
                        c("S16","S18","Phy","Zoo", "S18zoo"))
      
      mid_mid_netw <- cbind(c(sign_groups[[1]][[3]][[variable]], sign_groups[[2]][[3]][[variable]], sign_groups[[3]][[3]][[variable]], 
                              sign_groups[[4]][[3]][[variable]], sign_groups[[5]][[3]][[variable]]),
                            rep("100 km", 5), 
                            c("S16","S18","Phy","Zoo", "S18zoo"))
      
      small_netw <- cbind(c(sign_groups[[1]][[4]][[variable]], sign_groups[[2]][[4]][[variable]], sign_groups[[3]][[4]][[variable]], 
                            sign_groups[[4]][[4]][[variable]], sign_groups[[5]][[4]][[variable]]),
                          rep("65 km", 5), 
                          c("S16","S18","Phy","Zoo", "S18zoo"))
      
      min_netw <- cbind(c(sign_groups[[1]][[5]][[variable]], sign_groups[[2]][[5]][[variable]], sign_groups[[3]][[5]][[variable]], 
                          sign_groups[[4]][[5]][[variable]], sign_groups[[5]][[5]][[variable]]),
                        rep("6.5 km", 5), 
                        c("S16","S18","Phy","Zoo", "S18zoo"))
      
      fluv_netw <- cbind(c(flu_sign_groups[[1]][[variable]], flu_sign_groups[[2]][[variable]], flu_sign_groups[[3]][[variable]], 
                           flu_sign_groups[[4]][[variable]], flu_sign_groups[[5]][[variable]]),
                         rep("Fluvial", 5), 
                         c("S16","S18","Phy","Zoo", "S18zoo"))
      
      
      dataset_pval <- as.data.frame(rbind(max_netw,mid_netw,mid_mid_netw,small_netw,min_netw,fluv_netw))
      colnames(dataset_pval) <- c("val","Network","Group")
      dataset_pval$pval <-as.numeric(dataset_pval$pval)
      dataset_pval$Network <- factor(dataset_pval$Network,
                                     levels = c("650 km", "325 km","100 km","65 km","6.5 km", "Fluvial"))
      dataset_pval$Group <- factor(dataset_pval$Group,
                                   levels = c("S16","S18","Phy","Zoo", "S18zoo"))
      significants <- rep("NoSign",5*6)
      significants[which(dataset_pval$pval<0.05)]<- "Sign"
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
    Names_Networks <- c("650 km", "325 km","100 km","65 km","6.5 km")
    for (groups in 1:5) {
      color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")
      for (net in 1:5) {
        coin <- Network_results[[net]][c(length(Network_results[[net]])-54):length(Network_results[[net]])]
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
    coin <- Network_results[[net]][c(length(Network_results[[net]])-54):length(Network_results[[net]])]
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
    
    png(filename =paste("C:/Users/Cunilleramontcusi/","All_Divers",biod_names[[groups]],"_",names(Network_results)[[net]],".png"),
        width =582*2 ,height =629*2 ,units ="px",res = 200)
    grid.arrange(plots_grups[[1]],plots_grups[[2]],
                 plots_grups[[3]],plots_grups[[4]],
                 plots_grups[[5]],
                 ncol=2,nrow=3, top=names(Network_results)[[net]])
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
  
  coin <- Fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
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
  
  png(filename =paste("C:/Users/Cunilleramontcusi/","All_Fluvial_Divers",biod_names[[groups]],"_",names(Network_results)[[net]],".png"),
      width =582*2 ,height =629*2 ,units ="px",res = 200)
  grid.arrange(plots_grups[[1]],plots_grups[[2]],
               plots_grups[[3]],plots_grups[[4]],
               plots_grups[[5]],
               ncol=2,nrow=3, top=names(Network_results)[[net]])
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
    coin <- Network_results[[net]][c(length(Network_results[[net]])-54):length(Network_results[[net]])]
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




















