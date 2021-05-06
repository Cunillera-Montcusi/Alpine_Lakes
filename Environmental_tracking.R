
##########################################
##########                      ##########
##########  PACKAGE UPLOADING   ##########
##########                      ##########
##########################################
#Anlisys
library(vegan)
library(mgcv)
library(MuMIn)
#Plot and data management
library(ggplot2)
library(tidyverse)
library(gtable)    
library(grid)
library(gridExtra) 

# Define the workind directory
setwd("C:/Users/Cunilleramontcusi/Dropbox/DAVID DOC/LLAM al DIA/1. Lunz al DIA/ALPINE_Lakes/Alpine_Lakes/upload_Envtrack")

source("CUNILLERA_palette.R")

##########################################
# LOAD DATASET from ALPINE LAKES _______________________________####
##### 16S bacterial plankton
load("com.16S.Rdata")
load("env.16S.Rdata")
com.16S <- com
env.16S <- env 
# select OTUS with min occurence of 10 sites (=exclude rare OTUs)
seltax<-which(colSums(com.16S>0)>9)
com_S16<-com.16S[,seltax]

##### 18S bacterial plankton
load("com.18S.Rdata")
load("env.18S.Rdata")
com.18S <- com
env.18S <- env 
# select OTUS with min occurence of 10 sites (=exclude rare OTUs)
seltax<-which(colSums(com.18S>0)>9)
com.18S<-com.18S[,seltax]

##### Phytoplankton bacterial plankton
load("com.phyt.Rdata")
load("env.phyt.Rdata")
com.phy <- com
env.phy <- env 
# select OTUS with min occurence of 2 sites (=exclude rare OTUs)
seltax<-which(colSums(com.phy>0)>1)
com.phy<-com.phy[,seltax]

##### Zooplankton bacterial plankton
load("com.zoo.Rdata")
load("env.zoo.Rdata")
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

load("Norw.com.fish.Rdata")
load("Norw.com.zoo.Rdata")
load("Norw.com.phyto.Rdata")

Norw.fish <- fish
apply(Norw.fish,1,sum)

Norw.zoo <- ifelse(zoo>0,1,0)
apply(Norw.zoo,1,sum)
Norw.zoo <- Norw.zoo[-5,] 

Norw.phyto <- ifelse(com>0,1,0)
apply(Norw.phyto,1,sum)

load("Norw.env.phyto.Rdata")
Norw.env <- env[,c(7,8,14,16:23)]
Norw.env[is.na(Norw.env)] <- 0
Norw.env <- as.matrix(log(Norw.env+1))

Norw.env.zoo <- Norw.env[-5,] 

# LOAD DATASET from SODA PANS _______________________________####
load("Seewinkel_predictors.Rdata")
load("Seewinkel_zoopl.Rdata")

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
                             as.matrix(vegdist(communnity_dataset, method = distance))[e,-e], method = "pearson") # Same for the fitteds
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
                             as.numeric(communnity_dataset[e,]), method = "pearson")
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
Results_treatment <- function(Env_Track, lon, lat,number_colors=30, title, subtitle){
Best_Models <- list()
Plots <- list()

value.results <- as.data.frame(cbind(Env_Track,lon,lat))
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

# Generate dataset to plot all the values
z <- as.data.frame(matrix(pred2, nrow=30))  
plot_data <- data.frame(x=xy_layout[,1],y=xy_layout[,2], z=pivot_longer(z,cols=1:30))
colors <- CUNILLERA_pal("LGTBI_pale",reverse = F)(number_colors)[cut(value.results$EnvTrack, number_colors)]

Plots[[1]] <- ggplot(plot_data)+
  # Raster with all the data, hjust and vjust hide the bakground and INTERPOLATE= T dissolve the rasters 
  geom_raster(aes(x = x,y=y,fill=z.value),hjust=1,vjust=1,interpolate = T)+
  # Scale y axis and set the levels
  scale_y_continuous(limits = c(min(plot_data$y, na.rm=TRUE),max(plot_data$y, na.rm=TRUE)),
                     expand = c(0,0),
                     breaks = seq(from=min(plot_data$y, na.rm=TRUE), to=max(plot_data$y, na.rm=TRUE), by=0.2))+
  # Scale x axis and set the levels--> WARNING: The maximum distance is 21km!! IS A LIE!
  scale_x_continuous(limits = c(min(plot_data$x, na.rm=TRUE),max(plot_data$x, na.rm=TRUE)),
                     expand = c(0,0),
                     breaks = seq(from=min(plot_data$x, na.rm=TRUE), to=max(plot_data$x, na.rm=TRUE), by=0.2))+
  # Scalate thte colours, spectral means going from red to blue LIMITS force the range to be btween 0 & 100
  scale_fill_CUNILLERA(palette = "LGTBI_pale",discrete = F, reverse = F)+
  geom_contour(aes(x = x,y=y,z=z.value), col="grey45", size=0.1, linetype=2)+
  # Labs 
  labs(x="Longitude",y="Latitude",title = title,subtitle = subtitle)+
  geom_jitter(data =value.results, aes(y=lat, x=lon),fill=colors ,col="black", size=2, stroke=2,shape=21)+
  # Theme stuff
  theme_classic()+
  theme(axis.line = element_line(size = 0.05, linetype = "solid"), 
        axis.ticks = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(0.4,1.5,0.4,1.5, "cm"),
        axis.title = element_text(size = 10,face = "bold"), 
        axis.text = element_text(size = 10,face = "bold"), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.position = "right",
        panel.background = element_rect(colour = "black", size = 2.5, linetype = "solid"),
        plot.background = element_rect(colour = "black", size = 1, linetype = "solid"), 
        plot.subtitle = element_text(size = 12,face = "bold", colour = "gray0"),
        plot.title = element_text(size = 14, face = "bold"))

list(Best_Models[[1]], Plots[[1]])
}

# ALPINE RESULTS ####
# S16 ________________________________________________________________________________________#
#dbRDA
s16_dbRDA <- Results_treatment(Env_Track =output_dbRDA[[1]] ,lon = env.16S$lon,lat =env.16S$lat, number_colors=20, title="Alp_S16", subtitle="dbRDA")
#CCA
s16_CCA <- Results_treatment(Env_Track =output_CCA[[1]] ,lon = env.16S$lon,lat =env.16S$lat, number_colors=20, title="Alp_S16", subtitle="CCA")

# S18 ________________________________________________________________________________________#
#dbRDA
s18_dbRDA <- Results_treatment(Env_Track =output_dbRDA[[2]] ,lon = env.18S$lon,lat =env.18S$lat, number_colors=35, title="Alp_S18", subtitle="dbRDA")
#CCA
s18_CCA <- Results_treatment(Env_Track =output_CCA[[2]] ,lon = env.18S$lon,lat =env.18S$lat, number_colors=35, title="Alp_S18", subtitle="CCA")

# PHYT ________________________________________________________________________________________#
#dbRDA
phyt_dbRDA <- Results_treatment(Env_Track =output_dbRDA[[3]] ,lon = env.phy$lon,lat =env.phy$lat, number_colors=35, title="Alp_PHYT", subtitle="dbRDA")
#CCA
phyt_CCA <- Results_treatment(Env_Track =output_CCA[[3]] ,lon = env.phy$lon,lat =env.phy$lat, number_colors=20, title="Alp_PHYT", subtitle="CCA")

# ZOO ________________________________________________________________________________________#
#dbRDA
zoo_dbRDA <- Results_treatment(Env_Track =output_dbRDA[[4]] ,lon = env.zoo$lon,lat =env.zoo$lat, number_colors=35, title="Alp_ZOO", subtitle="dbRDA")
#CCA
zoo_CCA <- Results_treatment(Env_Track =output_CCA[[4]] ,lon = env.zoo$lon,lat =env.zoo$lat, number_colors=25, title="Alp_ZOO", subtitle="CCA")

# NORW RESULTS ####
# Fish ________________________________________________________________________________________#
#dbRDA
norw.fish_dbRDA <- Results_treatment(Env_Track =output_dbRDA[[5]] ,lon = env$Longitude,lat =env$Latitude, number_colors=40,title="Norw_Fish", subtitle="dbRDA")
#CCA
norw.fish_CCA <- Results_treatment(Env_Track =output_CCA[[5]] ,lon = env$Longitude,lat =env$Latitude, number_colors=40,title="Norw_Fish", subtitle="CCA")

# Zoo ________________________________________________________________________________________#
#dbRDA
norw.zoo_dbRDA <- Results_treatment(Env_Track =output_dbRDA[[6]] ,lon = env$Longitude[-5],lat =env$Latitude[-5], number_colors=45,title="Norw_ZOO", subtitle="dbRDA")
#CCA
norw.zoo_CCA <- Results_treatment(Env_Track =output_CCA[[6]] ,lon = env$Longitude[-5],lat =env$Latitude[-5], number_colors=35,title="Norw_ZOO", subtitle="CCA")

# Phy ________________________________________________________________________________________#
#dbRDA
norw.phy_dbRDA <- Results_treatment(Env_Track =output_dbRDA[[7]],lon = env$Longitude,lat =env$Latitude, number_colors=45,title="Norw_Phy", subtitle="dbRDA")
#CCA
norw.phy_CCA <- Results_treatment(Env_Track =output_CCA[[7]],lon = env$Longitude,lat =env$Latitude, number_colors=30,title="Norw_Phy", subtitle="CCA")

# SODA PANS RESULTS ####
# dbRDA __________________________________________________________________________________#
#dbRDA
Soda_dbRDA <- Results_treatment(Env_Track =output_dbRDA[[8]],lon = pred$lon,lat =pred$lat, number_colors=40,title="Soda", subtitle="dbRDA")

#CCA
Soda_CCA <- Results_treatment(Env_Track =output_CCA[[8]],lon = pred$lon,lat =pred$lat, number_colors=45,title="Soda", subtitle="CCA")

# Arrange the two charts
# The legend boxes are centered
grid.newpage()
png(filename = "Env_Track.png" ,width=3500,height=8000,units="px",res=300)
grid.arrange(s16_dbRDA[[2]],s16_CCA[[2]],
             s18_dbRDA[[2]],s18_CCA[[2]],
             phyt_dbRDA[[2]],phyt_CCA[[2]],
             zoo_dbRDA[[2]],zoo_CCA[[2]],
             norw.fish_dbRDA[[2]],norw.fish_CCA[[2]],
             norw.zoo_dbRDA[[2]],norw.zoo_CCA[[2]],
             norw.phy_dbRDA[[2]],norw.phy_CCA[[2]],
             Soda_dbRDA[[2]],Soda_CCA[[2]],
             ncol = 2, nrow=8)
dev.off()


All.longitudes <- list(env.16S$lon, env.18S$lon, env.phy$lon, env.zoo$lon,
                        env$Longitude, env$Longitude[-5], env$Longitude, pred$lon)
All.latitudes <- list(env.16S$lat, env.18S$lat, env.phy$lat, env.zoo$lat,
                       env$Latitude, env$Latitude[-5], env$Latitude, pred$lat)
titles.list <- c("Alp_S16","Alp_S18","Alp_PHYT","Alp_ZOO","Norw_Fish","Norw_ZOO","Norw_Phy","Soda")

library(corrmorant)

corplots_output <- list()
for (w in 1:8) {
a <- cbind(All.longitudes[[w]],All.latitudes[[w]], output_dbRDA[[w]],output_CCA[[w]])
colnames(a) <- c("Long","Lat","dbRDA","CCA")
corplots_output[[w]] <- ggcorrm(a, 
        mapping = aes(col = .corr, fill = .corr),
        bg_dia = "grey20", 
        rescale = "by_sd") +
  lotri(geom_smooth(method = "lm", size = .3)) +
  lotri(geom_point(alpha = 0.5)) +
  utri_corrtext(nrow = 2, squeeze = 0.6) +
  dia_names(y_pos = 0.15, size = 3, color = "white") +
  dia_density(lower = 0.3, color = "grey80", fill = "grey60", size = .3) +
  scale_color_corr(aesthetics = c("fill", "color"))+
  labs(title=titles.list[w])
}
grid.newpage()
png(filename = "Corplots_Env_Track.png" ,width=8000,height=4000,units="px",res=400)
grid.arrange(corplots_output[[1]],corplots_output[[2]],
             corplots_output[[3]],corplots_output[[4]],
             corplots_output[[5]],corplots_output[[6]],
             corplots_output[[7]],corplots_output[[8]],
             ncol = 4, nrow=2)
dev.off()


library(betareg)
summary(lm(output_CCA[[8]]~All.longitudes[[8]]))
summary(lm(output_CCA[[8]]~All.latitudes[[8]]))

summary(betareg(output_CCA[[8]]~All.longitudes[[8]]))
summary(betareg(output_CCA[[8]]~All.latitudes[[8]]))

par(mfrow=c(1,2))
plot(All.longitudes[[8]],output_CCA[[8]], main = titles.list[8], xlab = "Longitude")
abline(a=summary(betareg(output_CCA[[8]]~All.longitudes[[8]]))[1]$coefficients$mean[1,1],
       b=summary(betareg(output_CCA[[8]]~All.longitudes[[8]]))[1]$coefficients$mean[2,1])
plot(All.latitudes[[8]],output_CCA[[8]], main = titles.list[w], xlab = "Latitude")
abline(a=summary(betareg(output_CCA[[8]]~All.latitudes[[8]]))[1]$coefficients$mean[1,1],
       b=summary(betareg(output_CCA[[8]]~All.latitudes[[8]]))[1]$coefficients$mean[2,1],
       col="red",cex=3)
text(x = 47.71,y = 0.65,
  labels = paste("Pse-R2=",
   round(betareg(output_CCA[[8]]~All.latitudes[[8]])$pseudo.r.squared,2)))
   


# Check 

##########################################
##########                      ##########
##########  PACKAGE UPLOADING   ##########
##########                      ##########
##########################################
#Anlisys
library(vegan)
library(mgcv)
library(MuMIn)
#Plot and data management
library(ggplot2)
library(tidyverse)
library(gtable)    
library(grid)
library(gridExtra) 

# Define the workind directory
setwd("C:/Users/Cunilleramontcusi/Dropbox/DAVID DOC/LLAM al DIA/1. Lunz al DIA/ALPINE_Lakes/Alpine_Lakes/upload_Envtrack")

source("CUNILLERA_palette.R")

##########################################
# LOAD DATASET from ALPINE LAKES _______________________________####
##### 16S bacterial plankton
load("com.16S.Rdata")
load("env.16S.Rdata")
com.16S <- com
env.16S <- env 
# select OTUS with min occurence of 10 sites (=exclude rare OTUs)
seltax<-which(colSums(com.16S>0)>9)
com_S16<-com.16S[,seltax]

##### 18S bacterial plankton
load("com.18S.Rdata")
load("env.18S.Rdata")
com.18S <- com
env.18S <- env 
# select OTUS with min occurence of 10 sites (=exclude rare OTUs)
seltax<-which(colSums(com.18S>0)>9)
com.18S<-com.18S[,seltax]

##### Phytoplankton bacterial plankton
load("com.phyt.Rdata")
load("env.phyt.Rdata")
com.phy <- com
env.phy <- env 
# select OTUS with min occurence of 2 sites (=exclude rare OTUs)
seltax<-which(colSums(com.phy>0)>1)
com.phy<-com.phy[,seltax]

##### Zooplankton bacterial plankton
load("com.zoo.Rdata")
load("env.zoo.Rdata")
com.zoo <- com
env.zoo <- env 

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

ALP_list_env <- list(ALP_S16_env,ALP_S18_env,ALP_phyto_env,ALP_Zoo_env)
#_______________________________________________________________###

ALP_list_com_PA <- list(com.16S, com.18S, com.phy, com.zoo)

com.16S <- ifelse(com.16S>0,1,0)
com.18S <- ifelse(com.18S>0,1,0)
com.phy <- ifelse(com.phy>0,1,0)
com.zoo <- ifelse(com.zoo>0,1,0)

ALP_list_com <- list(com.16S, com.18S, com.phy, com.zoo)

CCA_Env_Track <- function(communnity_dataset,environmental_dataset){
  require(pscl)
  # Calculate the CCA with the community data and the environmental and the corresponding distance
  # P/A - Jaccard distance
  cca1<-cca(ALP_list_com_PA[[4]]~ALP_list_env[[4]])
  
  # Fitteds - we predict according to the model where our samples should be
  fitteds <- fitted(cca1, model="CCA", type= "response")
  
  tst_coeficient <- c()# Output to drop the coefficients
  # Correlation between predicted and observed values
  for (e in 1:nrow(ALP_list_com_PA[[4]])) {
    
    #tst_coeficient[e] <- cor(fitteds[e,],as.numeric(ALP_list_com_PA[[4]][e,]),method = "pearson")
    
    model <- glm(as.numeric(ALP_list_com_PA[[4]][e,])~fitteds[e,],family = "binomial")
    tst_coeficient[e] <- pR2(model)["McFadden"]
  }
  tst_coeficient
}


Output_PA <- list()
Output <- list()
for (i in 1:4) {
Output_PA[[i]]<- CCA_Env_Track(communnity_dataset =ALP_list_com_PA[[i]] ,environmental_dataset =ALP_list_env[[i]])
Output[[i]] <- CCA_Env_Track(communnity_dataset =ALP_list_com[[i]] ,environmental_dataset =ALP_list_env[[i]])  
}

corplots_output <- list()
for (w in 1:4) {
      a <- cbind(Output_PA[[w]],Output[[w]])
      colnames(a) <- c("CCA_PA","CCA")
      corplots_output[[w]] <- ggcorrm(a, 
                                      mapping = aes(col = .corr, fill = .corr),
                                      bg_dia = "grey20", 
                                      rescale = "by_sd") +
        lotri(geom_smooth(method = "lm", size = .3)) +
        lotri(geom_point(alpha = 0.5)) +
        utri_corrtext(nrow = 2, squeeze = 0.6) +
        dia_names(y_pos = 0.15, size = 3, color = "white") +
        dia_density(lower = 0.3, color = "grey80", fill = "grey60", size = .3) +
        scale_color_corr(aesthetics = c("fill", "color"))
}
    
  grid.newpage()
    png(filename = "Corplots_Env_Track.png" ,width=8000,height=4000,units="px",res=400)
    grid.arrange(corplots_output[[1]],corplots_output[[2]],
                 corplots_output[[3]],corplots_output[[4]],
                 ncol = 2, nrow=2)
    dev.off()



