
##########################################
##########                      ##########
##########  PACKAGE UPLOADING   ##########
##########                      ##########
##########################################

library(vegan)

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
#_______________________________________________________________###

# LOAD DATASET from NORWAY _______________________________####




# LOAD DATASET from SODA PANS _______________________________####




##########################################

a <- comm_data$S16
b <- as.matrix(env_data$env.S16[,4:6])

comm_data$S18
as.matrix(env_data$env.S18[,4:6])

comm_data$phyto
as.matrix(env_data$env.phyto[,4:6])

comm_data$zoo
as.matrix(env_data$env.zoo[,4:6])


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
                             as.matrix(vegdist(a, method = distance))[e,-e]) # Same for the fitteds
  }
  print(tst_coeficient)
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
                             as.numeric(communnity_dataset[e,]))
  }
  print(tst_coeficient)
}


### WARNING: CHANGE TO JACCARD IF WORKING WITH PA
dbRDA_Env_Track(communnity_dataset = a, environmental_dataset = b,distance = "jaccard")

CCA_Env_Track(communnity_dataset = a, environmental_dataset = b)

