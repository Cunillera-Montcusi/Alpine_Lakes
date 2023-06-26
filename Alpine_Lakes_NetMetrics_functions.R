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


# The function "max.comp_gradiente" calculates for all distances (form the min to the max distance)
# in given distance matrix a specific network. Then, it assess how many nodes are connected in it.
# Then, from the list generated we can obtain which is the first distance where only one component
# is obtained (a network where all nodes are connected). This distance and the corresponding network are
# called percolation distance and network. 

# max.comp_gradiente_NORW is a variation where you can set the "max" and "min" distances where to search the percolation 
#distance. It can be used for very big networks where computation times can become longer. 

max.comp_gradiente<-function(distancias,br){
  T<-seq(min(distancias),max(distancias),,br)
  out<-NULL 
  for(t in T){
    L<-ifelse(distancias>t,0,1)
    C<-component.dist(L)
    C.max<-max(C$csize)
    out<-rbind(out,c(t,C.max))
  }
  return(out)
  
}


max.comp_gradiente_NORW<-function(min_distancias, max_distancias,distancias,br){
  T<-seq(min_distancias,max_distancias,,br)
  out<-NULL 
  for(t in T){
    L<-ifelse(distancias>t,0,1)
    C<-component.dist(L)
    C.max<-max(C$csize)
    out<-rbind(out,c(t,C.max))
  }
  return(out)
}

# Function to extract x,y,z from a NMDS plot 
extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("x", "y", "z")
  return(xyz)
}
