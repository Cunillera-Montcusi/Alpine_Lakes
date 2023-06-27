
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

save.image("Database.RData")