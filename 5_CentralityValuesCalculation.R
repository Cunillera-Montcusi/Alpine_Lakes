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
#detach("package:sna", unload = TRUE)
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
for (e in 1:(length(cordenades_xarxes)-1)) {
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
Fluvial_network_results[[1]] <-log(fluvial_network_data[[1]][,4],10)
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
Col <- rbPal(length(Fluvial_network_results[[1]] ))[
  as.numeric(cut(Fluvial_network_results[[1]] ,
                 breaks = length(Fluvial_network_results[[1]] )))]
Col[which(cols=="red")] <- "red"
detach("package:sna", unload = TRUE)
library(igraph)
png(filename = "Figures/Alpine_fluvial_CLOS.png", width = 20000, height = 20000, res=500)
par(mar=c(0,0,0,0))
plot(GRAPH_xarxes_fluvial[[1]], vertex.label = NA, vertex.size = sizes, vertex.size2 = sizes, vertex.color=Col,
     edge.width=0.01, edge.color="grey50")
dev.off()


# Print network structure results_______________________________________________
#detach("package:sna", unload = TRUE)
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

save.image("Database.RData")

