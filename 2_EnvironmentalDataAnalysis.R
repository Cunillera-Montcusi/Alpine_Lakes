
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

save.image("Database.RData")

