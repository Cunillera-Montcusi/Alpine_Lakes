
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
