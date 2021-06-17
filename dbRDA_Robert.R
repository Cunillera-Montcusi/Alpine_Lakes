
#________________________________________#
#________________________________________#
##########  PACKAGE UPLOADING   ##########
#________________________________________#
#________________________________________#

# Data treatment and analysis
library(nlme)
library(vegan)
library(car)
library(mgcv)
library(Rmisc)
library(gamm4)
library(MuMIn)
library(usdm) # Collinearity
library(randomForestSRC) # RF
library(ggRandomForests) # RF
library(gbm) # BRT
library(dismo) # BRT
library(caret)


# Ploting and image 
library(ggplot2)
library(ggforce)
library(gtable)    
library(grid)
library(gridExtra) 
library(png)
library(viridis)
library(cowplot)
library(magick)

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

# 55 lakes closeness values 
closeness_values_600km
closeness_values_300km 
closeness_values_100km 
closeness_values_60km 
closeness_values_6km 

#Values corresponding to the lakes sampled for ZOO18S
coincidence_values_Zoo18S

# Extract the values of closeness that correspond to the sampled lakes
coin <- closeness_values_600km[c(length(closeness_values_600km)-54):length(closeness_values_600km)]
centr_iso <- cbind(coin[coincidence_values_Zoo18S[which(coincidence_values_Zoo18S>0)]])

# Community data for zooplankton at 18S
community_Zoo18S
# Environmental data for zooplankton at 18S
envir_Zoo18S 

# dbRDA
dbRDA <- capscale(community_Zoo18S~envir_Zoo18S[,1]+envir_Zoo18S[,2]+envir_Zoo18S[,3]+envir_Zoo18S[,4]+envir_Zoo18S[,5],
                  distance = "jaccard",add = "lingoes")
# Extract x and y axis of the first two unconstrained dimensions.
x<- dbRDA$CA$u[,1]
y <- dbRDA$CA$u[,2]

# Build a dataset with all the values
dataset <- data.frame(x,y,centr_iso)
# Build the model with the ordisurf function
NMDS_model <- ordisurf(dbRDA$CA$u ~ dataset[,3], plot = F,)
# Extract model results
NMDS_model_results <- summary(NMDS_model)

# Extract coordinates of the surfaces to be ploted 
extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("x", "y", "z")
  return(xyz)
}
contour.vals <- extract.xyz(obj = NMDS_model)

#Plot based on all this informacion
ggplot(dataset, aes(x=x,y=y))+
  geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+
  geom_jitter(shape=21, size=5, alpha=0.8, aes(fill=centr_iso))+
  geom_contour(data=contour.vals, aes(x, y, z = z, colour = ..level..))+
  scale_colour_viridis(alpha = 1,begin = 1,end = 0)+
  scale_fill_viridis(alpha = 1,begin = 1,end = 0)+
  labs(title = paste("Zoo18S", "Network xxx km"),
       subtitle = paste("Unc. Inert.=",round(dbRDA$CA$tot.chi*1/dbRDA$tot.chi,2), 
                        "R2=", round(NMDS_model_results$r.sq,2),
                        "ED=",round(NMDS_model_results$dev.expl,2)))+
  xlab("NMDS1")+ylab("NMDS2")+
  theme_classic()+
  theme(legend.position = "none",
        panel.background=element_rect(colour="black", fill=alpha("blue", 0.1)))








