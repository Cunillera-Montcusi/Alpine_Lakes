
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
##########  NET vs DIV VALUES   ##########
#________________________________________#
#________________________________________#

# Extracting the values of which lakes have been sampled for each taxonomic group
# S16 = 52 lakes
# S18 = 48 lakes
# Phy = 50 lakes
# Zoo = 52 lakes
# 18S.Zoo = 48 lakes

biod_names <- c("S16","S18","Phy","Zoo", "18S Zoo")

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
histo_list[[6]] <- ggplot(data.frame(y=Fluvial_network_results[[1]]))+
  geom_histogram(aes(x=y), bins=10,colour="black",fill="grey60")+
  geom_vline(xintercept = coin, colour="black")+
  scale_y_continuous(expand = c(0,0))+theme_classic()+labs(x="Out-closeness", title = "Fluvial")


png(filename ="Figures/Clos_Histo.png",
    width =582*4 ,height =629*4 ,units ="px",res = 300)
grid.arrange(histo_list[[1]],histo_list[[2]],histo_list[[3]],histo_list[[4]],histo_list[[5]],histo_list[[6]])
dev.off()


# GAM models_______________________________________________________________________________ ####
Names_Networks <- c("650 km", "325 km","100 km","65 km","6.5 km")
Metr_names <- c("Rich.","LCBD","Repl.","Rich.Diff.")

GAMmodel_resutls_total <- list()
GAM_model_resutls <- list()
GAM_direct_model <- list()
GAM_direct_model_total <- list()

# For Euclidean networks with 600, 300 100, 60 and 6 km
for (groups in 1:5) {
  color_groups <- as.character(CUNILLERA_cols("yellow","blue","green","red","cyan"))
  for (net in 1:5) {
    
    if(net==4 | net==5){biod <- biod_group_packs[[2]]}else{biod <- biod_group_packs[[1]]} # We select which "subgraph" we use
    
    coin <- Network_results[[net]][c(length(Network_results[[net]])-54):length(Network_results[[net]])]
    dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                     biod[[groups]][,1:4])
    colnames(dataset)[1] <-c("Network")
    colnames(dataset)[2:5] <- paste(biod_names[groups],Metr_names)
    
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
  
  biod <- biod_group_packs[[2]]
  
  color_groups <- as.character(CUNILLERA_cols("yellow","blue","green","red","cyan"))
  coin <- Fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
  dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                   biod[[groups]][,1:4])
  colnames(dataset)[1] <-c("Network")
  colnames(dataset)[2:5] <- colnames(dataset)[2:5] <- paste(biod_names[groups],Metr_names)
  
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
biod_names <- c("S16","S18","Phy","Zoo", "18S Zoo")
Names_Networks <- c("650 km", "325 km","100 km","65 km","6.5 km", "Fluvial")
Metr_names <- c("Rich.","LCBD","Repl.","Rich.Diff.")

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
      
      if(net==4 | net==5){biod <- biod_group_packs[[2]]}else{biod <- biod_group_packs[[1]]} # We select which "subgraph" we use
      
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
    if(net==4 | net==5){biod <- biod_group_packs[[2]]}else{biod <- biod_group_packs[[1]]} # We select which "subgraph" we use
    
    coin <- Network_results[[net]][c(length(Network_results[[net]])-54):length(Network_results[[net]])]
    dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                     biod[[groups]][,1:4])
    colnames(dataset)[1] <-c("Network")
    colnames(dataset)[2:5] <- colnames(dataset)[2:5] <- paste(biod_names[groups],Metr_names)
    
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
  
  biod <- biod_group_packs[[2]]
  
  color_groups <- CUNILLERA_cols("yellow","blue","green","red","cyan")
  coin <- Fluvial_network_results[[1]][all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]]]
  dataset <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]],
                   biod[[groups]][,1:4])
  colnames(dataset)[1] <-c("Network")
  colnames(dataset)[2:5] <- colnames(dataset)[2:5] <- paste(biod_names[groups],Metr_names)
  
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
    width =500*12 ,height =650*10 ,units ="px",res = 300)
grid.arrange(
arrangeGrob(
arrangeGrob(arrangeGrob(ggplot() + geom_blank()),nrow=1),
arrangeGrob(arrangeGrob(ggplot() + geom_blank()),nrow=1),
arrangeGrob(
arrangeGrob(plot_plot_sign_plot[[2]],top =textGrob("A)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
arrangeGrob(plot_plot_sign_plot[[11]],top =textGrob("E)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),nrow=2), # 100 km
arrangeGrob(
arrangeGrob(plot_plot_sign_plot[[1]],top =textGrob("B)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
arrangeGrob(plot_plot_sign_plot[[8]],top =textGrob("F)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),nrow=2), # 325 km
arrangeGrob(
arrangeGrob(plot_plot_sign_plot[[3]],top =textGrob("C)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
arrangeGrob(plot_plot_sign_plot[[6]],top =textGrob("G)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),nrow=2),# 650km
arrangeGrob(
arrangeGrob(plot_plot_sign_plot[[16]],top =textGrob("D)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
arrangeGrob(plot_plot_sign_plot[[17]],top =textGrob("H)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),nrow=2), #Fluvial

left="Richness", nrow=1,ncol=6),

arrangeGrob(
arrangeGrob(arrangeGrob(ggplot() + geom_blank()),nrow=1),
arrangeGrob(arrangeGrob(ggplot() + geom_blank()),nrow=1),
arrangeGrob(
arrangeGrob(plot_plot_sign_plot[[5]],top =textGrob("I)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
arrangeGrob(plot_plot_sign_plot[[13]],top =textGrob("L)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),nrow=2), # 100 km
arrangeGrob(
arrangeGrob(plot_plot_sign_plot[[4]],top =textGrob("J)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
arrangeGrob(plot_plot_sign_plot[[10]],top =textGrob("M)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),nrow=2), # 325 km
arrangeGrob(arrangeGrob(ggplot() + geom_blank()),nrow=1),
arrangeGrob(
arrangeGrob(plot_plot_sign_plot[[14]],top =textGrob("K)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
arrangeGrob(plot_plot_sign_plot[[15]],top =textGrob("N)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),nrow=2), # FluviaL
left="Replacement",nrow=1,ncol=6),

arrangeGrob(
arrangeGrob(arrangeGrob(ggplot() + geom_blank()),nrow=1),
arrangeGrob(arrangeGrob(ggplot() + geom_blank()),nrow=1),
arrangeGrob(arrangeGrob(ggplot() + geom_blank()),nrow=1),
arrangeGrob(arrangeGrob(ggplot() + geom_blank()),nrow=1),
arrangeGrob(arrangeGrob(ggplot() + geom_blank()),nrow=1),
arrangeGrob(
arrangeGrob(plot_plot_sign_plot[[18]],top =textGrob("O)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
arrangeGrob(ggplot() + geom_blank()),nrow=2),
left="Rich.Diff.",nrow=1,ncol=6),

arrangeGrob(
arrangeGrob(arrangeGrob(ggplot() + geom_blank()),nrow=1),
arrangeGrob(arrangeGrob(ggplot() + geom_blank()),nrow=1),
arrangeGrob(arrangeGrob(plot_plot_sign_plot[[12]],top =textGrob("P)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
            arrangeGrob(ggplot() + geom_blank()),nrow=2),# 100 km  
arrangeGrob(arrangeGrob(plot_plot_sign_plot[[9]],top =textGrob("Q)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
            arrangeGrob(ggplot() + geom_blank()),nrow=2),# 325 km  
arrangeGrob(arrangeGrob(plot_plot_sign_plot[[7]],top =textGrob("R)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
            arrangeGrob(ggplot() + geom_blank()),nrow=2),# 650 km  
arrangeGrob(arrangeGrob(ggplot() + geom_blank()),nrow=1),# Fluvial km  
left = "LCBD",nrow=1,ncol=6),

ncol=1,nrow=4, top="")
dev.off()


# NMDS plots_______________________________________________________________________________ ####
biod_names <- c("S16","S18","Phy","Zoo", "18S Zoo")
net_names <- c("650 km", "325 km","100 km","65 km","6.5 km", "Fluvial")
color_groups <- as.character(CUNILLERA_cols("yellow","blue","green","red","cyan"))

# For Euclidean network
plots_NMDS <- list()
plots_NMDS_total <- list()

plots_NMDS_model_result <- list()
plots_NMDS_total_model_result <- list()

for (groups in 1:5) {
  for (net in 1:5) {
    
    comun_data <- list()
    envir_data <- list()
    netw_centr <- list()
    if(net==4 | net==5){
    #cordenades_xarxes[[net]][which(components(graph.adjacency(MAPS_xarxes[[net]], mode = "undirected",diag = F))$membership==2),]
    
    Gral_coin <- which(1:length(Network_results[[net]])>=length(Network_results[[net]])-54)[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]]  
      
    Swiss <- as.vector(the_plot_db %>% filter(Longitude<10) %>% select(Lake_name))$Lake_name
    comun_data[[1]] <- comm_data[[groups]][rownames(comm_data[[groups]])%in%Swiss,]
    envir_data[[1]] <- env_data[[groups]][rownames(env_data[[groups]])%in%Swiss,]
    
    filter_subgraph<- which(which(components(graph.adjacency(MAPS_xarxes[[net]],mode = "undirected",diag = F))$membership==2)>=c(length(Network_results[[net]])-54))
    filter_subgraph_position <-  which(components(graph.adjacency(MAPS_xarxes[[net]],mode = "undirected",diag = F))$membership==2)[filter_subgraph]
    netw_centr[[1]] <- filter_subgraph_position[which(filter_subgraph_position%in%Gral_coin)]
    
    
    Aust <- as.vector(the_plot_db %>% filter(Longitude>10) %>% select(Lake_name))$Lake_name
    comun_data[[2]] <- comm_data[[groups]][rownames(comm_data[[groups]])%in%Aust,]
    envir_data[[2]] <- env_data[[groups]][rownames(env_data[[groups]])%in%Aust,]
    
    filter_subgraph<- which(which(components(graph.adjacency(MAPS_xarxes[[net]],mode = "undirected",diag = F))$membership==1)>=c(length(Network_results[[net]])-54))
    filter_subgraph_position <-  which(components(graph.adjacency(MAPS_xarxes[[net]],mode = "undirected",diag = F))$membership==1)[filter_subgraph]
    netw_centr[[2]] <- filter_subgraph_position[which(filter_subgraph_position%in%Gral_coin)]
    
    }else{
    comun_data[[1]] <- comm_data[[groups]]
    envir_data[[1]] <- env_data[[groups]]
    netw_centr[[1]] <-which(which(components(graph.adjacency(MAPS_xarxes[[net]], 
                                                             mode = "undirected",diag = F))$membership==1)>=c(length(Network_results[[net]])-54))
    } # We select which "subgraph" we use

    group_pack_NMDS_model_result <- list()
    group_pack_plots_NMDS <- list()
    SubGraphName <- c("- Switzerland","- Austria")
    for (group_pack in 1:length(comun_data)) {
    if(length(comun_data)==1){SubGraphName <- ""}  
      
    coin <- Network_results[[net]][netw_centr[[group_pack]]] # OLD way [c(length(Network_results[[net]])-54):length(Network_results[[net]])]
    if (length(comun_data)==1) {centr_iso <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]])}else{centr_iso <- coin}
    
    dbRDA <- capscale(comun_data[[group_pack]]~
                      envir_data[[group_pack]][,4]+envir_data[[group_pack]][,5]+envir_data[[group_pack]][,6]+envir_data[[group_pack]][,7]+envir_data[[group_pack]][,8],
                      distance = "jaccard",add = "lingoes")
    
    x<- dbRDA$CA$u[,1]
    y <- dbRDA$CA$u[,2]
    
    dataset <- data.frame(x,y,centr_iso)
    NMDS_model <- ordisurf(dbRDA$CA$u ~ dataset[,3], plot = F,)
    NMDS_model_results <- summary(NMDS_model)
    
    group_pack_NMDS_model_result[[group_pack]] <- NMDS_model_results
    
    contour.vals <- extract.xyz(obj = NMDS_model)
    
    the_plot <-ggplot(dataset, aes(x=x,y=y))+
      geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+
      geom_jitter(shape=21, size=5, alpha=0.8, aes(fill=centr_iso))+
      geom_contour(data=contour.vals, aes(x, y, z = z, colour = ..level..))+
      #geom_path(data = df_ellipse, aes(x=x, y=y, colour=Group), size=2, show.legend = FALSE)+
      scale_colour_viridis(alpha = 1,begin = 1,end = 0)+
      #manual(values = viridis_pal(0.9,1,0,direction = -1)(length(unique(df_ellipse$Group))))+
      scale_fill_viridis(alpha = 1,begin = 1,end = 0)+
      labs(title = paste(biod_names[groups], net_names[net],SubGraphName[group_pack]),
           subtitle = paste("Unc. Inert.=",round(dbRDA$CA$tot.chi*1/dbRDA$tot.chi,2), 
                            "R2=", round(NMDS_model_results$r.sq,2),
                            "ED=",round(NMDS_model_results$dev.expl,2)))+
      xlab("dbRDA1")+ylab("dbRDA2")+
      theme_classic()+
      theme(legend.position = "none",
            panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)))
    
    group_pack_plots_NMDS[[group_pack]] <- ggdraw() +  draw_plot(the_plot)+
      draw_image(magick::image_read(image_list[[groups]]),
                 scale = 0.2,x = 0.4,y = 0.38) 
    

    group_pack_NMDS_model_result[[group_pack]]
    } # Group_pack
    if (length(group_pack_plots_NMDS)==1){plots_NMDS[[net]] <- group_pack_plots_NMDS[[1]]}else{plots_NMDS[[net]] <- group_pack_plots_NMDS}
    if (length(group_pack_plots_NMDS)==1){plots_NMDS_model_result[[net]] <- group_pack_NMDS_model_result[[1]]}else{plots_NMDS_model_result[[net]] <- group_pack_NMDS_model_result}
  }# Net
  plots_NMDS_total[[groups]] <- plots_NMDS
  plots_NMDS_total_model_result[[groups]] <- plots_NMDS_model_result
}

# For Fluvial network
plots_NMDS_fluvial_model_result <- list()
plots_NMDS_fluvial_total_model_result <- list()
plots_NMDS_total_fluvial<- list()
plots_NMDS_fluvial_total_model_result <- list()
for (groups in 1:5) {
  comun_data <- list()
  envir_data <- list()
  netw_centr <- list()
 
  #cordenades_xarxes[[net]][which(components(graph.adjacency(MAPS_xarxes[[net]], mode = "undirected",diag = F))$membership==2),]
  
  Gral_coin <- V(GRAPH_xarxes_fluvial[[1]])[c(all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]])[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]]]
  
  Swiss <- as.vector(the_plot_db %>% filter(Longitude<10) %>% select(Lake_name))$Lake_name
  comun_data[[1]] <- comm_data[[groups]][rownames(comm_data[[groups]])%in%Swiss,]
  envir_data[[1]] <- env_data[[groups]][rownames(env_data[[groups]])%in%Swiss,]
  out_site <- c()           
  for (ds in 1:length(Gral_coin)) {
    out_temp<- which(which(components(GRAPH_xarxes_fluvial[[1]])$membership==1)==Gral_coin[ds])
    if (length(out_temp)==0) {out_temp <- 0}
    out_site[ds] <- out_temp}
  netw_centr[[1]] <-   which(components(GRAPH_xarxes_fluvial[[1]])$membership==1)[out_site[which(out_site>0)]]
    
  Aust <- as.vector(the_plot_db %>% filter(Longitude>10) %>% select(Lake_name))$Lake_name
  comun_data[[2]] <- comm_data[[groups]][rownames(comm_data[[groups]])%in%Aust,]
  envir_data[[2]] <- env_data[[groups]][rownames(env_data[[groups]])%in%Aust,]
  out_site <- c()           
  for (ds in 1:length(Gral_coin)) {
    out_temp<- which(which(components(GRAPH_xarxes_fluvial[[1]])$membership==2)==Gral_coin[ds])
    if (length(out_temp)==0) {out_temp <- 0}
    out_site[ds] <- out_temp}
  netw_centr[[2]] <-   which(components(GRAPH_xarxes_fluvial[[1]])$membership==2)[out_site[which(out_site>0)]]

  group_pack_NMDS_model_result <- list()
  group_pack_plots_NMDS <- list()
  
  SubGraphName <- c("- Switzerland","- Austria")
  for (group_pack in 1:length(comun_data)) {
  if(length(comun_data)==1){SubGraphName <- ""}  
    
  coin <- Fluvial_network_results[[1]][netw_centr[[group_pack]]]
  centr_iso <- coin 
    
  dbRDA <- capscale(comun_data[[group_pack]]~
                    envir_data[[group_pack]][,4]+envir_data[[group_pack]][,5]+envir_data[[group_pack]][,6]+envir_data[[group_pack]][,7]+envir_data[[group_pack]][,8],
                    distance = "jaccard",add = "lingoes")
  
  x<- dbRDA$CA$u[,1]
  y <- dbRDA$CA$u[,2]
  
  dataset <- data.frame(x,y,centr_iso)
  NMDS_model <- ordisurf(dbRDA$CA$u ~ dataset[,3], plot = F)
  NMDS_model_results <- summary(NMDS_model)
  
  group_pack_NMDS_model_result[[group_pack]] <- NMDS_model_results
  
  contour.vals <- extract.xyz(obj = NMDS_model)
  
  the_plot<- ggplot(dataset, aes(x=x,y=y))+
    geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+
    geom_jitter(shape=21, size=5, alpha=0.8, aes(fill=centr_iso))+
    geom_contour(data=contour.vals, aes(x, y, z = z, colour = ..level..))+
    scale_colour_viridis(alpha = 1,begin = 1,end = 0)+
    scale_fill_viridis(alpha = 1,begin =1 ,end = 0)+
    labs(title = paste(biod_names[groups], net_names[6],SubGraphName[group_pack]),
         subtitle = paste("Unc. Inert.=",round(dbRDA$CA$tot.chi*1/dbRDA$tot.chi,2), 
                          "R2=", round(NMDS_model_results$r.sq,2),
                          "ED=",round(NMDS_model_results$dev.expl,2)))+
    xlab("dbRDA1")+ylab("dbRDA2")+
    theme_classic()+
    theme(legend.position = "none",
          panel.background=element_rect(colour="black", fill=alpha(color_groups[groups], 0.1)))
  
  group_pack_plots_NMDS[[group_pack]]<-ggdraw() +  draw_plot(the_plot)+
    draw_image(magick::image_read(image_list[[groups]]),
               scale = 0.2,x = 0.4,y = 0.38) 
  }# group_pack
  plots_NMDS_total_fluvial[[groups]] <- group_pack_plots_NMDS
  plots_NMDS_fluvial_total_model_result[[groups]] <- group_pack_NMDS_model_result
}


# Print NMDS
png(filename ="Figures/NMDS_Diverse.png",
    width =729*12, height =629*10 ,units ="px",res = 300)
grid.arrange(plots_NMDS_total[[1]][[1]],plots_NMDS_total[[1]][[2]],plots_NMDS_total[[1]][[3]],
             plots_NMDS_total[[1]][[4]][[1]],plots_NMDS_total[[1]][[4]][[2]],
             plots_NMDS_total[[1]][[5]][[1]], plots_NMDS_total[[1]][[5]][[2]],
             plots_NMDS_total_fluvial[[1]][[1]],plots_NMDS_total_fluvial[[1]][[2]],
             
             plots_NMDS_total[[2]][[1]],plots_NMDS_total[[2]][[2]],plots_NMDS_total[[2]][[3]],
             plots_NMDS_total[[2]][[4]][[1]],plots_NMDS_total[[2]][[4]][[2]],
             plots_NMDS_total[[2]][[5]][[1]], plots_NMDS_total[[2]][[5]][[2]],
             plots_NMDS_total_fluvial[[2]][[1]],plots_NMDS_total_fluvial[[2]][[2]],
             
             plots_NMDS_total[[3]][[1]],plots_NMDS_total[[3]][[2]],plots_NMDS_total[[3]][[3]],
             plots_NMDS_total[[3]][[4]][[1]],plots_NMDS_total[[3]][[4]][[2]],
             plots_NMDS_total[[3]][[5]][[1]], plots_NMDS_total[[3]][[5]][[2]],
             plots_NMDS_total_fluvial[[3]][[1]],plots_NMDS_total_fluvial[[3]][[2]],
             
             plots_NMDS_total[[4]][[1]],plots_NMDS_total[[4]][[2]],plots_NMDS_total[[4]][[3]],
             plots_NMDS_total[[4]][[4]][[1]],plots_NMDS_total[[4]][[4]][[2]],
             plots_NMDS_total[[4]][[5]][[1]], plots_NMDS_total[[4]][[5]][[2]],
             plots_NMDS_total_fluvial[[4]][[1]],plots_NMDS_total_fluvial[[4]][[2]],
             
             plots_NMDS_total[[5]][[1]],plots_NMDS_total[[5]][[2]],plots_NMDS_total[[5]][[3]],
             plots_NMDS_total[[5]][[4]][[1]],plots_NMDS_total[[5]][[4]][[2]],
             plots_NMDS_total[[5]][[5]][[1]], plots_NMDS_total[[5]][[5]][[2]],
             plots_NMDS_total_fluvial[[5]][[1]],plots_NMDS_total_fluvial[[5]][[2]],
             
             ncol=9,nrow=5, top="NMDS")
dev.off()


# NMDS plots significant Ordisurfs_________________________________________________________ ####

# For Euclidean network
plots_NMDS_sign <- list()
ref_value <- 0

for (groups in 1:5) {
  for (net in 1:5) {
    
    comun_data <- list()
    envir_data <- list()
    netw_centr <- list()
    if(net==4 | net==5){
      #cordenades_xarxes[[net]][which(components(graph.adjacency(MAPS_xarxes[[net]], mode = "undirected",diag = F))$membership==2),]
      
      Gral_coin <- which(1:length(Network_results[[net]])>=length(Network_results[[net]])-54)[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]]  
      
      Swiss <- as.vector(the_plot_db %>% filter(Longitude<10) %>% select(Lake_name))$Lake_name
      comun_data[[1]] <- comm_data[[groups]][rownames(comm_data[[groups]])%in%Swiss,]
      envir_data[[1]] <- env_data[[groups]][rownames(env_data[[groups]])%in%Swiss,]
      
      filter_subgraph<- which(which(components(graph.adjacency(MAPS_xarxes[[net]],mode = "undirected",diag = F))$membership==2)>=c(length(Network_results[[net]])-54))
      filter_subgraph_position <-  which(components(graph.adjacency(MAPS_xarxes[[net]],mode = "undirected",diag = F))$membership==2)[filter_subgraph]
      netw_centr[[1]] <- filter_subgraph_position[which(filter_subgraph_position%in%Gral_coin)]
      
      
      Aust <- as.vector(the_plot_db %>% filter(Longitude>10) %>% select(Lake_name))$Lake_name
      comun_data[[2]] <- comm_data[[groups]][rownames(comm_data[[groups]])%in%Aust,]
      envir_data[[2]] <- env_data[[groups]][rownames(env_data[[groups]])%in%Aust,]
      
      filter_subgraph<- which(which(components(graph.adjacency(MAPS_xarxes[[net]],mode = "undirected",diag = F))$membership==1)>=c(length(Network_results[[net]])-54))
      filter_subgraph_position <-  which(components(graph.adjacency(MAPS_xarxes[[net]],mode = "undirected",diag = F))$membership==1)[filter_subgraph]
      netw_centr[[2]] <- filter_subgraph_position[which(filter_subgraph_position%in%Gral_coin)]
      
    }else{
      comun_data[[1]] <- comm_data[[groups]]
      envir_data[[1]] <- env_data[[groups]]
      netw_centr[[1]] <-which(which(components(graph.adjacency(MAPS_xarxes[[net]], 
                                                               mode = "undirected",diag = F))$membership==1)>=c(length(Network_results[[net]])-54))
    } # We select which "subgraph" we use
    
    SubGraphName <- c("- Switzerland","- Austria")
    for (group_pack in 1:length(comun_data)) {
    if(length(comun_data)==1){SubGraphName <- ""}  
    
    coin <- Network_results[[net]][netw_centr[[group_pack]]] # OLD way [c(length(Network_results[[net]])-54):length(Network_results[[net]])]
    if (length(comun_data)==1) {centr_iso <- cbind(coin[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]])}else{centr_iso <- coin}
    
    dbRDA <- capscale(comun_data[[group_pack]]~
                      envir_data[[group_pack]][,4]+envir_data[[group_pack]][,5]+envir_data[[group_pack]][,6]+envir_data[[group_pack]][,7]+envir_data[[group_pack]][,8],
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
        labs(title = paste(biod_names[groups], net_names[net],SubGraphName[group_pack]),
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
   }# Group packs
  }
  plots_NMDS_sign
}

# For Fluvial network
plots_NMDS_fluvial_sign<- list()
ref_value <- 0
for (groups in 1:5) {
  
  comun_data <- list()
  envir_data <- list()
  netw_centr <- list()
  
  #cordenades_xarxes[[net]][which(components(graph.adjacency(MAPS_xarxes[[net]], mode = "undirected",diag = F))$membership==2),]
  
  Gral_coin <- V(GRAPH_xarxes_fluvial[[1]])[c(all_lakes_BASINS_fluvial[[1]][correspondence_BASINS_fluvial[[1]]])[coincidence_values[[groups]][which(coincidence_values[[groups]]>0)]]]
  
  Swiss <- as.vector(the_plot_db %>% filter(Longitude<10) %>% select(Lake_name))$Lake_name
  comun_data[[1]] <- comm_data[[groups]][rownames(comm_data[[groups]])%in%Swiss,]
  envir_data[[1]] <- env_data[[groups]][rownames(env_data[[groups]])%in%Swiss,]
  out_site <- c()           
  for (ds in 1:length(Gral_coin)) {
    out_temp<- which(which(components(GRAPH_xarxes_fluvial[[1]])$membership==1)==Gral_coin[ds])
    if (length(out_temp)==0) {out_temp <- 0}
    out_site[ds] <- out_temp}
  netw_centr[[1]] <-   which(components(GRAPH_xarxes_fluvial[[1]])$membership==1)[out_site[which(out_site>0)]]
  
  Aust <- as.vector(the_plot_db %>% filter(Longitude>10) %>% select(Lake_name))$Lake_name
  comun_data[[2]] <- comm_data[[groups]][rownames(comm_data[[groups]])%in%Aust,]
  envir_data[[2]] <- env_data[[groups]][rownames(env_data[[groups]])%in%Aust,]
  out_site <- c()           
  for (ds in 1:length(Gral_coin)) {
    out_temp<- which(which(components(GRAPH_xarxes_fluvial[[1]])$membership==2)==Gral_coin[ds])
    if (length(out_temp)==0) {out_temp <- 0}
    out_site[ds] <- out_temp}
  netw_centr[[2]] <-   which(components(GRAPH_xarxes_fluvial[[1]])$membership==2)[out_site[which(out_site>0)]]
  
  SubGraphName <- c("- Switzerland","- Austria")
  for (group_pack in 1:length(comun_data)) {
  if(length(comun_data)==1){SubGraphName <- ""}    
    
  coin <- Fluvial_network_results[[1]][netw_centr[[group_pack]]]
  centr_iso <- coin 
  
  dbRDA <- capscale(comun_data[[group_pack]]~
                    envir_data[[group_pack]][,4]+envir_data[[group_pack]][,5]+envir_data[[group_pack]][,6]+envir_data[[group_pack]][,7]+envir_data[[group_pack]][,8],
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
      labs(title = paste(biod_names[groups], "Fluvial", SubGraphName[group_pack]),
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
  }# group_pack
  plots_NMDS_fluvial_sign
}

png(filename ="Figures/NMDS_Diverse_Sign.png",
    width =750*8 ,height =629*5 ,units ="px",res = 250)
grid.arrange(
  arrangeGrob(
    arrangeGrob(plots_NMDS_sign[[5]],top=textGrob("A)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
    arrangeGrob(plots_NMDS_sign[[6]],top=textGrob("G)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
    arrangeGrob(plots_NMDS_sign[[10]],top=textGrob("H)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
    arrangeGrob(plots_NMDS_sign[[15]],top=textGrob("M)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
    nrow=4),# 6.5
  arrangeGrob(
    arrangeGrob(plots_NMDS_sign[[4]],top=textGrob("B)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
    arrangeGrob(ggplot() + geom_blank()),
    arrangeGrob(ggplot() + geom_blank()),
    arrangeGrob(plots_NMDS_sign[[14]],top=textGrob("N)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
    nrow=4), # 65
  arrangeGrob(
    arrangeGrob(plots_NMDS_sign[[3]],top=textGrob("C)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
    arrangeGrob(ggplot() + geom_blank()),
    arrangeGrob(plots_NMDS_sign[[9]],top=textGrob("I)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
    arrangeGrob(plots_NMDS_sign[[13]],top=textGrob("O)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
    nrow=4), # 100
  arrangeGrob(
    arrangeGrob(plots_NMDS_sign[[2]],top=textGrob("D)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
    arrangeGrob(ggplot() + geom_blank()),
    arrangeGrob(plots_NMDS_sign[[8]],top=textGrob("J)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
    arrangeGrob(plots_NMDS_sign[[12]],top=textGrob("P)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
    nrow=4), # 325
  arrangeGrob(
    arrangeGrob(plots_NMDS_sign[[1]],top=textGrob("E)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
    arrangeGrob(ggplot() + geom_blank()),
    arrangeGrob(plots_NMDS_sign[[7]],top=textGrob("K)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
    arrangeGrob(plots_NMDS_sign[[11]],top=textGrob("Q)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
    nrow=4), # 650
  arrangeGrob(
    arrangeGrob(plots_NMDS_fluvial_sign[[1]],top=textGrob("F)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
    arrangeGrob(ggplot() + geom_blank()),
    arrangeGrob(plots_NMDS_fluvial_sign[[2]],top=textGrob("L)",x = 0.2, hjust = 0,gp=gpar(fontsize=30,font=2))),
    arrangeGrob(ggplot() + geom_blank()),
    nrow=4), # Fluvial
  ncol=6,nrow=1, top="")
dev.off()

# GAM NMDS models result in table format - Supplementary like______________________________ ####
biod_names <- c("S16","S18","Phy","Zoo", "18S Zoo")
Names_Networks <- c("650 km", "325 km","100 km","65 km","6.5 km", "Fluvial")

# Add the fluvial as a sixth network
for (t in 1:5) {
  plots_NMDS_total_model_result[[t]][[6]] <- plots_NMDS_fluvial_total_model_result[[t]]
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
    if(netw==4 | netw==5){sign[5] <-plots_NMDS_total_model_result[[group]][[netw]][[1]][[8]]}else{sign[5] <-plots_NMDS_total_model_result[[group]][[netw]][[8]]}
    if(netw==4 | netw==5){sign[6] <-plots_NMDS_total_model_result[[group]][[netw]][[2]][[8]]}else{sign[6] <-1}
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
  sign[5] <-plots_NMDS_fluvial_total_model_result[[group]][[1]][[8]]
  sign[6] <-plots_NMDS_fluvial_total_model_result[[group]][[2]][[8]]
  flu_sign_groups[[group]] <- sign
}


Names_Variab <- c("Species richness", "LCBD", "Replacement", "Richness difference","dbRDA","dbRDA")
plots_significance <- list()
for (variable in 1:6) {
  
  max_netw <- cbind(c(sign_groups[[1]][[1]][[variable]], sign_groups[[2]][[1]][[variable]], sign_groups[[3]][[1]][[variable]], 
                      sign_groups[[4]][[1]][[variable]], sign_groups[[5]][[1]][[variable]]),
                    rep("650 km", 5), 
                    c("S16","S18","Phy","Zoo", "S18 Zoo"))
  
  mid_netw <- cbind(c(sign_groups[[1]][[2]][[variable]], sign_groups[[2]][[2]][[variable]], sign_groups[[3]][[2]][[variable]], 
                      sign_groups[[4]][[2]][[variable]], sign_groups[[5]][[2]][[variable]]),
                    rep("325 km", 5), 
                    c("S16","S18","Phy","Zoo", "S18 Zoo"))
  
  mid_mid_netw <- cbind(c(sign_groups[[1]][[3]][[variable]], sign_groups[[2]][[3]][[variable]], sign_groups[[3]][[3]][[variable]], 
                          sign_groups[[4]][[3]][[variable]], sign_groups[[5]][[3]][[variable]]),
                        rep("100 km", 5), 
                        c("S16","S18","Phy","Zoo", "S18 Zoo"))
  
  small_netw <- cbind(c(sign_groups[[1]][[4]][[variable]], sign_groups[[2]][[4]][[variable]], sign_groups[[3]][[4]][[variable]], 
                        sign_groups[[4]][[4]][[variable]], sign_groups[[5]][[4]][[variable]]),
                      rep("65 km", 5), 
                      c("S16","S18","Phy","Zoo", "S18 Zoo"))
  
  min_netw <- cbind(c(sign_groups[[1]][[5]][[variable]], sign_groups[[2]][[5]][[variable]], sign_groups[[3]][[5]][[variable]], 
                      sign_groups[[4]][[5]][[variable]], sign_groups[[5]][[5]][[variable]]),
                    rep("6.5 km", 5), 
                    c("S16","S18","Phy","Zoo", "S18 Zoo"))
  
  fluv_netw <- cbind(c(flu_sign_groups[[1]][[variable]], flu_sign_groups[[2]][[variable]], flu_sign_groups[[3]][[variable]], 
                       flu_sign_groups[[4]][[variable]], flu_sign_groups[[5]][[variable]]),
                     rep("Fluvial", 5), 
                     c("S16","S18","Phy","Zoo", "S18 Zoo"))
  
  dataset_pval <- as.data.frame(rbind(max_netw,mid_netw,mid_mid_netw,small_netw,min_netw,fluv_netw))
  
  colnames(dataset_pval) <- c("pval","Network","Group")
  dataset_pval$pval <-as.numeric(dataset_pval$pval)
  dataset_pval$Network <- factor(dataset_pval$Network,
                                 levels = c("650 km", "325 km","100 km","65 km","6.5 km", "Fluvial"))
  dataset_pval$Group <- factor(dataset_pval$Group,
                               levels = c("S16","S18","Phy","Zoo", "S18 Zoo"))
  significants <- rep(">0.05",5*6)
  significants[which(dataset_pval$pval<0.05)]<- "<0.05"
  dataset_pval$Sign <- factor(significants)
  if(variable==5){dataset_pval_Temp <- dataset_pval}
  if(variable==6){dataset_pval <- bind_rows(dataset_pval_Temp,dataset_pval)}

  color_groups <- as.character(CUNILLERA_cols("yellow","blue","green","red","cyan"))
  
#  if (variable==4) {
#plots_significance[[variable]] <-  ggplot(dataset_pval, aes(x=Network, y=as.numeric(pval)))+
#  geom_abline(slope = 0,intercept = 0.05, colour="black", linetype=2,size=1)+
#  geom_jitter(aes(fill=Group, alpha=Sign, size=Sign),shape=21,width = 0.5)+
#  scale_x_discrete(limits=c("6.5 km","65 km", "100 km","325 km","650 km", "Fluvial"))+
#  scale_alpha_manual(values = c(0.3))+
#  scale_size_manual(values = c(2))+
#  scale_fill_manual(values=c(color_groups[1],color_groups[2],color_groups[3],
#                             color_groups[4],color_groups[5],color_groups[6]))+
#  scale_y_continuous(expand = c(0.2,0.01),
#                     breaks =c(0.05,0.2,0.4,0.6,0.8,1) )+
#  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5), size=1, colour="grey70")+
#  labs(title="", subtitle=Names_Variab[variable])+ylab("p-values")+xlab("")+
#  theme_classic()+theme(legend.position = "none")    
#  }else{
plots_significance[[variable]] <-  ggplot(dataset_pval, aes(x=Network, y=as.numeric(pval)))+
  geom_abline(slope = 0,intercept = 0.05, colour="black", linetype=2,size=1)+
  geom_jitter(aes(fill=Group, alpha=Sign, size=Sign),shape=21,width = 0.3)+
  scale_x_discrete(limits=c("6.5 km","65 km", "100 km","325 km","650 km", "Fluvial"))+
  scale_alpha_manual(values = c(0.9,0.3))+
  scale_size_manual(values = c(5,2))+
  scale_fill_manual(values=c(color_groups[1],color_groups[2],color_groups[3],
                             color_groups[4],color_groups[5],color_groups[6]))+
  scale_y_continuous(expand = c(0.2,0.01),
                     breaks =c(0.05,0.2,0.4,0.6,0.8,1) )+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5), size=1, colour="grey30")+
  labs(title="",subtitles=Names_Variab[variable])+ylab("p-values")+xlab("")+
  theme_classic()+theme(legend.position = "none")
#  }
  
  
}


plot_for_legend <-  ggplot(dataset_pval, aes(x=Network, y=as.numeric(pval)))+
  geom_abline(slope = 0,intercept = 0.05, colour="black", linetype=2,size=1)+
  geom_jitter(aes(fill=Group, alpha=Sign, size=Sign),shape=21,width = 0.5)+
  scale_x_discrete(limits=c("6.5 km","65 km", "100 km","325 km","650 km", "Fluvial"))+
  scale_alpha_manual(values = c(0.9,0.3))+
  scale_size_manual(values = c(7,2))+
  scale_fill_manual(values=c(color_groups[1],color_groups[2],color_groups[3],
                             color_groups[4],color_groups[5],color_groups[6]))+
  scale_y_continuous(expand = c(0.2,0.01),
                     breaks =c(0.05,0.2,0.4,0.6,0.8,1) )+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5), size=1, colour="grey70")+
  guides(fill = guide_legend(override.aes = list(size = 7)))+
  labs(title="", subtitle = "dbRDA ordisurf")+ylab("p-values")+xlab("")+
  theme_classic()+theme(legend.position = "bottom",
                        legend.key.size = unit(.25, "cm"),
                        legend.box="vertical")


legend_try <- cowplot::get_legend(plot_for_legend)


png(filename =paste("Figures/All_Significance_Groups.png"),
    width =400*3 ,height =650*6 ,units ="px",res = 300)
grid.arrange(plots_significance[[1]],
             plots_significance[[3]],plots_significance[[4]],plots_significance[[2]],
             plots_significance[[6]],
             legend_try,
             ncol=1,nrow=6, top="Significance values")
dev.off()      

save.image("Database.RData")




# OLD STUFF 
#Summary plots NMDS _________________________________________________________ ##

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

Group_names <- c("S16","S18","Phy","Zoo", "S18 Zoo")
output <- data.frame()
for (group in 1:length(Group_names)) {
  out <- data.frame("pval"=c(unlist(TOTAL_pval_summary_NMDS[[group]]),unlist(TOTAL_FLU_pval_summary_NMDS[[group]])),
                    "Network"=c("650 km", "325 km","100 km","65 km","6.5 km", "Fluvial"),
                    "Group"=rep(Group_names[group],6)) %>% 
    mutate(Sign=ifelse(pval<0.05,"<0.05",">0.05"))
  output <- bind_rows(output,out)
}
dataset_pval <- output
dataset_pval$Network <- factor(dataset_pval$Network, levels = c("650 km", "325 km","100 km","65 km","6.5 km", "Fluvial"))
dataset_pval$Group <- factor(dataset_pval$Group,levels = c("S16","S18","Phy","Zoo", "S18 Zoo"))

color_groups <- as.character(CUNILLERA_cols("yellow","blue","green","red","cyan"))


