### ------------------------------------------------------------------------ ### 
###  Ensemble - Main script
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list = ls()); graphics.off(); gc()

# Set data path
path_data <- paste0('C:/Users/claer14/OneDrive - University of Cambridge/III. Research/Data/')
path_git <- paste0('C:/Users/claer14/Documents/GitHub/hisp-climate-health-forecast/')
setwd(path_git)

#### DEPENDENCIES ####

# Install packages if not already in library & load
pkgs <- c("tidyverse", 'plyr', "dplyr", "ggplot2", "zoo", 'Rcpp', 'dfoptim', 'forecast', 'MASS',
          'msm', 'mgcv', 'RColorBrewer', 'gridExtra', 'cowplot')
install.packages(pkgs[!pkgs %in% installed.packages()[, "Package"]])
invisible(lapply(pkgs, library, character.only = TRUE))

# Load utilities and set output paths
source('3_Ensemble/Ensemble - Utilities.R')
file_import <- 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_ensemble_smooth'

####  MAIN CODE  #####

# 1. Load data and set parameters ####

# Load smoothed data
load(paste0(path_project, file_import, '.RData'))
data_rt <- df; rm(df)
data_rt <- data_rt %>% dplyr::rename('region'='location',
                              'iHosp_smooth'='total_patients_sm',
                              'iHosp_smooth_lwr'='total_patients_lwr',
                              'iHosp_smooth_upr'='total_patients_upr')
min(data_rt$date); max(data_rt$date)
length(unique(data_rt$region))
length(unique(data_rt$date))
length(unique(data_rt$date_analysis))

# Number of removed data points (not consolidated)
remove_last_n <- 2

# Set maximum prediction horizon
max_prediction_horizon<-14+remove_last_n

# List of regions (hospitals in this case)
# NB: Something weird with Hamirpur
regions_metro <- unique(data_rt$region)[1:2]
regions <-c("metropolitan", regions_metro)

# List of models
models <- c("baseline", 'ARIMA1', 'GAM1')  #"Const2","Const7","PL","GAM1", "ARDL","MLR","ARIMA2","GAM2","RF","BRT")

# List of models to include in the ensemble
# models_ensemble<- c("ARDL","MLR","GAM2","ARIMA2","RF","BRT")

# 2. Evaluation period ####

period<-1
dates_analysis<-seq(as.Date("2017-01-01"),as.Date("2017-06-01"),by="day") 

# Run component models
for (model in models) {

  print(model)
  source(paste0("3_Ensemble/Component_models/", model,".R"))

}

# Plot results (F1 and SF4)
source("3_Ensemble/Ensemble - Figures_period1.R")
