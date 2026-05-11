### ------------------------------------------------------------------------ ### 
###  Ensemble - Main script
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list = ls()); graphics.off(); gc()

# Set data path
path_data <- paste0('C:/Users/claer14/OneDrive - University of Cambridge/III. Research/Data/')
path_git <- paste0('C:/Users/claer14/Documents/GitHub/hisp-climate-health-forecast/')
setwd(path_git)

start <- Sys.time()

#### DEPENDENCIES ####

# Install packages if not already in library & load
pkgs <- c('tidyverse', 'plyr', 'dplyr', 'ggplot2', 'zoo', 'Rcpp', 'dfoptim', 'forecast', 'MASS',
          'msm', 'mgcv', 'RColorBrewer', 'gridExtra', 'cowplot')
# install.packages(pkgs[!pkgs %in% installed.packages()[, 'Package']])
invisible(lapply(pkgs, library, character.only = TRUE))

# Load utilities and set output paths
source('3_Ensemble/3_Utilities.R')
# file_import <- 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_ensemble_smooth'

####  MAIN CODE  #####

# 1. Load data and set parameters ####

# Load smoothed data
load(paste0(path_project, file_import, '.RData'))
data_rt <- df; rm(df)

# Set primary outcome
outcome <- 'total_patients'
data_rt <- data_rt %>% rename_with(~ gsub(outcome, 'outcome', .x), contains(outcome))

str(data_rt)
min(data_rt$date); max(data_rt$date)
length(unique(data_rt$location))
length(unique(data_rt$date))
length(unique(data_rt$date_analysis))

# Number of removed data points (not consolidated)
remove_last_n <- 2

# Set maximum prediction horizon
max_prediction_horizon<-14+remove_last_n

# List of locations (hospitals in this case)
# NB: Something weird with Hamirpur
locations_hosp <- unique(data_rt$location)[1:2]#[-length(unique(data_rt$location))]
locations <-c('metropolitan', locations_hosp)

# List of models
models <- c('baseline', 'ARIMA1', 'GAM1', 'ARIMA2', 'GAM2', 'TSGC')

# List of models to include in the ensemble
models_ensemble<- c('ARIMA1', 'GAM2', 'TSGC')


# 2. Evaluation period ####

period<-1
dates_analysis<-seq(as.Date('2017-01-01'),as.Date('2017-03-01'),by='day') #'2017-06-01'

# Run component models
for (dates_analysis in models) {

  print(model)
  source(paste0('3_Ensemble/Component_models/', model,'.R'))

}

# Plot results (F1 and SF4)
source('3_Ensemble/3B_Figures_period1.R', echo=TRUE)


# 3. Test period ####

# period<-2
# dates_analysis<-seq(as.Date('2017-06-02'),as.Date('2017-09-30'),by='day') #'2017-09-30'
# 
# # Run component models
# for (model in models) {
# 
#   print(model)
#   source(paste0('3_Ensemble/Component_models/', model,'.R'))
# 
# }


# Compute ensemble forecasts 
# for hospital admissions and project the 3 other targets (ICU admissions, general ward beds and ICU beds) 

source('3_Ensemble/3C_Ensemble_model.R', echo=TRUE)


# Plots  results

source('3_Ensemble/3D_Figures_period2.R')


# End timer: 1 hours 57 minutes 58 seconds (31/03/2026) ####
end <- Sys.time()
dt <- round(as.numeric(end - start, units = 'secs'))
hours   <- dt %/% 3600
minutes <- (dt %% 3600) %/% 60
seconds <- dt %% 60
sprintf('%d hours %d minutes %.d seconds', hours, minutes, seconds)


