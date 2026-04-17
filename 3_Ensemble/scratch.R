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
file_import <- 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_ensemble_smooth'

####  MAIN CODE  #####

# 1. Load data and set parameters ####

# Load smoothed data
load(paste0(path_project, file_import, '.RData'))
data_rt <- df; rm(df)

data_rt <- data_rt %>% dplyr::rename('region'='location') %>%
  rename_with(~ gsub('total_patients', 'iHosp_smooth', .x), contains('total_patients')) %>%
  dplyr::rename('iHosp_r' = 'iHosp_smooth_r') %>%
  rename_with(~ gsub('F', 'iICU_smooth', .x), contains('F')) %>%
  rename_with(~ gsub('M', 'inHosp_smooth', .x), contains('M')) %>%
  rename_with(~ gsub('age25_64', 'inICU_smooth', .x), contains('age25_64')) %>%
  mutate(iICU_smooth = as.integer(iICU_smooth))
dates_analysis<-seq(as.Date('2017-01-01'),as.Date('2017-06-01'),by='day') #'2017-06-01'

# Example data frame: df with columns growth_rate, residential, transit_stations,
# abs_humidity, gr_nb_positive, gr_prop_positive, etc.
iDate <- 1
dat_train <- data_rt |>
  filter(date_analysis == dates_analysis[iDate]) |>
  mutate(iHosp_r_covar = iHosp_r)

# GAM
fit_gam <- mgcv::gam(iHosp_r ~ s(t2m) + s(tp), data = dat_train)


library(gratia)
draw(fit_gam)  # quick version


# MLR
fit_mlr <- lm(iHosp_smooth ~ t2m + tp, data = data_rt)