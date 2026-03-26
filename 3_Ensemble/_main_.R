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
source('3_Ensemble/3_Utilities.R')
file_import <- 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_ensemble_smooth'

####  MAIN CODE  #####

# 1. Load data and set parameters ####

# Load smoothed data
load(paste0(path_project, file_import, '.RData'))
data_rt <- df; rm(df)
data_rt <- data_rt %>% dplyr::rename('region'='location') %>%
  rename_with(~ gsub("total_patients", "iHosp", .x), contains("total_patients")) %>%
  rename_with(~ gsub("F", "iICU", .x), contains("F")) %>%
  rename_with(~ gsub("M", "inHosp", .x), contains("M")) %>%
  rename_with(~ gsub("age25_64", "inICU", .x), contains("age25_64")) %>%
  mutate(inICU = as.integer(inICU))

str(data_rt)
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
models <- c("baseline", 'ARIMA1', 'GAM1', 'ARIMA2', 'GAM2')

# List of models to include in the ensemble
models_ensemble<- c('ARIMA2', 'GAM2')

# Checks ####

data_rt %>%
  filter(date_analysis == max(date_analysis),
         region=='metropolitan') %>%
  arrange(region, date) %>%
  ggplot(aes(x = date)) +
  # Reference line at r = 0 (stable)
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#888780", linewidth = 0.4) +
  # Growth rate line
  geom_line(aes(y = iHosp), colour = "#378ADD", linewidth = 0.7, na.rm = TRUE) +
  facet_wrap(~ region, scales = "free_y") +
  labs(
    title    = "Log growth rate of smoothed patients by location",
    subtitle = paste("Snapshot:", format(max(data_rt$date_analysis), "%d %b %Y")),
    x        = NULL,
    y        = "Log growth rate"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.subtitle      = element_text(colour = "grey50"),
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(face = "bold")
  )

# 2. Evaluation period ####

period<-1
dates_analysis<-seq(as.Date("2017-01-01"),as.Date("2017-06-01"),by="day") 

# # Run component models
# for (model in models) {
# 
#   print(model)
#   source(paste0("3_Ensemble/Component_models/", model,".R"))
# 
# }
# 
# Plot results (F1 and SF4)
source("3_Ensemble/3B_Figures_period1.R")


# 3. Test period ####

period<-2
dates_analysis<-seq(as.Date("2017-06-02"),as.Date("2017-09-30"),by="day") 

# # Run component models
# for (model in models) {
#   
#   print(model)
#   source(paste0("3_Ensemble/Component_models/", model,".R"))
# 
# }


# Compute ensemble forecasts 
# for hospital admissions and project the 3 other targets (ICU admissions, general ward beds and ICU beds) 

# source("3_Ensemble/3C_Ensemble_model.R", echo=TRUE)


# Plots  results  ----------------------------------------

source("3_Ensemble/3D_Figures_period2.R")

