### ------------------------------------------------------------------------ ### 
###  TSGC - Sensitivity analysis (HISP: Dharamshala)
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list = ls()); graphics.off(); gc()

# Set data path
setwd('C:/Users/claer14/Documents/GitHub/hisp-climate-health-forecast/')

#### DEPENDENCIES ####
source('3_TSGC/TSGC - Dependencies.R')
setwd(path_git)
file_import <- 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide'

# Function to iteratively run models and get results
run_models <- function(Y, y, exposures, controls, new.data, estimation.date.end, 
                       estimation.window, n.ahead, q, confidence.level){
  
  # Set parameters
  estimation.date.start <- estimation.date.end %m-% months(estimation.window)
  plt.length <- 30
  plt.start.date <- estimation.date.end - plt.length
  
  # Set sample period
  y.cum <- Y[estimation.date.start <= index(Y) & index(Y) <= estimation.date.end]
  # Get growth rate
  y.lg <- tsgc::df2ldl(y.cum[,outcome])
  # Get data post-sample period & convert to growth rates (for accuracy evaluation)
  y.eval <- Y[,outcome] %>%
    subset(index(.) > estimation.date.end) %>%
    tsgc::df2ldl() 
  # Get results
  results <- get_results(Y, y, exposures, controls, y.lg, y.cum, y.eval, new.data, 
                         n.ahead, q, confidence.level, sea.period, sea.type, 
                         plt.start.date)
  
  return(results)
}

####  MAIN CODE  #####

# Load data
load(paste0('3_TSGC/', file_import, '.RData'))


# 1. Set base model parameters ####

df <- df %>% filter(location=='Dharamshala')

# Variables
outcome <- 'total_patients_c' 
controls <- c('closed')
exposures <- c()
plot_components <- c(paste(controls, collapse = ' + '), paste(exposures, collapse = ' + '))

# Create xts object and subset relevant variables as numeric
Y <- xts(subset(df, select=c(outcome, exposures, controls)),  order.by=df$date)

# Parameters to test
q <- 0.005
estimation.window <- 3
estimation.date.end <- as.Date('2016-08-15')
n.ahead <- 10
confidence.level <- 0.68

# NOT TO CHANGE (global variables required)
# Set seasonal parameters
sea.type <- 'trigonometric'
sea.period <- 7

# 2. Sensitivity test: signal-to-noise ratio ####

# Test q 0.000-0.1 (increments of 0.005)
name <- 'q: '
test <- seq(0.000, 0.1, by=0.005) # 0.005 is most popular
signal <- data.frame()

for (t in test){
  
  # plot_title <- paste0(name, t, ' (', outcome ,' + ', paste(plot_components[plot_components != ''], collapse = ' + '), ')')
  plot_title <- paste0(name, t)
  
  # TEST
  q <- t
  # New data (to forecast)
  new.period <- paste0(estimation.date.end+1, '/', estimation.date.end+n.ahead)
  new.data <- Y[new.period,c(exposures, controls)]
  # Get results
  results <- run_models(Y, y, exposures, controls, new.data, 
                        estimation.date.end, estimation.window, n.ahead, 
                        q, confidence.level)
  results$plots
  signal <- rbind(signal, results$criteria)
  
}

# SET Q
signal
q <- 0.01

# 3. Sensitivity test: estimation window ####

# Test estimation window
name <- 'Estimation window: '
test <- seq(1, 12, by=1)
estimation <- data.frame()

for (t in test){
  
  # plot_title <- paste0(name, t, ' (', outcome ,' + ', paste(plot_components[plot_components != ''], collapse = ' + '), ')')
  plot_title <- paste0(name, t)
  
  # TEST
  estimation.window <- t
  # New data (to forecast)
  new.period <- paste0(estimation.date.end+1, '/', estimation.date.end+n.ahead)
  new.data <- Y[new.period,c(exposures, controls)]
  # Get results
  results <- run_models(Y, y, exposures, controls, new.data, 
                        estimation.date.end, estimation.window, n.ahead, 
                        q, confidence.level)
  results$plots
  estimation <- rbind(estimation, results$criteria)
  
}

# SET ESTIMATION WINDOW
estimation
estimation.window <- 2


# 4. Sensitivity test: forecast horizon ####

# Test n.ahead 1-20 (increments of 1) 
name <- 'Forecast horizon: '
test <- seq(1, 20, by=1)
horizon <- data.frame()

for (t in test){
  
  # plot_title <- paste0(name, t, ' (', outcome ,' + ', paste(plot_components[plot_components != ''], collapse = ' + '), ')')
  plot_title <- paste0(name, t)
  
  # TEST
  n.ahead <- t
  # New data (to forecast)
  new.period <- paste0(estimation.date.end+1, '/', estimation.date.end+n.ahead)
  new.data <- Y[new.period,c(exposures, controls)]
  # Get results
  results <- run_models(Y, y, exposures, controls, new.data, 
                        estimation.date.end, estimation.window, n.ahead, 
                        q, confidence.level)
  results$plots
  horizon <- rbind(horizon, results$criteria)
  
}

# SET N.AHEAD
horizon
n.ahead <- 14

# 5. Sensitivity test: exposure threshold ####

# IMD THRESHOLD Variables
exposures <- c('mx2t_IMD_b', 'mn2t_IMD_b', 'tp_IMD_b', 'fg10_IMD_b')
plot_components <- c(paste(controls, collapse = ' + '), paste(exposures, collapse = ' + '))
plot_title <- 'IMD thresholds'
# Get results
Y <- xts(subset(df, select=c(outcome, exposures, controls)),  order.by=df$date)
new.period <- paste0(estimation.date.end+1, '/', estimation.date.end+n.ahead)
new.data <- Y[new.period,c(exposures, controls)]
results <- run_models(Y, y, exposures, controls, new.data, 
                      estimation.date.end, estimation.window, n.ahead, 
                      q, confidence.level)
IMD <- results$criteria

# CLIMATE PERCENTILE Variables
exposures <- c('mx2t_p90', 'mn2t_p10', 'tp_p90', 'fg10_p90')
plot_components <- c(paste(controls, collapse = ' + '), paste(exposures, collapse = ' + '))
plot_title <- 'Climatology thresholds'
# Get results
Y <- xts(subset(df, select=c(outcome, exposures, controls)),  order.by=df$date)
new.period <- paste0(estimation.date.end+1, '/', estimation.date.end+n.ahead)
new.data <- Y[new.period,c(exposures, controls)]
results <- run_models(Y, y, exposures, controls, new.data, 
                      estimation.date.end, estimation.window, n.ahead, 
                      q, confidence.level)
clim <- results$criteria

exposures <- c('cold20_lag0', 'cold20_lag3', 'prec80_lag0', 'prec80_lag10')
plot_components <- c(paste(controls, collapse = ' + '), paste(exposures, collapse = ' + '))
plot_title <- 'DLNM thresholds'
# Get results
Y <- xts(subset(df, select=c(outcome, exposures, controls)),  order.by=df$date)
new.period <- paste0(estimation.date.end+1, '/', estimation.date.end+n.ahead)
new.data <- Y[new.period,c(exposures, controls)]
results <- run_models(Y, y, exposures, controls, new.data, 
                      estimation.date.end, estimation.window, n.ahead, 
                      q, confidence.level)
dlnm <- results$criteria

thresh <- rbind(IMD, clim, dlnm)
thresh <- thresh[order(thresh$BIC), ]
thresh

# 6. Sensitivity test: exposure lags ####

# NB: Lower AIC/BIC is better; Higher MLE is better
# Select exposures a priori based on Indian weather forecast system
w_sel_col <- c('cold20_lag0', 'cold20_lag3', 'prec80_lag0', 'prec80_lag10')


# Check ACF & PACF
#### ACF: Most have a strong sinusoidal pattern (indicating seasonal trend) and 
#### tapers off slowly (indicating non-stationarity).
#### PACF: Sharp spikes at lags 1-7, which may be important for modelling
# lag_plot <- list()
# for (c in w_sel_col) {
#   vec <- na.omit(df[[c]])
#   
#   # Skip if vector is empty after NA removal
#   if (length(vec) == 0) next
#   
#   # Skip if all values are the same
#   if (var(vec) == 0) {
#     message(paste("Skipping", c, "- no variation in data"))
#     next
#   }
#   
#   ts_data <- ts(na.omit(df[[c]]))
#   
#   par(mfrow = c(1, 2))  # Side-by-side ACF/PACF
#   acf(ts_data, lag.max = 50, main = "ACF")
#   pacf(ts_data, lag.max = 50, main = "PACF")
#   par(mfrow = c(1, 1))  # Reset layout
#   
#   lag_plot[[c]] <- recordPlot()
#   
# }
# 
# # Check model comparison criteria
# lags <- data.frame()
# for (w in w_sel_col){
#   
#   # cat(paste0(toupper(w), ':\n'))
#   # Variables
#   exp <- c(w, grep(paste0('^', w, '_lag(\\d{1,2})$'), names(df), value = TRUE))
#   res1 <- data.frame()
#   res2 <- data.frame()
#   for (e in exp){
#     
#     # !!! Ensure that predict_all_MULTI correctly reads exposure Q matrix length 
#     exposures <- e
#     
#     plot_title <- e
#     # Create xts object and subset relevant variables as numeric
#     Y <- xts(subset(df, select=c(outcome, e, controls)),  order.by=df$date)
#     
#     # New data (to forecast)
#     new.period <- paste0(estimation.date.end+1, '/', estimation.date.end+n.ahead)
#     new.data <- Y[new.period,c(e, controls)]
#     # Get results
#     results <- run_models(Y, y, e, controls, new.data, estimation.date.end, 
#                           estimation.window, n.ahead, q, confidence.level)
#     res1 <- rbind(res1, results$criteria)
#     res2 <- rbind(res2, results$estimates[1,])
#     # print(results$criteria)
#     # print(results$estimates[1:(nrow(results$estimates) - 9), ])
#     
#   }
#   # print(res1)
#   lags <- rbind(lags, res1)
#   # cat('CRITERIA:\n')
#   # print(res1[order(res1$BIC),])
#   # cat('ESTIMATES:\n')
#   # print(res2[order(res1$BIC),])
# }
# 
# lags 
# # # # Best exposures for IMD and other based on tool (using total_patient_c outcome)
# # # exposures <- c('mx2t_IMD_b_lag1', 'mn2t_IMD_b_lag6', 'tp_IMD_b_lag4', 'fg10_IMD_b_lag4') # tp_IMD_lag21
# # # exposures <- c('mx2t_p90_lag6', 'mn2t_p10_lag5', 'tp_p90_lag1', 'fg10_p90_lag3') # mx2t_p90_lag24; mx2t_p90_lag9; mn2t_p10_lag1; tp_p90_lag25; fg10_p90_lag21  

#### FORMAT & EXPORT ####

# Create tables for export
# NB: Run for each outcome of interest, then cbind and export
total_patients_all <- rbind(signal,estimation,horizon,IMD,clim, dlnm)
# total_patients_lag <- lags

general_patients_all <- rbind(signal,estimation,horizon,IMD,clim, dlnm)[-1]
# general_patients_lags <- lags[-1]

result_all <- round_numeric(cbind(total_patients_all, general_patients_all), 2)
# result_lags <- round_numeric(cbind(total_patients_lag, general_patients_lags), 2)

# Fix columns to your data frame
# All params
split_cols <- do.call(rbind, strsplit(result_all$Model, ": "))
result_all$Parameter <- split_cols[, 1]
result_all$Test <- split_cols[, 2]
result_all$Model <- NULL
names(result_all) <- make.unique(names(result_all))
result_all <- result_all %>% dplyr::select(Parameter, Test, everything())
names(result_all) <- gsub(".1", "", names(result_all))
result_all$Test[54:56] <- NA  

# Save tables
write.xlsx(result_all, file = paste0(path_report, "III. tsgc/tab_dharam_sens_param.xlsx"))
# write.xlsx(result_lags, file = paste0(path_report, "III. tsgc/tab_dharam_sens_lag.xlsx"))


