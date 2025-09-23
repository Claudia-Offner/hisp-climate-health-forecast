### ------------------------------------------------------------------------ ### 
###  TSGC Multivariate Model - Results                                       ###
### ------------------------------------------------------------------------ ### 

####  DEPENDENCIES  #####

# Packages
library(openxlsx); library(stringr)

# Turn off scientific notation
options(scipen=999)

# Set data path
user <- 'claer14'
file_import <- 'ind_5km_BigThaw-HISP-chamba_2015-2023_hosp_wide'
path1 <- paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/')
path2 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
path3 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Chapter 3 - Real time forecasts/')
path4 <- paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/3_TSGC/- Reporting/- figures/')
setwd(path1)

### Get model functions (make sure wd is in repo folder)
source(paste0(path1, '3_TSGC/TSGC_Model - Functions.R'))

run_models <- function(Y, y, exposures, controls, new.data,
                       estimation.date.end, estimation.window, 
                       n.ahead, q, confidence.level){
  
  # Set parameters
  estimation.date.start <- estimation.date.end %m-% months(estimation.window)
  plt.length <- 31
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

# Set identified parameters
q <- 0.005
estimation.window <- 3 # months
estimation.date.end <- as.Date('2016-08-24')
n.ahead <- 10
confidence.level <- 0.68
sea.type <- 'trigonometric'
sea.period <- 7


# 1. Model Results (exposures): Coefficients ####

outcomes <- c('total_patients_c', 
              'F_c','M_c','under_5_c','age5_14_c','age15_24_c','age25_64_c','over_65_c',
              'general_c','ophthalmology_c','pediatrics_c','ent_c')
controls <- c('closed')
exposures <- c('mx2t_IMD_b_lag1', 'tp_IMD_b_lag1', 'fg10_IMD_b_lag1')

# Outputs
esti_df <- data.frame()
comp_df <- data.frame()
plot_list <- list()

for (i in outcomes) {
  
  # Variables
  outcome <- i
  plot_components <- c(paste(controls, collapse = ' + '), paste(exposures, collapse = ' + '))
  plot_title <- outcome
  
  # Create xts object and subset relevant variables as numeric
  Y <- xts(subset(df, select=c(outcome, exposures, controls)),  order.by=df$date)
  
  # New data (to forecast)
  new.period <- paste0(estimation.date.end+1, '/', estimation.date.end+n.ahead)
  new.data <- Y[new.period,c(exposures, controls)]
  
  # Get results
  results <- run_models(Y, y, exposures, controls, new.data, 
                        estimation.date.end, estimation.window, n.ahead, 
                        q, confidence.level)
  
  if (i == 'total_patients_c') {
    log_growth_plot <- results$p1
  }
  
  # Get plots
  n <- sub('_c$', '', outcome)
  n <- replace_exact_matches(n, dict = nm)
  p <- results$p3$figure
  p <- p + labs(title = paste(n)) + theme(plot.title = element_text(hjust = 0.5))
  plot_list[[outcome]] <- p
  
  # Get comparison stats
  results$criteria$Model <- sub('_c$', '', results$criteria$Model)
  comp_df <- rbind(comp_df, replace_exact_matches(results$criteria, dict = nm))
  
  # Transform coeficient outputs
  df_transformed <- results$estimates[1:8, ] %>%
    mutate(across(where(is.numeric), ~ (exp(.) - 1) * 100))
  
  # Get estimates
  df_formatted <- df_transformed %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
    mutate(
      Estimate = ifelse(abs(Estimate) < 1e-6, 0, Estimate),
      CI_Lower = ifelse(abs(CI_Lower) < 1e-6, 0, CI_Lower),
      CI_Upper = ifelse(abs(CI_Upper) < 1e-6, 0, CI_Upper),
      Lag = str_extract(State, 'lag\\d+'),
      State = str_replace(State, "_lag\\d+", ""),           # Remove lag part from State column
      Value = sprintf('%.2f (%.2f, %.2f)', Estimate, CI_Lower, CI_Upper)
    ) %>%
    select(State, Lag, Value) %>%
    pivot_wider(names_from = c(State, Lag), values_from = Value) %>%
    rename_with(~ gsub('_lag', ' Lag ', .x)) %>%
    mutate(
      MAPE = paste0(results$p3$MAPE, '%'),
      Outcome = n
    ) %>%
    relocate(Outcome)
  
  esti_df <- rbind(esti_df, df_formatted)
}


# Reset names
nm_m <- nm[!(names(nm) == "M")]
names(esti_df) <- sub('_b', '', names(esti_df))
colnames(esti_df) <- names(esti_df) %>%
  sapply(function(x) {
    for(key in names(nm_m)) {
      if(str_detect(x, fixed(key))) {
        x <- str_replace(x, fixed(key), nm_m[[key]])
      }
    }
    x
  })

# Update header
lag_row <- str_extract(colnames(esti_df), "\\Lag \\d+")
lag_row[is.na(lag_row)] <- " "
new_colnames <- str_trim(str_remove(colnames(esti_df), "Lag \\d+"))
colnames(esti_df) <- new_colnames
esti_df <- rbind(as.list(lag_row), esti_df)


# 2. Model Results (exposures): Plots ####

ncol <- 4
nplots <- length(plot_list)
nrow <- ceiling(nplots / ncol)
bottom_row_indices <- ((nrow - 1) * ncol + 1):nplots
left_col_indices <- seq(1, nplots, by = ncol)

plot_list <- lapply(seq_along(plot_list), function(i) {
  p <- plot_list[[i]]
  
  # Remove legend and minimize title size
  p <- p + theme(
    legend.position = 'none',
    plot.title = element_text(size = 10)
  )
  
  # Remove x-axis for plots not in bottom row
  if (!(i %in% bottom_row_indices)) {
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  } else{
    p <- p + theme(axis.title.x = element_blank())
  }
  
  # Remove y-axis for plots not in bottom row
  if (!(i %in% left_col_indices)) {
    p <- p + theme(axis.title.y = element_blank())
  }
  p
})

# Add legend back
plot_list[[1]] <- plot_list[[1]] + theme(legend.position = "bottom")

# Arrange plots
forecast_plots <- wrap_plots(plot_list, ncol = ncol) + 
  plot_layout(guides = "collect") & theme(legend.position = "bottom")  

#### EXPORT RESULTS (exposures) ####

# Round df if needed
comp_df <- round_numeric(comp_df, 2)
# Remove irrelevant cols
names(esti_df) <- make.unique(names(esti_df))
esti_df <- esti_df %>% select(matches("Outcome|Mx2t|Tp|Fg10|MAPE"))
# Drop header
esti_df <- esti_df[-1, ]

# Save tables
write.xlsx(esti_df, file = paste0(path3, "/III. Tables/MT_weat_coefs_Aug16_ALL.xlsx"))
write.xlsx(comp_df, file = paste0(path3, "/III. Tables/ST_comp_model_Aug16_ALL.xlsx"))

# Save figures
ggsave(paste0(path3, '/IV. Figures/MF_forecast_Aug16_ALL.png'), forecast_plots, width = 10, height = 6, dpi = 300)
ggsave(paste0(path4, '/MF_forecast_Aug16_ALL.png'), forecast_plots, width = 10, height = 6, dpi = 300)
# ggsave(paste0(path3, '/IV. Figures/MF_log_Aug16_ALL.png'), log_growth_plot, width = 7, height = 5, dpi = 300)


# 3. Model Results (no exposures): Coefficients ####

outcomes <- c('total_patients_c', 
              'F_c','M_c','under_5_c','age5_14_c','age15_24_c','age25_64_c','over_65_c',
              'general_c','ophthalmology_c','pediatrics_c','ent_c') # 'orthopedics_c','gynecology_c'
controls <- c('closed')
exposures <- c()

# Outputs
esti_df <- data.frame()
comp_df <- data.frame()
plot_list <- list()

for (i in outcomes) {
  
  # Variables
  outcome <- i
  plot_components <- c(paste(controls, collapse = ' + '), paste(exposures, collapse = ' + '))
  plot_title <- outcome
  
  # Create xts object and subset relevant variables as numeric
  Y <- xts(subset(df, select=c(outcome, exposures, controls)),  order.by=df$date)
  
  # New data (to forecast)
  new.period <- paste0(estimation.date.end+1, '/', estimation.date.end+n.ahead)
  new.data <- Y[new.period,c(exposures, controls)]
  
  # Get results
  results <- run_models(Y, y, exposures, controls, new.data, 
                        estimation.date.end, estimation.window, n.ahead, 
                        q, confidence.level)
  
  if (i == 'total_patients_c') {
    log_growth_plot <- results$p1
  }
  
  # Get plots
  n <- sub('_c$', '', outcome)
  n <- replace_exact_matches(n, dict = nm)
  p <- results$p3$figure
  p <- p + labs(title = paste(n)) + theme(plot.title = element_text(hjust = 0.5))
  plot_list[[outcome]] <- p
  
  # Get comparison stats
  results$criteria$Model <- sub('_c$', '', results$criteria$Model)
  comp_df <- rbind(comp_df, replace_exact_matches(results$criteria, dict = nm))
  
  # Transform coeficient outputs
  df_transformed <- results$estimates[1:8, ] %>%
    mutate(across(where(is.numeric), ~ (exp(.) - 1) * 100))
  
  # Get estimates
  df_formatted <- df_transformed %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
    mutate(
      Estimate = ifelse(abs(Estimate) < 1e-6, 0, Estimate),
      CI_Lower = ifelse(abs(CI_Lower) < 1e-6, 0, CI_Lower),
      CI_Upper = ifelse(abs(CI_Upper) < 1e-6, 0, CI_Upper),
      Lag = str_extract(State, 'lag\\d+'),
      State = str_replace(State, "_lag\\d+", ""),           # Remove lag part from State column
      Value = sprintf('%.2f (%.2f, %.2f)', Estimate, CI_Lower, CI_Upper)
    ) %>%
    select(State, Lag, Value) %>%
    pivot_wider(names_from = c(State, Lag), values_from = Value) %>%
    rename_with(~ gsub('_lag', ' Lag ', .x)) %>%
    mutate(
      MAPE = paste0(results$p3$MAPE, '%'),
      Outcome = n
    ) %>%
    relocate(Outcome)
  
  esti_df <- rbind(esti_df, df_formatted)
}


# Reset names
nm_m <- nm[!(names(nm) == "M")]
names(esti_df) <- sub('_b', '', names(esti_df))
colnames(esti_df) <- names(esti_df) %>%
  sapply(function(x) {
    for(key in names(nm_m)) {
      if(str_detect(x, fixed(key))) {
        x <- str_replace(x, fixed(key), nm_m[[key]])
      }
    }
    x
  })

# Update header
lag_row <- str_extract(colnames(esti_df), "\\Lag \\d+")
lag_row[is.na(lag_row)] <- " "
new_colnames <- str_trim(str_remove(colnames(esti_df), "Lag \\d+"))
colnames(esti_df) <- new_colnames
esti_df <- rbind(as.list(lag_row), esti_df)


# 4. Model Results (no exposures): Plots ####

ncol <- 4
nplots <- length(plot_list)
nrow <- ceiling(nplots / ncol)
bottom_row_indices <- ((nrow - 1) * ncol + 1):nplots
left_col_indices <- seq(1, nplots, by = ncol)

plot_list <- lapply(seq_along(plot_list), function(i) {
  p <- plot_list[[i]]
  
  # Remove legend and minimize title size
  p <- p + theme(
    legend.position = 'none',
    plot.title = element_text(size = 10)
  )
  
  # Remove x-axis for plots not in bottom row
  if (!(i %in% bottom_row_indices)) {
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  } else{
    p <- p + theme(axis.title.x = element_blank())
  }
  
  # Remove y-axis for plots not in bottom row
  if (!(i %in% left_col_indices)) {
    p <- p + theme(axis.title.y = element_blank())
  }
  p
})

# Add legend back
plot_list[[1]] <- plot_list[[1]] + theme(legend.position = "bottom")

# Arrange plots
forecast_plots <- wrap_plots(plot_list, ncol = ncol) + 
  plot_layout(guides = "collect") & theme(legend.position = "bottom")  

#### EXPORT RESULTS (no exposures) ####

# Round df if needed
comp_df <- round_numeric(comp_df, 2)
# Remove irrelevant cols
names(esti_df) <- make.unique(names(esti_df))
# esti_df <- esti_df %>% select(-matches("Mx2t|Tp"))
names(esti_df) <- gsub(".1", "", names(esti_df))
esti_df[is.na(esti_df)] <- ""
# Drop header
names(esti_df) <- esti_df[1,]
names(esti_df)[1] <- "Outcome"
names(esti_df)[ncol(esti_df)] <- "MAPE"
esti_df <- esti_df[-1,]

# Save tables
write.xlsx(comp_df, file = paste0(path3, "/III. Tables/ST_comp_model_Aug16_noExp.xlsx"))
# Save figures
ggsave(paste0(path3, '/IV. Figures/SF_forecast_Aug16_noExp.png'), forecast_plots, width = 10, height = 6, dpi = 300)
ggsave(paste0(path4, 'SF_forecast_Aug16_noExp.png'), forecast_plots, width = 10, height = 6, dpi = 300)


