### ------------------------------------------------------------------------ ### 
###  TSGC - Forecast Nov 2025 (HISP: Shimla)
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list = ls()); graphics.off(); gc()

# Set data path
setwd('C:/Users/offne/Documents/GitHub/hisp-climate-health-forecast/')

#### DEPENDENCIES ####
source('3_TSGC/TSGC - Dependencies.R')
setwd(path_git)
file_import <- 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide'
file_export <- 'shimla-2025-11'

### Get model functions (make sure wd is in repo folder)
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
estimation.window <- 4 # months
estimation.date.end <- as.Date('2025-11-14')
n.ahead <- 10
confidence.level <- 0.68
sea.type <- 'trigonometric'
sea.period <- 7

df <- df %>% filter(location == 'Shimla')
startdate <- as.Date('2025-01-01')
enddate <- as.Date('2025-12-31')
df <- df[df$date >= startdate & df$date <= enddate,]


# 1. Model Results (no exposures): Coefficients ####

outcomes <- c('total_patients_c', 
              'F_c','M_c')
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
      State = str_replace(State, '_lag\\d+', ''),           # Remove lag part from State column
      Value = sprintf('%.2f (%.2f, %.2f)', Estimate, CI_Lower, CI_Upper)
    ) %>%
    dplyr::select(State, Lag, Value) %>%
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
nm_m <- nm[!(names(nm) == 'M')]
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
lag_row <- str_extract(colnames(esti_df), '\\Lag \\d+')
lag_row[is.na(lag_row)] <- ' '
new_colnames <- str_trim(str_remove(colnames(esti_df), 'Lag \\d+'))
colnames(esti_df) <- new_colnames
esti_df <- rbind(as.list(lag_row), esti_df)


# 2. Model Results (no exposures): Plots ####

ncol <- 3
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
plot_list[[1]] <- plot_list[[1]] + theme(legend.position = 'bottom')

# Arrange plots
(forecast_plots <- wrap_plots(plot_list, ncol = ncol) + 
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Shimla without exposures (Nov 2025)",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    ) &
    theme(legend.position = "bottom"))
#### EXPORT RESULTS (no exposures) ####

# Round df if needed
comp_df <- round_numeric(comp_df, 2)
# Remove irrelevant cols
names(esti_df) <- make.unique(names(esti_df))
# esti_df <- esti_df %>% select(-matches('Mx2t|Tp'))
names(esti_df) <- gsub('.1', '', names(esti_df))
# Drop header
names(esti_df) <- esti_df[1,]
names(esti_df)[1] <- 'Outcome'
names(esti_df)[ncol(esti_df)] <- 'MAPE'
esti_df <- esti_df[-1,]

# Save tables
# write.xlsx(comp_df, file = paste0(path_report, 'III. tsgc/tab_comp_model_', file_export, '_noExp.xlsx'))
# Save figures
ggsave(paste0(path_report, 'III. tsgc/fig_forecast_', file_export, '_noExp.png'), forecast_plots, width = 10, height = 5, dpi = 300)


