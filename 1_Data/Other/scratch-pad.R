### ------------------------------------------------------------------------ ### 
### HISP Data cleaning - SCRATCH PAD                                        ###
### ------------------------------------------------------------------------ ### 

# SCRATCH PAD

getwd()
y <- read.csv('2_Active_Data/201_heal/- Archive/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba_clean_wide.csv')

y <- read.csv('2_Active_Data/201_heal/- Archive/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba_clean.csv')
x <- load('2_Active_Data/201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_ALL_clean_all.RData')
x <- x[x$location=='Chamba',]

any(duplicated(y))
y[duplicated(y), ]


y <- y[y$date >= startdate & y$date <= enddate,]
# x <- df_wide[df_wide$location=='Chamba',]
# x <- df_all[df_all$location=='Chamba',]

min(x$total_patients, na.rm=TRUE)
max(x$total_patients, na.rm=TRUE)
min(y$total_patients, na.rm=TRUE)
max(y$total_patients, na.rm=TRUE)


# Time series plots
plot_ts <- function(df_all, column_name) {
  ggplot(df_all, aes(x = date, y = .data[[column_name]])) + 
    geom_line(color = mypalette[1]) +
    ggtitle(column_name) +
    xlab('Date') +                
    ylab('Patients') +      
    theme_minimal() +            
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
    scale_x_date(breaks = seq.Date(min(df_all$date), max(df_all$date), by = '6 months'), 
                 labels = scales::date_format('%b %Y'))
}

# ACF plots with custom ggplots
plot_acf <- function(series, title) {
  series <- na.omit(series)
  acf_data <- acf(as.numeric(series), lag.max = length(series), plot = FALSE)
  n <- length(series)
  ci <- 1.96 / sqrt(n)  # 95% confidence interval
  
  df_all_acf <- data.frame(
    Lag = as.numeric(acf_data$lag),
    ACF = as.numeric(acf_data$acf)
  )
  
  ggplot(df_all_acf, aes(x = Lag, y = ACF)) +
    geom_segment(aes(xend = Lag, yend = 0), color = 'gray20') +
    geom_hline(yintercept = 0, color = 'black') +
    geom_hline(yintercept = c(-ci, ci), linetype = 'dashed', color = 'blue') +
    labs(title = title, x = 'Lag', y = 'ACF') +
    theme_minimal()
}

# Compare stationary vs non stationary data
plot_stationarity <- function(df_all, base_col) {
  v <- grep(base_col, names(df_all), value = TRUE)
  
  if (length(v) < 2) {
    warning('Need at least two columns: original and detrended/stationary version.')
    return(invisible(NULL))
  }
  
  # Time series plots
  p1 <- ggplot(df_all, aes(x = date, y = .data[[v[1]]])) +
    geom_line(color = mypalette[1]) +
    labs(title = 'Non-Stationary Time Series', x = 'Time (t)', y = base_col) +
    theme_classic()
  
  p3 <- ggplot(df_all, aes(x = date, y = .data[[v[2]]])) +
    geom_line(color = mypalette[2]) +
    labs(title = 'Stationary Time Series', x = 'Time (t)', y = base_col)+
    theme_classic()
  
  p2 <- plot_acf(df_all[[v[1]]], 'ACF - Non-Stationary')+theme_classic()
  p4 <- plot_acf(df_all[[v[2]]], 'ACF - Stationary') +  theme_classic()
  
  # Combine into a 2x2 layout
  layout <- (p1 | p2) / (p3 | p4) +
    patchwork::plot_annotation(
      title = base_col,
      theme = theme(plot.title = element_text(size = 16, face = 'bold', hjust = 0.5))
    )
  
  
  return(layout)
}

# Get and visualize anomalies with percentiles
plot_anomalies <- function(ts_data){
  
  # Function to calculate the 5th and 95th percentiles
  calculate_percentiles <- function(df_all, quant=0.95) {
    bottom_p <- quantile(df_all, 1 - quant, na.rm = TRUE)
    top_p <- quantile(df_all, quant, na.rm = TRUE)
    return(c(bottom_p, top_p))
  }
  
  # Function to detect anomalies using Z-score
  detect_anomalies <- function(df_all) {
    z_scores <- scale(df_all)  # Compute Z-scores
    anomalies <- which(abs(z_scores) > 3)  # Identify anomalies (Z-score > 3)
    return(anomalies)
  }
  
  # Create a function to plot each variable with 5th and 95th percentiles
  plot_detrend <- function(var_name, df_all, anomalies, percentiles) {
    ggplot(data = data.frame(date = ts_data$date, value = df_all), aes(x = date, y = value)) +
      geom_line() +
      geom_point(data = data.frame(date = ts_data$date[anomalies], value = df_all[anomalies]), 
                 aes(x = date, y = value), color = mypalette[2], size = 2) +
      geom_hline(yintercept = percentiles[1], color = mypalette[1], linetype = 'dashed') +  # 5th percentile line
      geom_hline(yintercept = percentiles[2], color = mypalette[1], linetype = 'dashed') +  # 95th percentile line
      ggtitle(paste(var_name)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(y = 'Detrended Value')
  }
  
  # Calculate percentiles for each detrended variable
  percentiles_list <- setNames(
    apply(ts_data[,-1], 2, calculate_percentiles),
    colnames(ts_data)[-1]  # Name the list elements
  )
  # Apply anomaly detection to each detrended variable
  anomalies_list <- lapply(ts_data[,-1], detect_anomalies)
  
  # List of plots for each variable
  plots <- list()
  for (i in 2:ncol(ts_data)) {
    var_name <- colnames(ts_data)[i]
    anomalies <- anomalies_list[[i - 1]]  # Get anomalies for this variable
    percentiles <- percentiles_list[, var_name]  # Get percentiles for this variable
    plots[[var_name]] <- plot_detrend(var_name, ts_data[[i]], anomalies, percentiles)
  }
  
  return(list(anomalies = anomalies_list, plots = plots))
  
}
