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

# Function to summarize ANOVA
anova_summary <- function(dep_var, factor_var, data){
  formula <- as.formula(paste(dep_var, '~', factor_var))
  fit <- aov(formula, data = data)
  sum_fit <- summary(fit)[[1]]
  df1 <- sum_fit[1, 'Df']
  df2 <- sum_fit[2, 'Df']
  F_val <- sum_fit[1, 'F value']
  p_val <- sum_fit[1, 'Pr(>F)']
  
  data.frame(
    Variable = dep_var,
    Factor = factor_var,
    Df_between = df1,
    Df_within = df2,
    F = round(F_val, 3),
    p = round(p_val, 3)
  )
}

# Function to summarize ANOVA interaction tests
anova_interaction_summary <- function(dep_var, factor1, factor2, data){
  
  if (factor1 == factor2) {
    return(NULL)
  }
  
  formula <- as.formula(paste(dep_var, '~', factor1, '*', factor2))
  fit <- aov(formula, data = data)
  sum_fit <- summary(fit)[[1]]
  
  results <- data.frame(
    Variable = dep_var,
    Term = rownames(sum_fit),
    Df = sum_fit[,'Df'],
    F = round(sum_fit[,'F value'], 3),
    p = round(sum_fit[,'Pr(>F)'], 3))
  
  # results <- results[-nrow(results), ]
  # results <- results[3, ]
  results <- results[grepl(':', results$Term), ]
  
  
  return(results)
}


# INTERACTION TESTS ####

# Non-parametric spearman correlation tests (raw data)
sub <- df_wide %>% select(total_patients, all_of(w_cols[!grepl('_p|IMD', w_cols)]))
cor_matrix <- cor(sub, use = 'complete.obs', method = 'spearman')
melted_corr <- melt(cor_matrix)
cor_matrix <- ggplot(data = melted_corr, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = 0, limit = c(-1,1)) +
  theme_minimal() +
  scale_x_discrete(labels = function(x) {
    # If x is in nm, replace with nm value, else keep original
    sapply(x, function(lbl) {
      if (lbl %in% names(nm)) {
        nm[[lbl]]
      } else {
        lbl
      }})}) +
  scale_y_discrete(labels = function(x) {
    # If x is in nm, replace with nm value, else keep original
    sapply(x, function(lbl) {
      if (lbl %in% names(nm)) {
        nm[[lbl]]
      } else {
        lbl
      }})}) +
  geom_text(aes(label=round(value,2)), color='black', size=4) +
  labs(title = 'Spearman Correlation Matrix', x='', y='') +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # rotate x-axis labels
    axis.text.y = element_text(angle = 0)                # rotate y-axis labels if needed
  )

# IMD ANOVA tests
sub_IMD <- df_wide %>% select(total_patients, season, all_of(w_cols[grepl('IMD', w_cols)]))
test_cols <- names(sub_IMD[,-c(1, 2)])  # everything except total_patients
(IMD_anova <- do.call(rbind, lapply(c('season', test_cols), function(f) anova_summary('total_patients', f, sub_IMD))))

IMD_anova$Factor <- gsub('_b', '', IMD_anova$Factor)
IMD_anova <- replace_exact_matches(IMD_anova, nm)
names(IMD_anova) <- c('Variable', 'Factor', 'Df Between','Df Within','F','p')

# Interaction tests
# Season
(IMD_inter_seas <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'season', f, sub_IMD))))
# Humidity
(IMD_inter_hum_h <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'rh_hig_IMD_b', f, sub_IMD))))
(IMD_inter_hum_l <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'rh_low_IMD_b', f, sub_IMD))))
# Pollution
(IMD_inter_pm25 <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'pm25_IMD_b', f, sub_IMD))))
# Wind
(IMD_inter_wind <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'fg10_IMD_b', f, sub_IMD))))

IMD_interaction_tests <- rbind(IMD_inter_seas, IMD_inter_hum_h, IMD_inter_hum_l, IMD_inter_pm25, IMD_inter_wind)
row.names(IMD_interaction_tests) <- NULL


# CLIM ANOVA tests
sub_clim <- df_wide %>% select(total_patients, season, all_of(w_cols[grepl('_p', w_cols)]))
test_cols <- names(sub_clim[,-c(1, 2)])  # everything except total_patients
(CLIM_anova <- do.call(rbind, lapply(c('season', test_cols), function(f) anova_summary('total_patients', f, sub_clim))))
test_cols <- c('season', 'tp_p90', 'fg10_p90', 'rh_mean_p10', 'rh_mean_p90')

CLIM_anova$Factor <- gsub('_b', '', CLIM_anova$Factor)
CLIM_anova <- replace_exact_matches(CLIM_anova, nm)
names(CLIM_anova) <- c('Variable', 'Factor', 'Df Between','Df Within','F','p')


# Interaction tests
# Season
(clim_inter_seas <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'season', f, sub_clim))))
# Humidity
(clim_inter_hum_h <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'rh_mean_p90', f, sub_clim))))
(clim_inter_hum_l <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'rh_mean_p10', f, sub_clim))))
# Pollution
(clim_inter_pm25 <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'pm25_p90', f, sub_clim))))
# Wind
(clim_inter_wind <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'pm25_p90', f, sub_clim))))

CLIM_interaction_tests <- rbind(clim_inter_seas, clim_inter_hum_h, clim_inter_hum_l, clim_inter_pm25, clim_inter_wind)
row.names(CLIM_interaction_tests) <- NULL
