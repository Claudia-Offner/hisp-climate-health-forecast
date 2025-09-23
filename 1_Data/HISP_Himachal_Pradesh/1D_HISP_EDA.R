### ------------------------------------------------------------------------ ### 
### HISP Exploratory Data Analysis                                           ###
### ------------------------------------------------------------------------ ### 


####  DEPENDENCIES ####

# Turn off scientific notation
options(scipen=999)

# Set data path
user <- 'claer14'
path1 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
path2 <- paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/1_Data/HISP_Himachal_Pradesh/- Reporting/- figures/')
setwd(path1)

# Packages
library(tidyverse)
library(dplyr)
library(skimr)
library(DT)
library(gridExtra)
library(ggplot2)
library(zoo)
library(forecast)
library(patchwork)
library(purrr)

# Function to create a histogram for a column
create_histogram <- function(df, column_name, stat=NULL) {
  # df <- sub
  # stat <- skm
  # column_name <- 'mn2t'
  hist_data <- hist(df[[column_name]], breaks=40, plot = FALSE)
  y_max_idx <- which.max(hist_data$density)
  y_med_idx <- which.min(abs(hist_data$density - median(hist_data$density)))
  position <- y_max_idx > y_med_idx
  
  if(position == TRUE){
    x_pos <- min(df[[column_name]], na.rm = TRUE)
    h_pos <- 0
  } else {
    x_pos <- max(df[[column_name]], na.rm = TRUE)
    h_pos <- 1
  }
  
  plot <- ggplot(df, aes_string(x = column_name)) +
    geom_histogram(fill = mypalette[1], color = 'black') +
    theme_minimal() +
    theme(axis.title.y = element_blank()) +
    ggtitle(paste(column_name)) 
  
  if (!is.null(stat)){
    plot <- plot +
      annotate("text",
               x = x_pos, y = Inf,
               label = paste("\nNA: ", stat[stat$skim_variable == column_name, "n_missing"],
                             "\nMin: ", stat[stat$skim_variable == column_name, "min"],
                             "\nMedian: ", stat[stat$skim_variable == column_name, "median"],
                             "\nMax: ", stat[stat$skim_variable == column_name, "max"]),
               hjust = h_pos, vjust = 1, size = 3, color = "black", fontface = "italic")
  }
  
  return(plot)

}

# Time series plots
plot_ts <- function(df, column_name) {
  ggplot(df, aes(x = date, y = .data[[column_name]])) + 
    geom_line(color = mypalette[1]) +
    ggtitle(column_name) +
    xlab('Date') +                
    ylab('Patients') +      
    theme_minimal() +            
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
    scale_x_date(breaks = seq.Date(min(df$date), max(df$date), by = '6 months'), 
                 labels = scales::date_format('%b %Y'))
}

# ACF plots with custom ggplots
plot_acf <- function(series, title) {
  series <- na.omit(series)
  acf_data <- acf(as.numeric(series), lag.max = length(series), plot = FALSE)
  n <- length(series)
  ci <- 1.96 / sqrt(n)  # 95% confidence interval
  
  df_acf <- data.frame(
    Lag = as.numeric(acf_data$lag),
    ACF = as.numeric(acf_data$acf)
  )
  
  ggplot(df_acf, aes(x = Lag, y = ACF)) +
    geom_segment(aes(xend = Lag, yend = 0), color = "gray20") +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = c(-ci, ci), linetype = "dashed", color = "blue") +
    labs(title = title, x = "Lag", y = "ACF") +
    theme_minimal()
}

# Compare stationary vs non stationary data
plot_stationarity <- function(df, base_col) {
  v <- grep(base_col, names(df), value = TRUE)
  
  if (length(v) < 2) {
    warning("Need at least two columns: original and detrended/stationary version.")
    return(invisible(NULL))
  }
  
  # Time series plots
  p1 <- ggplot(df, aes(x = date, y = .data[[v[1]]])) +
    geom_line(color = mypalette[1]) +
    labs(title = "Non-Stationary Time Series", x = "Time (t)", y = base_col) +
    theme_classic()
  
  p3 <- ggplot(df, aes(x = date, y = .data[[v[2]]])) +
    geom_line(color = mypalette[2]) +
    labs(title = "Stationary Time Series", x = "Time (t)", y = base_col)+
    theme_classic()
  
  p2 <- plot_acf(df[[v[1]]], "ACF - Non-Stationary")+theme_classic()
  p4 <- plot_acf(df[[v[2]]], "ACF - Stationary") +  theme_classic()
  
  # Combine into a 2x2 layout
  layout <- (p1 | p2) / (p3 | p4) +
  patchwork::plot_annotation(
    title = base_col,
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )
  
  
  return(layout)
}

# Get and visualize anomalies with percentiles
plot_anomalies <- function(ts_data){
  
  # Function to calculate the 5th and 95th percentiles
  calculate_percentiles <- function(df, quant=0.95) {
    bottom_p <- quantile(df, 1 - quant, na.rm = TRUE)
    top_p <- quantile(df, quant, na.rm = TRUE)
    return(c(bottom_p, top_p))
  }
  
  # Function to detect anomalies using Z-score
  detect_anomalies <- function(df) {
    z_scores <- scale(df)  # Compute Z-scores
    anomalies <- which(abs(z_scores) > 3)  # Identify anomalies (Z-score > 3)
    return(anomalies)
  }
  
  # Create a function to plot each variable with 5th and 95th percentiles
  plot_detrend <- function(var_name, df, anomalies, percentiles) {
    ggplot(data = data.frame(date = ts_data$date, value = df), aes(x = date, y = value)) +
      geom_line() +
      geom_point(data = data.frame(date = ts_data$date[anomalies], value = df[anomalies]), 
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

mypalette <- c('#1D98EB','#D74B3C','#A6B727','#FEC306','#EC800A','#E569BF','#8968CD',
               '#A1C8F2','#F1A6A6','#C2D77D','#F6D974','#F9B34D','#F1C2D1','#B8A8D9')

# weekday_levels <- c('Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday')
# month_levels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
# season_levels <- c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)', 'Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)')
# w_cols <- c('t2m','d2m','mn2t','mx2t','tp','fg10','v10','u10','sp')
# age_levels <- c('TOTAL', '<1y', '1-17y', '18-60y', '>60y', 'Unknown')
# gender_levels <- c('TOTAL', 'F', 'M', 'O')
# ward_levels <- c('TOTAL', 'General Medicine', 'Emergency', 'Ophthalmology', 
#                  'Pediatrics', 'Orthopedics', 'ENT', 'Gynecology', 'Dermatology',
#                  'Surgery', 'Dentistry', 'Psychiatry', 'Oncology')


####  MAIN CODE ####
# Load data ####
setwd(path1)

# Patient data format (!! not processed for NAs !!)
df <- read.csv('203_comb/ind_5km_BigThaw-HISP-chamba_2015-2023_hosp_all.csv', check.names = FALSE, fileEncoding = 'UTF-8')
(column_types <- sapply(df, class))
df$date <- as.Date(df$date)

# Wide data format
df_w <- read.csv('203_comb/ind_5km_BigThaw-HISP-chamba_2015-2023_hosp_wide.csv', check.names = FALSE, fileEncoding = 'UTF-8')
(column_types <- sapply(df_w, class))
df_w$date <- as.Date(df_w$date)

# Remove COVID years for now:
startdate <- as.Date('2015-02-14')
enddate <- as.Date('2020-02-14')  
df_w <- df_w[df_w$date >= startdate & df_w$date <= enddate,]

# Back transform
df_w$tp <- exp(df_w$tp)

# 1. Numeric: Total Patients & Weather ####
# Use wide format to avoid issues with date duplicates


sub <- df_w %>% select(total_patients,mn2t,mx2t,tp,fg10)
skm <- skim(sub) %>% rename(min = numeric.p0, Q1 = numeric.p25, 
                            median = numeric.p50, Q3 = numeric.p75, 
                            max = numeric.p100) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

colnames(skm) <- gsub('numeric.', '', colnames(skm))
skm <- as.data.frame(skm)


# Histograms
cols_to_plot <- colnames(sub[-1])
plots <- lapply(cols_to_plot, function(col) create_histogram(sub, col, skm))
hist_weat <- grid.arrange(grobs=plots, nrow=2)
hist_tota <- create_histogram(sub, 'total_patients', skm)
hist_tota <- hist_tota +
  ggtitle("Chamba total patient distribution") +
  xlab("Total patients") +           # Change Y-axis label
  theme(
    legend.position = "none",        # Remove legend
    plot.title = element_text(hjust = 0.5)        # Center the title
  )


# # Density curves by season
# labels_map <- c("Winter (Jan-Feb)" = "Jan-Feb",
#                 "Pre-Monsoon (Mar-May)" = "Mar-May",
#                 "Monsoon (Jun-Sep)" = "Jun-Sep",
#                 "Post-Monsoon (Oct-Dec)" = "Oct-Dec")
# 
# plots <- lapply(cols_to_plot, function(c) {
#   ggplot(df_w, aes_string(x = c, color = "season", fill = "season")) +
#     geom_density(alpha = 0.3) +
#     theme_minimal() +
#     labs(
#       title = c,
#       x = c,
#       y = "Density",
#       fill = "Season",
#       color = "Season"
#     ) +
#     scale_fill_manual(values=c(mypalette[1], mypalette[4], mypalette[2], mypalette[3]), labels=labels_map) +
#     scale_color_manual(values=c(mypalette[1], mypalette[4], mypalette[2], mypalette[3]), labels=labels_map)
# })
# grid.arrange(grobs = c(plots[c(1:3)],list(ggplot() + theme_void()), plots[c(4:5)]), ncol=3, nrow=2)

# 2. Factor: Wards ####
# Use OG data to easily analyse groups and avoid col issues

# CHECK FREQUENCY
cat_ward <- ggplot(df, aes(x = ward_name, fill = ward_name)) +
            geom_bar(color = 'black') + 
            ggtitle('Chamba Ward Distribution') +  
            xlab('Ward') + 
            ylab('Count') + 
            scale_y_continuous(labels = scales::comma) + 
            scale_fill_manual(values = mypalette) +
            theme_minimal() +   
            theme(legend.position = 'none', 
                  axis.text.x = element_text(angle = 45, hjust = 1),
                  plot.title = element_text(hjust = 0.5))

# 3. Factor: Age ####
# Use OG data to easily analyse groups and avoid col issues

# CHECK FREQUENCY
df <- df %>% mutate(age_group = factor(age_group, levels = c('under_5', 'age5_14', 'age15_24', 'age25_64', 'over_65', 'Unknown')))
cat_age <- ggplot(df, aes(x = age_group, fill = age_group)) +
            geom_bar(color = 'black') +  # Bar color
            ggtitle('Chamba Age Distribution') +  # Title
            xlab('Age Groups') +  # X-axis label
            ylab('Count') +   # Y-axis label
            scale_fill_manual(values = mypalette) +
            scale_y_continuous(labels = scales::comma) +  # Remove scientific notation (comma format)
            theme_minimal() +   # Minimal theme for clean appearance
            theme(legend.position = 'none',
                  plot.title = element_text(hjust = 0.5))


# CHECK DISTRIBUTION
con_age <- ggplot(df, aes(x = age)) +
  geom_histogram(binwidth = 1, fill = mypalette[1], color = 'black') +  # Set bin width to 1
  ggtitle('Histogram of Chamba Age Distribution') +
  xlab('Age') +
  ylab('Frequency') +
  xlim(-25, 125) +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5))

# CHECK ANOMALY PATTERNS
# There are 286 anomalies in feasible age
table(df$age_group)
sort(unique(df$age))
# No clear pattern across wards
df %>% filter(age_group == 'Unknown') %>% count(age) # No clear pattern
df %>% filter(age_group == 'Unknown') %>% count(ward_name) # No clear pattern
df %>% filter(age_group == 'Unknown') %>% count(gender) 

# 4. Factor: Gender ####
# Use OG data to easily analyse groups and avoid col issues

# CHECK FREQUENCY
cat_sex <- ggplot(df, aes(x = gender, fill = gender)) +
            geom_bar(color = 'black') +  # Bar color
            ggtitle('Chamba Gender Distribution') +  # Title
            xlab('Gender') +  # X-axis label
            ylab('Count') +   # Y-axis label
            scale_fill_manual(values = mypalette) +  # Manual colors for Male and Female
            scale_y_continuous(labels = scales::comma) +  # Remove scientific notation (comma format)
            theme_minimal() +   # Minimal theme for clean appearance
            theme(legend.position = 'none',
                  plot.title = element_text(hjust = 0.5))

# CHECK ANOMALY PATTERNS
# Only 500 listed as 'Other' - 276 of which are from Pediatrics and aged 1-17y
table(df$gender)
df %>% filter(gender == 'O') %>% count(ward_name)
df %>% filter(gender == 'O') %>% count(age_group)



##### EXPORT ####
w <- 6
h <- 6

ggsave(paste0(path2, 'desc_hist_weat.png'), hist_weat, width = 8, height = h, dpi = 300)
ggsave(paste0(path2, 'desc_hist_tota.png'), hist_tota, width = w, height = h, dpi = 300)
ggsave(paste0(path2, 'desc_cat_ward.png'), cat_ward, width = 8, height = h, dpi = 300)
ggsave(paste0(path2, 'desc_cat_sex.png'), cat_sex, width = w, height = h, dpi = 300)
ggsave(paste0(path2, 'desc_cat_age.png'), cat_age, width = w, height = h, dpi = 300)
ggsave(paste0(path2, 'desc_con_age.png'), con_age, width = w, height = h, dpi = 300)

####  TEMPORAL  ####

# 2.1 Check weather time series ####

w_cols_to_plot <- c('mn2t','mx2t','tp','fg10')

# Check time series pattern
#### There is obvious seasonality in this data
ts_plots <- lapply(w_cols_to_plot, function(col) plot_ts(df_w, col))
grid.arrange(grobs = ts_plots[1:4], ncol = 2, nrow = 2)

# Check trend (yearly average)
#### For most variables, there is no strong change in trend BUT mn2t_mean and  
#### tp_mean have the most year to year change, though no clear pattern.

for (c in w_cols_to_plot){
  trends <- df_w %>% 
     group_by(year) %>%
     summarise(across(all_of(c),
                      list(mean = ~mean(. , na.rm = TRUE),
                           sd = ~sd(. , na.rm = TRUE),
                           min = ~min(. , na.rm = TRUE),
                           max = ~max(. , na.rm = TRUE)),
                      .names = "{.col}_{.fn}"))
  print(trends)
}

(trends <- df_w %>% 
    group_by(year) %>%
    summarise(across(all_of(w_cols_to_plot),
                     list(mean = ~mean(. , na.rm = TRUE)),
                     .names = "{.col}")))

trends_long <- trends %>%
  pivot_longer(cols = -year, names_to = "variable", values_to = "mean_value") %>%
  mutate(variable = factor(variable, levels = w_cols_to_plot))

ggplot(trends_long, aes(x = year, y = mean_value)) +
  geom_line(color = mypalette[2]) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "Yearly Mean Values", x = "Year", y = "Mean Value") +
  theme_minimal()

# Check ACF & PACF
#### ACF: Most have a strong sinusoidal pattern (indicating seasonal trend) and 
#### tapers off slowly (indicating non-stationarity).
#### PACF: Sharp spikes at lags 1-7, which may be important for modelling

for (c in w_cols_to_plot) {
  ts_data <- ts(df_w[[c]], frequency=365)
  par(mfrow = c(1, 2))  # Side-by-side
  acf(ts_data, lag.max=100, main = paste("ACF of", c))
  pacf(ts_data, lag.max=100, main = paste("PACF of", c))
  par(mfrow = c(1, 1))  # Reset layout
}

#### Data should be transformed to achieve stationary before applying to model, 
#### which may require a lag of 7 

# Check decomposition
#### Confirms previous tests, showing clear seasonality and trend.
#### Data will be trended for analysis
for (c in w_cols_to_plot){
  # Detrend and plot
  ts_data <- ts(df_w[[c]], frequency=365)
  stl_decomp <- stl(ts_data, s.window = 'periodic')
  plot(stl_decomp, main= paste('Decomposition of', c))
  # Create detrend variables
  # remainder <- stl_decomp$time.series[, 'remainder']
  # df_w[[paste(c, "_dt", sep = "")]] <- remainder
}

# Comparison of stationary vs non-stationary ACF tests
# lapply(w_cols_to_plot, function(col) plot_stationarity(df_w, col))
# But how would we interpret detrended data within the model??

# # Check lags: They have no effect.
# dt_w_cols_to_plot <- grep('_dt', names(df_w), value = TRUE)
# num <- 7
# for (c in dt_w_cols_to_plot){
#   lagged_var <- paste0(c, '_lag',num)
#   df_w <- df_w %>% mutate(!!lagged_var := lag(as.vector(df_w[[c]]), n = num))
# }
# lag_w_cols_to_plot <- grep(paste0('dt_lag', num), names(df_w), value = TRUE)
# for (c in lag_w_cols_to_plot) {
#   ts_data <- ts(df_w[[c]], frequency=1)
#   par(mfrow = c(1, 2))  # Side-by-side
#   acf(na.omit(ts_data), lag.max=100, main = paste("ACF of", c))
#   pacf(na.omit(ts_data), lag.max=100, main = paste("PACF of", c))
#   par(mfrow = c(1, 1))  # Reset layout
# }

# 2.2 Check health time series ####

h_cols_to_plot <- c('total_patients','F','M','under_5','age5_14','age15_24','age24_64','over_65',
                    'general','emergency','ophthalmology','pediatrics','orthopedics','ent',
                    'gynecology','dermatology','surgery','dentistry','psychiatry','oncology')

# Check time series pattern
#### There is a lot of noise. Obvious trends associated with COVID/new wards, 
#### but perhaps other phenomena (i.e. hospital management/operational changes) 
####and possible seasonality. 
ts_plots <- lapply(h_cols_to_plot, function(col) plot_ts(df_w, col))
grid.arrange(grobs = c(ts_plots[c(1,2,4,5)], list(ggplot() + theme_void()), ts_plots[c(3,6,7)]), ncol = 4, nrow = 2)
grid.arrange(grobs = ts_plots[8:19], ncol = 4, nrow = 3)

# Check trend (yearly average)
for (c in h_cols_to_plot){
  trends <- df_w %>% 
    group_by(year) %>%
    summarise(across(all_of(c),
                     list(mean = ~mean(. , na.rm = TRUE),
                          sd = ~sd(. , na.rm = TRUE),
                          min = ~min(. , na.rm = TRUE),
                          max = ~max(. , na.rm = TRUE)),
                     .names = "{.col}_{.fn}"))
  print(trends)
}

(trends <- df_w %>% 
    group_by(year) %>%
    summarise(across(all_of(h_cols_to_plot),
                     list(mean = ~mean(. , na.rm = TRUE)),
                     .names = "{.col}")))

trends_long <- trends %>%
  pivot_longer(cols = -year, names_to = "variable", values_to = "mean_value") %>%
  mutate(variable = factor(variable, levels = h_cols_to_plot))

ggplot(trends_long, aes(x = year, y = mean_value)) +
  geom_line(color = mypalette[2]) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "Yearly Mean Values", x = "Year", y = "Mean Value") +
  theme_minimal()

# Check ACF & PACF
#### ACF: Most have a strong sinusoidal pattern (indicating seasonal trend) and 
#### tapers off slowly (indicating non-stationarity).
#### PACF: Sharp spikes at lags 1-7, which may be important for modelling
for (c in h_cols_to_plot) {
  ts_data <- ts(na.omit(df_w[[c]]), frequency=365)
  par(mfrow = c(1, 2))  # Side-by-side
  acf(ts_data, lag.max=100, main = paste("ACF of", c))
  pacf(ts_data, lag.max=100, main = paste("PACF of", c))
  par(mfrow = c(1, 1))  # Reset layout
}

# Check decomposition
#### Confirms previous tests, showing clear seasonality and trend.
#### Data will be trended for analysis
for (c in h_cols_to_plot){
  # Detrend and plot (removing NAs)
  na_loc <- attr(na.omit(df_w[[c]]), "na.action")
  ts_data <- ts(na.omit(df_w[[c]]), frequency=365)
  stl_decomp <- stl(ts_data, s.window = 'periodic')
  plot(stl_decomp, main= paste('Decomposition of', c))
  # Create detrend variables
  remainder <- stl_decomp$time.series[, 'remainder']
  # Add NAs back in
  # remainder_vector <- as.vector(remainder)
  # for (i in sort(na_loc)) remainder_vector <- append(remainder_vector, NA, after = i - 1)
  # df_w[[paste(c, "_dt", sep = "")]] <- remainder_vector
}


# Comparison of stationary vs non-stationary ACF tests
# lapply(h_cols_to_plot, function(col) plot_stationarity(df_w, col))

#### WEATHER VS HEALTH ####
# library(corrplot)
# selected_df <- df_w[, c('total_patients','<1y', '1-17y', '18-60y', '>60y',
#                         'mn2t', 'mx2t', 'tp', 'fg10')]
# cor_matrix <- cor(selected_df, use = "complete.obs")
# corrplot(cor_matrix, method = "circle")


################################################################################

# # Long data format
# df_l <- read.csv('203_comb/ind_30km_HISP-HP_2015-2023_hosp-era5_chamba_long.csv', check.names = FALSE, fileEncoding = 'UTF-8')
# (column_types <- sapply(df_l, class))
# df_l$date <- as.Date(df_l$date)

# # Get and plot anomalies from detrended data
# dt_cols_to_plot <- grep('_dt', names(df_w), value = TRUE)
# filtered_dt_cols <- dt_cols_to_plot[sub('_dt$', '', dt_cols_to_plot) %in% cols_to_plot]
# detrended_df <- df_w[, c('date', filtered_dt_cols)]
# out <- plot_anomalies(detrended_df)
# anomalies <- out$anomalies
# plots <- out$plots
# 
# grid.arrange(grobs = c(plots[c(1:9)]), ncol=3, nrow=3)
# grid.arrange(grobs = c(plots[c(6,7,10)],list(ggplot() + theme_void()), plots[c(8:9)]), ncol=3, nrow=2)
# 
# 
# for (chunk in split(plots, ceiling(seq_along(plots) / 10))) {
#   grid.arrange(grobs = chunk, ncol = 2)
# }
# 
# # Create binary anomalies variables
# for (var in names(anomalies)) {
#   anom_var_name <- paste0(gsub("_dt", "", var), "_anom")  # e.g., "t2m_anom"
#   df_w[[anom_var_name]] <- 0
#   df_w[[anom_var_name]][anomalies[[var]]] <- 1
# }

# Filter to relevant data & organise
# df_w <- df_w %>% 
#   select(-contains('_'), total_patients) %>%
#   select(location, date, total_patients, everything()) %>%
#   select(-lat, -lon, everything(), lat, lon)

# # Get groups as lists
# wards <- c('total_patients', 'Dentistry','Dermatology','Emergency','ENT','General Medicine',
#            'Gynecology','Oncology','Ophthalmology','Orthopedics','Pediatrics',
#            'Psychiatry','Surgery')
# age <- c('total_patients', '<1y','1-17y','18-60y','>60y' )
# gender <- c('total_patients', 'M', 'F')
#
# # Calculate the moving averages for WARDS
# window_size <- 7
# ma_df <- df_w %>% select(date, Dentistry, ENT, `General Medicine`, Gynecology,
#                          Ophthalmology, Orthopedics, Pediatrics, Surgery,         
#                          Dermatology, Emergency, Psychiatry, Oncology)
# 
# ma_df[,-1] <- apply(ma_df[,-1], 2, function(x) {
#   rollapply(x, width = window_size, FUN = function(z) mean(z, na.rm = TRUE), fill = NA, align = "center")
# })
# long_ma_df <- pivot_longer(ma_df, cols = -date, names_to = 'variable', values_to = 'moving_average')
# ggplot(long_ma_df, aes(x = date, y = moving_average, color = variable)) +
#   geom_line() +
#   labs(title = 'Moving Averages of Variables', x = 'Date', y = 'Moving Average') +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# # Calculate the moving averages for AGES
# window_size <- 7
# ma_df <- df_w %>% select(date,`<1y`, `1-17y`, `18-60y`, `>60y`)
# ma_df[,-1] <- apply(ma_df[,-1], 2, function(x) {
#   rollapply(x, width = window_size, FUN = function(y) mean(y, na.rm = TRUE), 
#             fill = NA, align = 'center')
# })
# long_ma_df <- pivot_longer(ma_df, cols = -date, names_to = 'variable', values_to = 'moving_average')
# ggplot(long_ma_df, aes(x = date, y = moving_average, color = variable)) +
#   geom_line() +
#   labs(title = 'Moving Averages of Variables', x = 'Date', y = 'Moving Average') +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))



