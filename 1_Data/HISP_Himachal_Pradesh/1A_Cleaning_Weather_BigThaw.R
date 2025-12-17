### ------------------------------------------------------------------------ ### 
### HP Weather cleaning - BigThaw Data                                       ###
### ------------------------------------------------------------------------ ### 

setwd(here())
source('1_Data/HISP_Himachal_Pradesh/HISP_Cleaning - Dependencies.R')

# Set data paths
setwd(path_data)
pm25_import <- 'ind_adm4_LongPMInd_PM25_2000-2022'
weat_import <- 'ind_5km_JASMIN-BigThaw_2015-2023_ALL'
weat_export <- 'ind_5km_JASMIN-BigThaw_2015-2023_anom'

####  FUNCTIONS  #####

### Function to calculate SD anomalies dynamically with col names
calculate_anomalies <- function(df, variables, window_size = 15, days = 366) {
  # Dependencies: uses dplyr and rlang (no need to attach them if you prefer the :: style)
  if (!'day' %in% names(df)) stop('Data must contain a "day" column.')
  if (!all(variables %in% names(df))) stop('Some variables are not columns in df.')
  if (!is.numeric(window_size) || window_size < 0) stop('window_size must be a non-negative integer.')
  
  # Work from a clean copy of the original data (so smoothing is stable)
  base_df   <- df
  # Normalize day to 1..days just in case
  base_df$day <- ((as.integer(base_df$day) - 1) %% days) + 1
  
  # Start result with original rows
  df_result <- base_df
  
  for (variable in variables) {
    # Compute smooth DOY mean for this variable using the centered window
    day_vec <- seq_len(days)
    day_stats <- vapply(day_vec, function(d) {
      window_days <- ((d - window_size):(d + window_size) - 1) %% days + 1
      vals <- base_df[[variable]][ base_df$day %in% window_days ]
      mean(vals, na.rm = TRUE)
    }, numeric(1))
    
    # Names for new columns
    smooth_col   <- paste0(variable, '_sm')
    anom_col     <- paste0(variable, '_anom')
    p01_col      <- paste0(variable, '_p01')
    p05_col      <- paste0(variable, '_p05')
    p10_col      <- paste0(variable, '_p10')
    p90_col      <- paste0(variable, '_p90')
    p95_col      <- paste0(variable, '_p95')
    p99_col      <- paste0(variable, '_p99')
    
    # Build smooth dataframe and join
    smooth_df <- data.frame(day = day_vec, day_stats, stringsAsFactors = FALSE)
    names(smooth_df)[2] <- smooth_col
    df_result <- dplyr::left_join(df_result, smooth_df, by = 'day')
    
    # Calculate anomaly
    df_result <- dplyr::mutate(
      df_result,
      !!rlang::sym(anom_col) := .data[[variable]] - .data[[smooth_col]]
    )
    
    # Compute thresholds on the anomaly vector (global across all rows for this variable)
    anom_vals <- df_result[[anom_col]]
    thr01 <- stats::quantile(anom_vals, probs = 0.01, na.rm = TRUE)
    thr05 <- stats::quantile(anom_vals, probs = 0.05, na.rm = TRUE)
    thr10 <- stats::quantile(anom_vals, probs = 0.10, na.rm = TRUE)
    thr90 <- stats::quantile(anom_vals, probs = 0.90, na.rm = TRUE)
    thr95 <- stats::quantile(anom_vals, probs = 0.95, na.rm = TRUE)
    thr99 <- stats::quantile(anom_vals, probs = 0.99, na.rm = TRUE)
    
    # Create integer flags (1 if anomaly > threshold, else 0). NA stays NA.
    df_result <- dplyr::mutate(
      df_result,
      !!rlang::sym(p01_col) := ifelse(is.na(.data[[anom_col]]), NA_integer_, as.integer(.data[[anom_col]] < thr01)),
      !!rlang::sym(p05_col) := ifelse(is.na(.data[[anom_col]]), NA_integer_, as.integer(.data[[anom_col]] < thr05)),
      !!rlang::sym(p10_col) := ifelse(is.na(.data[[anom_col]]), NA_integer_, as.integer(.data[[anom_col]] < thr10)),      !!rlang::sym(p90_col) := ifelse(is.na(.data[[anom_col]]), NA_integer_, as.integer(.data[[anom_col]] > thr90)),
      !!rlang::sym(p90_col) := ifelse(is.na(.data[[anom_col]]), NA_integer_, as.integer(.data[[anom_col]] > thr90)),
      !!rlang::sym(p95_col) := ifelse(is.na(.data[[anom_col]]), NA_integer_, as.integer(.data[[anom_col]] > thr95)),
      !!rlang::sym(p99_col) := ifelse(is.na(.data[[anom_col]]), NA_integer_, as.integer(.data[[anom_col]] > thr99))
    )
  }
  
  return(df_result)
}


####  MAIN CODE ####

# Load data
df_w <- read.csv(paste0('1_Original_Data/118_BigThaw/', weat_import, '.csv'))
df_p <- get(load(paste0('2_Active_Data/204_misc/', pm25_import, '.RData')))

# Remove duplicates
df_w <- df_w %>% distinct()

# Set hourly time parameters
df_p$date <- as.Date(df_p$date)
df_w <- df_w %>%
  mutate(
    timestamp = ymd_hms(time),
    date = as.Date(timestamp),
    year = year(timestamp),
    month = lubridate::month(timestamp, label = TRUE, abbr = TRUE),
    day = day(timestamp),
    hour = hour(timestamp)
  ) %>%
  dplyr::select(date, year, month, day, hour, location, everything(), -time)

# 1. Clean data ####

# Calculate the overall day and night mean, min and max (convert from K to C)
temp_all <- df_w %>%
  group_by(date,location) %>%
  summarise(
    t2m = mean(temp_mean, na.rm = TRUE) - 273.15,
    mn2t = min(temp_min, na.rm = TRUE) - 273.15,
    mx2t = max(temp_max, na.rm = TRUE) - 273.15
  )

# Calculate the day (06:00 to 18:00) mean, min and max (convert from K to C)
temp_day <- df_w %>%
  dplyr::filter(hour >= 6 & hour < 18) %>%
  group_by(date,location) %>%
  summarise(
    t2m_day = mean(temp_mean, na.rm = TRUE) - 273.15,
    mn2t_day = min(temp_min, na.rm = TRUE) - 273.15,
    mx2t_day = max(temp_max, na.rm = TRUE) - 273.15
  )

# Calculate the night (18:00 to 06:00) mean, min and max (convert from K to C)
temp_night <- df_w %>%
  dplyr::filter(hour < 6 | hour >= 18) %>%
  group_by(date,location) %>%
  summarise(
    t2m_night = mean(temp_mean, na.rm = TRUE) - 273.15,
    mn2t_night = min(temp_min, na.rm = TRUE) - 273.15,
    mx2t_night = max(temp_max, na.rm = TRUE) - 273.15
  )

# Calculate the overall daily precipitation, including rain and snowfall (convert from m to mm)
precip <- df_w %>%
  group_by(date,location) %>%
  summarise(
    tp = sum(rain_max+snow_max, na.rm = TRUE)*1000
  )

# Aggregate wind to daily levels (convert to km per hour)
wind <- df_w %>%
  group_by(date,location) %>%
  summarise(
    fg10 = max(wind_max, na.rm = TRUE) * 3.6
  )

# Select all humidity variables
humid <- df_w %>%
  dplyr::select(date, location, contains('humid', ignore.case = TRUE)) %>%
  rename_with(~ str_replace_all(., regex('rel_humid', ignore_case = TRUE), 'rh')) %>%
  rename_with(~ str_replace_all(., regex('spec_humid', ignore_case = TRUE), 'sh')) %>%
  group_by(date,location) %>%
  summarise(
    rh_mean = mean(rh_mean, na.rm = TRUE),
    rh_min = min(rh_min, na.rm = TRUE),
    rh_max = max(rh_max, na.rm = TRUE),
    sh_mean = mean(sh_mean, na.rm = TRUE),
    sh_min = min(sh_min, na.rm = TRUE),
    sh_max = max(sh_max, na.rm = TRUE),
  )
  
# Merge processed temp data
df_w <- temp_all %>%
  full_join(temp_day, by = c('date', 'location')) %>%
  full_join(temp_night, by = c('date', 'location')) %>%
  full_join(precip, by = c('date', 'location')) %>%
  full_join(wind, by = c('date', 'location')) %>%
  full_join(humid, by = c('date', 'location')) %>%
  full_join(df_p, by = c('date', 'location'))
  
# Set daily time parameters
# NB: Seasons were set based on HP crop calendar, pg. 10 
# (https://himachal.nic.in/WriteReadData/l892s/13_l892s/1722326458.pdf)
df_w <- df_w %>%
  mutate(
    year = year(date),
    month = lubridate::month(date, label = TRUE, abbr = TRUE),
    week = isoweek(date),
    weekday = lubridate::wday(date, label = TRUE, abbr = TRUE),
    day = yday(date),
    season = case_when(
      month %in% c('Jan', 'Feb') ~ 'Winter (Jan-Feb)',
      month %in% c('Mar', 'Apr', 'May') ~ 'Pre-Monsoon (Mar-May)',
      month %in% c('Jun', 'Jul', 'Aug', 'Sep') ~ 'Monsoon (Jun-Sep)',
      month %in% c('Oct', 'Nov', 'Dec') ~ 'Post-Monsoon (Oct-Dec)'
    )) %>%
  dplyr::select(location, date, year, season, month, week, weekday, day, everything())

# *** Quick extreme temperature checks *** ####
# summary(df_w$t2m)
# summary(df_w$mx2t)
# summary(df_w$mn2t)
# filtered_df <- df_w[df_w$location == 'Chamba', ]
# summary(filtered_df$t2m)
# summary(filtered_df$mx2t)
# summary(filtered_df$mn2t)

# 2. Set thresholds provided by IMD ####
# Source: https://mausam.imd.gov.in/imd_latest/contents/pdf/forecasting_sop.pdf

# TEMPERATURE 
mx2t_IMD <- 30
mx2t_IMD_sd <- c(4.5, 6.4)
mn2t_IMD <- 0
mn2t_IMD_sd <- c(-4.5, -6.4)
df_w <- df_w %>%
  mutate(
    mx2t_IMD = case_when(
      mx2t >= (mx2t_IMD + mx2t_IMD_sd[2]) ~ 2, # severe_heatwave
      mx2t >= mx2t_IMD ~ 1, # heatwave
      TRUE ~ 0 # normal
    ),
    mx2t_IMD_b = if_else(mx2t_IMD >=1, 1, 0),
    mn2t_IMD = case_when(
      mn2t <= mn2t_IMD + mn2t_IMD_sd[2] ~ 2, # severe_coldwave
      mn2t <= mn2t_IMD ~ 1, # coldwave
      TRUE ~ 0 # normal
    ),
    mn2t_IMD_b = if_else(mn2t_IMD >=1, 1, 0)
  )

# RAINFALL
tp_IMD <- list(
  very_light = 2.4,
  light = c(2.4, 15.5),
  moderate = c(15.5, 64.4),
  heavy = c(64.4, 115.5),
  very_heavy = c(115.5, 204.4),
  ex_heavy = 204.4
)
df_w <- df_w %>%
  mutate(tp_IMD = case_when(
    tp <= 0 ~ 0,  # no rain or zero
    tp > 0 & tp <= tp_IMD$very_light ~ 1,
    tp > tp_IMD$light[1] & tp <= tp_IMD$light[2] ~ 2,
    tp > tp_IMD$moderate[1] & tp <= tp_IMD$moderate[2] ~ 3,
    tp > tp_IMD$heavy[1] & tp <= tp_IMD$heavy[2] ~ 4,
    tp > tp_IMD$very_heavy[1] & tp <= tp_IMD$very_heavy[2] ~ 5,
    tp > tp_IMD$ex_heavy ~ 6
  ),
  tp_IMD_b = if_else(tp_IMD >=3, 1, 0))

# STORM
fg10_IMD <- list(light=40,
                 moderate=c(41,61),
                 severe=c(62,87),
                 very_severe=87)
df_w <- df_w %>%
  mutate(
    fg10_IMD = case_when(
      fg10 <= fg10_IMD$light ~ 1,                            # light: ≤ 40
      fg10 > fg10_IMD$light & fg10 <= fg10_IMD$moderate[2] ~ 2,  # moderate: >40 to 61
      fg10 > fg10_IMD$moderate[2] & fg10 <= fg10_IMD$severe[2] ~ 3, # severe: >61 to 87
      fg10 > fg10_IMD$very_severe ~ 4,                        # very_severe: >87
      TRUE ~ 0                                               
    ),
    fg10_IMD_b = if_else(fg10_IMD >=2, 1, 0)
    )


# NOT IMD: Humidity
df_w <- df_w %>%
  mutate(rh_IMD = case_when(
    rh_mean < 30 ~ -1, # Low (30%)
    rh_mean >= 30 & rh_mean <= 60 ~ 0, # Optimal (30-60%)
    rh_mean > 60 ~ 1, # High (>60%)
    TRUE ~ NA_real_),
    rh_low_IMD_b = if_else(rh_IMD < 0, 1, 0),
    rh_hig_IMD_b = if_else(rh_IMD > 0, 1, 0)
    ) %>%
  mutate(sh_IMD = case_when(
    sh_mean < 5 ~ -1, #'Low (<5g/kg)',
    sh_mean >= 5 & sh_mean <= 15 ~ 0, #'Optimal (5-15g/kg)',
    sh_mean > 15 ~ 1, #'High (>15g/kg)',
    TRUE ~ NA_real_),
    sh_low_IMD_b = if_else(sh_IMD < 0, 1, 0),
    sh_hig_IMD_b = if_else(sh_IMD > 0, 1, 0)
    )

# NOT IMD: POLLUTION
# https://safar.tropmet.res.in/AQI-47-12-Details
df_w <- df_w %>%
  mutate(pm25_IMD = case_when(
    pm25 < 30 ~ 0, # 'Good',
    pm25 >= 30 & pm25 <= 60 ~ 1, #'Moderate',
    pm25 >= 60 & pm25 <= 90 ~ 2, #'Poor',
    pm25 >= 90 & pm25 <= 120 ~ 3, #'Unhealthy',
    pm25 >= 120 & pm25 <= 250 ~ 4, #'Severe',
    pm25 > 250 ~ 5, #'Hazardous',
    TRUE ~ NA_real_),
    pm25_IMD_b = if_else(pm25_IMD >= 2, 1, 0),
    )


# 3. Calculate climatology anomalies ####

result <- data.frame()
for (l in location_levels){
  t <- df_w %>% dplyr::filter(location==l)
  t <- calculate_anomalies(t, w_cols)
  result <- rbind(result, t)
}

# Check extremes across locations and years (+ NA check)
(na_counts <- sapply(result, function(x) sum(is.na(x))))
df_na <- result %>% dplyr::filter(is.na(pm25)) # check cols
df_w <- result

# # Check distribution across time and location
# table(lubridate::year(result$date)[result$t2m_p90 == 1], result$location[result$t2m_p90 == 1])
# dplyr::filtered_result <- result %>% dplyr::filter(location == 'Chamba')
# table(dplyr::filtered_result$year, dplyr::filtered_result$t2m_p90)

# 4. Create lag variables ####

# Add lags 1-21 for all w_anom_cols (remove the date with NA for certain variables, across locations)
w_anom_cols <- grep('IMD|_p', names(df_w), value = TRUE)
df_weat <- df_w %>% group_by(location) %>% arrange(date) %>% ungroup()
for (i in 1:21) {
  df_weat <- df_weat %>%
    group_by(location) %>%
    arrange(date) %>%
    mutate(across(all_of(w_anom_cols), 
                  ~lag(.x, i), 
                  .names = '{.col}_lag{i}')) %>%
    ungroup()
}


#### EXPORT ####
# write.csv(result, file=paste0('2_Active_Data/202_weat/', weat_export, '.csv'), row.names = FALSE)
save(df_weat, file=paste0('2_Active_Data/202_weat/', weat_export, '.RData'))

