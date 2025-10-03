### ------------------------------------------------------------------------ ### 
### Weather cleaning - ERA5 Data                                          ###
### ------------------------------------------------------------------------ ### 

# Turn off scientific notation
options(scipen=999)

# Set data path
user <- 'offne'
path <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
setwd(path)

####  PACKAGES & FUNCTIONS ####
library(dplyr)
library(tidyverse)
library(lubridate)
library(data.table)

# Calculate SD anomalies dynamically with col names
calculate_anomalies <- function(df_w, variables, window_size = 15) {
  
  # Ensure 'day' column exists
  if (!'day' %in% names(df_w)) {
    stop('Data must contain a "day" column.')
  }
  
  df_result <- df_w
  
  for (variable in variables) {
    
    # Calculate mean and SD for each DOY window in one go
    day_stats <- sapply(1:366, function(d) {
      # Define 31-day window (15 before and 15 after)
      window_days <- ((d - window_size):(d + window_size) - 1) %% 366 + 1
      
      # Filter data matching these days
      window_data <- df_w %>% filter(day %in% window_days)
      
      # Return both mean and SD for the specified variable
      c(mean = mean(window_data[[variable]], na.rm = TRUE)        )
    })
    
    # Build smooth curve data frame
    smooth_curve <- data.frame(
      day = 1:366,
      smooth_value = day_stats
    )
    
    # Merge smooth values and calculate anomalies
    df_result  <- df_result  %>%
      left_join(smooth_curve, by = 'day') %>%
      dplyr::mutate( # Calculate the anomaly first
        !!paste0(variable, '_sm') := smooth_value,
        !!paste0(variable, '_anom') := .[[variable]] - smooth_value) %>%
      dplyr::mutate(# Now calculate the other columns using the anomaly
        !!paste0(variable, '_u1SD') := as.integer(.[[paste0(variable, '_anom')]] > (1 * sd(.[[paste0(variable, '_anom')]], na.rm = TRUE))),
        !!paste0(variable, '_l1SD') := as.integer(.[[paste0(variable, '_anom')]] < (-1 * sd(.[[paste0(variable, '_anom')]], na.rm = TRUE))),
        !!paste0(variable, '_u2SD') := as.integer(.[[paste0(variable, '_anom')]] > (2 * sd(.[[paste0(variable, '_anom')]], na.rm = TRUE))),
        !!paste0(variable, '_l2SD') := as.integer(.[[paste0(variable, '_anom')]] < (-2 * sd(.[[paste0(variable, '_anom')]], na.rm = TRUE))),
        !!paste0(variable, '_p01') := as.integer(.[[paste0(variable, '_anom')]] < quantile(.[[paste0(variable, '_anom')]], 0.01, na.rm = TRUE)),
        !!paste0(variable, '_p05') := as.integer(.[[paste0(variable, '_anom')]] < quantile(.[[paste0(variable, '_anom')]], 0.05, na.rm = TRUE)),
        !!paste0(variable, '_p10') := as.integer(.[[paste0(variable, '_anom')]] < quantile(.[[paste0(variable, '_anom')]], 0.10, na.rm = TRUE)),
        !!paste0(variable, '_p90') := as.integer(.[[paste0(variable, '_anom')]] > quantile(.[[paste0(variable, '_anom')]], 0.90, na.rm = TRUE)),
        !!paste0(variable, '_p95') := as.integer(.[[paste0(variable, '_anom')]] > quantile(.[[paste0(variable, '_anom')]], 0.95, na.rm = TRUE)),
        !!paste0(variable, '_p99') := as.integer(.[[paste0(variable, '_anom')]] > quantile(.[[paste0(variable, '_anom')]], 0.99, na.rm = TRUE))
      ) %>%
      dplyr::select(-c(smooth_value))
    
    }
  return(df_result)
}

####  MAIN CODE ####

# 1. Clean data ####

# Load data
df_w <- read.csv('202_weat/ind_30km_JASMIN-ERA5_1979-2023_HISP.csv')

# Remove duplicates (quickly)
df_w <- unique(as.data.table(df_w))
rownames(df_w) <- NULL

# Set time parameters
# NB: Seasons were set based on HP crop calendar, pg. 10 (https://himachal.nic.in/WriteReadData/l892s/13_l892s/1722326458.pdf)
(column_types <- sapply(df_w, class))
df_w <- df_w %>% rename(date = time)
df_w$date <- as.Date(df_w$date)
df_w <- df_w %>%
  dplyr::mutate(
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

# Fix col locations & set temperature values to Celsius
df_w <- df_w %>% relocate(fg10 , .before = v10 ) %>% relocate(d2m , .after = t2m ) %>% relocate(fg10 , .after = tp )
temp <- c('t2m', 'd2m', 'mn2t', 'mx2t')
for (t in temp) {
  df_w[[t]] <- df_w[[t]]-273.15
}

# Convert rainfall from meters to mm
df_w$tp <- df_w$tp*1000

# Convert max wind speed to km per hour
df_w$fg10 <- df_w$fg10 * 3.6

# Check missing data
(na_counts <- sapply(df_w, function(x) sum(is.na(x))))
df_with_na <- df_w %>% filter(if_any(everything(), ~ is.na(.)))
df_w <- na.omit(df_w)

################

# summary(df_w$t2m)
# summary(df_w$mx2t)
# summary(df_w$mn2t)
# filtered_df <- df_w[df_w$location == "Chamba", ]
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
    fg10_IMD_b = if_else(fg10_IMD >=2, 1, 0))


# 3. Calculate climatology anomalies ####

# Transform tp and sp
df_w$tp <- log(df_w$tp+ 1e-6)

# For each location, calculate anomalies from 50-year moving average (+/- 15 days).
w_cols <- c('t2m','d2m','mn2t','mx2t','tp','fg10','v10','u10','sp')
locations <- unique(df_w$location)
result <- data.frame()
for (l in locations){
  t <- df_w %>% filter(location==l)
  t <- calculate_anomalies(t, w_cols)
  result <- rbind(result, t)
}

# # Check extremes across locations and years (+ NA check)
# (na_counts <- sapply(result, function(x) sum(is.na(x))))

# # Check distribution across time and location
# table(lubridate::year(result$date)[result$t2m_p10 == 1], result$location[result$t2m_p10 == 1])
# filtered_result <- t %>% filter(location == 'Chamba')
# table(filtered_result$year, filtered_result$t2m_p10)

# 4. Create lag variables ####

# Add lags 1-30 for all w_anom_cols (remove the date with NA for certain variables, across locations)
w_anom_cols <- grep('_', names(result), value = TRUE)
for (i in 1:30) {
  result <- result %>%
    group_by(location) %>%
    arrange(date) %>%
    mutate(across(all_of(w_anom_cols), 
                  ~lag(.x, i), 
                  .names = '{.col}_lag{i}'))%>%
    ungroup() %>%
    arrange(location, date) 
}


# Export
write.csv(result, file='202_weat/ind_30km_JASMIN-ERA5_1979-2023_HISP_anom.csv', row.names = FALSE)

