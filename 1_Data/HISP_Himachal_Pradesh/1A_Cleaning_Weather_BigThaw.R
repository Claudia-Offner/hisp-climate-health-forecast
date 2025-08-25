### ------------------------------------------------------------------------ ### 
### Weather cleaning - BigThaw Data                                          ###
### ------------------------------------------------------------------------ ### 

# Turn off scientific notation
options(scipen=999)

# Set data path
user <- 'offne'
path <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
setwd(path)

####  PACKAGES & FUNCTIONS ####
library(lubridate)
library(dplyr)

####  MAIN CODE ####

# 1. Clean data ####
df_w <- read.csv('202_weat/ind_10km_JASMIN-BigThaw_2015-2016.csv')

# Set hourly time parameters
df_w <- df_w %>%
  mutate(
    timestamp = ymd_hms(time),
    date = as.Date(timestamp),
    year = year(timestamp),
    month = month(timestamp, label = TRUE, abbr = TRUE),
    day = day(timestamp),
    hour = hour(timestamp)
  ) %>%
  select(date, year, month, day, hour, location, everything(), -time)

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
  filter(hour >= 6 & hour < 18) %>%
  group_by(date,location) %>%
  summarise(
    t2m_day = mean(temp_mean, na.rm = TRUE) - 273.15,
    mn2t_day = min(temp_min, na.rm = TRUE) - 273.15,
    mx2t_day = max(temp_max, na.rm = TRUE) - 273.15
  )

# Calculate the night (18:00 to 06:00) mean, min and max (convert from K to C)
temp_night <- df_w %>%
  filter(hour < 6 | hour >= 18) %>%
  group_by(date,location) %>%
  summarise(
    t2m_night = mean(temp_mean, na.rm = TRUE) - 273.15,
    mn2t_night = min(temp_min, na.rm = TRUE) - 273.15,
    mx2t_night = max(temp_max, na.rm = TRUE) - 273.15
  )

# Calculate the overall daily precipitation, including rain and snowfall 
# (convert from m to mm)
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

# Merge processed temp data
df_w <- temp_all %>%
  full_join(temp_day, by = c("date", "location")) %>%
  full_join(temp_night, by = c("date", "location")) %>%
  full_join(precip, by = c("date", "location")) %>%
  full_join(wind, by = c("date", "location"))

# Set daily time parameters
# NB: Seasons were set based on HP crop calendar, pg. 10 (https://himachal.nic.in/WriteReadData/l892s/13_l892s/1722326458.pdf)
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
  select(location, date, year, season, month, week, weekday, day, everything())

###############
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


# take log of tp
df_w$tp <- log(df_w$tp+ 1e-6)

# 3. Calculate climatology anomalies ####
# Not calculated for now because not enough data from BigThaw (2015-2016)

# 
# 
# # Calculate SD anomalies dynamically with col names
# calculate_anomalies <- function(df_w, variables, window_size = 15) {
#   
#   # Ensure 'day' column exists
#   if (!'day' %in% names(df_w)) {
#     stop('Data must contain a "day" column.')
#   }
#   
#   df_result <- df_w
#   
#   for (variable in variables) {
#     # variable <- 'mx2t'
#     # Calculate mean and SD for each DOY window in one go
#     day_stats <- sapply(1:366, function(d) {
#       # Define 31-day window (15 before and 15 after)
#       window_days <- ((d - window_size):(d + window_size) - 1) %% 366 + 1
#       
#       # Filter data matching these days
#       window_data <- df_w %>% filter(day %in% window_days)
#       
#       # Return both mean and SD for the specified variable
#       c(mean = mean(window_data[[variable]], na.rm = TRUE)        )
#     })
#     
#     # Build smooth curve data frame
#     smooth_curve <- data.frame(
#       day = 1:366,
#       smooth_value = day_stats
#     )
#     
#     # Merge smooth values and calculate anomalies
#     x  <- df_result  %>%
#       dplyr::left_join(smooth_curve, by = 'day') %>%
#       dplyr::mutate( # Calculate the anomaly first
#         !!paste0(variable, '_sm') := smooth_value,
#         !!paste0(variable, '_anom') := .data[[variable]] - smooth_value) %>%
#       dplyr::mutate(# Now calculate the other columns using the anomaly
#         !!paste0(variable, '_u1SD') := as.integer(.data[[paste0(variable, '_anom')]] > (1 * sd(.data[[paste0(variable, '_anom')]], na.rm = TRUE))),
#         !!paste0(variable, '_l1SD') := as.integer(.data[[paste0(variable, '_anom')]] < (-1 * sd(.data[[paste0(variable, '_anom')]], na.rm = TRUE))),
#         !!paste0(variable, '_u2SD') := as.integer(.data[[paste0(variable, '_anom')]] > (2 * sd(.data[[paste0(variable, '_anom')]], na.rm = TRUE))),
#         !!paste0(variable, '_l2SD') := as.integer(.data[[paste0(variable, '_anom')]] < (-2 * sd(.data[[paste0(variable, '_anom')]], na.rm = TRUE))),
#         !!paste0(variable, '_p01') := as.integer(.data[[paste0(variable, '_anom')]] < quantile(.data[[paste0(variable, '_anom')]], 0.01, na.rm = TRUE)),
#         !!paste0(variable, '_p05') := as.integer(.data[[paste0(variable, '_anom')]] < quantile(.data[[paste0(variable, '_anom')]], 0.05, na.rm = TRUE)),
#         !!paste0(variable, '_p10') := as.integer(.data[[paste0(variable, '_anom')]] < quantile(.data[[paste0(variable, '_anom')]], 0.10, na.rm = TRUE)),
#         !!paste0(variable, '_p90') := as.integer(.data[[paste0(variable, '_anom')]] > quantile(.data[[paste0(variable, '_anom')]], 0.90, na.rm = TRUE)),
#         !!paste0(variable, '_p95') := as.integer(.data[[paste0(variable, '_anom')]] > quantile(.data[[paste0(variable, '_anom')]], 0.95, na.rm = TRUE)),
#         !!paste0(variable, '_p99') := as.integer(.data[[paste0(variable, '_anom')]] > quantile(.data[[paste0(variable, '_anom')]], 0.99, na.rm = TRUE))
#       ) %>%
#       dplyr::select(-c(smooth_value))
#     
#   }
#   return(df_result)
# }
# 
# 
# w_cols <- c('mn2t','mx2t','tp','fg10')
# 
# locations <- unique(df_w$location)
# result <- data.frame()
# for (l in locations){
#   t <- df_w %>% filter(location==l)
#   t <- calculate_anomalies(t, w_cols)
#   result <- rbind(result, t)
# }

# 4. Create lag variables ####

# Add lags 1-30 for all w_anom_cols (remove the date with NA for certain variables, across locations)
w_anom_cols <- grep('IMD', names(df_w), value = TRUE)
result <- df_w %>% group_by(location) %>% arrange(date) %>% ungroup()
for (i in 1:30) {
  result <- result %>%
    group_by(location) %>%
    arrange(date) %>%
    mutate(across(all_of(w_anom_cols), 
                  ~lag(.x, i), 
                  .names = "{.col}_lag{i}")) %>%
    ungroup()
}


# Export
write.csv(result, file='202_weat/ind_10km_BigThaw_2015-2016_anom.csv', row.names = FALSE)


