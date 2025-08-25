### ------------------------------------------------------------------------ ### 
### HISP weather and hospital data merge                                     ###
### ------------------------------------------------------------------------ ### 
# NB COLOURS: https://r-charts.com/colors/

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

####  MAIN CODE ####

# 1. Merge weather-health data sets ####
# df1_daily_weat <- read.csv('202_weat/ind_30km_JASMIN-ERA5_1979-2023_HISP_anom.csv')
df1_daily_weat <- read.csv('202_weat/ind_10km_BigThaw_2015-2016_anom.csv')

# df_heal <- read.csv('201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba_clean.csv', check.names = FALSE, fileEncoding = 'UTF-8')
df_heal_w <- read.csv('201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba_clean_wide.csv', check.names = FALSE, fileEncoding = 'UTF-8')
# df_heal_l <- read.csv('201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba_clean_long.csv', check.names = FALSE, fileEncoding = 'UTF-8')

# Fix dtypes
# (column_types <- sapply(df_heal, class))
(column_types <- sapply(df_heal_w, class))
# (column_types <- sapply(df_heal_l, class))
(column_types <- sapply(df1_daily_weat, class))
df1_daily_weat$date <- as.Date(df1_daily_weat$date)
# df_heal$date <- as.Date(df_heal$date)
df_heal_w$date <- as.Date(df_heal_w$date)
# df_heal_l$date <- as.Date(df_heal_l$date)

# # Join dfs
# df_all <- df_heal %>%
#   left_join(df1_daily_weat %>% distinct(location, date, .keep_all = TRUE), by = c('location', 'date')) %>%
#   select(-lat, -lon, everything(), lat, lon)  %>%
#   select(location, date, year, season, month, week, weekday, day, everything())
# 
# df_long <- df_heal_l %>%
#   select(-weekday) %>%
#   left_join(df1_daily_weat %>% distinct(location, date, .keep_all = TRUE), by = c('location', 'date')) %>%
#   select(-lat, -lon, everything(), lat, lon) %>%
#   select(location, date, year, season, month, week, weekday, day, everything())

df_wide <- df_heal_w %>%
  dplyr::select(-weekday) %>%
  left_join(df1_daily_weat %>% distinct(location, date, .keep_all = TRUE), by = c('location', 'date')) %>%
  # dplyr::select(-lat, -lon, everything(), lat, lon) %>%
  dplyr::select(location, date, year, season, month, week, weekday, day, everything())


# Export
# write.csv(df_all, file='203_comb/ind_30km_HISP-HP_2015-2023_hosp-era5_chamba_all.csv', row.names = FALSE)
# write.csv(df_long, file='203_comb/ind_30km_HISP-HP_2015-2023_hosp-era5_chamba_long.csv', row.names = FALSE)
# write.csv(df_wide, file='203_comb/ind_30km_HISP-HP_2015-2023_hosp-era5_chamba_wide.csv', row.names = FALSE)
write.csv(df_wide, file='203_comb/ind_10km_HISP-HP_2015-2027_hosp-bigthaw_chamba_wide.csv', row.names = FALSE)


# 2. Subset data for model####

# # Load data
# df_wide <- read.csv('203_comb/ind_30km_HISP-HP_2015-2023_hosp-era5_chamba_wide.csv', check.names = FALSE, fileEncoding = 'UTF-8')
# (column_types <- sapply(df_wide, class))
# df_wide$date <- as.Date(df_wide$date)

# Select variables of interest
w_cols <- c('t2m', 'mn2t','mx2t','tp','fg10')
h_cols <- c('total_patients','F','M','under_5','age5_14','age15_24','age25_64','over_65',
            'general','emergency','ophthalmology','pediatrics','orthopedics','ent',
            'gynecology','dermatology','surgery','dentistry','psychiatry','oncology')

# Identify column names that match appropriate exposure suffixes
exposure_pattern <- paste(w_cols, collapse = '|')
suffix_pattern   <- paste(c('_IMD', '_p10', '_p90'), collapse = '|')
matched_cols <- names(df_wide)[str_detect(names(df_wide), exposure_pattern) & str_detect(names(df_wide), suffix_pattern)]

# Create cumulative health outcomes
### NB: Set NA values to 1 in the cumulative variable - otherwise the growth
### curve model gets confused!!
df_wide <- df_wide %>% mutate(across(all_of(h_cols), ~ cumsum(replace_na(., 1)), .names = '{.col}_c'))

# Select the desired columns
res <- df_wide %>%
  dplyr::select(location, date, year, season, month, week, weekday, closed, holiday,
         any_of(c(h_cols, grep('_c', names(df_wide), value = TRUE),
                  w_cols, matched_cols)))

# Export model df
# write.csv(res, file='203_comb/ind_30km_HISP-HP_2015-2023_hosp-era5_chamba_wide_model.csv', row.names = FALSE)
write.csv(res, file='203_comb/ind_10km_HISP-HP_2015-2017_hosp-bigthaw_chamba_wide_model.csv', row.names = FALSE)
