### ------------------------------------------------------------------------ ### 
### HISP weather and hospital data merge                                     ###
### ------------------------------------------------------------------------ ### 
# NB COLOURS: https://r-charts.com/colors/

#### DEPENDENCIES ####

# Turn off scientific notation
options(scipen=999)

# Set data path
user <- 'claer14'
weat_import <- 'ind_5km_JASMIN-BigThaw_2015-2023_anom.csv'
heal_import <- 'ind_adm4_HISP-HP_2015-2023_patient_visits_chamba_clean'
file_export <- 'ind_5km_BigThaw-HISP-chamba_2015-2023_hosp'
path <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
setwd(path)

# Packages
library(dplyr); library(tidyverse); library(lubridate)

####  MAIN CODE ####

# 1. Merge weather-health data sets ####
df1_daily_weat <- read.csv(paste0('202_weat/', weat_import))
df_heal_w <- read.csv(paste0('201_heal/', heal_import, '_wide.csv'), check.names = FALSE, fileEncoding = 'UTF-8')
df_heal_l <- read.csv(paste0('201_heal/', heal_import, '_long.csv'), check.names = FALSE, fileEncoding = 'UTF-8')
df_heal <- read.csv(paste0('201_heal/', heal_import, '.csv'), check.names = FALSE, fileEncoding = 'UTF-8')

# Fix dtypes
(column_types <- sapply(df1_daily_weat, class))
(column_types <- sapply(df_heal_w, class))
(column_types <- sapply(df_heal_l, class))
(column_types <- sapply(df_heal, class))
df1_daily_weat$date <- as.Date(df1_daily_weat$date)
df_heal_w$date <- as.Date(df_heal_w$date)
df_heal_l$date <- as.Date(df_heal_l$date)
df_heal$date <- as.Date(df_heal$date)

df_wide <- df_heal_w %>%
  dplyr::select(-weekday) %>%
  left_join(df1_daily_weat %>% distinct(location, date, .keep_all = TRUE), by = c('location', 'date')) %>%
  # dplyr::select(-lat, -lon, everything(), lat, lon) %>%
  dplyr::select(location, date, year, season, month, week, weekday, day, everything())

df_long <- df_heal_l %>%
  select(-weekday) %>%
  left_join(df1_daily_weat %>% distinct(location, date, .keep_all = TRUE), by = c('location', 'date')) %>%
  # select(-lat, -lon, everything(), lat, lon) %>%
  select(location, date, year, season, month, week, weekday, day, everything())

# Join dfs
df_all <- df_heal %>%
  left_join(df1_daily_weat %>% distinct(location, date, .keep_all = TRUE), by = c('location', 'date')) %>%
  # select(-lat, -lon, everything(), lat, lon)  %>%
  select(location, date, year, season, month, week, weekday, day, everything())


# 2. Subset data for model ####

# Select variables of interest
w_cols <- c('t2m', 'mn2t','mx2t','tp','fg10')
h_cols <- c('total_patients','F','M','under_5','age5_14','age15_24','age25_64','over_65',
            'general','emergency','ophthalmology','pediatrics','orthopedics','ent',
            'gynecology','dermatology','surgery','dentistry','psychiatry','oncology')

# Identify column names that match appropriate exposure suffixes
exposure_pattern <- paste(w_cols, collapse = '|')
suffix_pattern   <- paste(c('_IMD', '_p'), collapse = '|')
matched_cols <- names(df_wide)[str_detect(names(df_wide), exposure_pattern) & str_detect(names(df_wide), suffix_pattern)]

# Create cumulative health outcomes
### NB: Set NA values to 1 in the cumulative variable - otherwise the growth
### curve model gets confused!!
df_wide <- df_wide %>% mutate(across(all_of(h_cols), ~ cumsum(replace_na(., 1)), .names = '{.col}_c'))

# Select the desired columns
df_wide <- df_wide %>%
  dplyr::select(location, date, year, season, month, week, weekday, closed, holiday,
         any_of(c(h_cols, grep('_c', names(df_wide), value = TRUE),
                  w_cols, matched_cols)))

# Export
write.csv(df_wide, file=paste0('203_comb/', file_export, '_wide.csv'), row.names = FALSE)
# write.csv(df_long, file=paste0('203_comb/', file_export, '_long.csv'), row.names = FALSE)
# write.csv(df_all, file=paste0('203_comb/', file_export, '_all.csv'), row.names = FALSE)
