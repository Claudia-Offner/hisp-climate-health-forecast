### ------------------------------------------------------------------------ ### 
### HISP weather and hospital data merge                                     ###
### ------------------------------------------------------------------------ ### 
# NB COLOURS: https://r-charts.com/colors/

#### DEPENDENCIES ####

# Turn off scientific notation
options(scipen=999)

# Set data path
user <- 'claer14'
weat_import <- 'ind_5km_JASMIN-BigThaw_2015-2023_anom'
heal_import <- 'ind_adm4_HISP-HP_2015-2023_patient_visits_ALL_clean'
file_export <- 'ind_5km_BigThaw-HISP-ALL_2015-2023_daily_hosp'
path <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
setwd(path)

# Packages
library(dplyr); library(tidyverse); library(lubridate); library(data.table);
library(duckdb)


####  MAIN CODE ####

load(paste0('202_weat/', weat_import, '.RData'))
load(paste0('201_heal/', heal_import, '.RData'))

# Prep data for joining
df_heal_w$weekday <- NULL
df_heal_l$weekday <- NULL

# Find columns to exclude from weat for large dfs (names ending with 14-21)
df_weat_sub <- df_weat %>% select(location, date, year, season, month, week, weekday, day,
                                  t2m, mn2t, mx2t, t2m_day, mn2t_day, mx2t_day, 
                                  t2m_night, mn2t_night, mx2t_night, 
                                  tp, fg10, rh_mean, rh_min, rh_max, 
                                  sh_mean, sh_min, sh_max, pm25)

# 1. Merge weather-health data sets ####

con <- dbConnect(duckdb())
duckdb::duckdb_register(con, "weat_sub", df_weat_sub)
duckdb::duckdb_register(con, "weat", df_weat)
duckdb::duckdb_register(con, "heal_w", df_heal_w)
duckdb::duckdb_register(con, "heal_l", df_heal_l)
duckdb::duckdb_register(con, "heal", df_heal)

# WIDE FORMAT
df_wide <- dbGetQuery(con, "
  SELECT *
  FROM heal_w
  LEFT JOIN (
    FROM weat
  ) weat USING (location, date)
")

# LONG FORMAT
df_long <- dbGetQuery(con, "
  SELECT *
  FROM heal_l
  LEFT JOIN (
    FROM weat_sub
  ) weat USING (location, date)
")


# ALL PATIENTS
df_all <- dbGetQuery(con, "
  SELECT *
  FROM heal
  LEFT JOIN (
    FROM weat_sub
  ) weat USING (location, date)
")


# Set column order
desired_order <- c("location", "date", "year", "season", "month", "week", "weekday", "day")
setcolorder(df_wide,  c(desired_order, setdiff(names(df_wide),  desired_order)))
setcolorder(df_long,  c(desired_order, setdiff(names(df_long),  desired_order)))
setcolorder(df_all,   c(desired_order, setdiff(names(df_all),   desired_order)))


# 2. Process wide data for modelling ####

# Select variables of interest
w_cols <- c('t2m', 'mn2t','mx2t','tp','fg10', 'rh_mean', 'sh_mean', 'pm25')
h_cols <- c('total_patients', 'general', 'emergency', 'ophthalmology', 'ent', 'paediatrics',
            'gynaecology', 'orthopedics', 'dermatology', 'surgery',
            'psychiatry', 'pulmonary', 'disability', 'dentistry', 'oncology',
            'ncd', 'antiretroviral', 'substance_counselling', 'physiotherapy')

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
         any_of(c(h_cols, grep('_c$', names(df_wide), value = TRUE),
                  w_cols, matched_cols)))%>%
  arrange(location, date)

# Export ####

save(df_wide, file = paste0('203_comb/', file_export, '_wide.RData'))
save(df_long, file = paste0('203_comb/', file_export, '_long.RData'))
save(df_all, file = paste0('203_comb/', file_export, '_all.RData'))

