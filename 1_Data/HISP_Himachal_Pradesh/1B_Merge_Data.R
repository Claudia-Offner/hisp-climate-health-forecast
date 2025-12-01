### ------------------------------------------------------------------------ ### 
### HISP weather and hospital data merge                                     ###
### ------------------------------------------------------------------------ ### 

library(here)
setwd(here())
source('1_Data/HISP_Himachal_Pradesh/HISP_Cleaning - Dependencies.R')

# Set data path
setwd(path_data)
weat_import <- 'ind_5km_JASMIN-BigThaw_2015-2023_anom'
heal_import <- 'ind_adm4_HISP-HP_2015-2023_patient_visits_ALL_clean'
cont_import <- 'ind_misc_cont_adm4_pop_elev_2000-2023'
file_export <- 'ind_5km_BigThaw-HISP-ALL_2015-2023_daily_hosp'

####  FUNCTIONS  #####

####  MAIN CODE ####

load(paste0('2_Active_Data/202_weat/', weat_import, '.RData'))
load(paste0('2_Active_Data/201_heal/', heal_import, '.RData'))
load(paste0('2_Active_Data/204_misc/', cont_import, '.RData'))

# Add control information to weat

df_weat <- df_weat %>%
  left_join(elev_res, by='location') %>%
  left_join(pop_res, b=c('location', 'year')) %>%
  arrange(date, location) %>%
  rename(popu=population_density, 
         elev=elevation) %>%
  dplyr::select(location, date, year, season, month, week, 
         weekday, day, elev, popu, everything())


# Prep data for heal-weat joining
df_heal_w$weekday <- NULL
df_heal_l$weekday <- NULL

# Remove lagged/cumulative/percentile/IMD weat variables from long dfs (names ending with 14-21)
check <- df_weat[, !grepl('lag|_c|_p|IMD', names(df_weat))]
names(check)
df_weat_sub <- df_weat %>% dplyr::select(location, date, year, season, month, week, 
                                  weekday, day, elev, popu,
                                  t2m_day, mn2t_day, mx2t_day, 
                                  t2m_night, mn2t_night, mx2t_night, 
                                  any_of(w_cols), 
                                  rh_min, rh_max, sh_min, sh_max)

# 1. Merge weather-health-control data sets ####

con <- dbConnect(duckdb())
duckdb::duckdb_register(con, 'weat_sub', df_weat_sub)
duckdb::duckdb_register(con, 'weat', df_weat)
duckdb::duckdb_register(con, 'heal_w', df_heal_w)
duckdb::duckdb_register(con, 'heal_l', df_heal_l)
duckdb::duckdb_register(con, 'heal', df_heal)

# WIDE FORMAT
df_wide <- dbGetQuery(con, '
  SELECT *
  FROM heal_w
  LEFT JOIN (
    FROM weat
  ) weat USING (location, date)
')

# LONG FORMAT
df_long <- dbGetQuery(con, '
  SELECT *
  FROM heal_l
  LEFT JOIN (
    FROM weat_sub
  ) weat USING (location, date)
')


# ALL PATIENTS
df_all <- dbGetQuery(con, '
  SELECT *
  FROM heal
  LEFT JOIN (
    FROM weat_sub
  ) weat USING (location, date)
')


# Set column order
order <- c('location', 'date', 'year', 'season', 'month', 'week', 'weekday', 'day', 'closed', 'holiday', 'elev', 'popu')
setcolorder(df_wide,  c(order, setdiff(names(df_wide),  order)))
setcolorder(df_long,  c(order, setdiff(names(df_long),  order)))
setcolorder(df_all,   c(order, setdiff(names(df_all),   order)))


# 2. Process wide data for modelling ####

# Set factor levels
df_wide$location <- factor(df_wide$location, levels = location_levels)
df_wide$weekday <- factor(df_wide$weekday, levels = weekday_levels)
df_wide$month <- factor(df_wide$month, levels = month_levels)
df_wide$season <- factor(df_wide$season, levels = season_levels)

# Identify column names that match appropriate exposure suffixes
w_cols <- c(w_cols, 'rh_', 'sh_')
exposure_pattern <- paste(w_cols, collapse = '|')
suffix_pattern   <- paste(c('_IMD', '_p'), collapse = '|')
matched_cols <- names(df_wide)[str_detect(names(df_wide), exposure_pattern) & str_detect(names(df_wide), suffix_pattern)]

# Create cumulative health outcomes
### NB: Set NA values to 1 in the cumulative variable - otherwise the growth curve model gets confused!!
df_wide <- df_wide %>%
  group_by(location) %>%
  arrange(location, date, .by_group = TRUE) %>%
  mutate(across(all_of(h_cols), ~ cumsum(replace_na(., 1)), .names = '{.col}_c'))

# Select the desired columns
df_wide <- df_wide %>%
  dplyr::select(location, date, year, season, month, week, weekday, closed, holiday, elev, popu,
         any_of(c(h_cols, grep('_c$', names(df_wide), value = TRUE),
                  w_cols, matched_cols)))%>%
  arrange(location, date)

# 3. Process all data for formatting ####
df_all$location <- factor(df_all$location, levels = location_levels)
df_all$weekday <- factor(df_all$weekday, levels = weekday_levels)
df_all$month <- factor(df_all$month, levels = month_levels)
df_all$season <- factor(df_all$season, levels = season_levels)
df_all$age_group <- factor(df_all$age_group, levels = age_levels)
df_all$gender <- factor(df_all$gender, levels = gender_levels)

#### EXPORT ####

save(df_wide, file = paste0('2_Active_Data/203_comb/', file_export, '_wide.RData'))
save(df_long, file = paste0('2_Active_Data/203_comb/', file_export, '_long.RData'))
save(df_all, file = paste0('2_Active_Data/203_comb/', file_export, '_all.RData'))

