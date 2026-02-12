### ------------------------------------------------------------------------ ### 
###  TSGC - Data Cleaning
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list = ls()); graphics.off(); gc()
setwd('C:/Users/claer14/Documents/GitHub/hisp-climate-health-forecast/')


#### DEPENDENCIES ####
source('3_TSGC/TSGC - Dependencies.R')
setwd(path_data)
file_export <- 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide'


#### MAIN CODE ####

# Load data
load('2_Active_Data/203_comb/ind_5km_BigThaw-HISP-ALL_2015-2023_daily_hosp_wide.Rdata')
df <- df_wide

# Set date type, 
df$date <- as.Date(df$date, format='%Y-%m-%d')  # Adjust format as needed

# Set factor data
df <- df %>%
  mutate(weekday=factor(weekday, levels=weekday_levels),
         month=factor(month, levels=month_levels), 
         season=factor(season, levels=season_levels),
         year=factor(year))

# Select year of interest
startdate <- as.Date('2015-02-14')
enddate <- as.Date('2020-02-13')
df <- df[df$date >= startdate & df$date <= enddate,]


# Remove Bilaspur
df <- df %>% filter(location != 'Bilaspur')

# Save as Rdata file
save(df, file = paste0(path_git, '3_TSGC/', file_export, '.RData'))

