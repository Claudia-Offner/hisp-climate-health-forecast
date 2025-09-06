### ------------------------------------------------------------------------ ### 
###  DNLM - Tests and cleaning                                               ###
### ------------------------------------------------------------------------ ### 

#### DEPENDENCIES ####

# Turn off scientific notation
options(scipen=999)

# Path
user <- 'offne'
file_name <- 'ind_5km_BigThaw-HISP-chamba_2015-2023_hosp_wide'
path1 <- paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/')
path2 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
path3 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Chapter 3 - Real time forecasts/')
setwd(path2)

# Packages
library(dplyr); library(tidyr); library(lubridate)

# Variables
weekday_levels <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
month_levels <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
season_levels <- c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)', 'Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)')

#### MAIN CODE ####

# Load data
df <- read.csv(paste0('203_comb/', file_name, '.csv'), check.names=FALSE, fileEncoding='UTF-8')

# Set date type and factors 
df$date <- as.Date(df$date, format='%Y-%m-%d')  # Adjust format as needed
df <- df %>%
  mutate(weekday=factor(weekday, levels=weekday_levels),
         month=factor(month, levels=month_levels), 
         season=factor(season, levels=season_levels))

# Select variables of interest
df <- df %>%
  mutate(doy = yday(date),
         time = row_number()) %>%
  select(
    location, time, date, year, season, month, week, weekday, doy, closed, holiday, total_patients, 
    F, M, under_5, age5_14, age15_24, age25_64, over_65,
    general, emergency, ophthalmology, pediatrics, orthopedics, ent, gynecology, dermatology, surgery, dentistry, psychiatry, oncology, 
    t2m, mn2t, mx2t, tp, fg10) #,
    # mx2t_IMD, mx2t_IMD_b, mn2t_IMD, mn2t_IMD_b, tp_IMD, tp_IMD_b, fg10_IMD, fg10_IMD_b, 
    # t2m_p10, t2m_p90, mn2t_p10, mn2t_p90, mx2t_p10, mx2t_p90, tp_p10, tp_p90, fg10_p10, fg10_p90

df$tp <- exp(df$tp)
summary(df)

# Save as Rdata file
save(df, file = paste0(path1, '2_DLNM/', file_name, '.RData'))
