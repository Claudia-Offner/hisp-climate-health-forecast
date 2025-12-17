### ------------------------------------------------------------------------ ### 
###  TSGC - Tests and cleaning                                               ###
### ------------------------------------------------------------------------ ### 

#### DEPENDENCIES ####

# Turn off scientific notation
options(scipen=999)

# Path
user <- 'claer14'
file_name <- 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide'
path1 <- paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/')
path2 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
path3 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Chapter 3 - Real time forecasts/')
setwd(path2)

# Packages
library(dplyr); library(tidyr); library(lubridate); library(zoo)

# Variables
weekday_levels <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
month_levels <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
season_levels <- c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)', 'Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)')

#### MAIN CODE ####

# Load data
# df <- read.csv(paste0('203_comb/', file_name, '.csv'), check.names=FALSE, fileEncoding='UTF-8')
load('203_comb/ind_5km_BigThaw-HISP-ALL_2015-2023_daily_hosp_wide.Rdata')
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

# Back transform tp variables
df$tp <- exp(df$tp)

# Save as Rdata file
save(df, file = paste0(path1, '3_TSGC/', file_name, '.RData'))
