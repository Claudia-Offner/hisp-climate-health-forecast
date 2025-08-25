### ------------------------------------------------------------------------ ### 
###  DNLM - Descriptive Analysis                                             ###
### ------------------------------------------------------------------------ ### 

#### DEPENDENCIES ####

# Turn off scientific notation
options(scipen=999)

# Path
user <- 'claer14'
path1 <- paste0('C:/Users/', user, '/Documents/GitHub/HISP_TSGC/')
path2 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
path3 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Chapter 3 - Real time forecasts/')
setwd(path2)

# Packages
library(dlnm)
library(splines)
library(dplyr)
library(lubridate)

# Functions



# Variables
weekday_levels <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
month_levels <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
season_levels <- c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)', 'Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)')


#### MAIN CODE ####

# 1. Load data ####

# 2. FIG: Study area ####

# 3. TAB: Data characteristics ####


