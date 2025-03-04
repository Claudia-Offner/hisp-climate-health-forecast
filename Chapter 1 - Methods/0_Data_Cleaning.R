### --------------------------------------- ### 
### Simulating Daily Health Data from Scotland ###
### --------------------------------------- ### 


# SET WORKING DIRECTORY
setwd('C:/Users/offne/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')

# 0. Load dependencies ####
# install.packages('ggmap')

library(readxl)
library(dplyr)
library(sf)
library(sp)
library(tmap)
# library(ggmap)
# library(ggplot2)

# Function to distribute the weekly total across the days with noise
distribute_with_noise <- function(weekly_total, base_weights, noise_SD=5) {
  
  # Normalize the base weights so that they sum to 1
  norm_weights <- base_weights / sum(base_weights)
  
  # Compute initial daily estimates (based on the normalized weights)
  daily_estimates <- norm_weights * weekly_total
  
  # Add noise (uniform noise between -5 and 5)
  noise <- rnorm(length(daily_estimates), mean = 0, sd = noise_SD)
  
  # Adjust the daily estimates by adding noise
  noisy_estimates <- daily_estimates + noise
  
  # Ensure the noisy estimates are non-negative and the total is the same
  noisy_estimates <- pmax(0, noisy_estimates)
  
  # Adjust the sum to match the weekly total (by scaling the noisy estimates)
  scale_factor <- weekly_total / sum(noisy_estimates)
  noisy_estimates <- round(noisy_estimates * scale_factor)
  
  return(noisy_estimates)
}

# Function to systematically simulate daily data, across data sets
get_daily <- function(df_weekly, date, out_list, weights, noise_SD=5) {
  
  # Update date format & select relevant outcomes
  df_weekly <- df_weekly %>%
    mutate(date = as.Date(!!sym(date),format = "%Y%m%d")) %>% 
    select(date, all_of(out_list))    
  
  # Create empty data frame with daily dates and week days 
  daily_dates <- seq(as.Date(min(df_weekly$date)-6), max(df_weekly$date), by = "day")
  df_daily <- data.frame(date = daily_dates)
  df_daily$weekday <- weekdays(df_daily$date) 
  
  # Define weekday weights (Monday to Sunday, sum should be 1) - baesd on article
  # https://www.independent.co.uk/news/health/nhs-waits-accident-emergency-data-b2171791.html
  if (is.null(weights)){
    weights <- c(Monday = 0.16, Tuesday = 0.155, Wednesday = 0.13, Thursday = 0.14, Friday = 0.13, Saturday = 0.14, Sunday = 0.145) 
    } else{
    weights <- weights
  }
  
  # Iteratively apply noisy weights function to get daily estimates

    for (var in out_list){
    
    w_out <- data.frame()
    
    for (week in 1:nrow(df_weekly)){
      
      w_variable <- as.data.frame(distribute_with_noise(df_weekly[[var]][week], weights, noise_SD))
      colnames(w_variable) <- paste0('w_', var)
      w_out <- rbind(w_out, w_variable)
      
    }
    
    # Create cumulative variable
    w_out <- cbind(w_out, cumsum(w_out))
    colnames(w_out)[ncol(w_out)] <- paste0("c_", var)
    
    # Append to df_daily
    df_daily <- cbind(df_daily, w_out)
    
  }
  
  return(df_daily)
}

# Set seed for reproducibility
set.seed(123)

# 1. ENGLAND: Simulate Daily Mortality - National Level ####

# Get file list
file_list <- list.files(path='201_heal/ONS/', pattern = '\\.xlsx$', full.names = TRUE)

# Read all files and combine them into a single dataframe without column names
df_weekly <- do.call(rbind, lapply(file_list, function(file) {
  read_excel(file)
}))

# Data exploration & filter (if necessary)
colnames(df_weekly)

# Set parameters
out_list <- c('total_deaths', 'resp_deaths')
date <- 'week_end'

# Get daily estimates
results <- get_daily(df_weekly=df_weekly, date=date, out_list=out_list, 
                     weights=NULL, noise_SD=10)

# Visualize and export
head(results, 10)
plot(df_weekly[[date]], df_weekly[[out_list[1]]])
plot(results$date, results[,3])

# Export the data frame to a CSV file
write.csv(results, file = "eng_adm0_ONS_2019-2010_daily-national-mortality.csv", row.names = FALSE)


# 2. SCOTLAND: Simulate Daily A&E - National Level ####

# Get file list
file_list <- list.files(path='201_heal/NHS Scotland/', pattern = '\\.csv$', full.names = TRUE)

# Read all files and combine them into a single data frame without column names
df_weekly <- do.call(rbind, lapply(file_list, function(file) {
  read.csv(file)  # Read without column names
}))

# Data exploration & cleaning (if necessary)
head(df_weekly)
colnames(df_weekly)
df_weekly <- df_weekly %>% 
  filter(AttendanceCategory == "All") %>%
  group_by(WeekEndingDate) %>%
  summarize(across(where(is.numeric), sum)) # NB: duplicates in "G405H" - aggregate via SUM
df_weekly$WeekEndingDate <- as.character(df_weekly$WeekEndingDate)

# Set parameters
out_list <- c('NumberOfAttendancesEpisode')
date <- 'WeekEndingDate'

# Get daily estimates
results <- get_daily(df_weekly=df_weekly, date=date, out_list=out_list, 
                     weights=NULL, noise_SD=10)


# Visualize and export
head(results, 10)
plot(df_weekly[[date]], df_weekly[[out_list[1]]])
plot(results$date, results[,3])

# Export the data frame to a CSV file
write.csv(results, file = "agy_adm0_NHS-Scotland_2015-2025_daily-national-ae.csv", row.names = FALSE)



# 3. SCOTLAND: Simulate Daily A&E - Hospital Level ####

# Get file list
file_list <- list.files(path='201_heal/NHS Scotland/', pattern = '\\.csv$', full.names = TRUE)

# Read all files and combine them into a single data frame without column names
df_weekly <- do.call(rbind, lapply(file_list, function(file) {
  read.csv(file)  # Read without column names
}))

# Data exploration & cleaning (if necessary)
head(df_weekly)
colnames(df_weekly)
df_weekly <- df_weekly %>% 
  filter(AttendanceCategory == "All") %>%
  group_by(WeekEndingDate, TreatmentLocation) %>%
  summarize(across(where(is.numeric), sum)) # NB: duplicates in "G405H" - aggregate via SUM
df_weekly$WeekEndingDate <- as.character(df_weekly$WeekEndingDate)

# Set parameters
out_list <- c('NumberOfAttendancesEpisode')
date <- 'WeekEndingDate'
hosp_list <- unique(df_weekly$TreatmentLocation)

# Iterate over each hospital (TreatmentLocation)
results <- data.frame()
for (hospital in hosp_list) {
  
  # Filter data for the current hospital
  df_hospital <- df_weekly %>% filter(TreatmentLocation == hospital)
  
  # Apply function to the filtered data & keep hosptial code
  df_hospital <- get_daily(df_weekly=df_hospital, date=date, out_list=out_list, 
                       weights=NULL, noise_SD=10)
  df_hospital <- cbind(sitecode = hospital, df_hospital)
  
  # Store the result in the list with the hospital_id as the name
  results <- rbind(results, df_hospital)
}

# Visualize and export
head(results, 10)
df_weekly <- df_weekly %>% filter(TreatmentLocation=='A111H') # relevant sites
df_weekly$date <- as.Date(df_weekly$WeekEndingDate,format = "%Y%m%d")
results <- results %>% filter(sitecode=='A111H') # relevant sites
par(mfrow=c(1,2))  
plot(df_weekly$date, df_weekly$NumberOfAttendancesEpisode, main='Crosshouse Hospital - Weekly')
plot(results$date, results$w_NumberOfAttendancesEpisode, main='Crosshouse Hospital - Daily')

# Export the data frame to a CSV file
write.csv(results, file = "agy_adm4_NHS-Scotland_2015-2025_daily-hosptial-ae.csv", row.names = FALSE)



# 4. SCOTLAND: Combine hospital & weather time series with locations ####

# Open shape files and csv
pnt <- st_read('201_heal/NHS Scotland/agy_adm4_NHS-Scotland_hosptials.shp')
poly <- st_read('201_heal/NHS Scotland/agy_adm2_NHS-Scotland_HealthBoards_2019.shp')
weather <- read.csv('202_temp/agy_adm3_NHS-Scotland_2015-2020_daily-temperature.csv')
health <- read.csv('201_heal/agy_adm4_NHS-Scotland_2015-2025_daily-hosptial-ae.csv')

# CLEAN DATA FOR MERGES & VISUALISING
### Transform from ONS British National Grid to WGS84 (for ggmap use)
poly <- st_transform(poly, crs = 4326 ) # transform
pnt <- st_transform(pnt, crs = 4326 ) 
pnt$long <- st_coordinates(pnt)[,1] # get lat/long
pnt$lat <- st_coordinates(pnt)[,2]
### Visualize the data & choose health catchment boundary (HealthBoards)
pnt <- pnt %>% filter(sitecode %in% unique(health$sitecode)) # relevant sites
# pnt <- distinct(pnt, hbname, .keep_all= TRUE) # get specific 
# st_bbox(pnt)
tmap_mode("view") # set to interactive
tm_shape(poly) + 
  tm_polygons(col = "blue", alpha = 0.3, border.col = "darkblue") + 
  tm_shape(pnt) +
  tm_symbols(col = "red", size = 0.5) +
  tm_layout(title = "Interactive Map of Hospitals & Health Authority Boundaries")

# SPATIAL JOIN (ensure HBCodes are correct)
pnt <- st_join(pnt, poly, join = st_within)
pnt <- as(pnt, 'Spatial') # make spatial object
poly <- as(poly, 'Spatial') 

# MERGE WEATHER & POLYGONS
l_weather <- poly@data %>%
  left_join(weather, by=c('HBCode', 'HBName')) %>%
  rename(hbcode=HBCode, hbname=HBName)

# MERGE HEATLH & POINTS
l_health <- health %>%
  left_join(pnt@data, by = 'sitecode') %>%
  select(-hbcode, -hbname) %>%
  rename(hbcode=HBCode, hbname=HBName, 
         w_Num_AE_Attend=w_NumberOfAttendancesEpisode, 
         c_Num_AE_Attend=c_NumberOfAttendancesEpisode)

# MERGE WEATHER & HEALTH
final <- l_health %>%
  left_join(l_weather, by=c('hbname', 'hbcode', 'date')) %>%
  select(hbcode, hbname, sitecode, postcode, long, lat, sitename, type, date, weekday,
         temp, w_Num_AE_Attend, c_Num_AE_Attend)


head(final, 10)

