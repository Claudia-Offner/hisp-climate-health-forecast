### ------------------------------------------------------------------------ ### 
### Hospital cleaning - HISP Data                                            ###
### ------------------------------------------------------------------------ ### 

# Turn off scientific notation
options(scipen=999)

# Set data path
user <- 'claer14'
path1 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
path2 <- paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/1_Data/NHS_Cambrigeshire/')
path3 <- paste0('C:/Users/', user, '/Downloads/data_public')
# setwd(path2)

####  DEPENDENCIES  ####

# Packages
library(dplyr); library(purrr); library(tidyverse); library(lubridate); 
library(tidyselect); library(sf); library(tmap)


####  MAIN CODE ####

# Load data
df <- read.csv(paste0(path1, '201_heal/gbr_adm4_NHS-CP_2023-2025_ae_raw.csv'))
nhs_icb_sf <- st_read(paste0(path3, '/spatial/NHS_2023_ICB_EN/NHS_2023_ICB_EN.shp'))
icb_sf <- nhs_icb_sf %>% filter(ICB23NM == 'NHS Cambridgeshire and Peterborough Integrated Care Board')
nhs_regions_sf <- st_read(paste0(path3, '/spatial/NHS_2024_Regions_EN/NHS_2024_Regions_EN.shp'))
region_sf <- nhs_regions_sf %>% filter(NHSER24NM == 'East of England')

# Check spatial data
provider_unique <- df %>% distinct(provider_name, .keep_all = TRUE)
provider_sf <- st_as_sf(provider_unique, coords = c("provider_lon", "provider_lat"), crs = 4326)
practice_unique <- df %>% distinct(practice_name, .keep_all = TRUE)
practice_sf <- st_as_sf(practice_unique, coords = c("practice_lon", "practice_lat"), crs = 4326)
tmap_mode("view")
tm_shape(region_sf) +
  tm_polygons(col = "darkseagreen2", alpha = 1) +
  tm_text("NHSER24NM", size = 1.5, col = "black", xmod = 1) + 
  tm_shape(icb_sf) +
  tm_polygons(col = "lightblue", alpha = 1) +
  tm_shape(practice_sf) +
  tm_dots(col = 'darkblue', border.col=NA, size = 0.5) +
  tm_shape(provider_sf) +
  tm_dots(col = 'red', border.col=NA, size = 0.5) +
  tm_layout(main.title = "ICB, LSOA, Practices, Providers") +
  tm_add_legend(
    type = "symbol",
    labels = c("ICB", "Practices", "Providers"),
    col = c("lightblue", "darkblue", "red"),
    size = 0.5
  )

# Address missing values and set time series range ####

# Count missing values per column
(missing_counts <- colSums(is.na(df)))
# Check rows with any missing values
(df_missing_any <- df[rowSums(is.na(df)) > 0, ])

# Identify and set data types ####

# Get wide format daily counts ####

# Total daily count of patients
total <- df %>%
  group_by(location, start_date) %>%
  summarise(total_patients = n())

# Daily count by ward
ward_wide <- df %>%
  group_by(location, date, ward_name) %>%
  summarise(patients = n(), .groups = 'drop') %>%
  pivot_wider(names_from = ward_name, values_from = patients, values_fill = 0) %>%
  select(location, date, !!!syms(ward_levels))

# Aggregate by age
age_wide <- df %>%
  group_by(location, date, age_group) %>%
  summarise(patients = n()) %>%
  pivot_wider(names_from = age_group, values_from = patients,
              values_fill = list(patients = 0))

# Aggregate by gender
gender_wide <- df %>%
  group_by(location, date, gender) %>%
  summarise(patients = n()) %>%
  pivot_wider(names_from = gender, values_from = patients,
              values_fill = list(patients = 0))

# Aggregate by ward & age
age_ward_wide <- df %>%
  group_by(location, date, ward_name, age_group) %>%
  summarise(patients = n()) %>%
  pivot_wider(names_from = c(ward_name, age_group),
              values_from = patients,
              values_fill = list(patients = 0))

# Aggregate by ward & gender
gender_ward_wide <- df %>%
  group_by(location, date, ward_name, gender) %>%
  summarise(patients = n()) %>%
  pivot_wider(names_from = c(ward_name, gender),
              values_from = patients,
              values_fill = list(patients = 0))


unique(df$provider_lat)
unique(df$provider_lon)
