### ------------------------------------------------------------------------ ### 
### Hospital cleaning - HISP Data                                            ###
### ------------------------------------------------------------------------ ### 

# Turn off scientific notation
options(scipen=999)

# Set data path
user <- 'claer14'
path1 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
path2 <- paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/1_Data/NHS_Cambrigeshire/')
path3 <- paste0('C:/Users/', user, '/Downloads/data_icb')
path4 <- paste0('C:/Users/', user, '/Downloads/data_public')
# setwd(path2)

####  DEPENDENCIES  ####

# Packages
library(dplyr); library(purrr); library(tidyverse); library(lubridate); 
library(tidyselect); library(sf); library(tmap); library(janitor)


####  MAIN CODE ####

# Load data
df <- read.csv(paste0(path1, '201_heal/gbr_adm4_NHS-CP_2023-2025_ae_raw.csv'))
nhs_icb_sf <- st_read(paste0(path4, '/spatial/NHS_2023_ICB_EN/NHS_2023_ICB_EN.shp'))
icb_sf <- nhs_icb_sf %>% filter(ICB23NM == 'NHS Cambridgeshire and Peterborough Integrated Care Board')
nhs_regions_sf <- st_read(paste0(path4, '/spatial/NHS_2024_Regions_EN/NHS_2024_Regions_EN.shp'))
region_sf <- nhs_regions_sf %>% filter(NHSER24NM == 'East of England')

# # 1. Process raw data ####
# 
# # A&E (unplanned)
# ae_23 <- read.delim(paste0(path3, '/A&E_2324 070825.txt'), sep = '|', header = TRUE, stringsAsFactors = FALSE)|> clean_names()
# ae_24 <- read.delim(paste0(path3, '/A&E_2425 070825.txt'), sep = '|', header = TRUE, stringsAsFactors = FALSE)|> clean_names()
# identical(names(ae_23), names(ae_24))
# ae <- rbind(ae_23, ae_24)
# ae$type <- 'ae'
# 
# # Get PROVIDER codes (from all datasets)
# provider_lookup_file <- paste0(path3, '/provider_lookup.csv')
# provider_codes <- unique(ae$site_code_of_treatment)
# provider_lookup <- read.csv(provider_lookup_file)
# na_rows <- provider_lookup[is.na(provider_lookup$lat), ] # Check NA
# ### Create spatial file
# all_providers_sf <- st_as_sf(provider_lookup %>% filter(!is.na(lat)), coords = c("lon", "lat"), crs = 4326, remove = FALSE)
# all_providers_sf <- st_transform(all_providers_sf, crs = 27700) # BNG
# st_write(all_providers_sf, paste0(path4, '/spatial/NHS_2025_Provider_UK/NHS_2025_Provider_UK.shp'), delete_layer = TRUE)
# 
# # Get PRACTICE codes (from all datasets)
# practice_lookup_file <- paste0(path3, '/practice_lookup.csv')
# provider_codes <- unique(ae$dmic_practice_code)
# practice_lookup <- read.csv(practice_lookup_file)
# na_rows <- practice_lookup[is.na(practice_lookup$lat), ] # Check NA
# ### Create spatial file
# all_practices_sf <- st_as_sf(practice_lookup %>% filter(!is.na(lat)), coords = c("lon", "lat"), crs = 4326, remove = FALSE)
# all_practices_sf <- st_transform(all_practices_sf, crs = 27700) # BNG
# st_write(all_practices_sf, paste0(path4, '/spatial/NHS_2025_Practice_UK/NHS_2025_Practice_UK.shp'), delete_layer = TRUE)
# 
# # Get ICB and ONS boundaries
# nhs_icb_sf <- st_read(paste0(path4, '/spatial/NHS_2023_ICB_EN/NHS_2023_ICB_EN.shp'))
# all_lsoa_sf <- st_read(paste0(path4, '/spatial/ONS_2011_LSOA_EN/ONS_2011_LSOA_EN.shp'))
# ru_lsoa_sf <- st_read(paste0(path4, '/spatial/ONS_2011_LSOA_EN/ONS_2011_LSOA_EN_RU.shp'))
# 
# # Add urban/rural information to LSOA
# ru_lsoa <- ru_lsoa_sf %>% st_drop_geometry() %>% select(LSOA11CD, RUC11)
# all_lsoa_sf <- all_lsoa_sf %>% left_join(ru_lsoa, by = "LSOA11CD")
# all_lsoa_sf$area_class <- ifelse(grepl("Urban", all_lsoa_sf$RUC11, ignore.case = TRUE), 1,
#                                  ifelse(grepl("Rural", all_lsoa_sf$RUC11, ignore.case = TRUE), 0, NA))
# # Filter spatial data to ICB
# icb_sf <- nhs_icb_sf %>% filter(ICB23NM == 'NHS Cambridgeshire and Peterborough Integrated Care Board')
# 
# # Add LSOA data to practices and format providers/practices for joining
# points_with_lsoa <- st_join(all_practices_sf, all_lsoa_sf, join = st_within)
# all_practices_df <- st_set_geometry(points_with_lsoa, NULL) %>%
#   select(code, name, postcode, lon, lat, LSOA11CD, LSOA11NM, RUC11, area_class) %>%
#   rename('lsoa_code'='LSOA11CD','lsoa_name'='LSOA11NM', 'area_type'='RUC11') %>%
#   rename_with(~ paste0("practice_", .x))
# all_provider_df <- st_set_geometry(all_providers_sf, NULL) %>%
#   rename_with(~ paste0("provider_", .x))
# 
# # FORMAT FOR PROCESSING
# ae <- ae %>%
#   # Add missing cols
#   mutate(uid = paste0("ae", row_number()))  %>%
#   # Add practice/provider/rural/urban info
#   left_join(all_practices_df, by=c('dmic_practice_code'='practice_code')) %>%
#   left_join(all_provider_df, by=c('site_code_of_treatment'='provider_code')) %>%
#   left_join(ru_lsoa, by = c('lsoa'="LSOA11CD")) %>%
#   # Identify patient urban/rural area
#   mutate(patient_area_class = ifelse(grepl("Urban", RUC11, ignore.case = TRUE), 1,
#                                      ifelse(grepl("Rural", RUC11, ignore.case = TRUE), 0, NA))) %>%
#   filter(!is.na(lsoa) & lsoa != "" & lsoa != "NULL") %>%
#   # Select and name relevant cols
#   select(uid, type, year=fin_year, start_date=a_and_e_departure_date, end_date=a_and_e_departure_date,
#          patient_lsoa_code=lsoa, patient_area_type=RUC11, patient_area_class,
#          practice_code=dmic_practice_code, practice_name, practice_postcode, practice_lon, practice_lat,practice_lsoa_code, practice_lsoa_name,
#          provider_code=site_code_of_treatment, provider_name, provider_postcode, provider_lon, provider_lat,
#          patient_age=dmic_age, patient_imd_decile=decile, patient_imd_value=imd_2019,
#          hrg_code=hrg, hrg_desc, treat_desc=dmic_ae_banding,
#          practice_area_type, practice_area_class)
# 
# write.csv(ae, paste0(path1, '201_heal/gbr_adm4_NHS-CP_2023-2025_ae_raw.csv'))





# 2. Check spatial data ####

provider_unique <- df %>% distinct(provider_name, .keep_all = TRUE)
provider_sf <- st_as_sf(provider_unique, coords = c("provider_lon", "provider_lat"), crs = 4326)
practice_unique <- df %>% distinct(practice_name, .keep_all = TRUE)  %>% filter(!is.na(practice_lon), !is.na(practice_lat))
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

# How many patients per provider? (remove STAMFORD & RUTLAND HOSPITAL / RGN49)
table(df$provider_name)

# Address missing values and set time series range ####

# Count missing values per column
(missing_counts <- colSums(is.na(df)))
# Check rows with any missing values
(df_missing_any <- df[rowSums(is.na(df)) > 0, ])

age_levels <- c('under_5', 'age5_14', 'age15_24', 'age25_64',  'over_65','Unknown')

# Classify age groups
df <- df %>%
  mutate(age_group = case_when(
    patient_age < 0 ~ 'Unknown',
    patient_age >= 0 & patient_age < 5 ~ 'under_5',
    patient_age >=5 & patient_age < 16 ~ 'age5_14',
    patient_age >= 15 & patient_age < 25 ~ 'age15_24',
    patient_age >= 25 & patient_age < 65 ~ 'age25_64',
    patient_age >= 65 & patient_age <= 110 ~ 'over_65',
    patient_age > 110 ~ 'Unknown',
    TRUE ~ 'Unknown'
  )) %>%
  mutate(age_group = factor(age_group, levels = age_levels))
sort(unique(df$age))
table(df$age_group)
df %>% filter(age_group =='Unknown')


# Identify and set data types ####

names(df)
unique(df$hrg_desc)

df_heal <- df %>% select(date=start_date, 
              location=provider_name, 
              condition=treat_desc,
              age_group, patient_imd_decile)

# Re-code conditions
unique(df_heal$condition)
condition_patterns <- c(
  'Category 1 - Minor Treatment' = 'Category 1',
  'Category 5 - Major Treatment' = 'Category 5'
)

df_heal$condition_group <- df_heal$condition
for (pattern in names(condition_patterns)) {
  df_heal$condition_group[grepl(pattern, df_heal$condition, ignore.case = TRUE)] <- condition_patterns[pattern]
}

sort(unique(df_heal$condition_group))
sort(table(df_heal$condition_group))
df_heal$condition <- df_heal$condition_group
df_heal$condition_group <- NULL


# Recode deciles
unique(df_heal$patient_imd_decile)
df_heal <- df_heal %>%
  mutate(patient_quintile_group = case_when(
    patient_imd_decile <= 2 ~ 'Q1',
    patient_imd_decile > 2 & patient_imd_decile <= 4 ~ 'Q2',
    patient_imd_decile > 4 & patient_imd_decile <= 6 ~ 'Q3',
    patient_imd_decile > 6 & patient_imd_decile <= 8 ~ 'Q4',
    patient_imd_decile > 8 & patient_imd_decile <= 10 ~ 'Q5',
    TRUE ~ 'Unknown'
  )) %>%
  mutate(patient_quintile_group = factor(patient_quintile_group, levels = c('Q1','Q2', 'Q3', 'Q4', 'Q5', 'Unknown')))


# Get wide format daily counts ####

# Total daily count of patients
total <- df1 %>%
  group_by(location, date) %>%
  summarise(total_patients = n())

# Daily count by ward
ward_wide <- df1 %>%
  group_by(location, date, condition) %>%
  summarise(patients = n(), .groups = 'drop') %>%
  pivot_wider(names_from = condition, values_from = patients, values_fill = 0) %>%
  select(location, date, everything())

# Aggregate by age
age_wide <- df1 %>%
  group_by(location, date, age_group) %>%
  summarise(patients = n()) %>%
  pivot_wider(names_from = age_group, values_from = patients,
              values_fill = list(patients = 0))

# Aggregate by ward & age
age_ward_wide <- df1 %>%
  group_by(location, date, condition, age_group) %>%
  summarise(patients = n()) %>%
  pivot_wider(names_from = c(condition, age_group),
              values_from = patients,
              values_fill = list(patients = 0))


# Aggregate total patients per date
daily_totals <- aggregate(total_patients ~ date, data = total, FUN = sum)
daily_totals$date <- as.Date(daily_totals$date)
library(dplyr)
library(zoo)

daily_totals <- daily_totals %>%
  arrange(date) %>%  # make sure it's ordered
  mutate(ma7 = rollmean(total_patients, k = 7, fill = NA, align = "right"))


plot(daily_totals$date, daily_totals$total_patients, type = "l", col = "grey",
     xlab = "Date", ylab = "Total Patients")
lines(daily_totals$date, daily_totals$ma7, lwd = 2)  # default col is fine


