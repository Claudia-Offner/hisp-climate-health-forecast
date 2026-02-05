### ------------------------------------------------------------------------ ### 
### Hospital cleaning - NHS  Data
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list = ls()); graphics.off(); gc()
options(scipen=999)

# Set data path
user <- 'claer14'
path1 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
path2 <- paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/1_Data/NHS_Cambrigeshire/')

####  DEPENDENCIES  ####
# source('1_Data/HISP_Himachal_Pradesh/HISP_Cleaning - Dependencies.R')
heal_export <- 'eng_adm4_NHS_2024_ae_visits_ALL_clean'

# Packages
library(dplyr); library(purrr); library(tidyverse); library(lubridate); 
library(tidyselect); library(sf); library(tmap); library(janitor)

# Function to get columns in wide format for time series
send_wide <- function(df, col) {
  col_sym <- sym(col)
  col_levels <- levels(df[[col]])  # get factor levels
  
  df %>%
    group_by(provider_code, end_date, !!col_sym) %>%
    summarise(patients = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = !!col_sym,
      values_from = patients,
      values_fill = 0
    ) %>%
    select(provider_code, end_date, all_of(col_levels))  # enforce column order
}


####  MAIN CODE ####

# Load data
df <- read.csv(paste0(path1, '201_heal/gbr_adm4_NHS-CP_2023-2025_ae_raw.csv'), row.names = 1)
icb_sf <- st_read(paste0(path1, '/200_admn/eng_admn_ad3_py_s1_NHS_ICB_2023.shp'))
regions_sf <- st_read(paste0(path1, '/200_admn/eng_admn_ad1_py_s1_NHS_Regions_2024.shp'))
lsoa_sf <- st_read(paste0(path1, '/200_admn/eng_admn_ad4_py_s1_ONS_LSOA_2011.shp'))

# Filter data
icb_sf <- icb_sf %>% filter(ICB23NM == 'NHS Cambridgeshire and Peterborough Integrated Care Board')
region_sf <- regions_sf %>% filter(NHSER24NM == 'East of England')

# # 1. Process raw data ####
# 
# path3 <- paste0('C:/Users/', user, '/Downloads/data_icb')
# path4 <- paste0('C:/Users/', user, '/Downloads/data_public')
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
#     mutate(uid = paste0("ae", sample.int(1e9, n()))) %>%
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
  # tm_shape(lsoa_sf) +
  # tm_polygons(col = NA, border.col = "black", lwd = 0.8) +
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

# Get patient LSOA centroids
lsoa_df <- lsoa_sf %>% st_drop_geometry() %>% 
  select(patient_lsoa_code=LSOA11CD, patient_lon=LONG_, patient_lat=LAT)
df <- df %>% left_join(lsoa_df, by='patient_lsoa_code')

# 3. Set data types and categorical variables ####

# Count missing values per column
(missing_counts <- colSums(is.na(df)))
# Check rows with any missing values
(df_missing_any <- df[rowSums(is.na(df)) > 0, ])

# Re-code providers
df$provider_name <- str_to_title(df$provider_name)

# Re-code age groups
age_levels <- c('age_under_5', 'age_5_14', 'age_15_24', 'age_25_64',  'age_over_65','age_NA')
df <- df %>%
  mutate(
    patient_age_group = case_when(
      patient_age < 0 | patient_age > 110 ~ "age_NA",
      patient_age < 5  ~ "age_under_5",
      patient_age < 15 ~ "age_5_14",
      patient_age < 25 ~ "age_15_24",
      patient_age < 65 ~ "age_25_64",
      TRUE             ~ "age_over_65"
    ),
    patient_age_group = factor(patient_age_group, levels = age_levels)
  )
sort(table(df$patient_age_group))

# Re-code deciles
unique(df$patient_imd_decile)
df <- df %>%
  mutate(patient_imd_quintile= case_when(
    patient_imd_decile <= 2 ~ 'imd_q1',
    patient_imd_decile > 2 & patient_imd_decile <= 4 ~ 'imd_q2',
    patient_imd_decile > 4 & patient_imd_decile <= 6 ~ 'imd_q3',
    patient_imd_decile > 6 & patient_imd_decile <= 8 ~ 'imd_q4',
    patient_imd_decile > 8 & patient_imd_decile <= 10 ~ 'imd_q5',
    TRUE ~ 'imd_NA'
  ), patient_imd_quintile = factor(patient_imd_quintile, levels = c('imd_q1','imd_q2', 'imd_q3', 'imd_q4', 'imd_q5', 'imd_NA'))
)
sort(table(df$patient_imd_quintile))

# Re-code treatment codes
unique(df$treat_desc)
df <- df %>%
  mutate(
    treat_desc = recode(
      treat_desc,
      "Category 1 - Minor Treatment" = "treat_1",
      "Category 2" = "treat_2",
      "Category 3" = "treat_3",
      "Category 4" = "treat_4",
      "Category 5 - Major Treatment" = "treat_5",
      "Other" = 'treat_NA',
      .default = treat_desc
    ),
    treat_desc = factor(treat_desc, levels=c('treat_1', 'treat_2', 'treat_3', 'treat_4', 'treat_5', 'treat_NA'))
  )
sort(table(df$treat_desc))

# Re-code healthcare resource Groups (HRGs)
unique(df$hrg_desc)
df <- df %>%
  mutate(
    hrg_desc = recode(
      hrg_desc,
      "Emergency Medicine, No Investigation with No Significant Treatment" = "Investigation 0 with Treatment 0",
      "Emergency Medicine, Category 1 Investigation with Category 1-2 Treatment" = "Investigation 1 with Treatment 1-2",
      "Emergency Medicine, Category 1 Investigation with Category 3-4 Treatment" = "Investigation 1 with Treatment 3-4",
      "Emergency Medicine, Category 2 Investigation with Category 1 Treatment" = "Investigation 2 with Treatment 1-2",
      "Emergency Medicine, Category 2 Investigation with Category 2 Treatment" = "Investigation 2 with Treatment 1-2",
      "Emergency Medicine, Category 2 Investigation with Category 3 Treatment" = "Investigation 2 with Treatment 3-4",
      "Emergency Medicine, Category 2 Investigation with Category 4 Treatment" = "Investigation 2 with Treatment 3-4",
      "Emergency Medicine, Category 3 Investigation with Category 1-3 Treatment" = "Investigation 3 with Treatment 1-3",
      "Emergency Medicine, Category 3 Investigation with Category 4 Treatment" = "Investigation 3 with Treatment 4",
      "Emergency Medicine, Any Investigation with Category 5 Treatment" = "Any Investigation with Treatment 5",
      "Emergency Medicine, Patient Dead On Arrival" = NA_character_,
      "Emergency Medicine, Dental Care" = NA_character_,
      .default = hrg_desc
    ),
    hrg_desc = factor(hrg_desc, levels = c('Investigation 0 with Treatment 0', 
                                           'Investigation 1 with Treatment 1-2',
                                           'Investigation 1 with Treatment 3-4',
                                           'Investigation 2 with Treatment 1-2',
                                           'Investigation 2 with Treatment 3-4',
                                           'Investigation 3 with Treatment 1-3',
                                           'Investigation 3 with Treatment 4',
                                           'Any Investigation with Treatment 5')),
    hrg_num = factor(paste0("hrg_", as.numeric(hrg_desc)),
                     levels = paste0("hrg_", seq_along(levels(hrg_desc))))
  )
sort(table(df$hrg_desc))
sort(table(df$hrg_num))

# 4. Set binary thdf_widehold where relevant ####

# Re-code age groups
df <- df %>%
  mutate(
    patient_age_cat = case_when(
      patient_age < 0 | patient_age > 110 ~ "age_NA",
      patient_age < 18  ~ "age<18",
      patient_age < 64 ~ "age18-64",
      TRUE             ~ "age>64"
    ),
    patient_age_cat = factor(patient_age_cat, levels = c('age<18', 'age18-64', 'age>64', 'age_NA'))
  )
table(df$patient_age_cat)

# Re-code deciles
unique(df$patient_imd_decile)
df <- df %>%
  mutate(patient_imd_cat= case_when(
    patient_imd_decile <= 2 ~ 'imd_low_q',
    patient_imd_decile > 2 & patient_imd_decile <= 8 ~ 'imd_mid_q',
    patient_imd_decile > 8 & patient_imd_decile <= 10 ~ 'imd_hig_q',
    TRUE ~ 'imd_NA'
  ), patient_imd_cat = factor(patient_imd_cat, levels = c('imd_low_q','imd_mid_q', 'imd_hig_q', 'imd_NA'))
  )
sort(table(df$patient_imd_cat))

# Re-code land groups
df <- df %>%
  mutate(
    patient_area_class = case_when(
      patient_area_class == 0 ~ "land_rural",
      patient_area_class ==1 ~ "land_urban",
      TRUE             ~ "land_NA"
    ),
    patient_area_class = factor(patient_area_class, levels = c('land_rural', 'land_urban', 'land_NA'))
  )
table(df$patient_area_class)


# Re-code treatment codes
unique(df$treat_desc)
df <- df %>%
  mutate(
    treat_cat = case_when(
      treat_desc %in% c("treat_1", "treat_2", "treat_3") ~ "treat_mild",
      treat_desc %in% c("treat_4", "treat_5") ~ "treat_severe",
      TRUE ~ "treat_NA"
    ),
    treat_cat = factor(treat_cat, levels = c("treat_mild", "treat_severe", "treat_NA"))
  )
sort(table(df$treat_cat))


# 5. Get wide format daily counts ####

# Total daily count of patients
wide_total <- df %>% group_by(provider_code, end_date) %>% summarise(total_patients = n())

# Wide format categorical data
wide_imd <- send_wide(df, 'patient_imd_quintile')
wide_age <- send_wide(df, 'patient_age_group')
wide_land <- send_wide(df, 'patient_area_class')
wide_treat <- send_wide(df, 'treat_desc')
wide_hrg <- send_wide(df, 'hrg_num')

# Wide format thdf_widehold data
wide_age_cat <- send_wide(df, 'patient_age_cat')
wide_imd_cat <- send_wide(df, 'patient_imd_cat')
wide_treat_cat <- send_wide(df, 'treat_cat')

# Merge datasets
lat_lon <- df %>% select(provider_code, provider_name, provider_lon, provider_lat) %>% distinct()
df_wide <- wide_total %>% 
  left_join(wide_age, by=c('provider_code', 'end_date')) %>%
  left_join(wide_age_cat, by=c('provider_code', 'end_date')) %>%
  left_join(wide_imd, by=c('provider_code', 'end_date')) %>%
  left_join(wide_imd_cat, by=c('provider_code', 'end_date')) %>%
  left_join(wide_land, by=c('provider_code', 'end_date')) %>%
  left_join(wide_hrg, by=c('provider_code', 'end_date')) %>%
  left_join(wide_treat, by=c('provider_code', 'end_date')) %>%
  left_join(wide_treat_cat, by=c('provider_code', 'end_date')) %>%
  mutate(
    end_date = as.Date(end_date),
    year = year(end_date),
    month = month(end_date, label = TRUE, abbr = TRUE),
    week = week(end_date),
    doy = yday(end_date),
    dow = wday(end_date, label = TRUE, abbr = TRUE)) %>%
  select(-contains("_NA")) %>%
  left_join(lat_lon, by=c('provider_code')) %>% 
  select(date=end_date, year, month, week, doy, dow, provider_code, provider_name, provider_lon, provider_lat, everything()) %>%
  arrange(provider_code, date)

# 6. Check temporal information ####
weekday_levels  <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')

# 479 days missing between '2023-04-01' and '2025-03-31'
(missing <- df_wide %>%
   group_by(provider_name) %>%
   summarise(
     missing_days = list({
       date_range <- seq(min(date), max(date), by = 'day')
       as.Date(setdiff(date_range, date))
     }),
     .groups = 'drop'
   ) %>%
   mutate(
     missing_weekday = lapply(missing_days, function(d) {
       factor(format(d, "%a"), levels = weekday_levels)
     }),
     n_missing = sapply(missing_days, length)
   ))
min(df_wide$date)
max(df_wide$date)
sum(missing$n_missing)
missing_days <- missing$missing_days
missing_weekday <- missing$missing_weekday

missing_tables <- lapply(seq_along(missing_days), function(i) {
  table(missing_weekday[[i]], format(missing_days[[i]], '%Y'))
})
names(missing_tables) <- missing$provider_name  # label by location
missing_tables

# Create full time series with NA
df_wide <- df_wide %>%
  group_by(provider_code, provider_name, provider_lon, provider_lat) %>%
  complete(date = seq(min(date), max(date), by = "day")) %>%
  arrange(provider_code, date)

#### EXPORT ####

save(df, file=paste0(path1, '/201_heal/', heal_export, '.RData'))
save(df_wide, file=paste0(path1, '/201_heal/', heal_export, '_wide.RData'))
