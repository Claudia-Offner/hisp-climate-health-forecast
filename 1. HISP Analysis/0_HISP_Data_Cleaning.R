### ------------------ ### 
### HISP Data Cleaning ###
### ------------------ ### 

# NB COLOURS: https://r-charts.com/colors/

# Set data path
user <- 'claer14'
path <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')

####  PACKAGES  ####
library(tidyverse)
library(dplyr)
library(lubridate)

####  FUNCTIONS  ####

# Unpack the tt separated HISP data
unpack_HISP <- function(file_path, hospital, download_path){
  
  # Load data
  data <- read.csv(file_path)
  names(data) <- 'raw_column'

  # Splitting the columns & organising
  df <- data %>%
    separate(raw_column, into = c('ward_name', 'patient_id', 'age', 'gender', 'category', 'date', 'visit_type'),
             sep = 'tt', fill = 'right') %>%
    mutate(location = hospital)

  # Reorganize col order & set type
  df <- df[, c('location', 'date', 'patient_id', 'ward_name', 'age', 'gender', 'category', 'visit_type')]

  # Export
  write.csv(df, file=download_path, row.names = FALSE)
  
}

####  DEPENDENCIES  ####

weekday_levels <- c('Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday')
age_levels <- c('TOTAL', '<1y', '1-17y', '18-60y', '>60y', 'Unknown')
gender_levels <- c('TOTAL', 'F', 'M', 'O')
ward_levels <- c('TOTAL', 'General Medicine', 'Emergency', 'Ophthalmology', 
                 'Pediatrics', 'Orthopedics', 'ENT', 'Gynecology', 'Dermatology',
                 'Surgery', 'Dentistry', 'Phychiatry', 'Oncology')

#### ------------------ ####
#### PROCESS HEALTH DATA ###
#### ------------------ ####
# 1.0. Load data ####

# # Unpack HISP datasets
# # setwd('C:/Users/claer14/OneDrive - University of Cambridge/III. Research/Data/')
# # unpack_HISP(hospital = 'Chamba',
# #             file_path='1_Original_Data/108_HISP/20250326_Chamba_Data.csv',
# #             download_path='2_Active_Data/201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba.csv')

setwd(path)

# Load data
# NB: WorldPop only has national holidays from 2010-19. I supplement with data from "Office Holidays"
df1 <- read.csv('201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba_raw.csv')
hols <- read.csv('204_demo/glo_Adm0_WorldPop_2010-2019_public-school_holidays_daily.csv', check.names = FALSE, fileEncoding = 'UTF-8')
hols_hp <- read.csv('204_demo/ind_Adm2_OfficeHolidays_2015-2023_HP_holidays_daily.csv', check.names = FALSE, fileEncoding = 'UTF-8')

# Set data types
(column_types <- sapply(df1, class))
df1$date <- as.Date(df1$date)
df1$age <- as.numeric(df1$age)
hols$Date <- as.Date(hols$Date)
hols$weekday <- weekdays(hols$Date)

# 1.1. Address missing values and set time series range ####

# Count missing values per column
(missing_counts <- colSums(is.na(df1)))
# Check rows with any missing values
(df1_missing_any <- df1[rowSums(is.na(df1)) > 0, ])
# Remove rows with NA dates
df1 <- df1 %>% filter(!is.na(date))

# Check for missing days: 810 days missing between '2015-01-20' and '2024-06-26'
date_range <- seq(from = min(df1$date), to = max(df1$date), by = 'day')
(missing_days <- as.Date(setdiff(date_range, df1$date)))

# Cut large gaps in days: 268 days missing between '2015-02-14' and '2023-01-19'
startdate <- as.Date('2015-02-14')
enddate <- as.Date('2023-01-19')
df1 <- df1[df1$date >= startdate & df1$date <= enddate,]
date_range <- seq(from = min(df1$date), to = max(df1$date), by = 'day')
(missing_days <- as.Date(setdiff(date_range, df1$date)))


# 1.2. Assess data types and values ####

# WARD

# Change names
df1 <- df1 %>%
  mutate(ward_name = case_when(
    ward_name == "CHILD OPD ROOM NO 311" ~ "Pediatrics",
    ward_name == "DENTAL OPD ROOM NO 316" ~ "Dentistry",
    ward_name == "EMERGENCY ROOM NO 1" ~ "Emergency",
    ward_name == "ENT OPD ROOM NO 418" ~ "ENT", # Ears nose throat
    ward_name == "EYE OPD ROOM NO 423" ~ "Ophthalmology", # Eye care
    ward_name == "GENERAL MEDICINE OPD ROOM NO 100 101 " ~ "General Medicine", # Internal medicine & Chronic disease
    ward_name == "GENERAL SURGERY OPD ROOM NO 209" ~ "Surgery",
    ward_name == "GYNAE OPD ROOM NO 313" ~ "Gynecology",
    ward_name == "OPTHMALMOLOGY OPD ROOM NO 424" ~ "Ophthalmology",  # # Eye care - Duplicate
    ward_name == "ORTHO OPD ROOM NO 207" ~ "Orthopedics", # Bone & joint
    ward_name == "PHYCHIATRY OPD ROOM NO 422" ~ "Phychiatry",
    ward_name == "RADIOTHERAPY AND ONCOLOGY OPD ROOM NO 421" ~ "Oncology", # Cancer
    ward_name == "SKIN OPD ROOM NO 317" ~ "Dermatology",
    TRUE ~ ward_name  # Keep original if no match
  ))
sort(unique(df1$ward_name))
table(df1$ward_name)

# AGE GROUPS
### NB: This is sectioned into infants, children, adults and over 60s
### - NOT a standardized classification but ONS provides some (for UK)
### https://www.ons.gov.uk/census/census2021dictionary/variablesbytopic/demographyvariablescensus2021/age/classifications

# Classify age groups
df1 <- df1 %>%
  mutate(age_group = case_when(
    age < 0 ~ 'Unknown',
    age == 0 ~ '<1y',
    age >=1 & age < 18 ~ '1-17y',
    age >= 18 & age <= 60 ~ '18-60y',
    age > 60 & age <= 100 ~ '>60y',
    age > 100 ~ 'Unknown',
    TRUE ~ 'Unknown'
  )) %>%
  mutate(age_group = factor(age_group, levels = c(
    '<1y', '1-17y', '18-60y', '>60y', 'Unknown'
  )))
sort(unique(df1$age))
table(df1$age_group)
df1 %>% filter(age_group == 'Unknown') %>% count(ward_name) # No clear pattern

# GENDER
sort(unique(df1$gender))
table(df1$gender)
df1 %>% filter(gender == 'O') %>% count(ward_name) # More than half are in Child OPD

# 1.3. Get daily counts ####

# NB: need to double check if there are multiple daily visits by individual
# patients and figure out how to visit_type fits in !!!

# GET wide FORMAT DAILY

# Total daily count of patients
total <- df1 %>%
  group_by(location, date) %>%
  summarise(total_patients = n())

# Daily count by ward
ward_columns <- c("Dentistry", "Dermatology", "Emergency", "ENT",
                  "General Medicine", "Gynecology", "Oncology", "Ophthalmology",
                  "Orthopedics", "Pediatrics", "Phychiatry", "Surgery")
ward_wide <- df1 %>%
  group_by(location, date, ward_name) %>%
  summarise(patients = n()) %>%
  pivot_wider(names_from = ward_name,
              values_from = patients,
              values_fill = list(patients = 0)) %>%
  select(date, starts_with(ward_columns[1]), starts_with(ward_columns[2]),
         starts_with(ward_columns[3]), starts_with(ward_columns[4]),
         starts_with(ward_columns[5]), starts_with(ward_columns[6]),
         starts_with(ward_columns[7]), starts_with(ward_columns[8]),
         starts_with(ward_columns[9]), starts_with(ward_columns[10]),
         starts_with(ward_columns[11]), starts_with(ward_columns[12]))

# Aggregate by age
age_wide <- df1 %>%
  group_by(location, date, age_group) %>%
  summarise(patients = n()) %>%
  pivot_wider(names_from = age_group,
              values_from = patients,
              values_fill = list(patients = 0))

# Aggregate by gender
gender_wide <- df1 %>%
  group_by(location, date, gender) %>%
  summarise(patients = n()) %>%
  pivot_wider(names_from = gender,
              values_from = patients,
              values_fill = list(patients = 0))

# Aggregate by ward & age
age_ward_wide <- df1 %>%
  group_by(location, date, ward_name, age_group) %>%
  summarise(patients = n()) %>%
  pivot_wider(names_from = c(ward_name, age_group),
              values_from = patients,
              values_fill = list(patients = 0))

# Aggregate by ward & gender
gender_ward_wide <- df1 %>%
  group_by(location, date, ward_name, gender) %>%
  summarise(patients = n()) %>%
  pivot_wider(names_from = c(ward_name, gender),
              values_from = patients,
              values_fill = list(patients = 0))

# Join aggregate ward data on location and date
ward_disag <- ward_wide %>%
  left_join(gender_ward_wide, c('location', "date")) %>%
  left_join(age_ward_wide, c('location', "date"))  %>%
  select(date, starts_with(ward_columns[1]), starts_with(ward_columns[2]),
         starts_with(ward_columns[3]), starts_with(ward_columns[4]),
         starts_with(ward_columns[5]), starts_with(ward_columns[6]),
         starts_with(ward_columns[7]), starts_with(ward_columns[8]),
         starts_with(ward_columns[9]), starts_with(ward_columns[10]),
         starts_with(ward_columns[11]), starts_with(ward_columns[12]))

df1_daily_wide <- total %>%
  left_join(gender_wide, by = c('location', "date")) %>%
  left_join(age_wide, by = c('location', "date")) %>%
  left_join(ward_disag, by = c('location', "date")) #%>%
  # select(-matches(c("_O", "Unknown"), ignore.case = TRUE))

# GET long FORMAT DAILY

total_long <- total %>%
  mutate(ward_name = 'TOTAL', age_group = 'TOTAL', gender = 'TOTAL', patients = total_patients) %>%
  select(date, ward_name, age_group, gender, patients)

ward_long <- df1 %>%
  group_by(location, date, ward_name) %>%
  summarise(patients = n()) %>%
  mutate(age_group = 'TOTAL', gender = 'TOTAL', patients = patients) %>%
  select(date, ward_name, age_group, gender, patients)

age_long <- df1 %>%
  group_by(location, date, age_group) %>%
  summarise(patients = n()) %>%
  mutate(ward_name = 'TOTAL', gender = 'TOTAL', patients = patients) %>%
  select(date, ward_name, age_group, gender, patients)

gender_long <- df1 %>%
  group_by(location, date, gender) %>%
  summarise(patients = n()) %>%
  mutate(ward_name = 'TOTAL', age_group = 'TOTAL', patients = patients) %>%
  select(date, ward_name, age_group, gender, patients)

age_gender_ward_long <- df1 %>%
  group_by(location, date, ward_name, age_group, gender) %>%
  summarise(patients = n())

df1_daily_long <- bind_rows(age_gender_ward_long, ward_long, age_long, gender_long, total_long)%>%
  arrange(date)

# 1.4. Near zero values and missing dates ####
# 1.4. A: Run checks ####

# 268 days missing between '2015-02-14' and '2023-01-19'
date_range <- seq(from = min(df1_daily_wide$date), to = max(df1_daily_wide$date), by = 'day')
(missing_days <- as.Date(setdiff(date_range, df1_daily_wide$date)))
missing_weekday <- factor(weekdays(missing_days), levels=weekday_levels)

# 319 days with total_patient values < 100
(zero_days <- df1_daily_wide[df1_daily_wide$total_patients < 100, ]$date)
zero_weekday <- factor(weekdays(zero_days), levels=weekday_levels)

# 587 days where low are missing in total
all_days <- append(missing_days, zero_days)
all_weekday <- factor(weekdays(all_days), levels=weekday_levels)

### Maybe the practice was closed on Sundays & bank holidays? Yes
### Why are there suddenly so many 0 patient values as of 2018, before up until mid 2020? Change in management style
table(missing_weekday, format(missing_days, '%Y')) # NB: none missing in 2020
table(zero_weekday, format(zero_days, '%Y')) # NB: none missing in 2015
table(all_weekday, format(all_days, '%Y'))

# Compare against 'Sundays' and bank holidays (get list of all holidays - sundays, bank/school)
hols_f <- hols %>% filter(ADM_name == "India" & (holiday == 1 | school == 1 | weekday == "Sunday"))
hols_s <- df1_daily_wide %>% filter(weekdays(date) == "Sunday") 
all_hols <- unique(c(hols_f$Date, hols_hp$Date, hols_s$date))

# 429 accounted for by holidays and 158 not
(account_hols <- as.Date(intersect(all_days, all_hols)))
(account_not <- as.Date(setdiff(all_days, all_hols)))

# Check unaccounted dates (mostly in 2020, which can be due to COVID)
account_n_weeks <- factor(weekdays(account_not), levels=weekday_levels)
table(account_n_weeks, format(account_not, '%Y'))

# 1.4. B: Address wide format ####

# Add in missing dates 
full_dates_df <- data.frame(date = date_range)
df1_daily_wide <- merge(full_dates_df, df1_daily_wide, by = "date", all.x = TRUE)
df1_daily_wide$location <- "Chamba"

# Add columns for weekdays and holiday (includes Sunday or bank or school holiday)
df1_daily_wide$weekday <- weekdays(df1_daily_wide$date)
df1_daily_wide$holiday <- ifelse(df1_daily_wide$date %in% all_hols, 1, 0)
df1_daily_wide <- df1_daily_wide[, c('location', "date", 'weekday', "holiday", setdiff(names(df1_daily_wide), c('location', "date", 'weekday', "holiday")))]

# Set NAs for holidays & values < 50 (partial closure days skews model)
df1_daily_wide <- df1_daily_wide %>% mutate(total_patients = ifelse(holiday==1 & total_patients < 50, NA, total_patients))

# Expand NA values across other health variables
cols_exclude <- c("location", "date", "weekday", "holiday", "total_patients")
df1_daily_wide <- df1_daily_wide %>%
  mutate(across(.cols = -all_of(cols_exclude),
                .fns = ~ ifelse(is.na(total_patients), NA, .),
                .names = "{.col}"))
(missing_counts <- colSums(is.na(df1_daily_wide)))
(na_counts <- df1_daily_wide %>% filter(is.na(total_patients)) %>% group_by(holiday) %>% summarise(count = n()))

# 1.4. C: Address long format ####

# Add in missing dates and full grid of category combinations
full_grid <- expand.grid(
  location = unique(df1_daily_long$location),
  date = seq(min(df1_daily_long$date), max(df1_daily_long$date), by = "day"),
  ward_name = unique(df1_daily_long$ward_name),
  age_group = unique(df1_daily_long$age_group),
  gender = unique(df1_daily_long$gender),
  stringsAsFactors = FALSE
)

# Add holidays
full_grid$holiday <- ifelse(full_grid$date %in% all_hols, 1, 0)

# Left join original data onto the full grid, sorting by factor levels and date
df_long <- full_grid %>%
  left_join(df1_daily_long, by = c("location", "date", "ward_name", "age_group", "gender")) %>%
  mutate(
    date = as.Date(date),
    weekday = factor(weekdays(date), levels = weekday_levels), 
    age_group = factor(age_group, levels = age_levels),
    gender = factor(gender, levels = gender_levels),
    ward_name = factor(ward_name, levels = ward_levels) 
  ) %>%
  arrange(date, ward_name, age_group, gender) %>%
  select(location, date, weekday, holiday, everything())

# Set all previously identified dates to NA & add holidays
na_dates <- df1_daily_wide %>% filter(is.na(total_patients)) %>% pull(date)
df1_daily_long <- df_long %>% mutate(patients = ifelse(date %in% na_dates, NA, patients))

# 1.5. Export initial health data sets ####

write.csv(df1, file='201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba_clean.csv', row.names = FALSE)
write.csv(df1_daily_wide, file='201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba_clean_wide.csv', row.names = FALSE)
write.csv(df1_daily_long, file='201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba_clean_long.csv', row.names = FALSE)

#### ------------------ ####
#### JOIN WEATHER DATA ###
#### ------------------ ####

setwd(path)
df1_daily_weat <- read.csv('202_weat/ind_30km_JASMIN-ERA5_2015-2023_chamba.csv')
df_heal <- read.csv('201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba_clean.csv', check.names = FALSE, fileEncoding = "UTF-8")
df_heal_w <- read.csv('201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba_clean_wide.csv', check.names = FALSE, fileEncoding = "UTF-8")
df_heal_l <- read.csv('201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba_clean_long.csv', check.names = FALSE, fileEncoding = "UTF-8")

# Fix dtypes
(column_types <- sapply(df_heal, class))
(column_types <- sapply(df_heal_w, class))
(column_types <- sapply(df_heal_l, class))
(column_types <- sapply(df1_daily_weat, class))
df1_daily_weat <- df1_daily_weat %>% rename(date = time)
df1_daily_weat$date <- as.Date(df1_daily_weat$date)
df_heal$date <- as.Date(df_heal$date)
df_heal_w$date <- as.Date(df_heal_w$date)
df_heal_l$date <- as.Date(df_heal_l$date)

# Fix weather data
df1_daily_weat <- df1_daily_weat %>% relocate(fg10 , .before = v10 ) %>% relocate(sp , .after = u10 )
temp <- c('t2m', 'd2m', 'mn2t', 'mx2t')
for (t in temp) {
  df1_daily_weat[[t]] <- df1_daily_weat[[t]]-272.15
}

# Join dfs
df_all <- df_heal %>%
  left_join(df1_daily_weat %>% distinct(location, date, .keep_all = TRUE), by = c('location', 'date')) %>%
  select(-lat, -lon, everything(), lat, lon)

df_long <- df_heal_l %>%
  left_join(df1_daily_weat %>% distinct(location, date, .keep_all = TRUE), by = c('location', 'date')) %>%
  select(-lat, -lon, everything(), lat, lon)

df_wide <- df_heal_w %>%
  left_join(df1_daily_weat %>% distinct(location, date, .keep_all = TRUE), by = c('location', 'date'))

# Add in year and week controls
df_all <- df_all %>%
  mutate(
    year = year(date),
    month = month(date, label = TRUE, abbr = TRUE),
    week = isoweek(date),
    weekday = weekdays(date)
  ) %>%
  select(location, date, year, month, week, weekday, everything())

df_long <- df_long %>%
  mutate(
    year = year(date),
    month = month(date, label = TRUE, abbr = TRUE),
    week = isoweek(date)
  ) %>%
  select(location, date, year, month, week, weekday, everything())

df_wide <- df_wide %>%
  mutate(
    year = year(date),
    month = month(date, label = TRUE, abbr = TRUE),
    week = isoweek(date)
  ) %>%
  select(location, date, year, month, week, weekday, everything())



# Export
write.csv(df_all, file='203_comb/ind_30km_HISP-HP_2015-2023_hosp-era5_chamba_all.csv', row.names = FALSE)
write.csv(df_long, file='203_comb/ind_30km_HISP-HP_2015-2023_hosp-era5_chamba_long.csv', row.names = FALSE)
write.csv(df_wide, file='203_comb/ind_30km_HISP-HP_2015-2023_hosp-era5_chamba_wide.csv', row.names = FALSE)

#####################################

# Questions: 
### What do the category values mean (e.g. General, RSBY, Child Less Than 1yr, Antenatal, Staff, BPL, Other Free)?
### Why are there weird numbers in age (e.g. -1121, 1951, REG-INITIAL, REG-REVISIT)
### There are 268 missing days between "2015-02-14" and "2023-01-19",  155 of which are on Sundays - is this because the practice is closed on sundays & bank holidays?
### Why are there suddenly so many 0 patient values as of 2018, before up until mid 2020? (see time series plot attached)


# Join data sets
# df1_daily_wide <- total %>%
#   left_join(gender_wide, by = c('location', "date")) %>%
#   left_join(age_wide, by = c('location', "date")) %>%
#   left_join(ward_wide, by = c('location', "date")) %>%
#   select(-matches(c("_O", "Unknown"), ignore.case = TRUE))

# library(imputeTS)
# # Apply kalman filter on other rows (n=26)
# #### Remove holidays and <50 NAs
# no_hols <- df1_daily_wide %>% filter(!(holiday==1 & is.na(total_patients)))
# (missing_counts <- colSums(is.na(no_hols))) # n=56
# #### Apply Kalman filter only to the non-holiday NAs
# no_hols$imputed_value <- na_kalman(no_hols$total_patients, model = "StructTS")
# #### Merge the imputed values back to the main dataset
# df1_daily_wide <- df1_daily_wide %>%
#   left_join(no_hols %>% select(date, imputed_value), by = "date") %>%
#   mutate(total_patients = ifelse(is.na(total_patients) & !holiday & total_patients >= 50, imputed_value, total_patients)) %>%
#   select(-imputed_value)  # Drop the imputed_value column after merging
# (missing_counts <- colSums(is.na(df1_daily_wide)))
