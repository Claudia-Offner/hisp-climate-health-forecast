### ------------------------------------------------------------------------ ### 
### Hospital cleaning - HISP Data                                            ###
### ------------------------------------------------------------------------ ### 

# Turn off scientific notation
options(scipen=999)

# Set data path
user <- 'claer14'
path <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
setwd(path)

####  PACKAGES & FUNCTIONS ####
library(dplyr)
library(purrr)
library(tidyverse)
library(lubridate)
library(tidyselect)

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
age_levels <- c('under_5', 'age5_14', 'age15_24', 'age25_64',  'over_65','Unknown')
gender_levels <- c('F', 'M', 'O')
ward_levels <- c('general', 'emergency', 'ophthalmology', 'pediatrics', 
                 'orthopedics', 'ent', 'gynecology', 'dermatology', 'surgery', 
                 'dentistry', 'psychiatry', 'oncology')

####  MAIN CODE ####

# 1. Load data ####

# # Unpack HISP datasets
# # setwd('C:/Users/claer14/OneDrive - University of Cambridge/III. Research/Data/')
# # unpack_HISP(hospital = 'Chamba',
# #             file_path='1_Original_Data/108_HISP/20250326_Chamba_Data.csv',
# #             download_path='2_Active_Data/201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba.csv')

# Load data
# NB: WorldPop only has national holidays from 2010-19, which I supplement with 
# regional holidays from 'Office Holidays'.
df1 <- read.csv('201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba_raw.csv')
hols <- read.csv('204_demo/glo_Adm0_WorldPop_2010-2019_public-school_holidays_daily.csv', check.names = FALSE, fileEncoding = 'UTF-8')
hols_hp <- read.csv('204_demo/ind_Adm2_OfficeHolidays_2015-2023_HP_holidays_daily.csv', check.names = FALSE, fileEncoding = 'UTF-8')

# Set data types
(column_types <- sapply(df1, class))
df1$date <- as.Date(df1$date)
df1$age <- as.numeric(df1$age)
hols$Date <- as.Date(hols$Date)
hols$weekday <- weekdays(hols$Date)

# 2. Address missing values and set time series range ####

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


# 3. Assess data types and values ####

# WARD
df1 <- df1 %>%
  mutate(ward_name = case_when(
    ward_name == 'CHILD OPD ROOM NO 311' ~ 'pediatrics',
    ward_name == 'DENTAL OPD ROOM NO 316' ~ 'dentistry',
    ward_name == 'EMERGENCY ROOM NO 1' ~ 'emergency',
    ward_name == 'ENT OPD ROOM NO 418' ~ 'ent', # Ears nose throat
    ward_name == 'EYE OPD ROOM NO 423' ~ 'ophthalmology', # Eye care
    ward_name == 'GENERAL MEDICINE OPD ROOM NO 100 101 ' ~ 'general', # Internal medicine & Chronic disease
    ward_name == 'GENERAL SURGERY OPD ROOM NO 209' ~ 'surgery',
    ward_name == 'GYNAE OPD ROOM NO 313' ~ 'gynecology',
    ward_name == 'OPTHMALMOLOGY OPD ROOM NO 424' ~ 'ophthalmology',  # # Eye care - Duplicate
    ward_name == 'ORTHO OPD ROOM NO 207' ~ 'orthopedics', # Bone & joint
    ward_name == 'PHYCHIATRY OPD ROOM NO 422' ~ 'psychiatry',
    ward_name == 'RADIOTHERAPY AND ONCOLOGY OPD ROOM NO 421' ~ 'oncology', # Cancer
    ward_name == 'SKIN OPD ROOM NO 317' ~ 'dermatology',
    TRUE ~ ward_name  # Keep original if no match
  ))
sort(unique(df1$ward_name))
table(df1$ward_name)

# AGE GROUPS
### NB: This is sectioned into infants, children, adults and over 60s
### - NOT a standardized classification but ONS provides some (for UK)
### https://www.ons.gov.uk/census/census2021dictionary/variablesbytopic/demographyvariablescensus2021/age/classifications
### India  0-14 years (children), 15-24 years (early working age),  25-54 years 
### (prime working age), 55-64 years (mature working age), 65 years and over (elderly)

# Classify age groups
df1 <- df1 %>%
  mutate(age_group = case_when(
    age < 0 ~ 'Unknown',
    age >= 0 & age < 5 ~ 'under_5',
    age >=5 & age < 16 ~ 'age5_14',
    age >= 15 & age < 25 ~ 'age15_24',
    age >= 25 & age < 65 ~ 'age25_64',
    age >= 65 & age <= 100 ~ 'over_65',
    age > 100 ~ 'Unknown',
    TRUE ~ 'Unknown'
  )) %>%
  mutate(age_group = factor(age_group, levels = age_levels))
sort(unique(df1$age))
table(df1$age_group)
df1 %>% filter(age_group == 'Unknown') %>% count(ward_name) # No clear pattern

# GENDER
sort(unique(df1$gender))
table(df1$gender)
df1 %>% filter(gender == 'O') %>% count(ward_name) # More than half are in Child OPD

# 4. Get daily counts ####

# NB: need to double check if there are multiple daily visits by individual
# patients and figure out how to visit_type fits in !!!

# GET wide FORMAT DAILY

# Total daily count of patients
total <- df1 %>%
  group_by(location, date) %>%
  summarise(total_patients = n())

# Daily count by ward
ward_wide <- df1 %>%
  group_by(location, date, ward_name) %>%
  summarise(patients = n(), .groups = 'drop') %>%
  pivot_wider(names_from = ward_name, values_from = patients, values_fill = 0) %>%
  select(location, date, !!!syms(ward_levels))

# Aggregate by age
age_wide <- df1 %>%
  group_by(location, date, age_group) %>%
  summarise(patients = n()) %>%
  pivot_wider(names_from = age_group, values_from = patients,
              values_fill = list(patients = 0))

# Aggregate by gender
gender_wide <- df1 %>%
  group_by(location, date, gender) %>%
  summarise(patients = n()) %>%
  pivot_wider(names_from = gender, values_from = patients,
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
ward_select_exprs <- rlang::exprs(location, date, !!!lapply(ward_levels, function(w) rlang::expr(starts_with(!!w))))
ward_disag <- list(ward_wide, gender_ward_wide, age_ward_wide) %>%
  reduce(left_join, by = c('location', 'date')) %>%
  select(!!!ward_select_exprs)

df1_daily_wide <- total %>%
  left_join(gender_wide, by = c('location', 'date')) %>%
  left_join(age_wide, by = c('location', 'date')) %>%
  left_join(ward_disag, by = c('location', 'date')) #%>%
  # select(-matches(c('_O', 'Unknown'), ignore.case = TRUE))

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

# 5. Near zero values and missing dates ####
# 5. A: Run checks on total_patients ####

# 268 days missing between '2015-02-14' and '2023-01-19'
date_range <- seq(from = min(df1_daily_wide$date), to = max(df1_daily_wide$date), by = 'day')
missing_days <- as.Date(setdiff(date_range, df1_daily_wide$date))
missing_weekday <- factor(weekdays(missing_days), levels=weekday_levels)
(length(missing_days))

# 319 days with total_patient values < 100
zero_days <- df1_daily_wide[df1_daily_wide$total_patients < 100, ]$date
zero_weekday <- factor(weekdays(zero_days), levels=weekday_levels)
(length(zero_days))

# 587 days were low or missing in total
all_days <- append(missing_days, zero_days)
all_weekday <- factor(weekdays(all_days), levels=weekday_levels)
(length(all_days))

### Maybe the practice was closed on Sundays & bank holidays? Yes
### Why are there suddenly so many 0 patient values as of 2018, before up until mid 2020? Change in management style
table(missing_weekday, format(missing_days, '%Y')) # NB: none missing in 2020
table(zero_weekday, format(zero_days, '%Y')) # NB: none missing in 2015
table(all_weekday, format(all_days, '%Y'))

# Compare against 'Sundays' and bank holidays (get list of all holidays - sundays, public holiday)
hols_f <- hols %>% filter(ADM_name == 'India' & (holiday == 1 | weekday == 'Sunday')) 
hols_s <- data.frame(date = date_range) %>% filter(weekdays(date) == 'Sunday') 
all_hols <- unique(c(hols_f$Date, hols_hp$Date, hols_s$date))

# 438 accounted for by holidays (74.61%) and 149 not (25.38%)
account_hols <- as.Date(intersect(all_days, all_hols))
account_not <- as.Date(setdiff(all_days, all_hols))
length(account_hols)
length(account_hols)/length(all_days)*100
length(account_not)
length(account_not)/length(all_days)*100

# Check unaccounted dates (mostly in 2020, which can be due to COVID)
account_n_weeks <- factor(weekdays(account_not), levels=weekday_levels)
table(account_n_weeks, format(account_not, '%Y'))

# 5. B: Checks near zero days across wards ####

filtered_df <- df1_daily_wide %>%
  filter(date %in% zero_days) %>%
  select(location, date, total_patients, ward_levels)

count_zeros <- filtered_df %>%
  ungroup() %>% select(-c(location, date)) %>%
  summarise(across(everything(), ~sum(. == 0)))  # count 0s in each column

count_non_zeros <- filtered_df %>%
  ungroup() %>% select(-c(location, date)) %>%
  summarise(across(everything(), ~sum(. != 0)))  # count non-0s in each column

# Lowest proportion of near zero days were general and emergency
(x <- rbind(count_zeros, count_non_zeros) %>%
  mutate(count = c('Zeros', 'Non-zeros')) %>%
  select(c(count, everything())))


# 5. C: Address wide format ####

# Add in missing dates 
full_dates_df <- data.frame(date = date_range)
df1_daily_wide <- merge(full_dates_df, df1_daily_wide, by = 'date', all.x = TRUE)
df1_daily_wide$location <- 'Chamba'

# Add columns for weekdays, hospital closures (NA or total patients <50)
# holiday (includes Sunday or bank or school holiday)
df1_daily_wide$weekday <- weekdays(df1_daily_wide$date)
df1_daily_wide$closed <- ifelse(is.na(df1_daily_wide$total_patients) | df1_daily_wide$total_patients < 50, 1, 0)
df1_daily_wide$holiday <- ifelse(df1_daily_wide$date %in% all_hols, 1, 0)
df1_daily_wide <- df1_daily_wide[, c('location', 'date', 'weekday', 'closed', 'holiday', 
                                     setdiff(names(df1_daily_wide), c('location', 'date', 'weekday', 'closed', 'holiday')))]

# Set NAs for holidays & values < 50 (partial closure days skews model)
# 2025-05-12: We will change this to remove all days < 50 because this is messing with model interpretability
df1_daily_wide <- df1_daily_wide %>% mutate(total_patients = ifelse(total_patients < 50, NA, total_patients))

# Expand NA values across other health variables
cols_exclude <- c('location', 'date', 'weekday', 'closed', 'holiday', 'total_patients')
df1_daily_wide <- df1_daily_wide %>%
  mutate(across(.cols = -all_of(cols_exclude),
                .fns = ~ ifelse(is.na(total_patients), NA, .),
                .names = '{.col}'))
(missing_counts <- colSums(is.na(df1_daily_wide)))
(na_counts <- df1_daily_wide %>% filter(is.na(total_patients)) %>% group_by(holiday) %>% summarise(count = n()))


# 5. D: Address long format ####

# Add in missing dates and full grid of category combinations
full_grid <- expand.grid(
  location = unique(df1_daily_long$location),
  date = seq(min(df1_daily_long$date), max(df1_daily_long$date), by = 'day'),
  ward_name = unique(df1_daily_long$ward_name),
  age_group = unique(df1_daily_long$age_group),
  gender = unique(df1_daily_long$gender),
  stringsAsFactors = FALSE
)

# Add holidays & closures
closed_dates <- df1_daily_wide$date[df1_daily_wide$closed == 1]
full_grid$holiday <- ifelse(full_grid$date %in% all_hols, 1, 0)
full_grid$closed <- ifelse(full_grid$date %in% closed_dates, 1, 0)

# Left join original data onto the full grid, sorting by factor levels and date
df_long <- full_grid %>%
  left_join(df1_daily_long, by = c('location', 'date', 'ward_name', 'age_group', 'gender')) %>%
  mutate(
    date = as.Date(date),
    weekday = factor(weekdays(date), levels = weekday_levels), 
    age_group = factor(age_group, levels = c('TOTAL', age_levels)),
    gender = factor(gender, levels = c('TOTAL', gender_levels)),
    ward_name = factor(ward_name, levels = c('TOTAL', ward_levels)) 
  ) %>%
  arrange(date, ward_name, age_group, gender) %>%
  select(location, date, weekday, closed, holiday, everything())

# Set all previously identified dates to NA & add holidays
na_dates <- df1_daily_wide$date[is.na(df1_daily_wide$total_patients)]
df1_daily_long <- df_long %>% mutate(patients = ifelse(date %in% na_dates, NA, patients))

# 6. Export initial health data sets ####

write.csv(df1, file='201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba_clean.csv', row.names = FALSE)
write.csv(df1_daily_wide, file='201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba_clean_wide.csv', row.names = FALSE)
write.csv(df1_daily_long, file='201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_chamba_clean_long.csv', row.names = FALSE)


################################################################################

# Questions: 
### What do the category values mean (e.g. General, RSBY, Child Less Than 1yr, Antenatal, Staff, BPL, Other Free)?
### Why are there weird numbers in age (e.g. -1121, 1951, REG-INITIAL, REG-REVISIT)
### There are 268 missing days between '2015-02-14' and '2023-01-19',  155 of which are on Sundays - is this because the practice is closed on sundays & bank holidays?
### Why are there suddenly so many 0 patient values as of 2018, before up until mid 2020? (see time series plot attached)


# Join data sets
# df1_daily_wide <- total %>%
#   left_join(gender_wide, by = c('location', 'date')) %>%
#   left_join(age_wide, by = c('location', 'date')) %>%
#   left_join(ward_wide, by = c('location', 'date')) %>%
#   select(-matches(c('_O', 'Unknown'), ignore.case = TRUE))

# library(imputeTS)
# # Apply kalman filter on other rows (n=26)
# #### Remove holidays and <50 NAs
# no_hols <- df1_daily_wide %>% filter(!(holiday==1 & is.na(total_patients)))
# (missing_counts <- colSums(is.na(no_hols))) # n=56
# #### Apply Kalman filter only to the non-holiday NAs
# no_hols$imputed_value <- na_kalman(no_hols$total_patients, model = 'StructTS')
# #### Merge the imputed values back to the main dataset
# df1_daily_wide <- df1_daily_wide %>%
#   left_join(no_hols %>% select(date, imputed_value), by = 'date') %>%
#   mutate(total_patients = ifelse(is.na(total_patients) & !holiday & total_patients >= 50, imputed_value, total_patients)) %>%
#   select(-imputed_value)  # Drop the imputed_value column after merging
# (missing_counts <- colSums(is.na(df1_daily_wide)))
