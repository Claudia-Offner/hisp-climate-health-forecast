### ------------------------------------------------------------------------ ### 
### HP Hospital cleaning - HISP Data                                         ###
### ------------------------------------------------------------------------ ### 

setwd(here())
source('1_Data/HISP_Himachal_Pradesh/HISP_Cleaning - Dependencies.R')

# Set data paths
setwd(path_data)
heal_export <- 'ind_adm4_HISP-HP_2015-2023_patient_visits_ALL_clean'

####  FUNCTIONS  ####

# Unpack the tt separated HISP data
unpack_HISP <- function(file_path, hospital, download_path){
  

  # Replace commas with a different character
  txt <- readLines(file_path)
  txt <- gsub(',', ':', txt)
  writeLines(txt, file_path)
  
  # Load data
  data <- read.csv(file_path, header = FALSE, col.names = 'raw_column')

  # Splitting the columns & organising
  df_heal <- data %>%
    separate(raw_column, into = c('ward_name', 'patient_id', 'age', 'gender', 'category', 'date', 'visit_type'),
             sep = 'tt', fill = 'right') %>%
    mutate(location = hospital)

  # Reorganize col order & set type
  df_heal <- df_heal[, c('location', 'date', 'patient_id', 'ward_name', 'age', 'gender', 'category', 'visit_type')]

  # Export
  if (!file.exists(download_path)) {
    # File does not exist → create new CSV
    write.csv(df_heal, file = download_path, row.names = FALSE)
  } else {
    # File exists → append without column headers
    write.table(df_heal, file = download_path, sep = ',', 
                col.names = FALSE, row.names = FALSE, append = TRUE)
  }

}

####  MAIN CODE ####

# # Unpack HISP datasets
# datasets <- list.files('1_Original_Data/108_HISP/')
# for (s in location_levels){
#   file <- datasets[grepl(s, datasets)]
#   if (length(file) == 0) {
#     message(paste0('No files found for ', s, ' — skipping.'))
#   } else {
#     message(paste0('Extracting files for ', s))
#     unpack_HISP(hospital = s,
#               file_path=paste0('1_Original_Data/108_HISP/', file),
#               download_path='2_Active_Data/201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_ALL_raw.csv')
#   }
# }

# Load data
# NB: WorldPop only has national holidays from 2010-19, which I supplement with 
# regional holidays from 'Office Holidays'.
df_heal <- read.csv('2_Active_Data/201_heal/ind_adm4_HISP-HP_2015-2023_patient_visits_ALL_raw.csv')
hols <- read.csv('2_Active_Data/204_demo/glo_Adm0_WorldPop_2010-2019_public-school_holidays_daily.csv', check.names = FALSE, fileEncoding = 'UTF-8')
hols_hp <- read.csv('2_Active_Data/204_demo/ind_Adm2_OfficeHolidays_2015-2023_HP_holidays_daily.csv', check.names = FALSE, fileEncoding = 'UTF-8')

# Set data types
(column_types <- sapply(df_heal, class))
df_heal$date <- as.Date(df_heal$date)
df_heal$age <- as.numeric(gsub('[^0-9.]', '', df_heal$age))
hols$Date <- as.Date(hols$Date)
hols$weekday <- weekdays(hols$Date)

# 1. Address missing values and set time series range ####

# Count missing values per column
(missing_counts <- colSums(is.na(df_heal)))
# Check rows with any missing values
(df_heal_missing_any <- df_heal[rowSums(is.na(df_heal)) > 0, ])
# Remove rows with NA dates
df_heal <- df_heal %>% dplyr::filter(!is.na(date))

# Check for missing days for EACH site: 810+141+880 days missing between '2012-04-04' and '2024-06-26'
(missing_days <- df_heal %>% group_by(location) %>%
  summarise(min_date = min(date, na.rm = TRUE),
            max_date = max(date, na.rm = TRUE),
            missing_days = list(setdiff(seq(min_date, max_date, by = 'day'), date))))
# Check for pattern
for (l in unique(missing_days$location)){
  print(l)
  print(missing_days$missing_days[missing_days$location == l])
}

# Cut large gaps in days: 268+41+557 days missing between '2015-02-14' and '2023-01-19'
startdate <- as.Date('2015-02-14')
enddate <- as.Date('2023-01-19')
df_heal <- df_heal[df_heal$date >= startdate & df_heal$date <= enddate,]
date_range <- seq(from = min(df_heal$date), to = max(df_heal$date), by = 'day')
(missing_days <- df_heal %>% group_by(location) %>%
    summarise(min_date = min(date, na.rm = TRUE),
              max_date = max(date, na.rm = TRUE),
              missing_days = list(setdiff(seq(min_date, max_date, by = 'day'), date))))

# 2. Assess data types and values ####

# WARD
sort(table(df_heal$ward_name))
sort(unique(df_heal$ward_name))

ward_patterns <- c(
  'art|ictc' = 'antiretroviral',
  'skin'  = 'dermatology',
  'disability' = 'disability',
  'DEIC' = 'paediatrics',
  'emergency'='emergency',
  'ENT'= 'ent',
  'surgery|surgical' = 'surgery',
  'medicine|medical' = 'general',
  'general|new' = 'general',
  'gynae' = 'gynaecology',
  'mch' = 'mch',
  'NAI DISHA KENDRA' = 'substance_counselling',
  'ncd' = 'ncd',
  'ortho' = 'orthopedics',
  'cancer' = 'oncology',
  'radio|oncolo' = 'oncology',
  'eye|opth' = 'ophthalmology',
  'child|10 to 19 years' = 'paediatrics',
  'physio' = 'physiotherapy', 
  'pulmon' = 'pulmonary',
  'psychiatry|phych' = 'psychiatry',
  'dental' = 'dentistry',
  'FEMALE|MALE' = 'general'
  )

df_heal$ward_group <- df_heal$ward_name
for (pattern in names(ward_patterns)) {
  df_heal$ward_group[grepl(pattern, df_heal$ward_name, ignore.case = TRUE)] <- ward_patterns[pattern]
}
sort(unique(df_heal$ward_group))
sort(table(df_heal$ward_group))
df_heal$ward_name <- df_heal$ward_group
df_heal$ward_group <- NULL


# AGE GROUPS
### NB: This is sectioned into infants, children, adults and over 60s
### - NOT a standardized classification but ONS provides some (for UK)
### https://www.ons.gov.uk/census/census2021dictionary/variablesbytopic/demographyvariablescensus2021/age/classifications
### India  0-14 years (children), 15-24 years (early working age),  25-54 years 
### (prime working age), 55-64 years (mature working age), 65 years and over (elderly)

# Classify age groups
df_heal <- df_heal %>%
  mutate(age_group = case_when(
    age < 0 ~ 'Unknown',
    age >= 0 & age < 5 ~ 'under_5',
    age >=5 & age < 16 ~ 'age5_14',
    age >= 15 & age < 25 ~ 'age15_24',
    age >= 25 & age < 65 ~ 'age25_64',
    age >= 65 & age <= 110 ~ 'over_65',
    age > 110 ~ 'Unknown',
    TRUE ~ 'Unknown'
  )) %>%
  mutate(age_group = factor(age_group, levels = age_levels))
sort(unique(df_heal$age))
table(df_heal$age_group)
df_heal %>% dplyr::filter(age_group == 'Unknown') %>% count(ward_name) # No clear pattern

# GENDER
sort(unique(df_heal$gender))
table(df_heal$gender)
df_heal %>% dplyr::filter(gender == 'O') %>% count(ward_name) # More than half are in Child OPD

# 3. Get daily counts ####

# NB: need to double check if there are multiple daily visits by individual
# patients and figure out how to visit_type fits in !!!

# GET wide FORMAT DAILY
# Total daily count of patients
total <- df_heal %>%
  group_by(location, date) %>%
  summarise(total_patients = n())

# Daily count by ward
ward_wide <- df_heal %>%
  group_by(location, date, ward_name) %>%
  summarise(patients = n(), .groups = 'drop') %>%
  pivot_wider(names_from = ward_name, values_from = patients, values_fill = 0) %>%
  select(location, date, !!!syms(ward_levels))

# Aggregate by age
age_wide <- df_heal %>%
  group_by(location, date, age_group) %>%
  summarise(patients = n()) %>%
  pivot_wider(names_from = age_group, values_from = patients,
              values_fill = list(patients = 0)) %>%
  select(-matches('Unknown'))

# Aggregate by gender
gender_wide <- df_heal %>%
  group_by(location, date, gender) %>%
  summarise(patients = n()) %>%
  pivot_wider(names_from = gender, values_from = patients,
              values_fill = list(patients = 0)) %>%
  select(-matches('O'))

# Aggregate by ward & age
age_gender_wide <- df_heal %>%
  group_by(location, date, gender, age_group) %>%
  summarise(patients = n()) %>%
  pivot_wider(names_from = c(gender, age_group),
              values_from = patients,
              values_fill = list(patients = 0)) %>%
  select(-matches('Unknown|O_'))

# # Aggregate by ward & age
# age_ward_wide <- df_heal %>%
#   group_by(location, date, ward_name, age_group) %>%
#   summarise(patients = n()) %>%
#   pivot_wider(names_from = c(ward_name, age_group),
#               values_from = patients,
#               values_fill = list(patients = 0))
# 
# # Aggregate by ward & gender
# gender_ward_wide <- df_heal %>%
#   group_by(location, date, ward_name, gender) %>%
#   summarise(patients = n()) %>%
#   pivot_wider(names_from = c(ward_name, gender),
#               values_from = patients,
#               values_fill = list(patients = 0))

# Join aggregate ward data on location and date
ward_select_exprs <- rlang::exprs(location, date, !!!lapply(ward_levels, function(w) rlang::expr(starts_with(!!w))))
ward_disag <- list(ward_wide) %>%
  reduce(left_join, by = c('location', 'date')) %>%
  select(!!!ward_select_exprs)

df_heal_w <- total %>%
  left_join(gender_wide, by = c('location', 'date')) %>%
  left_join(age_wide, by = c('location', 'date')) %>%
  left_join(ward_disag, by = c('location', 'date')) %>%
  left_join(age_gender_wide, by = c('location', 'date'))

# GET long FORMAT DAILY
total_long <- total %>%
  mutate(ward_name = 'TOTAL', age_group = 'TOTAL', gender = 'TOTAL', patients = total_patients) %>%
  select(date, ward_name, age_group, gender, patients)

ward_long <- df_heal %>%
  group_by(location, date, ward_name) %>%
  summarise(patients = n()) %>%
  mutate(age_group = 'TOTAL', gender = 'TOTAL', patients = patients) %>%
  select(date, ward_name, age_group, gender, patients)

age_long <- df_heal %>%
  group_by(location, date, age_group) %>%
  summarise(patients = n()) %>%
  mutate(ward_name = 'TOTAL', gender = 'TOTAL', patients = patients) %>%
  select(date, ward_name, age_group, gender, patients)

gender_long <- df_heal %>%
  group_by(location, date, gender) %>%
  summarise(patients = n()) %>%
  mutate(ward_name = 'TOTAL', age_group = 'TOTAL', patients = patients) %>%
  select(date, ward_name, age_group, gender, patients)

age_gender_long <- df_heal %>%
  group_by(location, date, age_group, gender) %>%
  summarise(patients = n()) %>% # ward_name
  mutate(ward_name = 'TOTAL') %>%
  select(location, date, ward_name, everything())

df_heal_l <- bind_rows(age_gender_long, ward_long, age_long, gender_long, total_long)%>%
  arrange(location, date)

# 4. Near zero values and missing dates ####
# 4. A: Checks on total_patients ####

# 2065 days missing between '2015-02-14' and '2023-01-19'
(missing <- df_heal_w %>%
  group_by(location) %>%
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
sum(missing$n_missing)
missing_days <- missing$missing_days
missing_weekday <- missing$missing_weekday

# 407 days with total_patient values < 100
(zero <- df_heal_w %>%
  group_by(location) %>%
  summarise(
    zero_days = list(date[total_patients < 100]),
    .groups = 'drop'
  ) %>%
  mutate(
    zero_weekday = lapply(zero_days, function(d) {
      factor(format(d, "%a"), levels = weekday_levels)
    }),
    n_zero_days = sapply(zero_days, length)
  ))
sum(zero$n_zero_days)
zero_days <- zero$zero_days
zero_weekday <- zero$zero_weekday

# 2472 days were low or missing in total
(all <- missing %>%
  inner_join(zero, by = 'location') %>%
  mutate(
    all_days = Map(c, missing_days, zero_days),   # append lists
    all_weekday = lapply(all_days, function(d) {
      factor(format(d, "%a"), levels = weekday_levels)
    }),
    n_all_days = sapply(all_days, length)
  ))
sum(all$n_all_days)
all_days <- all$all_days
all_weekday <- all$all_weekday

### Maybe the practice was closed on Sundays & bank holidays? Yes
### Why are there suddenly so many 0 patient values as of 2018, before up until mid 2020? Change in management style
missing_tables <- lapply(seq_along(missing_days), function(i) {
  table(missing_weekday[[i]], format(missing_days[[i]], '%Y'))
})
names(missing_tables) <- missing$location  # label by location
missing_tables

zero_tables <- lapply(seq_along(zero_days), function(i) {
  table(zero_weekday[[i]], format(zero_days[[i]], '%Y'))
})
names(zero_tables) <- zero$location
zero_tables

all_tables <- lapply(seq_along(all_days), function(i) {
  table(all_weekday[[i]],format(all_days[[i]], '%Y'))
})
names(all_tables) <- all_days$location
all_tables

# Compare against 'Sundays' and bank holidays (get list of all holidays - Sundays, public holiday)
hols_f <- hols %>% dplyr::filter(ADM_name == 'India' & (holiday == 1 | weekday == 'Sunday')) 
hols_s <- data.frame(date = date_range) %>% dplyr::filter(weekdays(date) == 'Sunday') 
all_hols <- unique(c(hols_f$Date, hols_hp$Date, hols_s$date))

# 2156 accounted for by holidays (87%) and 316 not (13%)
account_hols_list <- lapply(all_days, function(d) {intersect(d, all_hols)})
account_not_list <- lapply(all_days, function(d) {setdiff(d, all_hols)})

# Counts per location
(account_summary <- data.frame(
  location = all$location,
  n_all_days = sapply(all$all_days, length),
  n_account_hols = sapply(account_hols_list, length),
  pct_account_hols = sapply(account_hols_list, length) / sapply(all$all_days, length) * 100,
  n_account_not = sapply(account_not_list, length),
  pct_account_not = sapply(account_not_list, length) / sapply(all$all_days, length) * 100
))
sum(account_summary$n_account_hols)/(sum(account_summary$n_account_hols)+sum(account_summary$n_account_not))*100
sum(account_summary$n_account_not)/(sum(account_summary$n_account_hols)+sum(account_summary$n_account_not))*100

# Check unaccounted dates (mostly in 2020, which can be due to COVID, but some random weird years)
account_not_tables <- lapply(account_not_list, function(d) {
  table(factor(format(d, "%a"), levels = weekday_levels),format(d, '%Y'))
})
names(account_not_tables) <- all$location
account_not_tables

# 4. B: Checks across wards ####

# Combine location names and corresponding zero_days
filtered_list <- mapply(function(loc, days) {
  df_heal_w %>%
    dplyr::filter(location == loc & date %in% days) %>%
    select(location, date, total_patients, ward_levels)
}, loc = zero$location, days = zero$zero_days, SIMPLIFY = FALSE)

# Combine into a single data frame
filtered_df_heal <- bind_rows(filtered_list)
filtered_df_heal$location <- as.character(filtered_df_heal$location)

# Compare zeros vs non-zeros across wards
count_zeros <- filtered_df_heal %>%
  group_by(location) %>%
  summarise(across(-c(date), ~sum(. == 0)), .groups = 'drop')%>%
  mutate(count = 'Zeros')
count_non_zeros <- filtered_df_heal %>%
  group_by(location) %>%
  summarise(across(-c(date), ~sum(. != 0)), .groups = 'drop')%>%
  mutate(count = 'Non-zeros')
(x <- bind_rows(count_zeros, count_non_zeros) %>%
  select(count, location, everything()))


# 4. C: Address wide format ####

# Add in missing dates 
full_dates_df_heal <- expand.grid(location = unique(df_heal_w$location),
                             date = date_range)
df_heal_w <- full_dates_df_heal %>%
  left_join(df_heal_w, by = c('location', 'date')) %>%
  arrange(location, date)

# Add columns for weekdays, hospital closures (NA or total patients <50)
# holiday (includes Sunday or bank or school holiday)
df_heal_w$weekday <- weekdays(df_heal_w$date)
df_heal_w$closed <- ifelse(is.na(df_heal_w$total_patients) | df_heal_w$total_patients < 50, 1, 0)
df_heal_w$holiday <- ifelse(df_heal_w$date %in% all_hols, 1, 0)
df_heal_w <- df_heal_w[, c('location', 'date', 'weekday', 'closed', 'holiday', 
                           setdiff(names(df_heal_w), c('location', 'date', 'weekday', 'closed', 'holiday')))]

# Set NAs for holidays & values < 50 (partial closure days skews model)
# 2025-05-12: We will change this to remove all days < 50 because this is messing with model interpretability
df_heal_w <- df_heal_w %>% mutate(total_patients = ifelse(total_patients < 50, NA, total_patients))

# Expand NA values across other health variables
cols_exclude <- c('location', 'date', 'weekday', 'closed', 'holiday', 'total_patients')
df_heal_w <- df_heal_w %>%
  mutate(across(.cols = -all_of(cols_exclude),
                .fns = ~ ifelse(is.na(total_patients), NA, .),
                .names = '{.col}'))
(missing_counts <- colSums(is.na(df_heal_w)))
(na_counts <- df_heal_w %>% dplyr::filter(is.na(total_patients)) %>% group_by(holiday) %>% summarise(count = n()))


# 4. D: Address long format ####

# Add in missing dates and full grid of category combinations
full_grid <- expand.grid(
  location = unique(df_heal_l$location),
  date = seq(min(df_heal_l$date), max(df_heal_l$date), by = 'day'),
  ward_name = unique(df_heal_l$ward_name),
  age_group = unique(df_heal_l$age_group),
  gender = unique(df_heal_l$gender),
  stringsAsFactors = FALSE
)

# Add holidays & closures
closed_dates <- df_heal_w$date[df_heal_w$closed == 1]
full_grid$holiday <- ifelse(full_grid$date %in% all_hols, 1, 0)
full_grid$closed <- ifelse(full_grid$date %in% closed_dates, 1, 0)

# Left join original data onto the full grid, sorting by factor levels and date
df_heal_long <- full_grid %>%
  left_join(df_heal_l, by = c('location', 'date', 'ward_name', 'age_group', 'gender')) %>%
  mutate(
    date = as.Date(date),
    weekday = factor(weekdays(date), levels = weekday_levels), 
    age_group = factor(age_group, levels = c('TOTAL', age_levels)),
    gender = factor(gender, levels = c('TOTAL', gender_levels)),
    ward_name = factor(ward_name, levels = c('TOTAL', ward_levels)) 
  ) %>%
  arrange(location, date, ward_name, age_group, gender) %>%
  select(location, date, weekday, closed, holiday, everything())


# Set all previously identified dates to NA & add holidays
na_dates <- df_heal_w %>% 
  dplyr::filter(is.na(total_patients)) %>% 
  select(location, date)

df_heal_l <- df_heal_long %>%
  left_join(na_dates %>% mutate(na_flag = TRUE), by = c('location', 'date')) %>%
  mutate(patients = ifelse(!is.na(na_flag), NA, patients)) %>%
  select(-na_flag)


#### EXPORT ####

save(df_heal_w, file=paste0('2_Active_Data/201_heal/', heal_export, '_wide.RData'))
save(df_heal_l, file=paste0('2_Active_Data/201_heal/', heal_export, '_long.RData'))
save(df_heal,  file=paste0('2_Active_Data/201_heal/', heal_export, '_all.RData'))
save(df_heal_w, df_heal_l, df_heal, file=paste0('2_Active_Data/201_heal/', heal_export, '.RData'))
