
# Load packages
library(tidyverse)


# Load packages and dependencies 
setwd('C:/Users/claer14/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/201_heal/HISP/')


# 0. Format data ####

# # Load packages and dependencies 
# setwd('C:/Users/claer14/OneDrive - University of Cambridge/III. Research/Data/')
# hospital <- 'Chamba'
#
# # Load data
# data <- read.csv('1_Original_Data/108_HISP/20250326_Chamba_Data.csv')
# names(data) <- 'raw_column'
# 
# # Splitting the columns & organising
# df <- data %>%
#   separate(raw_column, into = c('ward_name', 'patient_id', 'age', 'gender', 'category', 'date', 'visit_type'), 
#            sep = "tt", fill = "right") %>%
#   mutate(hospital_name = hospital) 
# 
# # Reorganize col order & set type
# df <- df[, c('hospital_name', 'date', 'patient_id', 'ward_name', 'age', 'gender', 'category', 'visit_type')]
# # table(sort(df$age))  # Counts unique values
# # df$age <- as.numeric(gsub("[^0-9]", "", df$age))
# df$date <- as.Date(df$date)
# 
# # Export 
# write.csv(df, file = "2_Active_Data/201_heal/HISP/ind_adm4_HISP_2015-2023_patient_visits_chamba.csv", row.names = FALSE)

################################################################################

data <- read.csv('ind_adm4_HISP_2015-2023_patient_visits_chamba.csv')




### NB: Start 2015-02-14; End 2023-01-19

# How many missing data?
### Remove rows with NA dates

# How many missing time steps? (how continuous)

# Categories age groups

# Aggregate DAILY counts as: totals; ward specific; age specific; gender specific (+ ward & age/gender specific)

# Question: 
### What do the visit_type and category values mean? 
### Why are there unviable numbers in age (e.g. -1121, 1951, REG-INITIAL, REG-REVISIT)


# Future query: Get date with the time stamp