### ------------------------------------------------------------------------ ### 
###  DNLM - Tests and cleaning                                               ###
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
library(dplyr); library(tidyr); library(lubridate)

# Functions

# Variables
weekday_levels <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
month_levels <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
season_levels <- c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)', 'Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)')

#### MAIN CODE ####

# Load data
# df <- read.csv('203_comb/ind_10km_HISP-HP_2015-2017_hosp-bigthaw_chamba_wide_model.csv', check.names=FALSE, fileEncoding='UTF-8')
df <- read.csv('203_comb/ind_30km_HISP-HP_2015-2023_hosp-era5_chamba_wide_model.csv', check.names=FALSE, fileEncoding='UTF-8')

# Set date type and factors 
df$date <- as.Date(df$date, format='%Y-%m-%d')  # Adjust format as needed
df <- df %>%
  mutate(weekday=factor(weekday, levels=weekday_levels),
         month=factor(month, levels=month_levels), 
         season=factor(season, levels=season_levels))

# Select variables of interest
df <- df %>%
  mutate(doy = yday(date),
         time = row_number()) %>%
  select(
    location, time, date, year, season, month, week, weekday, doy, closed, holiday, total_patients, 
    F, M, under_5, age5_14, age15_24, age25_64, over_65,
    general, emergency, ophthalmology, pediatrics, orthopedics, ent, gynecology, dermatology, surgery, dentistry, psychiatry, oncology, 
    t2m, mn2t, mx2t, tp, fg10 #,
    # mx2t_IMD, mx2t_IMD_b, mn2t_IMD, mn2t_IMD_b, tp_IMD, tp_IMD_b, fg10_IMD, fg10_IMD_b, 
    # t2m_p10, t2m_p90, mn2t_p10, mn2t_p90, mx2t_p10, mx2t_p90, tp_p10, tp_p90, fg10_p10, fg10_p90
  )

df$tp <- exp(df$tp)
summary(df)

# Test for interactions (p<0.05 indicates the interaction improves model)
test_interaction <- function(data, response, var1, var2, var3 = NULL, family = gaussian) {
  
  if (is.null(var3)) {
    # --- 2-way interaction ---
    f_main  <- as.formula(paste(response, "~", var1, "+", var2))
    f_inter <- as.formula(paste(response, "~", var1, "*", var2))
  } else {
    # --- 3-way interaction ---
    f_main  <- as.formula(paste(response, "~", var1, "+", var2, "+", var3, 
                                "+", var1, ":", var2, "+", var1, ":", var3, "+", var2, ":", var3))
    f_inter <- as.formula(paste(response, "~", var1, "*", var2, "*", var3))
  }
  
  # fit models
  m0 <- glm(f_main, data = data, family = family)
  m1 <- glm(f_inter, data = data, family = family)
  
  # likelihood ratio test
  result <- anova(m0, m1, test = "Chisq")
  return(result)
}

# Test for weather interactions
test_interaction(df, "total_patients", "tp", "fg10", family = gaussian) # p=0.01472 *
test_interaction(df, "total_patients", "t2m", "tp", family = gaussian) # p=0.15
test_interaction(df, "total_patients", "t2m", "fg10", family = gaussian) # p=0.0005258 ***
test_interaction(df, "total_patients", "mx2t", "tp", family = gaussian) # p=0.18
test_interaction(df, "total_patients", "mx2t", "fg10", family = gaussian) # p=0.00252 **
test_interaction(df, "total_patients", "mn2t", "tp", family = gaussian) # p=0.0756 .
test_interaction(df, "total_patients", "mn2t", "fg10", family = gaussian) # p=0.06709 .
# Check 3-way interaction (NULL)
test_interaction(df, "total_patients", "t2m", "tp", "fg10", family = gaussian) # p= 0.4012
# = Potential interactions: tp*fg10; t2m*fg10; mx2t*fg10


# Test for socio-dem interactions

test_interaction(df, "total_patients", 't2m', "under_5", family = gaussian) # p=0.003974 **
test_interaction(df, "total_patients", 't2m', "age5_14", family = gaussian) # p=0.003974 **
test_interaction(df, "total_patients", 't2m', "age15_24", family = gaussian) # p=0.003974 **
test_interaction(df, "total_patients", 't2m', "age25_64", family = gaussian) # p=0.003974 **
test_interaction(df, "total_patients", 't2m', "over_65", family = gaussian) # p=0.003974 **

test_interaction(df, "total_patients", 't2m', "F", family = gaussian) # p=0.003974 **
test_interaction(df, "total_patients", 't2m', "M", family = gaussian) # p=0.003974 **

test_interaction(df, "total_patients", 't2m', "under_5", "M", family = gaussian) # p=0.003974 **
test_interaction(df, "total_patients", 't2m', "under_5", "F", family = gaussian) # p=0.0000000000007516 ***
test_interaction(df, "total_patients", 't2m', "age5_14", "M", family = gaussian) # p=0.08904 .
test_interaction(df, "total_patients", 't2m', "age5_14", "F", family = gaussian) # p=0.00000000000643 ***
test_interaction(df, "total_patients", 't2m', "age15_24", "M", family = gaussian) # p=0.00000000000000022 ***
test_interaction(df, "total_patients", 't2m', "age15_24", "F", family = gaussian) # p=0.00000000009531 ***
test_interaction(df, "total_patients", 't2m', "age25_64", "M", family = gaussian) # p=0.00000000000003963 ***
test_interaction(df, "total_patients", 't2m', "age25_64", "F", family = gaussian) # p=0.00956 **
test_interaction(df, "total_patients", 't2m', "over_65", "M", family = gaussian) # p=0.00003949 ***
test_interaction(df, "total_patients", 't2m', "over_65", "F", family = gaussian) # p=0.0000000000000005428 ***

# Maybe interaction between age and gender across weather variables

# Save as Rdata file
save(df, file = paste0(path1, "DLNM/ind_30km_HISP-HP_2015-2023_hosp-era5_chamba_wide_model.RData"))

