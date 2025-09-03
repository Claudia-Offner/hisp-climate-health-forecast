### ------------------------------------------------------------------------ ### 
###  DNLM - Descriptive Analysis                                             ###
### ------------------------------------------------------------------------ ### 

#### DEPENDENCIES ####

# Turn off scientific notation
options(scipen=999)

# Path
user <- 'claer14'
path1 <- paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/')
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

# 4. Interaction Tests ####

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

