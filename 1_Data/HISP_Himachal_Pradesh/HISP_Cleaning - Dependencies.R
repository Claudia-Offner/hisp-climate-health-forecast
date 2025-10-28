### ------------------------------------------------------------------------ ### 
### HISP Data cleaning - Dependencies                                        ###
### ------------------------------------------------------------------------ ### 
# NB COLOURS: https://r-charts.com/colors/

# Turn off scientific notation
options(scipen=999)

#### PATHS ####
user <- 'claer14'
path_data <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/')
path_git <- paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/')
path_project <- paste0(path_git, '/1_Data/HISP_Himachal_Pradesh/')
path_report <- paste0(path_git, path_project, '- Reporting/')


#### PACKAGES ####

# List packages
packages <- c('tidyverse', 'terra', 'duckdb', 'tidyselect', 'sf',
              'ggplot2', 'data.table', 'tools', 'moments', 'grid',
              'patchwork', 'RColorBrewer', 'reshape2', 'zoo', 'ggnewscale',
              'skimr', 'DT', 'forecast', 'purrr', 'scales', 'xts', 'ggthemes',
              'elevatr', 'rnaturalearth', 'rnaturalearthdata')



# Install if needed
installed <- packages %in% rownames(installed.packages())
if(any(!installed)) {
  install.packages(packages[!installed])
}

# Load packages
lapply(packages, library, character.only = TRUE)

#### FUNCTIONS ####

# Function to replace exact matches from a list
replace_exact_matches <- function(df, dict) {
  df[] <- lapply(df, function(col) {
    if (is.factor(col)) col <- as.character(col)
    if (is.character(col)) {
      col <- ifelse(col %in% names(dict), dict[col], col)
    }
    col
  })
  df
}


#### DEPENDENCIES ####

mypalette <- c('#1D98EB','#D74B3C','#A6B727','#FEC306','#EC800A','#E569BF','#8968CD','#52C5B0FF', '#910030FF', '#127804FF',
               '#A1C8F2','#F1A6A6','#C2D77D','#F6D974','#F9B34D','#F1C2D1','#B8A8D9','#52C5B080', '#91003080', '#12780480') 

w_cols <- c('t2m', 'mn2t','mx2t','tp','fg10', 'rh_mean', 'sh_mean', 'pm25')
h_cols <- c('total_patients',
            'F','M','under_5','age5_14','age15_24','age25_64','over_65',
            'general', 'emergency', 'ophthalmology', 'ent', 'paediatrics',
            'gynaecology', 'orthopedics', 'dermatology', 'surgery',
            'psychiatry', 'pulmonary', 'disability', 'dentistry', 'oncology',
            'ncd', 'antiretroviral', 'substance_counselling', 'physiotherapy')

location_levels <- c('Bilaspur','Chamba', 'Dharamshala', 'Hamirpur', 'Kullu', 'Kyelang', 
                     'Mandi', 'Nahan', 'Reckong Peo', 'Shimla', 'Solan', 'Una')
season_levels   <- c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)', 'Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)')
month_levels    <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
weekday_levels  <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')
age_levels      <- c('under_5', 'age5_14', 'age15_24', 'age25_64',  'over_65','Unknown')
gender_levels   <- c('F', 'M', 'O')
ward_levels     <- c('general', 'emergency', 'ophthalmology', 'ent', 'paediatrics',
                     'gynaecology', 'orthopedics', 'dermatology', 'surgery',
                     'psychiatry', 'pulmonary', 'disability', 'dentistry', 'oncology',
                     'ncd', 'mch', 'antiretroviral', 'substance_counselling', 'physiotherapy')

# Dictionary of names for figures and tables
nm <- c(season='Season', closed='Closed', holiday='Holiday', total_patients='Total patients', 
        # POPULATION GROUPS
        M='Men', F='Women',
        under_5='Under 5', age5_14='Ages 5-14', age15_24='Ages15-24', age25_64='Ages 25-64', over_65='Ages 65+', 
        # WARDS
        general='General', emergency='Emergency', ophthalmology='Ophthalmology', 
        paediatrics='Paediatric', orthopedics='Orthopedics', ent='ENT', 
        gynaecology='Gynaecology', dermatology='Dermatology', surgery='Surgery', 
        dentistry='Dentistry', psychiatry='Psychiatry', oncology='Oncology', 
        pulmonary='Pulmonary', disability='Disability', ncd='NCD', mch='MCH',
        antiretroviral='HIV Treatment', substance_counselling='Substance counselling',
        physiotherapy='Physiotherapy', 
        # CONTINUOUS WEATHER
        t2m='Mean temp. (\u00B0C)', mx2t='Max temp. (\u00B0C)', mn2t='Min temp. (\u00B0C)',
        rh_mean='Relative humidity (%)', tp='Rainfall (mm)', fg10='Windspeed (km/hr)', pm25='Air pollution (µg/m³)',
        # IMD DEFINITIONS
        mx2t_IMD='IMD: Mx temp. >30\u00B0C', mn2t_IMD='IMD: Min temp. <0\u00B0C', 
        tp_IMD='IMD: Rainfall >15mm', fg10_IMD='IMD: Wind speed >40km/h', pm25_IMD='IMD: PM2.5 >60 µg/m³', 
        rh_low_IMD='IMD: Humidity <30%', rh_hig_IMD='IMD: Humidity >60%',
        # CLIM. DEFINITIONS
        mx2t_p90='Clim: Max temp. >90th', mn2t_p10='Clim: Min temp. <10th', 
        tp_p90='Clim: Rainfall >90th', fg10_p90='Clim: Wind speed >90th', pm25_p90='Clim: PM2.5 >90th', 
        rh_mean_p10 ='Clim: Humidity <10th', rh_mean_p90='Clim: Humidity >90th') 

