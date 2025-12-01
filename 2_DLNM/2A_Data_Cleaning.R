### ------------------------------------------------------------------------ ### 
###  DNLM - Data cleaning and tests                                          ###
### ------------------------------------------------------------------------ ### 

setwd(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/'))
source('2_DLNM/HISP_DLNM - Dependencies.R')

# Set data path
setwd(path_data)

#### MAIN CODE ####

load('2_Active_Data/203_comb/ind_5km_BigThaw-HISP-ALL_2015-2023_daily_hosp_wide.Rdata')
df <- df_wide

# Remove COVID years for now:
startdate <- as.Date('2015-02-14')
enddate <- as.Date('2020-02-14')  
df <- df_wide[df$date >= startdate & df$date <= enddate,]

# Select variables of interest
df <- df %>%
  mutate(
    doy = yday(date),
    time = row_number(),
    tp = exp(tp) # back transform tp
  ) %>%
  dplyr::select(
    location, time, date, year, season, month, week, weekday, doy, closed, holiday, elev, popu,
    all_of(h_cols),
    all_of(w_cols)
  )

# Rename rh
df$rh <- df$rh_mean
df$rh_mean <- NULL
df$sh_mean <- NULL

# Create a lag variable for outcome


# Add lat lon back into df

loc_points <- data.frame(
  location = c("Chamba", "Dharamshala", "Solan", "Hamirpur", "Kullu",
               "Reckong Peo", "Kyelang", "Nahan", "Bilaspur", "Mandi",
               "Shimla", "Una"),
  lon = c(76.127017, 76.321482, 77.099506, 76.522135, 77.107945,
          78.267181, 77.023559, 77.294546, 76.761165, 76.928977,
          77.170999, 76.276997),
  lat = c(32.554798, 32.216995, 30.906408, 31.684233, 31.956550,
          31.539809, 32.571883, 30.564781, 31.341244, 31.712559,
          31.104221, 31.479129)
)

df <- merge(df, loc_points, by='location')

# Save as Rdata file
save(df, file = paste0(path_project, 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))
save(df, file = paste0(path_report, 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))
