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
startdate <- as.Date('2015-01-01')
enddate <- as.Date('2020-02-14')  
df <- df_wide[df$date >= startdate & df$date <= enddate,]

# Select variables of interest
df <- df %>%
  mutate(
    doy = yday(date),
    time = row_number(),
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

# Filter out Bilaspur
df <- df %>% filter(location != 'Bilaspur')

# 7 day rolling sum
df <- df %>%
  group_by(location) %>%
  mutate(tp=rollapply(tp,   width = 7, FUN = sum,  align = "right", fill = NA)) %>%
  ungroup()

# Log transform tp
df$tp <- log(df$tp + 1e-6)

# Diurnal range
df$dt <- df$mx2t-df$mn2t



# WET BULB TEMP
# https://github.com/open-meteo/open-meteo/discussions/1091
# library(psychrolib)
pressure_from_elevation <- function(elevation_m) {
  101325 * (1 - 2.25577e-5 * elevation_m) ^ 5.25588
}
df$pre <- pressure_from_elevation(df$elev)
# (wb <- t2m*atan(0.151977*sqrt(rh+8.313659)) + atan(t2m+rh)-atan(rh-1.67631) + 0.00391838*rh^1.5*atan(0.023101*rh) - 4.686035)
library(psychrolib)
SetUnitSystem("SI")
# Safely compute wet bulb temperature
GetTw_safe <- Vectorize(function(T, RH, P) {
  if (is.na(T) | is.na(RH) | is.na(P)) return(NA)
  GetTWetBulbFromRelHum(T, RH, P)
})
df$wb <- GetTw_safe(df$t2m, df$rh / 100, df$pre)



df <- df %>% group_by(location) %>% arrange(date) %>% mutate(patients_lag = lag(total_patients, 1))


w_cols <- c(setdiff(w_cols, "rh_mean"), "rh", "dt", "wb")
sub <- df %>% ungroup() %>% dplyr::select(total_patients, all_of(w_cols[!grepl('_p|IMD', w_cols)]))
cor_matrix <- cor(sub, use = 'complete.obs', method = 'spearman')
melted_corr <- reshape2::melt(cor_matrix)
(cor_matrix <- ggplot(data = melted_corr, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = 0, limit = c(-1,1)) +
  theme_minimal() +
  scale_x_discrete(labels = function(x) {
    # If x is in nm, replace with nm value, else keep original
    sapply(x, function(lbl) {
      if (lbl %in% names(nm)) {
        nm[[lbl]]
      } else {
        lbl
      }})}) +
  scale_y_discrete(labels = function(x) {
    # If x is in nm, replace with nm value, else keep original
    sapply(x, function(lbl) {
      if (lbl %in% names(nm)) {
        nm[[lbl]]
      } else {
        lbl
      }})}) +
  geom_text(aes(label=round(value,2)), color='black', size=4) +
  labs(title = 'Spearman Correlation Matrix', x='', y='') +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # rotate x-axis labels
    axis.text.y = element_text(angle = 0)                # rotate y-axis labels if needed
  ))

# max_lag <- 14
# for (lag in 0:max_lag) {
#   df[[paste0('t2m_l', lag)]] <- dplyr::lag(df$t2m, lag)
#   df[[paste0('mn2t_l', lag)]]  <- dplyr::lag(df$mn2t, lag)
#   df[[paste0('mx2t_l', lag)]]  <- dplyr::lag(df$mx2t, lag)
#   df[[paste0('tp_l', lag)]]  <- dplyr::lag(df$tp, lag)
#   df[[paste0('fg10_l', lag)]]  <- dplyr::lag(df$fg10, lag)
#   df[[paste0('rh_l', lag)]]  <- dplyr::lag(df$rh, lag)
#   df[[paste0('pm25_l', lag)]] <- dplyr::lag(df$pm25, lag)
#   df[[paste0('dt_l', lag)]] <- dplyr::lag(df$dt, lag)
#   df[[paste0('wb_l', lag)]] <- dplyr::lag(df$wb, lag)
# }
# 
# 

# Elevation filter
# df <- df %>%
#   # filter(elev > 1000)
#   filter(elev < 1000)

# Seasonal filter
# df <- df %>%
#   filter(season %in% c('Winter (Jan-Feb)', 'Pre-Monsoon (Mar-May)'))
#   filter(season %in% c('Monsoon (Jun-Sep)', 'Post-Monsoon (Oct-Dec)'))

# # Categories for exp
# features <- names(df)[grepl('t2m|tp|rh|pm25|fg10', names(df))]
# df <- df %>%
#   mutate(across(
#     all_of(features),
#     ~ cut(.x, breaks = 20, labels = FALSE),  # divide into 20 bins
#     .names = "{.col}_cat"                     # new column names
#   ))

# # Test 3-way interactions for all exposures. 
# # NB: wind and humidity are not considered key exposures but potential exacerbaters
# m0 <- lm(total_patients ~ t2m + tp + fg10 + rh + pm25, data = df)
# m1 <- lm(total_patients ~ (t2m + tp + pm25 + fg10 + rh)^3, data = df)
# anova(m0, m1)
# summary(m0)

# Save as Rdata file
save(df, file = paste0(path_project, 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))
save(df, file = paste0(path_report, 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))
