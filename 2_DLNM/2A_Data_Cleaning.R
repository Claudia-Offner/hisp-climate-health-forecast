### ------------------------------------------------------------------------ ### 
###  DNLM - Data cleaning and tests                                          ###
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list = ls()); graphics.off(); gc()
setwd('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/')
source('2_DLNM/HISP_DLNM - Dependencies.R')

# Set data path
setwd(path_data)

#### MAIN CODE ####

df <- get(load('2_Active_Data/203_comb/ind_5km_BigThaw-HISP-ALL_2015-2023_daily_hosp_wide.Rdata'))
rm(df_wide)

# 1. Data formatting ####

# Add lat/lon back into df
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

# Filter out Bilaspur (incomplete series)
df <- df %>%   filter(!location %in% c("Bilaspur"))

# Filter out COVID years
startdate <- as.Date('2015-01-01')
enddate <- as.Date('2020-02-14')  
df <- df[df$date >= startdate & df$date <= enddate,]

# !!! Remove Sundays for now !!!
# df <- df %>% filter(weekday != 'Sun')

# Select variables of interest
names(df) <- gsub("rh_mean", "rh", names(df))
names(df) <- gsub("sh_mean", "sh", names(df))
df <- df %>%
  group_by(location) %>%
  arrange(date, .by_group = TRUE) %>% 
  # Compute time within each location
  mutate(doy = yday(date), time = row_number()) %>%
  ungroup() %>%
  # Select cols of interest
  dplyr::select(location, lat, lon, time, date, year, season, month, week, 
                weekday, doy, closed, holiday, elev, popu, all_of(h_cols), 
                all_of(w_cols))

# 2. Data transformations ####

# PRECIPITATION: Create 7-day rolling sum & log transform
df <- df %>%
  group_by(location) %>%
  mutate(tp_roll=rollapply(tp,   width = 7, FUN = sum,  align = "right", fill = NA)) %>%
  ungroup()
df$tp_roll_log <- log(df$tp_roll+10)

# Check
p1 <- ggplot(df, aes(x = tp)) + 
  geom_histogram(binwidth = 25, fill = "grey", color = "black") +
  labs(title = "Raw") + 
  theme_minimal()

p2 <- ggplot(df, aes(x = tp_roll)) +
  geom_histogram(binwidth = 75, fill = "grey", color = "black") +
  labs(title = "Rolling (7-day)") + 
  theme_minimal()

p3 <- ggplot(df, aes(x = tp_roll_log)) +
  geom_histogram(binwidth = 1, fill = "grey", color = "black") +
  labs(title = "Rolling + 10 with log") + 
  theme_minimal()

tp_dist <- wrap_plots(list(p1, p2, p3), nrow = 1) + plot_annotation(title = "Distribution of precipitation Variables")

df$tp <- df$tp_roll_log
df$tp_roll <- NULL
df$tp_roll_log <- NULL

# LAGS: Add lags for outcome & w_col
lag_cols <- c('total_patients', w_cols)
for (col in lag_cols) {
  for (i in 1:30) {
    df <- df %>%
      group_by(location) %>%
      mutate("{col}_lag{i}" := lag(.data[[col]], i)) %>%
      ungroup()
  }
}

# Forward-fill ONLY the lagged variables
lagged_vars <- grep("_lag\\d+$", names(df), value = TRUE)

df <- df %>%
  group_by(location) %>%
  fill(all_of(lagged_vars), .direction = "down") %>%
  ungroup()

# EXPORT ####
ggsave(paste0(path_report, '- figures/desc_tp_dist.png'), tp_dist, width=8, height=6, dpi = 300, bg='white')
save(df, file = paste0(path_project, 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))
save(df, file = paste0(path_report, 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_DLNM.RData'))
