### ------------------------------------------------------------------------ ### 
###  Ensemble - Data Cleaning (HISP)
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list = ls()); graphics.off(); gc()

# Set data path
setwd('C:/Users/claer14/Documents/GitHub/hisp-climate-health-forecast/')

#### DEPENDENCIES ####
source('1_Data/Data - Dependencies.R')
setwd(path_data)
file_export <- 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide'


#### MAIN CODE ####

# Load data
load('2_Active_Data/203_comb/ind_5km_BigThaw-HISP-ALL_2015-2023_daily_hosp_wide.Rdata')


# 1. DLNM analysis ####

controls <- c('elev', 'popu', "hosp_n", "hosp_teach", 'mcv1', "dtp1", "dtp3", "hsize", "anc", "insurance", "malnutrition", "pnc", "sba", "tt", "decision", "hwealth", "meducation")

# Add lat/lon back into df
loc_points <- data.frame(
    location = c('Chamba', 'Dharamshala', 'Solan',
                 'Hamirpur', 'Kullu', 'Kangra',
                 'Reckong Peo', 'Kyelang', 'Nahan',
                 'Bilaspur', 'Mandi', 'Shimla', 'Una'),
    lat = c(76.127017, 76.321482, 77.099506,
          76.522135, 77.107945, 76.274797,
          78.267181, 77.023559, 77.294546,
          76.761165, 76.928977, 77.170999, 76.276997),
    lon = c(32.554798, 32.216995, 30.906408,
          31.684233, 31.956550, 32.103220,
          31.539809, 32.571883, 30.564781,
          31.341244, 31.712559, 31.104221, 31.479129)
  )

df_dlnm <- merge(df_wide, loc_points, by='location')

# Filter out COVID years
startdate <- as.Date('2015-01-01')
enddate <- as.Date('2020-02-14')  
df_dlnm <- df_dlnm[df_dlnm$date >= startdate & df_dlnm$date <= enddate,]

# Select variables of interest
names(df_dlnm) <- gsub("rh_mean", "rh", names(df_dlnm))
names(df_dlnm) <- gsub("sh_mean", "sh", names(df_dlnm))
w_cols[w_cols %in% c("rh_mean", "sh_mean")] <- c("rh", "sh")

df_dlnm <- df_dlnm %>%
  group_by(location) %>%
  arrange(date, .by_group = TRUE) %>% 
  # Compute time within each location
  mutate(doy = yday(date), time = row_number()) %>%
  ungroup() %>%
  # Select cols of interest
  dplyr::select(location, lat, lon, time, date, year, season, month, week, 
                weekday, doy, closed, holiday, all_of(controls), all_of(h_cols), 
                all_of(w_cols))

# Create weekly df
df_weekly <- df_dlnm %>%
  group_by(location, year, week) %>%
  summarise(
    # --- identifiers / constants (take first value in the week) ---
    lat        = first(lat),
    lon        = first(lon),
    elev       = first(elev),
    week_start = min(date),
    season     = first(season),   # a week could span two seasons/months at the boundary;
    month      = first(month),    # `first()` just takes whichever comes first in the week
    
    # --- daily 0/1 indicators -> count of days in the week ---
    n_days     = n(),
    n_closed   = sum(closed, na.rm = TRUE),
    n_holiday  = sum(holiday, na.rm = TRUE),
    
    # --- slow-moving covariates (essentially constant within location/time) -> mean ---
    across(c(popu, mcv1, dtp1, dtp3, hsize, anc, insurance, malnutrition,
             pnc, sba, tt, decision, hwealth, meducation),
           ~ mean(.x, na.rm = TRUE)),
    
    # --- facility indicators -> first (assumed constant per location) ---
    hosp_n     = first(hosp_n),
    hosp_teach = first(hosp_teach),
    
    # --- patient visit counts -> sum across the week ---
    across(c(total_patients, F, M, under_5, age5_14, age15_24, age25_64,
             over_65, general, emergency, ophthalmology, ent, paediatrics,
             gynaecology, orthopaedics, dermatology, surgery, psychiatry,
             pulmonary, disability, dentistry, oncology, ncd,
             antiretroviral, substance_counselling, physiotherapy),
           ~ sum(.x, na.rm = TRUE)),
    
    # --- climate variables -> mean for temperature/humidity/wind types ---
    across(c(t2m, mn2t, mx2t, wb, dt, fg10, rh, sh, pm25),
           ~ mean(.x, na.rm = TRUE)),
    
    # --- precipitation -> sum (it's an accumulated quantity, not an instantaneous one) ---
    tp = sum(tp, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(location, week_start) %>%
  group_by(location) %>%
  mutate(time = row_number()) %>%
  ungroup()

# Create cumulative measures (windows 1-14) for w_cols
windows <- 1:14

df_dlnm <- df_dlnm %>%
  group_by(location) %>%
  mutate(
    across(
      all_of(w_cols),
      .fns = setNames(
        lapply(windows, function(w) {
          function(x) {
            rollapply(
              x,
              width = w,
              FUN = if (cur_column() %in% c("tp", "pm25")) sum else mean,
              align = "right",
              fill = NA
            )
          }
        }),
        paste0("roll", windows)
      )
    )
  ) %>%
  ungroup()

# LAGS: Add lags for outcome & w_col
lag_cols <- c(h_cols, w_cols)
for (col in lag_cols) {
  for (i in windows) {
    df_dlnm <- df_dlnm %>%
      group_by(location) %>%
      mutate("{col}_lag{i}" := lag(.data[[col]], i)) %>%
      ungroup()
  }
}
# LAGS: Add lags for outcome & w_col
for (col in lag_cols) {
  for (i in windows) {
    df_weekly <- df_weekly %>%
      group_by(location) %>%
      mutate("{col}_lag{i}" := lag(.data[[col]], i)) %>%
      ungroup()
  }
}


# Log transform tp variables
df_dlnm$tp_log <- log(df_dlnm$tp+2)
df_dlnm <- df_dlnm %>%
  mutate(
    across(
      matches("^tp.*roll"),
      ~ log(.x + 10),
      .names = "{.col}_log"
    ),
    across(
      matches("^tp.*lag"),
      ~ log(.x + 2),
      .names = "{.col}_log"
    )
  )


# Log transform tp variables
df_weekly$tp_log <- log(df_weekly$tp+2)
df_weekly <- df_weekly %>%
  mutate(
    across(
      matches("^tp.*roll"),
      ~ log(.x + 10),
      .names = "{.col}_log"
    ),
    across(
      matches("^tp.*lag"),
      ~ log(.x + 2),
      .names = "{.col}_log"
    )
  )

# Forward-fill ONLY the lagged variables
lagged_vars <- grep("_lag\\d+$", names(df_dlnm), value = TRUE)
lagged_vars <- grep("_lag\\d+$", names(df_weekly), value = TRUE)

df_dlnm <- df_dlnm %>%
  group_by(location) %>%
  fill(all_of(lagged_vars), .direction = "down") %>%
  ungroup()
df_weekly <- df_weekly %>%
  group_by(location) %>%
  fill(all_of(lagged_vars), .direction = "down") %>%
  ungroup()


# # Check
# p1 <- ggplot(df_dlnm, aes(x = tp)) + 
#   geom_histogram(binwidth = 25, fill = 'gray50', color = "black") +
#   labs(title = "Raw") + 
#   theme_minimal() +
#   theme(axis.title.x=element_blank())
# 
# p2 <- ggplot(df_dlnm, aes(x = tp_roll)) +
#   geom_histogram(binwidth = 75, fill = 'gray50', color = "black") +
#   labs(title = "Rolling (7-day)") +
#   theme_minimal() +
#   theme(axis.title.x=element_blank())
# 
# p3 <- ggplot(df_dlnm, aes(x = tp_log)) +
#   geom_histogram(binwidth = 1, fill = 'gray50', color = "black") +
#   labs(title = "Log + 1 (raw)") +
#   theme_minimal() +
#   theme(axis.title.x=element_blank())
# 
# p4 <- ggplot(df_dlnm, aes(x = tp_roll_log)) +
#   geom_histogram(binwidth =1, fill = 'gray50', color = "black") +
#   labs(title = "Log rolling + 1 (7-day)") + 
#   theme_minimal()+
#   theme(axis.title.x=element_blank())
# 
# tp_dist <- wrap_plots(list(p1, p2, p3, p4), nrow = 1) +
#   plot_annotation(
#     title = "Distribution of precipitation variables",
#     theme = theme(
#       plot.title = element_text(hjust = 0.5, face='bold')
#     )
#   )
# 
# # Select tp exposure
# df_dlnm$tp <- df_dlnm$tp_log
# df_dlnm$tp_roll <- df_dlnm$tp_roll_log

# 2. Ensemble analysis ####

# Remove Bilaspur
df_ensm <- df_wide %>% filter(location != 'Bilaspur')

# Identify locations
locations <- as.character(df_ensm %>% distinct(location, elev) %>% filter(!is.na(elev)) %>% arrange(elev) %>% pull(location))

# Select year of interest
startdate <- as.Date('2015-02-14')
enddate <- as.Date('2020-02-13')
df_ensm <- df_ensm [df_ensm $date >= startdate & df_ensm $date <= enddate,]

# Create 7-day rolling sum & log transform
df_ensm <- df_ensm %>%
  group_by(location) %>%
  mutate(tp_roll=rollapply(tp,   width = 7, FUN = sum,  align = "right", fill = NA)) %>%
  ungroup()
df_ensm $tp <- log(df_ensm $tp_roll+10)

# Create percentile look up
variables <- c('mx2t', 'tp')
cen_tables <- list()
datalist <- lapply(locations, function(region) df_ensm [df_ensm $location==region,])
for (v in variables) {
  cen <- sapply(datalist,function(x) {
    quantile(x[[v]],c(0:100)/100,na.rm=T)
  })
  if (v=='tp'){
    cen <- exp(cen)-10
  }
  colnames(cen) <- locations
  # cen <- cbind(cen, Average = rowMeans(cen))
  cen_tables[[v]] <- round(cen, 4)
}

# Set weather thresholds based on DLNM outputs
mx2t_20th <- as.data.frame(cen_tables$mx2t) %>%
  mutate(percentile = rownames(.)) %>%
  pivot_longer(cols = -percentile, names_to = "location", values_to = "thresh_mx2t") %>% 
  filter(percentile == "20%") %>% dplyr::select(location, thresh_mx2t)

tp_80th <- as.data.frame(cen_tables$tp) %>%
  mutate(percentile = rownames(.)) %>%
  pivot_longer(cols = -percentile, names_to = "location", values_to = "thresh_tp") %>% 
  filter(percentile == "80%") %>% dplyr::select(location, thresh_tp)

df_ensm <- df_ensm %>% 
  left_join(mx2t_20th , by = "location")  %>%
  left_join(tp_80th , by = "location")  %>%
  mutate(cold20_lag0 = ifelse(mx2t < thresh_mx2t, 1, 0),
         prec80_lag0 = ifelse(tp > thresh_tp, 1, 0)) %>%
  group_by(location) %>%
  mutate(
    cold20_lag3 = lag(cold20_lag0, 3),
    prec80_lag10 = lag(prec80_lag0, 10)
  ) %>%
  ungroup() %>%
  dplyr::select(-c(tp_roll, thresh_tp, thresh_mx2t), -matches("lag(2|[4-9]|1[1-4])"))

names(df_ensm )
table(df_ensm $cold20_lag0)
table(df_ensm $prec80_lag0)

#### EXPORT ####
# ggsave(paste0(path_report, 'I. desc/fig_4_tp_dist.png'), tp_dist, width=8, height=4, dpi = 300, bg='white')
save(df_dlnm, df_weekly, file = paste0(path_git, '2_DLNM/', file_export, '_dlnm.RData'))
save(df_ensm , file = paste0(path_git, '3_Ensemble/', file_export, '_ensemble.RData'))

