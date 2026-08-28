### ------------------------------------------------------------------------ ### 
###  DNLM – Tables & Figures (HISP)
### ------------------------------------------------------------------------ ### 

rm(list = ls()); graphics.off(); gc()
setwd('C:/Users/claer14/Documents/GitHub/hisp-climate-health-forecast/')

#### DEPENDENCIES ####

source('2_DLNM/2_Utilities.R')
path_data <- paste0('C:/Users/claer14/OneDrive - University of Cambridge/III. Research/Data/')
setwd(path_data)

# Packages
library(readxl); library(tidyverse); library(boot); library(deming); 
library(forecast); library(raster); library(sf)


# Functions
get_data <- function(data){
  
  data <- data[-1, ]
  data %>%
    # Fix col names
    rename_with(~ if_else(grepl("\\d{4}", .x), .x, paste0("19", .x))) %>%
    rename(date_raw = 1) %>%
    # Fix month and day cols
    mutate(
      date_raw = as.numeric(date_raw),
      is_month_start = date_raw > 1000,
      month_date = if_else(is_month_start, as.Date(date_raw, origin = "1899-12-30"), as.Date(NA)),
      month_date = zoo::na.locf(month_date, na.rm = FALSE),
      month = format(month_date, "%b"),
      day = if_else(is_month_start, 1, date_raw)
    ) %>%
    dplyr::select(-date_raw, -is_month_start, -month_date) %>%
    relocate(month, day) %>%
    # Pivot year data to long format
    mutate(across(-c(month, day), as.character)) %>%
    pivot_longer(
      cols = -c(month, day),
      names_to = "year",
      values_to = var
    ) %>%
    mutate(
      year = as.integer(year),
      across(all_of(var), as.numeric)
    ) %>%
    # Create date col
    mutate(
      date = as.Date(paste0(day, "-", month, "-", year), format = "%d-%b-%Y"),
      year = format(date, "%Y"),
      month = format(date, "%b"),
      location = loc
    ) %>%
    dplyr::select(date, year, month, location, all_of(var)) %>%
    arrange(date, location) %>%
    filter(!is.na(date))
  
}


####  MAIN CODE ####
# 
# loc_points <- st_as_sf(
#   data.frame(
#     location = c("Bhuntar", "Larji", "Manali", "Pandoh"),
#     lat = c(31.88318090301291, 31.725126418581784, 32.24346317173383, 31.667516894291232),
#     lon = c(77.14625339534699, 77.21866419672843, 77.19123937228814, 77.05351788783875)
#   ),
#   coords = c("lon", "lat"),
#   crs = 4326
# )
# loc_points_sf <- st_as_sf(loc_points, coords = c("lon", "lat"), crs = 4326) |>
#   st_transform(crs = 3857) |>
#   st_buffer(dist = 5000) |>
#   st_transform(crs = 4326)



# # A. Clean IMD data #####
# folder <- '1_Original_Data/124_IMD/rainfall_Temp/'
# files <- list.files(folder, pattern='xls')
# 
# mn2t_df <- data.frame()
# mx2t_df <- data.frame()
# 
# for(f in files) {
# 
#   # Get info
#   var <- sapply(strsplit(f, "_"), `[`, 1)
#   loc <- sapply(strsplit(f, "_"), `[`, 2)
# 
#   # Load data
#   if(loc=='Pandoh'){ s <- 1} else {s <- 2}
# 
#   data <- read_excel(paste0(folder, f), skip = s)
#   df <- get_data(data)
# 
#   if (var=='mn2t'){
#     mn2t_df <- rbind(mn2t_df, df)
#   } else {
#     mx2t_df <- rbind(mx2t_df, df)
#   }
# 
# }
# 
# imd_df <- mn2t_df %>% merge(mx2t_df, by=c('location', 'date', 'year', 'month'))
# 
# 
# # Check NA values over years
# table(imd_df$location[!is.na(imd_df$mx2t)], imd_df$year[!is.na(imd_df$mx2t)])
# table(imd_df$location[!is.na(imd_df$mn2t)], imd_df$year[!is.na(imd_df$mn2t)])
# # Check 2018 values
# test <- imd_df %>% filter(year== '2018')
# table(test$location[!is.na(test$mx2t)], test$month[!is.na(test$mx2t)])
# 
# # Set date range
# startdate <- as.Date('1986-01-01')
# enddate <- as.Date('2018-03-01')
# imd_df <- imd_df[imd_df$date >= startdate & imd_df$date <= enddate,]
# 
# # B. Clean BigThaw data #####
# 
# file <- 'ind_5km_JASMIN-BigThaw_1998-2018_4samples'
# bt_df <- read.csv(paste0('1_Original_Data/118_BigThaw/', file, '.csv'))
# 
# # Remove duplicates
# bt_df <- bt_df %>% distinct()
# 
# # Set hourly time parameters
# bt_df <- bt_df %>%
#   mutate(
#     timestamp = ymd_hms(time),
#     date = as.Date(timestamp),
#     year = year(timestamp),
#     month = lubridate::month(timestamp, label = TRUE, abbr = TRUE),
#     day = day(timestamp),
#     hour = hour(timestamp)
#   ) %>%
#   dplyr::select(date, year, month, day, hour, location, everything(), -time)
# 
# # Calculate the overall day and night mean, min and max (convert from K to C)
# temp_all <- bt_df %>%
#   group_by(date,location) %>%
#   summarise(
#     t2m = mean(temp_mean, na.rm = TRUE) - 273.15,
#     mn2t = min(temp_min, na.rm = TRUE) - 273.15,
#     mx2t = max(temp_max, na.rm = TRUE) - 273.15
#   )
# 
# # Calculate the day (06:00 to 18:00) mean, min and max (convert from K to C)
# temp_day <- bt_df %>%
#   dplyr::filter(hour >= 6 & hour < 18) %>%
#   group_by(date,location) %>%
#   summarise(
#     t2m_day = mean(temp_mean, na.rm = TRUE) - 273.15,
#     mn2t_day = min(temp_min, na.rm = TRUE) - 273.15,
#     mx2t_day = max(temp_max, na.rm = TRUE) - 273.15
#   )
# 
# # Calculate the night (18:00 to 06:00) mean, min and max (convert from K to C)
# temp_night <- bt_df %>%
#   dplyr::filter(hour < 6 | hour >= 18) %>%
#   group_by(date,location) %>%
#   summarise(
#     t2m_night = mean(temp_mean, na.rm = TRUE) - 273.15,
#     mn2t_night = min(temp_min, na.rm = TRUE) - 273.15,
#     mx2t_night = max(temp_max, na.rm = TRUE) - 273.15
#   )
# 
# # Merge processed temp data
# bt_df <- temp_all %>%
#   full_join(temp_day, by = c('date', 'location')) %>%
#   full_join(temp_night, by = c('date', 'location')) %>%
#   arrange(location, date) %>%
#   dplyr::select(location, date, mn2t, mx2t)
# 
# # C. Get elevation data (elevatr package & BigThaw orography) ####
# 
# # Load elevation rasters
# elev_raster <- raster("2_Active_Data/204_misc/ind_envi_s3_AWSOpenData_elevation_raster.tif")
# res(elev_raster) * 111000  # in metres
# loc_points$elev <- raster::extract(elev_raster, loc_points_sf, fun = mean, na.rm = TRUE)
# elev_res <- dplyr::select(loc_points, location, elev) %>% st_drop_geometry()
# elev_res$elev <- as.numeric(elev_res$elev[,1])
# 
# # Load BigThaw model orography
# oro_df <- read.csv(paste0('1_Original_Data/118_BigThaw/ind_5km_JASMIN-BigThaw_orography_17sites.csv'))
# elev <- merge(elev_res, oro_df, by = c("location"))
# names(elev) <- c('location', 'elev', 'orog')
# 
# # D. Merge datasets ####
# bt_df <- merge(bt_df, elev, by = c("location"))
# df <- merge(imd_df, bt_df, by = c("location", "date"), suffixes = c("_imd", "_bt"))
# save(df, file=paste0('2_Active_Data/202_weat/ind_city_obs_IMD-BigThaw_1998-2018.RData'))

# 0. Adjust for elevation differences ####

load('2_Active_Data/202_weat/ind_city_obs_IMD-BigThaw_1998-2018.Rdata')

# Get difference
df$elev_diff <- df$orog -df$elev
unique(df$elev_diff)

# Calculate lapse rate: Assuming 0.5°C decrease for every 100 metres of altitude gained
df$mn2t_bt <- df$mn2t_bt + df$elev_diff * 0.5 / 100
df$mx2t_bt <- df$mx2t_bt + df$elev_diff * 0.5 / 100

# NB: Differences are slight - outputs should not differ greatly from previous report

# 1. Basic deviation (bias, error magnitude, agreement) #####

compute_metrics <- function(obs, model) {
  err <- model - obs
  n   <- sum(!is.na(err))
  tibble(
    n     = n,
    ME    = mean(err, na.rm = TRUE),                       # mean bias (systematic error)
    MAE   = mean(abs(err), na.rm = TRUE),                  # mean absolute error
    RMSE  = sqrt(mean(err^2, na.rm = TRUE)),                # root mean sq error
    r     = cor(obs, model, use = "complete.obs"),          # linear association
    sd_IMD   = sd(obs, na.rm = TRUE),
    sd_BigThaw = sd(model, na.rm = TRUE)
  )
}

metrics_mn2t <- compute_metrics(df$mn2t_imd, df$mn2t_bt)
metrics_mx2t <- compute_metrics(df$mx2t_imd, df$mx2t_bt)

cat("=== Min temp (mn2t) deviation metrics, pooled ===\n"); print(metrics_mn2t)
cat("=== Max temp (mx2t) deviation metrics, pooled ===\n"); print(metrics_mx2t)




## Per-station breakdown — reanalysis error is rarely uniform across
## sites (elevation, terrain complexity, station siting all matter)
metrics_by_station <- df %>%
  group_by(location) %>%
  summarise(
    n     = sum(!is.na(mn2t_bt - mn2t_imd)),
    ME_mn    = mean(mn2t_bt - mn2t_imd, na.rm = TRUE),
    RMSE_mn  = sqrt(mean((mn2t_bt - mn2t_imd)^2, na.rm = TRUE)),
    ME_mx    = mean(mx2t_bt - mx2t_imd, na.rm = TRUE),
    RMSE_mx  = sqrt(mean((mx2t_bt - mx2t_imd)^2, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(desc(RMSE_mx))

print(metrics_by_station)

# 2. Error autocorrelation ####

## Daily bias in weather data is almost always autocorrelated
## (a warm-bias day is likely followed by another warm-bias day).
## This matters because it breaks the iid assumption behind a
## naive bootstrap or a naive SE = sd/sqrt(n) calculation.

one_station <- df %>% filter(location == df$location[1]) %>% arrange(date)
err_series  <- one_station$mx2t_bt - one_station$mx2t_imd
Acf(err_series, main = "ACF of BT-IMD error (mx2t), example station")
# If bars extend beyond the blue dashed lines for several lags,
# the error is autocorrelated -> use block bootstrap, not iid bootstrap.

# 3. Block bootstrap with CI for bias metrics ####

## Moving block bootstrap: resample contiguous blocks of days
## (preserves short-term autocorrelation structure) rather than
## resampling individual days independently.

block_bootstrap_ci <- function(obs, model, block_len = 10, R = 2000) {
  err <- model - obs
  err <- err[!is.na(err)]
  n <- length(err)
  n_blocks <- ceiling(n / block_len)
  
  stat_fun <- function(data, indices) {
    # indices select block START points; build resampled series from blocks
    resampled <- unlist(lapply(indices, function(i) {
      end <- min(i + block_len - 1, n)
      data[i:end]
    }))
    c(ME = mean(resampled), RMSE = sqrt(mean(resampled^2)))
  }
  
  starts <- 1:(n - block_len + 1)
  boot_out <- boot(
    data = err,
    statistic = function(data, i) stat_fun(data, sample(starts, n_blocks, replace = TRUE)),
    R = R
  )
  
  ci_ME   <- boot.ci(boot_out, type = "perc", index = 1)$percent[4:5]
  ci_RMSE <- boot.ci(boot_out, type = "perc", index = 2)$percent[4:5]
  
  list(
    ME_est   = mean(err), ME_CI   = ci_ME,
    RMSE_est = sqrt(mean(err^2)), RMSE_CI = ci_RMSE,
    boot_out = boot_out
  )
}

set.seed(123)
ci_mx2t <- block_bootstrap_ci(df$mx2t_imd, df$mx2t_bt, block_len = 10, R = 2000)
ci_mn2t <- block_bootstrap_ci(df$mn2t_imd, df$mn2t_bt, block_len = 10, R = 2000)

cat("\nmx2t bias:", round(ci_mx2t$ME_est,3),
    " 95% CI:", round(ci_mx2t$ME_CI,3), "\n")
cat("mx2t RMSE:", round(ci_mx2t$RMSE_est,3),
    " 95% CI:", round(ci_mx2t$RMSE_CI,3), "\n")

cat("\nmn2t bias:", round(ci_mn2t$ME_est,3),
    " 95% CI:", round(ci_mn2t$ME_CI,3), "\n")
cat("mn2t RMSE:", round(ci_mn2t$RMSE_est,3),
    " 95% CI:", round(ci_mn2t$RMSE_CI,3), "\n")

## block_len choice: pick roughly the decorrelation length seen in the
## ACF plot above (e.g. if ACF drops off after ~5-10 lags, block_len=10 is reasonable)

# 4. Bland-Altman agreement analysis ####

## Shows whether disagreement depends on the magnitude of temperature
## (e.g. reanalysis might be fine in mid-range temps but diverge at extremes)

bland_altman_plot <- function(obs, model, varname) {
  d <- tibble(
    mean_val = (obs + model) / 2,
    diff_val = model - obs
  ) %>% drop_na()
  
  bias <- mean(d$diff_val)
  sdiff <- sd(d$diff_val)
  loa_upper <- bias + 1.96 * sdiff
  loa_lower <- bias - 1.96 * sdiff
  
  ggplot(d, aes(mean_val, diff_val)) +
    geom_point(alpha = 0.2) +
    geom_hline(yintercept = bias, color = "blue") +
    geom_hline(yintercept = loa_upper, color = "red", linetype = "dashed") +
    geom_hline(yintercept = loa_lower, color = "red", linetype = "dashed") +
    labs(title = paste("Bland-Altman:", varname),
         subtitle = paste0("Bias=", round(bias,2),
                           "  Limits of agreement=[", round(loa_lower,2),
                           ", ", round(loa_upper,2), "]"),
         x = "Mean of IMD & BT (°C)", y = "BT - IMD (°C)") +
    theme_minimal()
}

print(bland_altman_plot(df$mx2t_imd, df$mx2t_bt, "mx2t"))
print(bland_altman_plot(df$mn2t_imd, df$mn2t_bt, "mn2t"))

# 5. Deming regression analysis (errors-in-variables) ####

## Ordinary least squares assumes the x-axis variable is error-free.
## Since BOTH IMD and BT carry uncertainty, Deming regression is more
## appropriate than OLS for assessing systematic (proportional) bias.

deming_mx2t <- deming(mx2t_bt ~ mx2t_imd, data = df)
print(summary(deming_mx2t))
# Slope != 1 or intercept != 0 indicates systematic (not just random) bias

# 6. Seasonal/Monthly breakdown of bias ####

seasonal_bias <- df %>%
  group_by(month) %>%
  summarise(
    ME_mx2t = mean(mx2t_bt - mx2t_imd, na.rm = TRUE),
    ME_mn2t = mean(mn2t_bt - mn2t_imd, na.rm = TRUE),
    RMSE_mx2t = sqrt(mean((mx2t_bt - mx2t_imd)^2, na.rm = TRUE)),
    .groups = "drop"
  )
print(seasonal_bias)

ggplot(seasonal_bias, aes(x = factor(month, levels = month.abb), y = ME_mx2t)) +
  geom_col(fill = "steelblue") +
  labs(title = "Monthly mean bias in mx2t (BT - IMD)", x = "Month", y = "Bias (°C)") +
  theme_minimal()


# 7. TREND CHECK: does reanalysis bias drift over the years? ####

yearly_bias <- df %>%
  group_by(year) %>%
  summarise(ME_mx2t = mean(mx2t_bt - mx2t_imd, na.rm = TRUE), .groups = "drop")

trend_model <- lm(ME_mx2t ~ year, data = yearly_bias)
print(summary(trend_model))
# Significant slope -> reanalysis bias is not stationary over time,
# which matters if you're using BT data for long-term trend analysis

ggplot(yearly_bias, aes(year, ME_mx2t)) +
  geom_line() + geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Yearly mean bias in mx2t (BT - IMD)", y = "Bias (°C)") +
  theme_minimal()

#### SUMMARY: what to report ####
## - Pooled and per-station ME/MAE/RMSE (Section 1)
## - Block-bootstrap 95% CI on ME and RMSE (Section 3) <- your main
##   "statistical uncertainty" numbers, honestly accounting for
##   autocorrelation, rather than pretending each day is independent
## - Bland-Altman limits of agreement (Section 4) <- practical range
##   within which BT deviates from IMD 95% of the time
## - Deming slope/intercept (Section 5) <- tells you if bias is
##   systematic/proportional vs random scatter
## - Seasonal and yearly bias patterns (Sections 6-7) <- tells you
##   whether uncertainty is stable or context-dependent```

