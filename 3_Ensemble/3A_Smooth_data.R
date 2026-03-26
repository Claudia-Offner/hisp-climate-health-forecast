### ------------------------------------------------------------------------ ### 
###  Ensemble - Smoothing functions
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list = ls()); graphics.off(); gc()

# Set data path
path_data <- paste0('C:/Users/claer14/OneDrive - University of Cambridge/III. Research/Data/')
path_git <- paste0('C:/Users/claer14/Documents/GitHub/hisp-climate-health-forecast/')
setwd(path_git)

#### DEPENDENCIES ####
source('3_Ensemble/Ensemble - Utilities.R')
file_import <- 'ind_5km_BigThaw-HISP-ALL_2015-2023_hosp_wide_ensemble'

# Install packages if not already in library & load
pkgs <- c('dplyr', 'ggplot2', 'zoo')
install.packages(pkgs[!pkgs %in% installed.packages()[, 'Package']])
invisible(lapply(pkgs, library, character.only = TRUE))

# Function applying 7-day moving average with CI's
smooth_with_bounds <- function(x, k = 14) {
  sm  <- rollmean(x, k = k, fill = NA, align = 'right', na.rm = TRUE)
  se  <- rollapply(x, width = k, FUN = function(w) sd(w, na.rm = TRUE) / sqrt(sum(!is.na(w))),
                   fill = NA, align = 'right')
  list(
    smooth  = sm,
    smooth_lwr = sm - 1.96 * se,
    smooth_upr = sm + 1.96 * se
  )
}

# Function to plot smoothed data
plot_smoothing <- function(data, location_name, col) {
  
  sm_col  <- paste0(col, '_smooth')
  lwr_col <- paste0(col, '_smooth_lwr')
  upr_col <- paste0(col, '_smooth_upr')
  
  # Use the most recent date_analysis to get the full series
  latest <- data %>%
    filter(location == location_name) %>%
    filter(date_analysis == max(date_analysis),
           between(date, as.Date('2017-01-01'), as.Date('2027-12-31')))
  
  ggplot(latest, aes(x = date)) +
    geom_ribbon(aes(ymin = .data[[lwr_col]], ymax = .data[[upr_col]]),
                fill = '#378ADD', alpha = 0.2) +
    geom_line(aes(y = .data[[col]]), colour = '#888780', linewidth = 0.5, na.rm = TRUE) +
    geom_line(aes(y = .data[[sm_col]]), colour = '#378ADD', linewidth = 1, na.rm = TRUE) +
    labs(
      title    = paste(location_name, '—', col),
      subtitle = 'Grey = raw  |  Blue = 7-day MA  |  Shaded = 95% CI',
      x        = NULL,
      y        = col
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.subtitle = element_text(colour = 'grey50'))
}

####  MAIN CODE  #####

# Load data
load(paste0('3_Ensemble/', file_import, '.RData'))
df <- df_ensm; rm(df_ensm)

# Subset relevant data
cols_to_smooth <- c('total_patients', 'F', 'M', 'age25_64', 't2m', 'tp')
sub <- df %>% 
  dplyr::select(location, date, all_of(cols_to_smooth)) %>%
  filter(between(date, as.Date('2016-05-15'), as.Date('2017-12-31')))
dates_analysis <- seq(as.Date('2017-01-01'), as.Date('2017-09-30'), by = 'day')

# Create regional estimates (i.e. metropolitan as per Paireau)
sub_regional <- sub %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(
    across(all_of(cols_to_smooth[1:3]), ~ sum(.x,  na.rm = TRUE)),
    across(all_of(cols_to_smooth[4:5]), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(location = "metropolitan")
sub <- bind_rows(sub, sub_regional)


# Create date analysis df, address NA's and smooth selected variables ####

# NB: For any given 'date_analysis' date, the dataset contains all the data 
# available at this date (the smoothing is done without using future data).

# Data analysis df
sub <- lapply(dates_analysis, function(d) {
  sub %>%
    filter(date <= d-2) %>%
    group_by(location) %>%
    arrange(date) %>%
    mutate(date_analysis = as.Date(d, origin = '1970-01-01'))
}) %>%
  bind_rows() %>% 
  dplyr::select(location, date, date_analysis, everything()) %>%
  arrange(location, date_analysis, date)

# Address NA's by assuming static position
sub <- sub %>%
  group_by(location, date_analysis) %>%
  tidyr::fill(where(is.numeric), .direction = "down") %>%
  ungroup()

# Calculate smoothed variables with growth rate
result <- sub %>%
  group_by(location, date_analysis) %>%
  mutate(across(
    all_of(cols_to_smooth),
    list(
      smooth     = ~ smooth_with_bounds(.x)$smooth,
      smooth_lwr = ~ smooth_with_bounds(.x)$smooth_lwr,
      smooth_upr = ~ smooth_with_bounds(.x)$smooth_upr
    ),
    .names = '{.col}_{.fn}'
  )) %>%
  mutate(total_patients_r = log(total_patients_smooth / lag(total_patients_smooth))) %>%
  ungroup() %>%
  filter(between(date, as.Date('2016-06-01'), as.Date('2017-12-31')))

# Single location / column
plot_smoothing(result, 'Una', 'total_patients')

# All locations at once with facets
result %>%
  filter(date_analysis == max(date_analysis),
         between(date, as.Date('2017-01-01'), as.Date('2027-12-31'))) %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = total_patients_smooth_lwr, ymax = total_patients_smooth_upr),
              fill = '#378ADD', alpha = 0.2) +
  geom_line(aes(y = total_patients),    colour = '#888780', linewidth = 0.4, na.rm = TRUE) +
  geom_line(aes(y = total_patients_smooth), colour = '#378ADD', linewidth = 0.8, na.rm = TRUE) +
  facet_wrap(~ location, scales = 'free_y') +
  labs(
    title    = 'Total patients — raw vs 7-day MA by location',
    subtitle = 'Grey = raw  |  Blue = 7-day MA  |  Shaded = 95% CI',
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.subtitle = element_text(colour = 'grey50'))

# Check growth rate
result %>%
  filter(date_analysis == max(date_analysis),
         location=='metropolitan') %>%
  arrange(location, date) %>%
  ggplot(aes(x = date)) +
  # Reference line at r = 0 (stable)
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#888780", linewidth = 0.4) +
  # Growth rate line
  geom_line(aes(y = total_patients_r), colour = "#378ADD", linewidth = 0.7, na.rm = TRUE) +
  facet_wrap(~ location, scales = "free_y") +
  labs(
    title    = "Log growth rate of smoothed patients by location",
    subtitle = paste("Snapshot:", format(max(result$date_analysis), "%d %b %Y")),
    x        = NULL,
    y        = "Log growth rate"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.subtitle      = element_text(colour = "grey50"),
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(face = "bold")
  )

df <- result

#### EXPORT
save(df , file = paste0(path_git, '3_Ensemble/', file_import, '_smooth.RData'))
