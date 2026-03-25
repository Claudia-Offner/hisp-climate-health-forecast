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
smooth_with_bounds <- function(x, k = 7) {
  sm  <- rollmean(x, k = k, fill = NA, align = 'right', na.rm = TRUE)
  se  <- rollapply(x, width = k, FUN = function(w) sd(w, na.rm = TRUE) / sqrt(sum(!is.na(w))),
                   fill = NA, align = 'right')
  list(
    sm  = sm,
    lwr = sm - 1.96 * se,
    upr = sm + 1.96 * se
  )
}

# Function to plot smoothed data
plot_smoothing <- function(data, location_name, col) {
  
  sm_col  <- paste0(col, '_sm')
  lwr_col <- paste0(col, '_lwr')
  upr_col <- paste0(col, '_upr')
  
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
cols_to_smooth <- c('total_patients', 'F', 'M', 't2m', 'tp')
sub <- df %>% dplyr::select(location, date, all_of(cols_to_smooth))
dates_analysis <- seq(as.Date('2017-01-01'), as.Date('2017-12-31'), by = 'day')


# Create regional estimates (i.e. metropolitan as per Paireau)
sub_regional <- sub %>%
  group_by(date) %>%
  summarise(
    across(all_of(cols_to_smooth[1:3]),  ~ sum(.x,  na.rm = TRUE)),
    across(all_of(cols_to_smooth[4:5]), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(location = "metropolitan")
sub <- bind_rows(sub, sub_regional)


# Create date_analysis and smooth select variables
# NB: For any given 'date_analysis' date, the dataset contains all the data 
# available at this date (the smoothing is done without using future data).

result <- lapply(dates_analysis, function(d) {
  sub %>%
    filter(date <= d) %>%
    arrange(location, date) %>%
    group_by(location) %>%
    mutate(across(
      all_of(cols_to_smooth),
      list(
        sm  = ~ smooth_with_bounds(.x)$sm,
        lwr = ~ smooth_with_bounds(.x)$lwr,
        upr = ~ smooth_with_bounds(.x)$upr
      ),
      .names = '{.col}_{.fn}'
    )) %>%
    ungroup() %>%
    mutate(date_analysis = as.Date(d, origin = '1970-01-01'))
}) %>%
  bind_rows() %>% 
  dplyr::select(location, date, date_analysis, everything())

# Single location / column
plot_smoothing(result, 'Una', 'total_patients')

# All locations at once with facets
result %>%
  filter(date_analysis == max(date_analysis),
         between(date, as.Date('2017-01-01'), as.Date('2027-12-31'))) %>%
  ggplot(aes(x = date)) +
  geom_ribbon(aes(ymin = total_patients_lwr, ymax = total_patients_upr),
              fill = '#378ADD', alpha = 0.2) +
  geom_line(aes(y = total_patients),    colour = '#888780', linewidth = 0.4, na.rm = TRUE) +
  geom_line(aes(y = total_patients_sm), colour = '#378ADD', linewidth = 0.8, na.rm = TRUE) +
  facet_wrap(~ location, scales = 'free_y') +
  labs(
    title    = 'Total patients — raw vs 7-day MA by location',
    subtitle = 'Grey = raw  |  Blue = 7-day MA  |  Shaded = 95% CI',
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.subtitle = element_text(colour = 'grey50'))

df <- result

#### EXPORT
save(df , file = paste0(path_git, '3_Ensemble/', file_import, '_smooth.RData'))
