### ------------------------------------------------------------------------ ### 
### Sensitivity analysis - BigThaw vs ERA5                                   ###
### ------------------------------------------------------------------------ ### 

# Turn off scientific notation
options(scipen=999)

# Set data path
user <- 'offne'
path <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
path1 <- paste0('C:/Users/', user, '/Documents/GitHub/HISP_TSGC/')
shapefile_path <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/1_Original_Data/100_GADM/gadm41_IND_shp/')
setwd(path)

####  PACKAGES & FUNCTIONS ####
library(tidyverse)
library(dplyr)
library(skimr)
library(DT)
library(gridExtra)
library(ggplot2)
library(zoo)
library(forecast)
library(patchwork)
library(purrr)
library(slider)
library(raster)
library(sf)
library(elevatr)
library(ggnewscale)
library(cowplot)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)

# Get functions (make sure wd is in repo folder)
source(paste0(path1, 'TSGC_Model - Functions.R'))

# Function to create a histogram for a column
create_histogram <- function(df, column_name, stat=NULL) {
  # df <- sub
  # stat <- skm
  # column_name <- 'mn2t'
  hist_data <- hist(df[[column_name]], breaks=40, plot = FALSE)
  y_max_idx <- which.max(hist_data$density)
  y_med_idx <- which.min(abs(hist_data$density - median(hist_data$density)))
  position <- y_max_idx > y_med_idx
  
  if(position == TRUE){
    x_pos <- min(df[[column_name]], na.rm = TRUE)
    h_pos <- 0
  } else {
    x_pos <- max(df[[column_name]], na.rm = TRUE)
    h_pos <- 1
  }
  
  plot <- ggplot(df, aes_string(x = column_name)) +
    geom_histogram(fill = mypalette[1], color = 'black') +
    theme_minimal() +
    theme(axis.title.y = element_blank()) +
    ggtitle(paste(column_name)) 
  
  if (!is.null(stat)){
    plot <- plot +
      annotate('text',
               x = x_pos, y = Inf,
               label = paste('\nNA: ', stat[stat$skim_variable == column_name, 'n_missing'],
                             '\nMin: ', stat[stat$skim_variable == column_name, 'min'],
                             '\nMedian: ', stat[stat$skim_variable == column_name, 'median'],
                             '\nMax: ', stat[stat$skim_variable == column_name, 'max']),
               hjust = h_pos, vjust = 1, size = 3, color = 'black', fontface = 'italic')
  }
  
  return(plot)
  
}

# Time series plots
plot_ts <- function(df, column_name) {
  ggplot(df, aes(x = date, y = .data[[column_name]])) + 
    geom_line(color = mypalette[1]) +
    ggtitle(column_name) +
    xlab('Date') +                
    ylab('Patients') +      
    theme_minimal() +            
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
    scale_x_date(breaks = seq.Date(min(df$date), max(df$date), by = '6 months'), 
                 labels = scales::date_format('%b %Y'))
}

# Elevation plot
create_elevation_map <- function(hp, ind, cities, locations, z = 7, buffer_dist = 10000) {
  
  # 1. Get bbox of state + buffer
  bbox_geom <- st_as_sfc(st_bbox(hp))
  bbox <- st_sf(st_buffer(bbox_geom, dist = buffer_dist))
  
  # 2. Get elevation data (zoom level z controls resolution)
  elev_raster <- get_elev_raster(locations = bbox, z = z, clip = 'bbox')
  
  # 3. Inverse mask out Admin1 outside HP
  bbox2 <- st_sf(st_buffer(bbox, dist = 10000))
  outside_hp <- st_difference(bbox2, st_union(hp))
  
  # 4. Compute terrain derivatives
  slope <- terrain(elev_raster, opt = 'slope')
  aspect <- terrain(elev_raster, opt = 'aspect')
  hillshade <- hillShade(slope, aspect, angle = 40, direction = 315)
  
  # 5. Convert for plotting
  elev_df <- as.data.frame(elev_raster, xy = TRUE)
  colnames(elev_df) <- c('x', 'y', 'elevation')
  hill_df <- as.data.frame(hillshade, xy = TRUE)
  colnames(hill_df) <- c('x', 'y', 'hill')
  
  # 6. Set colour scheme
  arcgis_topo_colors <- c(
    '#5A7D24',  # dark green (low elevation)
    '#9DBF68',  # light green
    '#F3E99F',  # pale yellow
    '#D2B48C',  # tan/light brown
    '#8B7A55',  # brown (higher elevation)
    '#FFFFFF'   # white (highest peaks/snow)
  )
  
  # 7. Create map
  main_map <- ggplot() +
    # Hillshade raster
    geom_raster(data = hill_df, aes(x = x, y = y, fill = hill), show.legend = FALSE) +
    scale_fill_gradient(low = 'grey20', high = 'grey90', guide = 'none') +
    new_scale_fill() +
    
    # Elevation raster
    geom_raster(data = elev_df, aes(x = x, y = y, fill = elevation), alpha = 0.6) +
    scale_fill_gradientn(colours = arcgis_topo_colors, name = 'Elevation (m)') +
    
    # Mask outside HP
    geom_sf(data = outside_hp, fill = 'grey50', alpha = 0.6) +
    
    # Boundaries and cities
    geom_sf(data = ind, fill = NA, color = 'grey50', size = 0.6) +
    geom_sf(data = hp, fill = NA, color = 'black', size = 0.6) +
    geom_sf(data = cities, color = 'darkred', size = 2) +
    
    # City labels
    geom_label(data = locations, aes(x = lon, y = lat, label = location),
               size = 3, color = 'darkred', fill = 'white', fontface = 'bold',
               label.padding = unit(0.15, 'lines'), label.r = unit(0.15, 'lines'),
               label.size = 0, nudge_y = 0.1) +
    # geom_label_repel(data = locations, aes(x = lon, y = lat,
    #     label = paste0(
    #       location, '\n',
    #       'mn2t: ', round(diff_mn2t, 1), '\n',
    #       'mx2t: ', round(diff_mx2t, 1), '\n',
    #       'tp: ', round(diff_tp, 1), '\n',
    #       'fg10: ', round(diff_fg10, 1)
    #     )),
    #   size = 3, color = 'darkred', fill = 'white', fontface = 'bold',
    #   label.padding = unit(0.15, 'lines'), label.r = unit(0.15, 'lines'),
    #   label.size = 0) +
    # 
    coord_sf(xlim = range(elev_df$x), ylim = range(elev_df$y), expand = FALSE) +
    labs(x = 'Longitude', y = 'Latitude') +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = 'white', color = NA),
      panel.background = element_rect(fill = 'white', color = NA)
    )
  
  return(main_map)
}


####  MAIN CODE ####

# Select year of interest
# era5 <- read.csv('202_weat/ind_30km_JASMIN-ERA5_1979-2023_HISP_anom.csv')
# era5 <- era5[era5$date >= '2015-01-01' & era5$date <= '2023-12-31',]
# write.csv(era5, file='202_weat/ind_30km_JASMIN-ERA5_2015-2023_HISP_anom.csv', row.names = FALSE)

# 1. Clean data ####
bigt <- read.csv('202_weat/ind_10km_BigThaw_2015-2016_anom.csv')
era5 <- read.csv('202_weat/ind_30km_JASMIN-ERA5_2015-2023_HISP_anom.csv')
heal <- read.csv('203_comb/ind_30km_HISP-HP_2015-2023_hosp-era5_chamba_wide_model.csv', check.names=FALSE, fileEncoding='UTF-8')

# Set date type, and time frame
bigt$date <- as.Date(bigt$date, format='%Y-%m-%d')  # Adjust format as needed
era5$date <- as.Date(era5$date, format='%Y-%m-%d')  # Adjust format as needed
heal$date <- as.Date(heal$date, format='%Y-%m-%d')  # Adjust format as needed
era5 <- era5[era5$date >= '2015-01-02' & era5$date <= '2017-01-02',]

# Back transform tp variables
era5$tp <- exp(era5$tp)
bigt$tp <- exp(bigt$tp)

# 2. Comparative overall distributions
# Bigt tails are longer because we can better capture extremes.


# 2. Comparative distributions, overall and by season ####
# Use wide format to avoid issues with date duplicates
df_w <- era5
sub <- df_w %>% dplyr::select(mn2t,mx2t,tp,fg10)
skm <- skim(sub) %>% rename(min = numeric.p0, Q1 = numeric.p25, 
                            median = numeric.p50, Q3 = numeric.p75, 
                            max = numeric.p100) %>% mutate(across(where(is.numeric), ~ round(., 2)))

colnames(skm) <- gsub('numeric.', '', colnames(skm))
skm <- as.data.frame(skm)

# Histograms
cols_to_plot <- colnames(sub)
plots <- lapply(cols_to_plot, function(col) create_histogram(sub, col, skm))
gridExtra::grid.arrange(grobs = plots[1:4], ncol = 4, nrow = 1)

# Density curves by season
labels_map <- c('Winter (Jan-Feb)' = 'Jan-Feb',
                'Pre-Monsoon (Mar-May)' = 'Mar-May',
                'Monsoon (Jun-Sep)' = 'Jun-Sep',
                'Post-Monsoon (Oct-Dec)' = 'Oct-Dec')
df_w$tp <- log(df_w$tp)
plots <- lapply(cols_to_plot, function(c) {
  p <- ggplot(df_w, aes_string(x = c, color = 'season', fill = 'season')) +
    geom_density(alpha = 0.3) +
    theme_minimal() +
    labs(
      title = c,
      x = c,
      y = 'Density',
      fill = 'Season',
      color = 'Season'
    ) +
    scale_fill_manual(values=c(mypalette[1], mypalette[4], mypalette[2], mypalette[3]), labels=labels_map) +
    scale_color_manual(values=c(mypalette[1], mypalette[4], mypalette[2], mypalette[3]), labels=labels_map)+ 
    theme(legend.position = 'none')
  # 
  #   if (c == 'tp') {
  #     p <- p + scale_x_continuous(
  #       breaks = log(c(0.0000000000001, 50, 100, 150, 200, 250)),
  #       labels = c(0, 50, 100, 150, 200, 250)
  #     )
  #   }
  p
})
gridExtra::grid.arrange(grobs = plots[1:4], ncol = 4, nrow = 1)

# 3. Comparative time series ####
# Temporal trends are the same, but spikes and dips more identifiable.

w_cols_to_plot <- c('mn2t','mx2t','tp','fg10')

# Check time series pattern
plot_ts <- function(df, column_name, window = 30, palette=1) {
  df <- df %>%
    arrange(date) %>%
    mutate(ma = slide_dbl(.data[[column_name]], mean, .before = window - 1, .complete = TRUE))
  
  ggplot(df, aes(x = date)) + 
    geom_line(aes(y = .data[[column_name]]), color = mypalette[palette], alpha = 0.5) +  # original line, lighter
    geom_line(aes(y = ma), color = mypalette[palette], size = 1) +                       # moving average
    ggtitle(column_name) +
    xlab('Date') +                
    ylab(column_name) +      
    theme_minimal() +            
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
    scale_x_date(breaks = seq.Date(min(df$date), max(df$date), by = '6 months'), 
                 labels = scales::date_format('%b %Y'))
}
df_w <- bigt
ts_plots <- lapply(w_cols_to_plot, function(col) plot_ts(df_w, col, palette=1))
grid.arrange(grobs = ts_plots[1:4], ncol = 4, nrow = 1)
df_w <- era5
ts_plots <- lapply(w_cols_to_plot, function(col) plot_ts(df_w, col, palette=1))
grid.arrange(grobs = ts_plots[1:4], ncol = 4, nrow = 1)

# Compare moving averages
plot_ts_compare <- function(df1, df2, column_name, window = 30, palette = c(1,2)) {
  
  # Prepare first df
  df1 <- df1 %>%
    arrange(date) %>%
    mutate(ma = slide_dbl(.data[[column_name]], mean, .before = window - 1, .complete = TRUE)) %>%
    mutate(source = 'BigThaw')
  
  # Prepare second df
  df2 <- df2 %>%
    arrange(date) %>%
    mutate(ma = slide_dbl(.data[[column_name]], mean, .before = window - 1, .complete = TRUE)) %>%
    mutate(source = 'ERA5')
  
  # Combine
  df_all <- bind_rows(df1, df2)
  
  ggplot(df_all, aes(x = date, color = source)) + 
    geom_line(aes(y = ma), size = 1) +                                # moving average lines
    ggtitle(column_name) +
    xlab('Date') +                
    ylab(column_name) +      
    theme_minimal() +            
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
    scale_color_manual(values = mypalette[palette], name = 'Dataset') +
    scale_x_date(breaks = seq.Date(min(df_all$date, na.rm=TRUE), 
                                   max(df_all$date, na.rm=TRUE), 
                                   by = '6 months'), 
                 labels = scales::date_format('%b %Y'))+
    theme(legend.position = 'none')
  
}
ts_plots <- lapply(w_cols_to_plot, function(col) plot_ts_compare(bigt, era5, col, palette = c(1,2)))
grid.arrange(grobs = ts_plots[1:4], ncol = 2, nrow = 2)


# # 4. Map differences across spaces ####
# # (i.e. heatmap of change from era5 to bigt across locations)
# # Illustrate that more mountainous areas are more susceptible to bigger differences
# 
# # Get admin data for india and filtered to Himachal Pradesh
# ind <- st_read(paste0(shapefile_path,'gadm41_IND_1.shp'))
# shp <- st_read(paste0(shapefile_path,'gadm41_IND_3.shp'))
# hp <- shp[shp$NAME_1 == 'Himachal Pradesh', ]
# 
# # Create a data frame with locations and coordinates
# locations <- data.frame(
#   location = c('Chamba', 'Dharamshala', 'Solan', 'Hamirpur', 'Kullu'),
#   lon = c(76.127017, 76.321482, 77.099506, 76.522135, 77.107945),
#   lat = c(32.554798, 32.216995, 30.906408, 31.684233, 31.956550)
# )
# cities <- st_as_sf(locations, coords = c('lon', 'lat'), crs = 4326)
# 
# # Average difference in mn2t, mx2t, tp and fg10 between bigt and era5, by location
# avg_diff <- bigt %>%
#   dplyr::select(date, location, mn2t, mx2t, tp, fg10) %>%
#   left_join(era5 %>% dplyr::select(date, location, mn2t, mx2t, tp, fg10), 
#             by = c('date', 'location'), suffix = c('_bigt', '_era5')) %>%
#   group_by(location) %>%
#   mutate(
#     diff_mn2t = mn2t_bigt - mn2t_era5,
#     diff_mx2t = mx2t_bigt - mx2t_era5,
#     diff_tp = tp_bigt - tp_era5,
#     diff_fg10 = fg10_bigt - fg10_era5
#   ) %>%
#   summarise(across(starts_with('diff_'), mean, na.rm = TRUE))
# 
# locations <- locations %>% left_join(avg_diff, by = 'location')
# main_map <- create_elevation_map(hp, ind, cities, locations)
# main_map
# 
# 5. Table comparing frequency of IMD cut offs between datasets ####

x <- bigt %>% 
  dplyr::select(date, season, location, contains('IMD_b') & !contains('lag')) %>% 
  left_join(era5 %>% dplyr::select(date, season, location, contains('IMD_b') & !contains('lag')), 
            by = c('date', 'location', 'season'), suffix = c('_bigt', '_era5')) 

table(x$tp_IMD_b_era5, x$tp_IMD_b_bigt)


# Get IMD columns excluding those with 'lag' from bigt and era5
imd_cols_bigt <- names(bigt)[grepl("IMD_b", names(bigt)) & !grepl("lag", names(bigt))]
imd_cols_era5 <- names(era5)[grepl("IMD_b", names(era5)) & !grepl("lag", names(era5))]
# Extract the base names of IMD variables (without suffix)
imd_base_names <- sub("_bigt$", "", imd_cols_bigt)

# Initialize an empty list to store tables
comparison_tables <- list()


# Loop through each IMD variable base name
imd_vars <- sub("_bigt$", "", imd_cols_bigt)
comparison_tables <- lapply(imd_vars, function(var) {
  var_bigt <- paste0(var, "_bigt")
  var_era5 <- paste0(var, "_era5")
  
  tbl <- table(x[[var_era5]], x[[var_bigt]])
  
  # Convert to data.frame for nicer presentation
  as.data.frame(tbl) %>%
    rename(era5 = Var1, bigt = Var2, count = Freq) %>%
    mutate(variable = var)
})

# Combine all into one big table
comparison_df <- bind_rows(comparison_tables) %>%
  arrange(variable, era5, bigt)

write.csv(comparison_df, paste0(path1,"IMD_comparison_table.csv"), row.names = FALSE)


# 6. Comparative quantile regression ####

# Format data for modelling
df <- bigt %>% 
  dplyr::select(date, season, location, contains('IMD_b') & !contains('lag')) %>% 
  left_join(era5 %>% dplyr::select(date, season, location, contains('IMD_b') & !contains('lag')), 
            by = c('date', 'location', 'season'), suffix = c('_bigt', '_era5')) %>%
  left_join(heal %>% dplyr::select(date, season, location, holiday, total_patients), 
            by = c('date', 'location', 'season')) %>%
  dplyr::select(location, date, season, holiday, total_patients, everything()) %>%
  filter(location== 'Chamba') %>%
  filter(date >= as.Date('2015-02-14'))




