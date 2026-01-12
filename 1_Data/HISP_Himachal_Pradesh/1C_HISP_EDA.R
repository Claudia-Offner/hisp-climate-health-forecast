### ------------------------------------------------------------------------ ### 
### HISP Exploratory Data Analysis                                           ###
### ------------------------------------------------------------------------ ### 

setwd(here())
source('1_Data/HISP_Himachal_Pradesh/HISP_Cleaning - Dependencies.R')

# Set data path
setwd(path_data)

#### FUNCTIONS ####

# Function to create a histogram for a column
create_histogram <- function(df_all, column_name, stat=NULL, color=mypalette[1]) {
  # df_all <- sub
  # stat <- skm
  # column_name <- 'fg10'
  hist_data <- hist(df_all[[column_name]], breaks=40, plot = FALSE)
  y_max_idx <- which.max(hist_data$counts)
  y_med_idx <- which.min(abs(hist_data$counts - median(hist_data$counts)))

  if(column_name %in% c('fg10')){
    position <- FALSE
  } else if (column_name %in% c('rh_mean')) {
    position <- TRUE
  } else{
    position <- y_max_idx > y_med_idx
  }
  
  if(position == TRUE){
    x_pos <- min(df_all[[column_name]], na.rm = TRUE)
    h_pos <- 0
  } else {
    x_pos <- max(df_all[[column_name]], na.rm = TRUE)
    h_pos <- 1
  }
  
  title_name <- ifelse(column_name %in% names(nm), nm[[column_name]], column_name)
  
  plot <- ggplot(df_all, aes_string(x = column_name)) +
    geom_histogram(fill = color, color = 'black') +
    theme_minimal() +
    theme(axis.title.y = element_blank()) +
    labs(
      title = title_name,
      # x = title_name
    )
  
  if (!is.null(stat)){
    plot <- plot +
      annotate('label',
               x = x_pos, y = Inf,
               label = paste(
                 'NA: ', stat[stat$skim_variable == column_name, 'n_missing'],
                 '\nMin: ', stat[stat$skim_variable == column_name, 'min'],
                 '\nMedian: ', stat[stat$skim_variable == column_name, 'median'],
                 '\nMax: ', stat[stat$skim_variable == column_name, 'max']
               ),
               hjust = h_pos, vjust = 1,
               size = 3, 
               color = 'black',           # text color
               fontface = 'italic',
               fill = 'white',            # background color
               label.size = 0.3,          # ✅ border thickness (0 = none)
               label.r = unit(0.15, 'lines')  # rounded corners (optional)
      )
  
  }
  
  return(plot)

}

# Function to get base stats of df
get_desc <- function(df, h_cols, w_cols) {
  
  # Handling for location specific tables
  if (length(unique(df$location)) == 1){
    h_cols <- h_cols[h_cols != "location"]
    check <- c('F', 'under_5', 'general', 't2m')
    labels <- c('Gender', 'Age', 'Ward', 'Weather Characteristics')
    } else{
    check <- c(sort(unique(as.character(df$location)))[1], 'F', 'under_5', 'general', 't2m')
    labels <- c('Location', 'Gender', 'Age', 'Ward', 'Weather Characteristics')
  }
  
  vars <- c(h_cols, w_cols)
  
  # Get descriptive stats 
  (desc <- vars %>%
      lapply(function(col) {
        x <- df[[col]]
        
        
        if(is.numeric(x)){
          
          # Check if binary
          is_binary <- all(x %in% c(0, 1, NA), na.rm = TRUE)
          
          perc_val <- if (is_binary) {
            # For binary: percent of times '1' appears relative to number of rows
            round((sum(x, na.rm = TRUE) / nrow(df)) * 100)
          } else {
            # For non-binary: sum of x relative to total patients, rounded to 1 decimal
            round((sum(x, na.rm = TRUE) / sum(df$total_patients, na.rm = TRUE)) * 100, 1)
          }
          
          mean_val <- if (is_binary) {
            NA
          } else {
            round(mean(x, na.rm = TRUE))        # round mean to whole number
          }
          
          sd_val <- if (is_binary) {
            NA
          } else {
            round(sd(x, na.rm = TRUE))          # round SD to whole number
          }
          
          min_val <- if (is_binary) {
            NA  
          } else {
            round(min(x,na.rm = TRUE))
          }
          
          max_val <- if (is_binary) {
            NA  
          } else {
            round(max(x,na.rm = TRUE))
          }
          
          skew_val <- if (is_binary) {
            NA  
          } else {
            round(skewness(x,na.rm = TRUE), 1)
          }
          
          data.frame(
            Variable = sub('_b$', '', col),
            total_visits = sum(x, na.rm = TRUE),
            perc_visits = perc_val,
            mean_daily_visits = mean_val,
            sd_daily_visits = sd_val,
            min_daily_visit = min_val,
            max_daily_visit = max_val,
            skew_daily_visit = skew_val
          )
          
          
        } else {
          
          df %>%
            group_by(location) %>%
            summarise(
              total_visits = sum(total_patients, na.rm = TRUE),
              mean_daily_visits = round(mean(total_patients, na.rm = TRUE)),
              sd_daily_visits = round(sd(total_patients, na.rm = TRUE)),
              min_daily_visit = round(min(total_patients, na.rm = TRUE)),
              max_daily_visit = round(max(total_patients, na.rm = TRUE)),
              skew_daily_visit = round(skewness(total_patients, na.rm = TRUE))
            ) %>%
            mutate(
              perc_visits = round((total_visits / sum(total_visits)) * 100, 1)
            ) %>%
            dplyr::select(Variable=location, total_visits, perc_visits, everything())
          
        }
        
      }) %>%
      bind_rows())
  
  
  (base_desc <- desc %>%
      mutate(
        Frequency = paste0(format(round(total_visits), big.mark = ','), ' (', perc_visits, '%)'),
        Average = ifelse(
          is.na(sd_daily_visits),
          NA,
          paste0(
            round(mean_daily_visits),
            ifelse(Variable == 'tp', ' + ', ' ± '),
            round(sd_daily_visits),
            ' (', min_daily_visit, ', ', max_daily_visit, ')'
          )

        ),
        # Range = ifelse(
        #   is.na(skew_daily_visit),
        #   NA,
        #   # paste0(min_daily_visit, ', ',max_daily_visit, ' (', skew_daily_visit, ')')
        #   paste0(sd_daily_visits, ' (', skew_daily_visit, ')')
        # )
      ) %>%
      dplyr::select(Variable, Frequency, Average))
  
  
  # Remove frequency values for raw weather data
  base_desc$Frequency[base_desc$Variable %in%  gsub('(_IMD_b|_p90|_p10)$', '', w_cols)] <- NA
  
  
  # FORMATING
  # Insert headings at fixed positions
  base_desc <- rbind(
    tibble(Variable = 'Hospital Characteristics', Frequency = NA, Average = NA),
    base_desc
  )
  # Find positions to insert before
  positions <- which(base_desc$Variable %in% check)
  new_rows <- tibble(
    Variable = labels,
    Frequency = NA,
    Average = NA,
  )
  # Insert rows in reverse order so positions don't shift
  for (i in rev(seq_along(positions))) {
    pos <- positions[i]
    base_desc <- bind_rows(
      base_desc[1:(pos-1), ],
      new_rows[i, ],
      base_desc[pos:nrow(base_desc), ]
    )
  }
  
  return(base_desc)
}

### Function for plotting data with moving average
plot_moving_avg <- function(variable, all.cases=TRUE, MA=7, y.lab='', colour='red', plt.start.date=NULL, est.date.index, estimation.date.end){
  
  # Determine plot length
  if (is.null(plt.start.date)) {plt.start.date <- head(est.date.index, 1)}
  tmp.date <- min(estimation.date.end - 4, plt.start.date)
  s <- sprintf('%s/', format(tmp.date, '%Y-%m-%d'))
  variable <- variable[s]
  
  # Calculate the 7 day centered moving average
  ma.cent <- zoo::rollmean(variable, MA, align='center') 
  
  # Create df for plot (Data data vs. MA data)
  d <- setNames(cbind(variable, ma.cent), c('New Cases', 'Centered 7-day MA'))
  df.f <- data.frame(Date=index(d), 
                     New.Cases=coredata(d[,1]),
                     Centered.7.day.MA=coredata(d[,2]))
  
  if(all.cases==TRUE){
    p1 <- ggplot(df.f, aes(x = Date)) + geom_line(aes(y = New.Cases, color = 'Real data'), linewidth = 0.1)
  } else{
    p1 <- ggplot(df.f, aes(x = Date))
  }
  
  # Plot Figure 1
  p1 <- p1 +
    geom_line(aes(y = Centered.7.day.MA, color = 'Centered 7 day MA'), linewidth = 1) +
    scale_color_manual(name = '', values = c('Real data' = 'black', 'Centered 7 day MA' = colour)) +
    scale_y_continuous(n.breaks = 6) +
    scale_x_date(labels = scales::date_format('%b %y'), 
                 breaks = seq(min(df.f$Date), max(df.f$Date), by = '1 year')) +
    labs(
      title = y.lab,
      y=y.lab, 
      x='Date' 
    ) +
    theme_economist_white(gray_bg = FALSE, base_size = 12) +
    theme(
      legend.position = 'none',
      axis.text.x = element_text(angle = 45, vjust = 0.75),
      axis.title = element_text(size = rel(1), margin = margin(t = 10, r = 10)),
      axis.text = element_text(size = rel(1)),
      plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5)),
      plot.caption = element_text(size = rel(1)),
      text = element_text(size = rel(1.1))
      ) +
    scale_linetype_manual(values = rep('solid', 3)) +
    scale_size_manual(values = rep(1, 3))
  
  return(p1) 
}

# Function to extract study area elevation map w/ point locations
plot_shp_elevation_map <- function(shp, loc_points, country_shp, zoom = 7, buffer_km = 20, stats=FALSE, new_points=NULL) {
  
  # 1. Convert sf_points to sf sf_points ----
  cities <- st_as_sf(loc_points, coords = c("lon", "lat"), crs = 4326)
  cities <- cbind(cities, st_coordinates(cities))
  names(cities)[names(cities) %in% c("X","Y")] <- c("lon","lat")
  
  
  # 2. Bounding box with buffer ----
  bbox <- st_as_sfc(st_bbox(shp)) |> st_buffer(buffer_km * 1000)
  
  # 3. Elevation and terrain ----
  elev_raster <- raster("2_Active_Data/204_misc/ind_envi_s3_AWSOpenData_elevation_raster.tif")
  slope <- terrain(elev_raster, opt = 'slope')
  aspect <- terrain(elev_raster, opt = 'aspect')
  hillshade <- hillShade(slope, aspect, angle = 40, direction = 315)
  
  # 4. Mask areas outside shp ----
  outside_shp <- st_difference(st_buffer(bbox, 10000), st_union(shp))
  
  # 5. Convert rasters to data frames ----
  elev_df <- as.data.frame(elev_raster, xy = TRUE)
  names(elev_df) <- c("x", "y", "elevation")
  hill_df <- as.data.frame(hillshade, xy = TRUE)
  names(hill_df) <- c("x", "y", "hill")
  
  # 6. Color palette ----
  arcgis_topo_colors <- c("#5A7D24","#9DBF68","#F3E99F","#D2B48C","#8B7A55","#FFFFFF")
  
  # 7. Main elevation map ----
  main_map <- ggplot() +
    geom_raster(data = hill_df, aes(x, y, fill = hill)) +
    scale_fill_gradient(low = "grey20", high = "grey90", guide = "none") +
    new_scale_fill() +
    geom_raster(data = elev_df, aes(x, y, fill = elevation), alpha = 0.6) +
    scale_fill_gradientn(colours = arcgis_topo_colors, name = "Elevation (m)") +
    geom_sf(data = outside_shp, fill = "grey70", alpha = 0.6) +
    geom_sf(data = country_shp, fill = NA, color = "grey50", size = 0.4) +
    geom_sf(data = shp, fill = NA, color = "black", size = 0.5) +
    geom_sf(data = cities, color = "darkred", size = 2) +
    coord_sf(xlim = range(elev_df$x), ylim = range(elev_df$y), expand = FALSE) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(plot.background = element_rect(fill = 'white', color = NA),
          panel.background = element_rect(fill = 'white', color = NA))
  
  if (stats==TRUE){
    main_map <- main_map +
      geom_label(data = cities[!complete.cases(st_set_geometry(cities, NULL)), ], aes(x = lon, y = lat, label = location),
                 size = 3, color = "darkred", fill = "white", fontface = "bold",
                 label.padding = unit(0.15, "lines"), label.r = unit(0.15, "lines"), 
                 label.size = 0,   nudge_y = 0.1) +
      geom_label_repel(
        data = cities[complete.cases(st_set_geometry(cities, NULL)), ],           # your locations + stats
        aes(x = lon, y = lat,
            label = paste0(location, 
                           "\nMean Temp.: ", t2m, 
                           "\nMin Temp.: ", mn2t, 
                           "\nMax Temp.: ", mx2t, 
                           "\nRainfall: ", tp, 
                           "\nWindspeed: ", fg10, 
                           "\nRel. Humidity: ", rh_mean, 
                           "\nPM2.5: ", pm25 
            )),
        box.padding = 0.5,               # spacing between boxes
        point.padding = 0.5,             # spacing from the point
        segment.color = "darkred",        # line connecting label to point
        segment.size = 0.5,
        force = 2,
        max.overlaps = 12,
        fill='white',
        color='darkred'
      )
  } else {
    
    main_map <- main_map +
      geom_label(data = cities, aes(x = lon, y = lat, label = location),
                 size = 3, color = "darkred", fill = "white", fontface = "bold",
                 label.padding = unit(0.15, "lines"), label.r = unit(0.15, "lines"), 
                 label.size = 0,   nudge_y = 0.1)
  }
  
  if (!is.null(new_points)){
    
    # Add points and labels
    main_map <- main_map + 
      geom_point(data=new_points, aes(x=lon, y=lat), color='black', shape=17, size=2.5) +
      geom_label(data=new_points, aes(x = lon, y = lat, label=name),
                 size = 3, color = "black", fill = "white", fontface = "bold",
                 label.padding = unit(0.15, "lines"), label.r = unit(0.15, "lines"), 
                 label.size = 0,   nudge_y = -0.1)
    
  }
  
  # 8. Inset map ----
  world <- ne_countries(scale = "medium", returnclass = "sf")
  oceans <- ne_download(scale = "medium", type = "ocean", category = "physical", returnclass = "sf")
  ind_bbox <- st_bbox(country_shp)
  inset_map <- ggplot() +
    geom_sf(data = oceans, fill = "lightblue") +
    geom_sf(data = world, fill = "white", color = "grey70") +
    geom_sf(data = country_shp, fill = "grey60", color = "grey60", size = 0.3) +
    geom_sf(data = shp, fill = "darkred", color = "darkred", size = 0.4) +
    coord_sf(xlim = c(ind_bbox$xmin - 2, ind_bbox$xmax + 2),
             ylim = c(ind_bbox$ymin - 2, ind_bbox$ymax + 2), expand = FALSE) +
    theme_void() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.6),
          plot.background = element_rect(fill = 'white', color = NA))
  
  # 9. Combine main map + inset ----
  main_map +
    inset_element(inset_map, left = 0.75, bottom = 0.75, right = 1.01, top = 1, align_to = "panel")
  
  
}

# Function to summarize ANOVA
anova_summary <- function(dep_var, factor_var, data){
  formula <- as.formula(paste(dep_var, '~', factor_var))
  fit <- aov(formula, data = data)
  sum_fit <- summary(fit)[[1]]
  df1 <- sum_fit[1, 'Df']
  df2 <- sum_fit[2, 'Df']
  F_val <- sum_fit[1, 'F value']
  p_val <- sum_fit[1, 'Pr(>F)']
  
  data.frame(
    Variable = dep_var,
    Factor = factor_var,
    Df_between = df1,
    Df_within = df2,
    F = round(F_val, 3),
    p = round(p_val, 3)
  )
}

# Function to summarize ANOVA interaction tests
anova_interaction_summary <- function(dep_var, factor1, factor2, data){
  
  if (factor1 == factor2) {
    return(NULL)
  }
  
  formula <- as.formula(paste(dep_var, '~', factor1, '*', factor2))
  fit <- aov(formula, data = data)
  sum_fit <- summary(fit)[[1]]
  
  results <- data.frame(
    Variable = dep_var,
    Term = rownames(sum_fit),
    Df = sum_fit[,'Df'],
    F = round(sum_fit[,'F value'], 3),
    p = round(sum_fit[,'Pr(>F)'], 3))
  
  # results <- results[-nrow(results), ]
  # results <- results[3, ]
  results <- results[grepl(':', results$Term), ]
  
  
  return(results)
}

####  MAIN CODE ####

# Load spatial data
ind <- st_read('2_Active_Data/201_admn/ind_admn_ad0_py_s1_GADM_pp_india.shp')
hp <- st_read('2_Active_Data/201_admn/ind_admn_ad3_py_s1_GADM_pp_Himachal_Pradesh.shp')
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

# Patient data format (!! not processed for NAs !!)
load('2_Active_Data/203_comb/ind_5km_BigThaw-HISP-ALL_2015-2023_daily_hosp_all.Rdata')
load('2_Active_Data/203_comb/ind_5km_BigThaw-HISP-ALL_2015-2023_daily_hosp_wide.Rdata')

# Remove COVID years for now:
startdate <- as.Date('2015-02-14')
enddate <- as.Date('2020-02-14')  
df_wide <- df_wide[df_wide$date >= startdate & df_wide$date <= enddate,]

# Checks
check <- df_wide[, !grepl('lag|_c|_p|IMD', names(df_wide))]
names(check)

# Expand cols of interest
h_cols <- c('closed', 'total_patients', 'location', h_cols[-1])
w_cols <- c(w_cols[-7],
            'mx2t_IMD_b', 'mn2t_IMD_b', 'tp_IMD_b', 'fg10_IMD_b', 
            'rh_low_IMD_b', 'rh_hig_IMD_b', 'pm25_IMD_b',
            'mx2t_p90', 'mn2t_p10', 'tp_p90', 'fg10_p90', 
            'rh_mean_p10', 'rh_mean_p90', 'pm25_p90')

location_levels <- unique(df_wide$location)

# GENERAL CHARACTERISTICS ####

# TAB: Base statistics (base_desc, base_desc_loc, base_weat_loc) ####

# OVERALL
base_desc <- get_desc(df_wide, h_cols, w_cols)
base_desc <- replace_exact_matches(base_desc, nm)

# SITE SPECIFIC
locations <- as.character(location_levels)
base_desc_loc <- NULL

for (l in locations) {
  sub <- df_wide %>% dplyr::filter(location == l)
  res <- get_desc(sub, h_cols, w_cols)
  names(res)[-1] <- paste0(names(res)[-1], ' (', l, ')')
  
  if (is.null(base_desc_loc)) {
    base_desc_loc <- res
  } else {
    base_desc_loc <- cbind(base_desc_loc, res[-1])  # drop first column if it's repeated
  }
}
base_desc_loc <- replace_exact_matches(base_desc_loc, nm)


# TAB: Table of weather averages & frequencies (weat_loc_tab) ####

weat_av <- base_desc_loc[33:39, grepl('Variable|Average', names(base_desc_loc))]
names(weat_av) <- c('Variable', as.character(location_levels))
weat_freq <- base_desc_loc[40:nrow(base_desc_loc), grepl('Variable|Frequency', names(base_desc_loc))]
names(weat_freq) <- c('Variable', as.character(location_levels))
weat_loc_tab <- rbind(weat_av, weat_freq)

# TAB: Table of population averages & frequencies (pop_loc_tab) ####
pop_av <- base_desc_loc[1:12, grepl('Variable|Average', names(base_desc_loc))]
names(pop_av) <- c('Variable', as.character(location_levels))
pop_av$Variable <- paste0(pop_av$Variable, ' (Average)')
pop_freq <- base_desc_loc[1:12, grepl('Variable|Frequency', names(base_desc_loc))]
names(pop_freq) <- c('Variable', as.character(location_levels))
pop_av$Variable <- paste0(pop_av$Variable, ' (Frequency)')
pop_loc_tab <- rbind(pop_freq, pop_av)

# TAB: Table of population averages & frequencies (ward_loc_tab) ####
ward_av <- base_desc_loc[13:31, grepl('Variable|Average', names(base_desc_loc))]
names(ward_av) <- c('Variable', as.character(location_levels))
ward_av$Variable <- paste0(ward_av$Variable, ' (Average)')
rownames(ward_av) <- NULL
ward_freq <- base_desc_loc[14:31, grepl('Variable|Frequency', names(base_desc_loc))]
names(ward_freq) <- c('Variable', as.character(location_levels))
ward_freq$Variable <- paste0(ward_freq$Variable, ' (Frequency)')
rownames(ward_freq) <- NULL
ward_loc_tab <- rbind(ward_freq, ward_av)

# FIG: Histograms total_patients, weather (hist_tota, hist_tota_loc, hist_weat) ####
sub1 <- df_wide %>% ungroup() %>% dplyr::select(total_patients, all_of(w_cols[!grepl('_p|IMD', w_cols)]))
skm1 <- skim(sub1) %>% rename(min = numeric.p0, Q1 = numeric.p25, 
                            median = numeric.p50, Q3 = numeric.p75, 
                            max = numeric.p100) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

colnames(skm1) <- gsub('numeric.', '', colnames(skm1))
skm1 <- as.data.frame(skm1)

# HISTOGRAMS
# Total patients only (individual)
hist_tota <- create_histogram(sub1, 'total_patients', skm1, color=mypalette[2])
(hist_tota <- hist_tota +
  ggtitle('All total patient distribution') +
  xlab('Total patients') +           # Change Y-axis label
  theme(
    legend.position = 'none',        # Remove legend
    plot.title = element_text(hjust = 0.5)        # Center the title
  ))



# Total patients only (by location)
sub2 <- df_wide %>% dplyr::select(total_patients, all_of(w_cols[!grepl('_p|IMD', w_cols)]))

plots_by_location <- lapply(location_levels, function(loc) {
  # Subset data for this location
  sub_loc <- subset(sub2, location == loc)
  skm_loc <- skim(sub_loc) %>%
    rename(min = numeric.p0, Q1 = numeric.p25, 
           median = numeric.p50, Q3 = numeric.p75, 
           max = numeric.p100) %>%
    mutate(across(where(is.numeric), ~ round(., 0)))
  colnames(skm_loc) <- gsub('numeric.', '', colnames(skm_loc))
  skm_loc <- as.data.frame(skm_loc)
  
  
  p <- create_histogram(sub_loc, 'total_patients', skm_loc, color=mypalette[2])
  p <- p +
    ggtitle(loc) +
    xlab('Total patients') +
    theme(
      legend.position = 'none',
      plot.title = element_text(hjust = 0.5)
    )
  
  return(p)
})
(hist_tota_loc <- wrap_plots(plots_by_location, ncol = 3))

# Weather only
cols_to_plot <- colnames(sub1[, -c(1, 2)])
plots <- lapply(cols_to_plot, function(col) create_histogram(sub1, col, skm1))
(hist_weat <- wrap_plots(plots, nrow=2))

# # Density curves by season
# labels_map <- c('Winter (Jan-Feb)' = 'Jan-Feb',
#                 'Pre-Monsoon (Mar-May)' = 'Mar-May',
#                 'Monsoon (Jun-Sep)' = 'Jun-Sep',
#                 'Post-Monsoon (Oct-Dec)' = 'Oct-Dec')
# 
# plots <- lapply(cols_to_plot, function(c) {
#   ggplot(df_wide, aes_string(x = c, color = 'season', fill = 'season')) +
#     geom_density(alpha = 0.3) +
#     theme_minimal() +
#     labs(
#       title = c,
#       x = c,
#       y = 'Density',
#       fill = 'Season',
#       color = 'Season'
#     ) +
#     scale_fill_manual(values=c(mypalette[1], mypalette[4], mypalette[2], mypalette[3]), labels=labels_map) +
#     scale_color_manual(values=c(mypalette[1], mypalette[4], mypalette[2], mypalette[3]), labels=labels_map)
# })
# (weat_dens <- wrap_plots(plots, ncol = 3))

# FIG: Bar charts Wards, age, gender (cat_ward, cat_ward_loc, cat_age, cat_age_loc, cat_sex, cat_sex_loc) ####
# Use OG data to easily analyse groups and avoid col issues

# FIG: WARD --
# CHECK OVERALL
df_all$ward_name <- factor(df_all$ward_name, levels = names(sort(table(df_all$ward_name), decreasing = TRUE)))
(cat_ward <- ggplot(df_all, aes(x = ward_name, fill = ward_name)) +
  geom_bar(color = 'black') + 
  scale_fill_manual(values = mypalette) +
  ggtitle('All Ward Distribution') +  
  xlab('Ward') + 
  ylab('Count') + 
  scale_y_continuous(labels = scales::comma) +
  scale_x_discrete(labels = function(x) {
    # If x is in nm, replace with nm value, else keep original
    sapply(x, function(lbl) {
      if (lbl %in% names(nm)) {
        nm[[lbl]]
      } else {
        lbl
      }})}) +
  theme_minimal() +   
  theme(
    legend.position = 'none', 
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ))

# CHECK BY LOCATION
colours <- mypalette
names(colours) <- ward_levels # keep colors static across wards
plots_by_location <- lapply(location_levels, function(loc) {
  df_loc <- subset(df_all, location == loc)
  
  # Order wards by frequency within this location
  df_loc$ward_name <- factor(df_loc$ward_name, 
                             levels = names(sort(table(df_loc$ward_name), decreasing = TRUE)))
  
  # Create plot
  p <- ggplot(df_loc, aes(x = ward_name, fill = ward_name)) +
    geom_bar(color = 'black') + 
    scale_fill_manual(values = colours) +
    ggtitle(loc) +  
    xlab('Ward') + 
    ylab('Count') + 
    scale_y_continuous(limits = c(0, 750000), labels = scales::comma) +  # Remove scientific notation (comma format)
    scale_x_discrete(labels = function(x) {
      sapply(x, function(lbl) if (lbl %in% names(nm)) nm[[lbl]] else lbl)
    }) +
    theme_minimal() +   
    theme(
      legend.position = 'none', 
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  return(p)
})
(cat_ward_loc <- wrap_plots(plots_by_location, ncol = 3))

# REMOVE wards with daily average < 50

# FIG: AGE--
# CHECK FREQUENCY
(cat_age <- ggplot(df_all, aes(x = age_group, fill = age_group)) +
            geom_bar(color = 'black') +  # Bar color
            ggtitle('All Age Distribution') +  # Title
            xlab('Age Groups') +  # X-axis label
            ylab('Count') +   # Y-axis label
            scale_x_discrete(labels = function(x) {
              # If x is in nm, replace with nm value, else keep original
              sapply(x, function(lbl) {
                if (lbl %in% names(nm)) {
                  nm[[lbl]]
                } else {
                  lbl
                }})}) +
            scale_fill_manual(values = mypalette) +
            scale_y_continuous(labels = scales::comma) +  # Remove scientific notation (comma format)
            theme_minimal() +   # Minimal theme for clean appearance
            theme(legend.position = 'none',
                  plot.title = element_text(hjust = 0.5)))

# CHECK BY LOCATION
colours <- mypalette
names(colours) <- age_levels # keep colors static across wards
plots_by_location <- lapply(location_levels, function(loc) {
  df_loc <- subset(df_all, location == loc)
  
  # Create plot
  p <- ggplot(df_loc, aes(x = age_group, fill = age_group)) +
        geom_bar(color = 'black') +  # Bar color
        ggtitle(loc) +  # Title
        xlab('Age Groups') +  # X-axis label
        ylab('Count') +   # Y-axis label
        scale_x_discrete(labels = function(x) {
          # If x is in nm, replace with nm value, else keep original
          sapply(x, function(lbl) {
            if (lbl %in% names(nm)) {
              nm[[lbl]]
            } else {
              lbl
            }})}) +
        scale_fill_manual(values = mypalette) +
        scale_y_continuous(limits = c(0, 1300000), labels = scales::comma) +  # Remove scientific notation (comma format)
        theme_minimal() +   
        theme(
          legend.position = 'none', 
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)
        )
  
  
  return(p)
})
(cat_age_loc <- wrap_plots(plots_by_location, ncol = 3))

# FIG: Gender--
# CHECK FREQUENCY
(cat_sex <- ggplot(df_all, aes(x = gender, fill = gender)) +
            geom_bar(color = 'black') +  # Bar color
            ggtitle('All Gender Distribution') +  # Title
            xlab('Gender') +  # X-axis label
            ylab('Count') +   # Y-axis label
            scale_fill_manual(values = mypalette) +  # Manual colors for Male and Female
            scale_y_continuous(labels = scales::comma) +  # Remove scientific notation (comma format)
            theme_minimal() +   # Minimal theme for clean appearance
            theme(legend.position = 'none',
                  plot.title = element_text(hjust = 0.5)))


# CHECK BY LOCATION
plots_by_location <- lapply(location_levels, function(loc) {
  df_loc <- subset(df_all, location == loc)
  
  # Create plot
  p <- ggplot(df_loc, aes(x = gender, fill = gender)) +
    geom_bar(color = 'black') +  # Bar color
    ggtitle(loc) +  # Title
    xlab('Gender') +  # X-axis label
    ylab('Count') +   # Y-axis label
    scale_fill_manual(values = mypalette) +  # Manual colors for Male and Female
    scale_y_continuous(limits = c(0, 1300000), labels = scales::comma) +  # Remove scientific notation (comma format)
    theme_minimal() +   # Minimal theme for clean appearance
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5))
  
  return(p)
})
(cat_sex_loc <- wrap_plots(plots_by_location, ncol = 3))

# TEMPORAL CHARACTERISTICS ####

clean <- df_wide[!is.na(df_wide$total_patients), ]

# FIG: 2-month moving average of  WEATHER across locations (ts_weat, ts_loc, ts_loc_weat) ####
w_vars <- w_cols[1:(which(w_cols == "pm25"))]
plots <- list()
plots_loc <- list()
for (v in w_vars){
  
  variable <- xts(clean[[v]],  order.by=clean$date)
  y_label <- nm[v]
  
  ts_plot <- plot_moving_avg(variable, all.cases=FALSE, MA=60, y.lab=y_label, colour=mypalette[1],
                             plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate)
  plots[[v]] <- ts_plot
  
  plots_by_location <- lapply(location_levels, function(loc) {
    df_loc <- subset(clean, location == loc)
    variable <- xts(df_loc[[v]],  order.by=df_loc$date)
    p <- plot_moving_avg(variable, all.cases=FALSE, MA=60, y.lab=y_label, colour=mypalette[1],
                         plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate)
    
    return(p)
  })
  names(plots_by_location) <- location_levels
  plots_loc[[v]] <- plots_by_location
}
(ts_weat <- wrap_plots(plots, ncol = 4) + 
    plot_annotation(title = 'Weather time series, across all locations', theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16))))

# # CHECK ALL BY LOCATION
# ts_weat_loc <- list()
# for (i in location_levels){
# 
#   (ts <- wrap_plots(lapply(plots_loc, `[[`, i), ncol = 4) +
#      plot_annotation(title = i, theme = theme(
#        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))))
#   ts_weat_loc[[i]] <- ts
# 
# }
# ts_weat_loc

# CHECK LOCATIONS BY WEATHER
ts_loc_weat <- list()
for (v in w_vars){
  
  weat_plots <- plots_loc[[v]]
  weat_plots_named <- Map(function(p, nm) {
    p + ggtitle(nm) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  }, weat_plots, names(weat_plots))
  title_label <- nm[v]
  (ts <- wrap_plots(weat_plots_named, ncol = 3)+ 
      plot_annotation(title = title_label, theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))))
  ts_loc_weat[[v]] <- ts
}
ts_loc_weat

# FIG: 2-month moving average of  POPULATION HEALTH across locations (ts_heal_pop, ts_pop_loc, ts_pop_loc_heal) #####
clean <- df_wide[!is.na(df_wide$total_patients), ]
vars <- setdiff(h_cols, c("closed", "location"))

# Identify variables of interest
pop_vars <- vars[1:(which(vars == "over_65"))]

plots <- list()
plots_loc <- list()
for (v in pop_vars){
  
  variable <- xts(clean[[v]],  order.by=clean$date)
  y_label <- nm[v]
  ts_plot <- plot_moving_avg(variable, all.cases=FALSE, MA=60, y.lab=y_label, colour=mypalette[2],
                            plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate)
  plots[[v]] <- ts_plot
  plots_by_location <- lapply(location_levels, function(loc) {
    df_loc <- subset(clean, location == loc)
    variable <- xts(df_loc[[v]],  order.by=df_loc$date)
    p <- plot_moving_avg(variable, all.cases=FALSE, MA=60, y.lab=y_label, colour=mypalette[2],
                         plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate)
    
    return(p)
  })
  names(plots_by_location) <- location_levels
  plots_loc[[v]] <- plots_by_location
}
(ts_heal_pop <- wrap_plots(plots, ncol = 4) + 
  plot_annotation(title = 'Population hospital time series, across all locations', theme = theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16))))

# CHECK ALL BY LOCATION
# ts_pop_loc <- list()
# for (i in location_levels){
#   
#   (ts <- wrap_plots(lapply(plots_loc, `[[`, i), ncol = 4) + 
#      plot_annotation(title = i, theme = theme(
#        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))))
#   ts_pop_loc[[i]] <- ts
#   
# }
# ts_pop_loc

# CHECK LOCATIONS BY POP GROUP
ts_pop_loc_heal <- list()
for (v in pop_vars){
  
  weat_plots <- plots_loc[[v]]
  weat_plots_named <- Map(function(p, nm) {
    p + ggtitle(nm) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  }, weat_plots, names(weat_plots))
  title_label <- nm[v]
  (ts <- wrap_plots(weat_plots_named, ncol = 3)+ 
      plot_annotation(title = title_label, theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))))
  ts_pop_loc_heal[[v]] <- ts
}
ts_pop_loc_heal



# FIG: 2-month moving average of  WARD HEALTH across locations (ts_heal_ward, ts_ward_loc, ts_ward_loc_heal) #####
ward_var <- vars[(which(vars == "general")):which(vars == "dermatology")]

# Identify variables of interest
pop_vars <- vars[1:(which(vars == "over_65"))]

plots <- list()
plots_loc <- list()
for (v in ward_var){
  
  variable <- xts(clean[[v]],  order.by=clean$date)
  y_label <- nm[v]
  ts_plot <- plot_moving_avg(variable, all.cases=FALSE, MA=60, y.lab=y_label, colour=mypalette[2],
                             plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate)
  plots[[v]] <- ts_plot
  plots_by_location <- lapply(location_levels, function(loc) {
    df_loc <- subset(clean, location == loc)
    variable <- xts(df_loc[[v]],  order.by=df_loc$date)
    p <- plot_moving_avg(variable, all.cases=FALSE, MA=60, y.lab=y_label, colour=mypalette[2],
                         plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate)
    
    return(p)
  })
  names(plots_by_location) <- location_levels
  plots_loc[[v]] <- plots_by_location
}
(ts_heal_ward <- wrap_plots(plots, ncol = 4) + 
    plot_annotation(title = 'Population hospital time series, across all locations', theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16))))

# # CHECK ALL BY LOCATION
# ts_ward_loc <- list()
# for (i in location_levels){
#   
#   (ts <- wrap_plots(lapply(plots_loc, `[[`, i), ncol = 4) + 
#      plot_annotation(title = i, theme = theme(
#        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))))
#   ts_ward_loc[[i]] <- ts
#   
# }
# ts_ward_loc

# CHECK LOCATIONS BY WARD
ts_ward_loc_heal <- list()
for (v in ward_var){
  
  weat_plots <- plots_loc[[v]]
  weat_plots_named <- Map(function(p, nm) {
    p + ggtitle(nm) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  }, weat_plots, names(weat_plots))
  title_label <- nm[v]
  (ts <- wrap_plots(weat_plots_named, ncol = 3)+ 
      plot_annotation(title = title_label, theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))))
  ts_ward_loc_heal[[v]] <- ts
}
ts_ward_loc_heal



# Test for temporal autocorrelation ####

# SPATIAL CHARACTERISTICS ####


# FIG: Spatial distribution of weather variables (study_area) ####

# !! Merge stats !!
map_stat <- base_desc_loc[33:39, !grepl('Frequency', names(base_desc_loc))]
rownames(map_stat) <- map_stat[,1]
map_stat <- map_stat[, !grepl("Range", names(map_stat))]
names(map_stat) <- gsub("Average|\\(|\\)|\\ ", "", names(map_stat))
map_stat <- as.data.frame(t(map_stat[,-1]))
map_stat$location <- rownames(map_stat)
rownames(map_stat) <- NULL
map_stat <- map_stat[, c(ncol(map_stat), 1:(ncol(map_stat)-1))]
rev_nm <- setNames(names(nm), nm)
names(map_stat)[-1] <- rev_nm[names(map_stat)[-1]]
map_stat$location[map_stat$location == "ReckongPeo"] <- "Reckong Peo"
loc_stats <- merge(loc_points, map_stat, by = "location", all.x = TRUE)
# Filter to exisitng locations
loc_stats <- loc_stats %>% filter(!is.na(t2m)) %>% filter(location != 'Bilaspur')


# New points
new_points <- data.frame(
  lon = c(76.92898, 76.27332772176327, 76.76117),
  lat = c(31.71256, 32.102957165868474, 31.34124),
  name = c("Mandi", "Kangra", 'Bilaspur')
)

# Map
# save_elev(hp, loc_stats, ind)
(study_area <- plot_shp_elevation_map(hp, loc_stats, ind, stats=FALSE, new_points=new_points))
# (study_area <- plot_shp_elevation_map(hp, loc_stats, ind, stats=TRUE))

# Test for spatial autocorrelation ####


# INTERACTION TESTS ####

# FIG: Non-parametric spearman correlation tests from raw data (cor_matrix) ####
sub <- df_wide %>% ungroup() %>% dplyr::select(total_patients, all_of(w_cols[!grepl('_p|IMD', w_cols)]))
cor_matrix <- cor(sub, use = 'complete.obs', method = 'spearman')
melted_corr <- melt(cor_matrix)
cor_matrix <- ggplot(data = melted_corr, aes(x=Var1, y=Var2, fill=value)) +
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
  )

# TAB: IMD ANOVA tests (IMD_anova, IMD_interaction_tests) ####
sub_IMD <- df_wide %>% dplyr::select(total_patients, season, all_of(w_cols[grepl('IMD', w_cols)]))
test_cols <- names(sub_IMD[,-c(1, 2)])  # everything except total_patients
(IMD_anova <- do.call(rbind, lapply(c('season', test_cols), function(f) anova_summary('total_patients', f, sub_IMD))))

IMD_anova$Factor <- gsub('_b', '', IMD_anova$Factor)
IMD_anova <- replace_exact_matches(IMD_anova, nm)
names(IMD_anova) <- c('Variable', 'Factor', 'Df Between','Df Within','F','p')

# Interaction tests
# Season
(IMD_inter_seas <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'season', f, sub_IMD))))
# Humidity
(IMD_inter_hum_h <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'rh_hig_IMD_b', f, sub_IMD))))
(IMD_inter_hum_l <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'rh_low_IMD_b', f, sub_IMD))))
# Pollution
(IMD_inter_pm25 <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'pm25_IMD_b', f, sub_IMD))))
# Wind
(IMD_inter_wind <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'fg10_IMD_b', f, sub_IMD))))

IMD_interaction_tests <- rbind(IMD_inter_seas, IMD_inter_hum_h, IMD_inter_hum_l, IMD_inter_pm25, IMD_inter_wind)
row.names(IMD_interaction_tests) <- NULL


# TAB: CLIM ANOVA tests (CLIM_anova, CLIM_interaction_tests) ####
sub_clim <- df_wide %>% dplyr::select(total_patients, season, all_of(w_cols[grepl('_p', w_cols)]))
test_cols <- names(sub_clim[,-c(1, 2)])  # everything except total_patients
(CLIM_anova <- do.call(rbind, lapply(c('season', test_cols), function(f) anova_summary('total_patients', f, sub_clim))))
test_cols <- c('season', 'tp_p90', 'fg10_p90', 'rh_mean_p10', 'rh_mean_p90')

CLIM_anova$Factor <- gsub('_b', '', CLIM_anova$Factor)
CLIM_anova <- replace_exact_matches(CLIM_anova, nm)
names(CLIM_anova) <- c('Variable', 'Factor', 'Df Between','Df Within','F','p')

# Interaction tests
# Season
(clim_inter_seas <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'season', f, sub_clim))))
# Humidity
(clim_inter_hum_h <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'rh_mean_p90', f, sub_clim))))
(clim_inter_hum_l <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'rh_mean_p10', f, sub_clim))))
# Pollution
(clim_inter_pm25 <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'pm25_p90', f, sub_clim))))
# Wind
(clim_inter_wind <- do.call(rbind, lapply(test_cols, function(f) anova_interaction_summary('total_patients', 'pm25_p90', f, sub_clim))))

CLIM_interaction_tests <- rbind(clim_inter_seas, clim_inter_hum_h, clim_inter_hum_l, clim_inter_pm25, clim_inter_wind)
row.names(CLIM_interaction_tests) <- NULL


##### EXPORT ####
w <- 6
h <- 6

# FIGURES
ggsave(paste0(path_report, '- figures/desc_hist_weat.png'), hist_weat, width = 12, height = h, dpi = 300)
ggsave(paste0(path_report, '- figures/desc_hist_tota.png'), hist_tota, width = w, height = h, dpi = 300)
ggsave(paste0(path_report, '- figures/desc_hist_tota_loc.png'), hist_tota_loc, width = 12, height = h, dpi = 300)
ggsave(paste0(path_report, '- figures/desc_cat_ward.png'), cat_ward, width = 8, height = h, dpi = 300)
ggsave(paste0(path_report, '- figures/desc_cat_ward_loc.png'), cat_ward_loc, width = 12, height = h, dpi = 300)
ggsave(paste0(path_report, '- figures/desc_cat_sex.png'), cat_sex, width = w, height = h, dpi = 300)
ggsave(paste0(path_report, '- figures/desc_cat_age.png'), cat_age, width = w, height = h, dpi = 300)
ggsave(paste0(path_report, '- figures/desc_con_age.png'), con_age, width = w, height = h, dpi = 300)
ggsave(paste0(path_report, '- figures/desc_cor_matrix.png'), cor_matrix, width = 9, height = h, dpi = 300)
ggsave(paste0(path_report, '- figures/desc_study_area.png'), study_area, width = 7, height = 7, dpi = 300)
ggsave(paste0('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts/2_DLNM/- Reporting/- figures/desc_study_area.png'), study_area, width = 7, height = 7, dpi = 300)



# TABLES
write.csv(base_desc, paste0(path_report, '- tables/desc_base.csv'), row.names = FALSE)
write.csv(base_weat_loc, paste0(path_report, '- tables/desc_base_weat_loc.csv'), row.names = FALSE)
write.csv(IMD_anova, paste0(path_report, '- tables/desc_IMD_anova.csv'), row.names = FALSE)
write.csv(IMD_interaction_tests, paste0(path_report, '- tables/desc_IMD_interaction_tests.csv'), row.names = FALSE)
write.csv(CLIM_anova, paste0(path_report, '- tables/desc_CLIM_anova.csv'), row.names = FALSE)
write.csv(CLIM_interaction_tests, paste0(path_report, '- tables/desc_CLIM_interaction_tests.csv'), row.names = FALSE)


