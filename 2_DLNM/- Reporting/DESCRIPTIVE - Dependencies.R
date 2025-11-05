
#### PACKAGES ####

# List packages
packages <- c('tidyverse', 'terra', 'duckdb', 'tidyselect', 'sf', 'knitr', 'kableExtra',
              'ggplot2', 'data.table', 'tools', 'moments', 'grid',
              'patchwork', 'RColorBrewer', 'reshape2', 'zoo', 'ggnewscale',
              'skimr', 'DT', 'forecast', 'purrr', 'scales', 'xts', 'ggthemes',
              'elevatr', 'rnaturalearth', 'rnaturalearthdata', 'dplyr', 'tidyr', 
              'lubridate', 'dlnm', 'splines', 'mvmeta')

# library(MASS); library(ggfortify); library(gridExtra)

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

# Function to get base stats of df
get_desc <- function(df, w_cols, h_cols=NULL) {
  
  # Handling for location specific tables
  if (length(unique(df$location)) == 1){
    h_cols <- h_cols[h_cols != "location"]
    check <- c('F', 'under_5', 'ophthalmology', 't2m')
    labels <- c('GENDER', 'AGE', 'WARD', 'WEATHER CHARACTERISTICS')
  } else{
    check <- c(sort(unique(as.character(df$location)))[1], 'F', 'under_5', 'ophthalmology', 't2m')
    labels <- c('LOCATION', 'GENDER', 'AGE', 'WARD', 'WEATHER CHARACTERISTICS')
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
            select(Variable=location, total_visits, perc_visits, everything())
          
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
            round(sd_daily_visits)
          )
        ),
        Range = ifelse(
          is.na(skew_daily_visit),
          NA,
          paste0(min_daily_visit, ', ',max_daily_visit, ' (', skew_daily_visit, ')')
        )
      ) %>%
      dplyr::select(Variable, Frequency, Average, Range))
  
  
  # Remove frequency values for raw weather data
  base_desc$Frequency[base_desc$Variable %in%  gsub('(_IMD_b|_p90|_p10)$', '', w_cols)] <- NA
  
  
  # FORMATING
  # Insert headings at fixed positions
  base_desc <- rbind(
    tibble(Variable = 'HOSPITAL CHARACTERISTICS', Frequency = NA, Average = NA, Range = NA),
    base_desc
  )
  # Find positions to insert before
  positions <- which(base_desc$Variable %in% check)
  new_rows <- tibble(
    Variable = labels,
    Frequency = NA,
    Average = NA,
    Range = NA
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


# Function to extract study area elevation map w/ point locations
plot_shp_elevation_map <- function(shp, loc_points, country_shp, zoom = 7, buffer_km = 20, stats=FALSE) {
  
  # 1. Convert sf_points to sf sf_points ----
  cities <- st_as_sf(loc_points, coords = c("lon", "lat"), crs = 4326)
  cities <- cbind(cities, st_coordinates(cities))
  names(cities)[names(cities) %in% c("X","Y")] <- c("lon","lat")
  
  
  # 2. Bounding box with buffer ----
  bbox <- st_as_sfc(st_bbox(shp)) |> st_buffer(buffer_km * 1000)
  
  # 3. Elevation and terrain ----
  elev_raster <- get_elev_raster(locations = st_sf(bbox), z = zoom, clip = "bbox")
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
                           "\nRel. Humidity: ", rh, 
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

#### DEPENDENCIES ####

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

mypalette <- c('#1D98EB','#D74B3C','#A6B727','#FEC306','#EC800A','#E569BF','#8968CD','#52C5B0FF', '#910030FF', '#127804FF',
               '#A1C8F2','#F1A6A6','#C2D77D','#F6D974','#F9B34D','#F1C2D1','#B8A8D9','#52C5B080', '#91003080', '#12780480') 

w_cols <- c('t2m', 'mn2t','mx2t','tp', 'rh', 'pm25')
h_cols <- c('total_patients', 'location',
            'F','M','under_5','age5_14','age15_24','age25_64','over_65',
            'ophthalmology', 'ent', 'paediatrics',
            'gynaecology', 'orthopedics')

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
        rh='Relative humidity (%)', tp='Rainfall (mm)', fg10='Windspeed (km/hr)', pm25='Air pollution (µg/m³)',
        # IMD DEFINITIONS
        mx2t_IMD='IMD: Mx temp. >30\u00B0C', mn2t_IMD='IMD: Min temp. <0\u00B0C', 
        tp_IMD='IMD: Rainfall >15mm', fg10_IMD='IMD: Wind speed >40km/h', pm25_IMD='IMD: PM2.5 >60 µg/m³', 
        rh_low_IMD='IMD: Humidity <30%', rh_hig_IMD='IMD: Humidity >60%',
        # CLIM. DEFINITIONS
        mx2t_p90='Clim: Max temp. >90th', mn2t_p10='Clim: Min temp. <10th', 
        tp_p90='Clim: Rainfall >90th', fg10_p90='Clim: Wind speed >90th', pm25_p90='Clim: PM2.5 >90th', 
        rh_p10 ='Clim: Humidity <10th', rh_p90='Clim: Humidity >90th') 


