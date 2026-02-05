### ------------------------------------------------------------------------ ### 
### NHS Exploratory Data Analysis
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list = ls()); graphics.off(); gc()
setwd('C:/Users/claer14/Documents/GitHub/PhD-Climate-Health-Forecasts')

# Set data path
user <- 'claer14'
path1 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
path2 <- paste0('C:/Users/', user, '/Documents/GitHub/PhD-Climate-Health-Forecasts/1_Data/NHS_Cambrigeshire/')


#### DEPENDENCIES ####

library(sf); library(tidyverse); library(e1071); library(patchwork); 
library(xts); library(ggthemes); library(skimr); library(tmap)

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
  
  # Handling for provider_name specific tables
  if (length(unique(df$provider_name)) == 1){
    h_cols <- h_cols[h_cols != "provider_name"]
    check <- c('age_under_5', 'imd_q1', 'land_rural', 'treat_1')
    labels <- c('Age', 'Deprivation Quintile', 'Patient Area Classification', 'Treatment Severity')
  } else{
    check <- c(sort(unique(as.character(df$provider_name)))[1], 'age_under_5', 'imd_q1', 'land_rural', 'treat_1')
    labels <- c('Provider', 'Age', 'Deprivation Quintile', 'Patient Area Classification', 'Treatment Severity')
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
            group_by(provider_name) %>%
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
            dplyr::select(Variable=provider_name, total_visits, perc_visits, everything())
          
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
  # base_desc$Frequency[base_desc$Variable %in%  gsub('(_IMD_b|_p90|_p10)$', '', w_cols)] <- NA
  
  
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
    geom_line(aes(y = Centered.7.day.MA, color = 'Centered 7 day MA'), linewidth = 0.5) +
    scale_color_manual(name = '', values = c('Real data' = 'black', 'Centered 7 day MA' = colour)) +
    scale_y_continuous(n.breaks = 6) +
    scale_x_date(labels = scales::date_format('%b %y'), 
                 breaks = seq(min(df.f$Date), max(df.f$Date), by = '3 months')) +
    labs(
      title = y.lab,
      y=y.lab, 
      x='Date' 
    ) +
    theme_economist_white(gray_bg = FALSE, base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 10, hjust = 0.5, face = "bold", margin = margin(b = 5)),
      axis.text  = element_text(size = 8),
      axis.title = element_text(size = 8),
      axis.text.x  = element_text(angle = 45, vjust =1, hjust=1, margin = margin(t = 2)),
      axis.title.x = element_text(margin = margin(t = 1)),
      axis.text.y  = element_text(margin = margin(r = 2)),
      axis.title.y = element_text(margin = margin(r = 3)),
      plot.margin = margin(5, 5, 5, 5, unit = "pt")
    ) +
    scale_linetype_manual(values = rep('solid', 3)) +
    scale_size_manual(values = rep(1, 3))
  
  return(p1) 
}

# Function to extract study area elevation map w/ point provider_names
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
  
  bbox_map <- st_bbox(shp)
  
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
          panel.background = element_rect(fill = 'white', color = NA),
          plot.margin = margin(0, 0, 0, 0))
  
  if (stats==TRUE){
    main_map <- main_map +
      geom_label(data = cities[!complete.cases(st_set_geometry(cities, NULL)), ], aes(x = lon, y = lat, label = provider_name),
                 size = 3, color = "darkred", fill = "white", fontface = "bold",
                 label.padding = unit(0.15, "lines"), label.r = unit(0.15, "lines"), 
                 label.size = 0,   nudge_y = 0.1) +
      geom_label_repel(
        data = cities[complete.cases(st_set_geometry(cities, NULL)), ],           # your provider_names + stats
        aes(x = lon, y = lat,
            label = paste0(provider_name, 
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
      geom_label(data = cities, aes(x = lon, y = lat, label = provider_name),
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
    inset_element(inset_map, left = 0.75, bottom = 0.75, right = 1.015, top = 1, align_to = "panel")
  
  
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

mypalette <- c('#1D98EB','#D74B3C','#A6B727','#FEC306','#EC800A','#E569BF','#8968CD','#52C5B0FF', '#910030FF', '#127804FF',
               '#A1C8F2','#F1A6A6','#C2D77D','#F6D974','#F9B34D','#F1C2D1','#B8A8D9','#52C5B080', '#91003080', '#12780480') 

h_cols <- c('total_patients', 'provider_name',
            'age_under_5', 'age_5_14', 'age_15_24', 'age_25_64', 'age_over_65',
            'imd_q1', 'imd_q2', 'imd_q3', 'imd_q4', 'imd_q5', 
            'land_rural', 'land_urban')
w_cols <- c('treat_1', 'treat_2', 'treat_3', 'treat_4', 'treat_5')
nm <- c(total_patients='Total patients', 
        # Age
        age_under_5='Under 5', age_5_14='Ages 5-14', age_15_24='Ages15-24', age_25_64='Ages 25-64', age_over_65='Ages 65+',
        # Index of Multiple Deprivation
        imd_q1='Q1 (High)',imd_q2='Q2', imd_q3='Q3', imd_q4='Q4', imd_q5='Q5 (Low)',
        land_rural='Rural', land_urban='Urban',
        # Treatment severity
        treat_1='Category 1 (minor)', treat_2='Category 2', treat_3='Category 3',  treat_4='Category 4',  treat_5='Category 5 (severe)'
        )
        
#### MAIN CODE ####

# Patient data format (!! not processed for NAs !!)
load(paste0(path1, '201_heal/eng_adm4_NHS_2024_ae_visits_ALL_clean.Rdata'))
load(paste0(path1, '201_heal/eng_adm4_NHS_2024_ae_visits_ALL_clean_wide.Rdata'))
df_all <- df

# Load spatial data
country_sf <- st_read(paste0(path1, '/200_admn/eng_admn_ad1_py_s1_NHS_Regions_2024.shp'))
icb_sf <- st_read(paste0(path1, '/200_admn/eng_admn_ad3_py_s1_NHS_ICB_2023.shp'))
icb_sf <- icb_sf %>% filter(ICB23NM == 'NHS Cambridgeshire and Peterborough Integrated Care Board')
locs <- df_wide %>% dplyr::select(contains('provider')) %>% distinct() %>%
  st_as_sf(coords = c("provider_lon", "provider_lat"), crs = 27700, remove = FALSE)

# Remove COVID years for now:
startdate <- as.Date(min(df_wide$date))
enddate <- as.Date(max(df_wide$date))  
df_wide <- df_wide[df_wide$date >= startdate & df_wide$date <= enddate,]

# Set location levels
location_levels <- unique(df_wide$provider_name)


#### GENERAL CHARACTERISTICS ####
# TAB: Base statistics (base, base_loc) ####

# Check missing data
(missing_counts <- colSums(is.na(df_wide)))
(df_missing_any <- df_wide[rowSums(is.na(df_wide)) > 0, ])

# OVERALL
base <- get_desc(df_wide, h_cols, w_cols)
base <- replace_exact_matches(base, nm)

# SITE SPECIFIC
base_loc <- NULL
for (l in unique(df_wide$provider_name)) {
  sub <- df_wide %>% dplyr::filter(provider_name == l)
  res <- get_desc(sub, h_cols, w_cols)
  names(res)[-1] <- paste0(names(res)[-1], ' (', l, ')')
  
  if (is.null(base_loc)) {
    base_loc <- res
  } else {
    base_loc <- cbind(base_loc, res[-1])  # drop first column if it's repeated
  }
}
base_loc <- replace_exact_matches(base_loc, nm)


# FIG: Histograms total_patients (hist_tota) ####
sub1 <- df_wide %>% ungroup() %>% dplyr::select(total_patients)
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
    xlab('Total patients') +
    theme(
      legend.position = 'none',  
      plot.title = element_text(hjust = 0.5)
    ))
# Total patients only (by location)
sub2 <- df_wide %>% dplyr::select(total_patients)
plots_by_location <- lapply(location_levels, function(loc) {
  # Subset data for this location
  sub_loc <- subset(sub2, provider_name == loc)
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

# FIG: Bar charts age, deprivation, land, treatment (cat_age, cat_imd, cat_lan, cat_treat) ####

# FIG: AGE--
# CHECK FREQUENCY
(cat_age <- ggplot(df_all, aes(x = patient_age_group, fill = patient_age_group)) +
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
names(colours) <- levels(df$patient_age_group) # keep colors static across wards
plots_by_location <- lapply(location_levels, function(loc) {
  df_loc <- subset(df_all, provider_name == loc)
  
  # Create plot
  p <- ggplot(df_loc, aes(x = patient_age_group, fill = patient_age_group)) +
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
    scale_y_continuous(limits = c(0, 150000), labels = scales::comma) +  # Remove scientific notation (comma format)
    theme_minimal() +   
    theme(
      legend.position = 'none', 
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  
  return(p)
})
(cat_age_loc <- wrap_plots(plots_by_location, ncol = 3))

# FIG: IMD --
# CHECK FREQUENCY
(cat_imd <- ggplot(df_all, aes(x = patient_imd_quintile, fill = patient_imd_quintile)) +
    geom_bar(color = 'black') +  # Bar color
    ggtitle('All Deprivation Distribution') +  # Title
    xlab('Index of Multiple Deprivation') +  # X-axis label
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
names(colours) <- levels(df$patient_imd_quintile) # keep colors static across wards
plots_by_location <- lapply(location_levels, function(loc) {
  df_loc <- subset(df_all, provider_name == loc)
  
  # Create plot
  p <- ggplot(df_loc, aes(x = patient_imd_quintile, fill = patient_imd_quintile)) +
    geom_bar(color = 'black') +  # Bar color
    ggtitle(loc) +  # Title
    xlab('Index of Multiple Deprivation') +  # X-axis label
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
    scale_y_continuous(limits = c(0, 125000), labels = scales::comma) +  # Remove scientific notation (comma format)
    theme_minimal() +   
    theme(
      legend.position = 'none', 
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  
  return(p)
})
(cat_imd_loc <- wrap_plots(plots_by_location, ncol = 3))

# FIG: LAND --
# CHECK FREQUENCY
(cat_lan <- ggplot(df_all, aes(x = patient_area_class, fill = patient_area_class)) +
    geom_bar(color = 'black') +  # Bar color
    ggtitle('All land classification') +  # Title
    xlab('Land Class') +  # X-axis label
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
names(colours) <- levels(df$patient_area_class) # keep colors static across wards
plots_by_location <- lapply(location_levels, function(loc) {
  df_loc <- subset(df_all, provider_name == loc)
  
  # Create plot
  p <- ggplot(df_loc, aes(x = patient_area_class, fill = patient_area_class)) +
    geom_bar(color = 'black') +  # Bar color
    ggtitle(loc) +  # Title
    xlab('Land Class') +  # X-axis label
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
    scale_y_continuous(limits = c(0, 250000), labels = scales::comma) +  # Remove scientific notation (comma format)
    theme_minimal() +   
    theme(
      legend.position = 'none', 
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  
  return(p)
})
(cat_lan_loc <- wrap_plots(plots_by_location, ncol = 3))

# FIG: TREAT --
# CHECK FREQUENCY
(cat_treat <- ggplot(df_all, aes(x = treat_desc, fill = treat_desc)) +
    geom_bar(color = 'black') +  # Bar color
    ggtitle('All treatment classification') +  # Title
    xlab('Treatment Class') +  # X-axis label
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
names(colours) <- levels(df$treat_desc) # keep colors static across wards
plots_by_location <- lapply(location_levels, function(loc) {
  df_loc <- subset(df_all, provider_name == loc)
  
  # Create plot
  p <- ggplot(df_loc, aes(x = treat_desc, fill = treat_desc)) +
    geom_bar(color = 'black') +  # Bar color
    ggtitle(loc) +  # Title
    xlab('Treatment Class') +  # X-axis label
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
    scale_y_continuous(limits = c(0, 150000), labels = scales::comma) +  # Remove scientific notation (comma format)
    theme_minimal() +   
    theme(
      legend.position = 'none', 
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  
  return(p)
})
(cat_treat_loc <- wrap_plots(plots_by_location, ncol = 3))

####  TEMPORAL CHARACTERISTICS ####

# FIG: 2-month moving average of  POPULATION HEALTH across locations (ts_heal_pop, ts_pop_loc, ts_pop_loc_heal) #####
clean <- df_wide[!is.na(df_wide$total_patients), ]
vars <- setdiff(h_cols, c("provider_name"))

# Identify variables of interest
pop_vars <- c(vars[1:(which(vars == "age_over_65"))], 'land_rural', 'land_urban')

plots <- list()
plots_loc <- list()
for (v in pop_vars){
  
  variable <- xts(clean[[v]],  order.by=clean$date)
  y_label <- nm[v]
  ts_plot <- plot_moving_avg(variable, all.cases=FALSE, MA=60, y.lab=y_label, colour=mypalette[2],
                             plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate)
  plots[[v]] <- ts_plot
  plots_by_location <- lapply(location_levels, function(loc) {
    df_loc <- subset(clean, provider_name == loc)
    variable <- xts(df_loc[[v]],  order.by=df_loc$date)
    p <- plot_moving_avg(variable, all.cases=FALSE, MA=60, y.lab=y_label, colour=mypalette[2],
                         plt.start.date=startdate, est.date.index=index, estimation.date.end=enddate)
    
    return(p)
  })
  names(plots_by_location) <- location_levels
  plots_loc[[v]] <- plots_by_location
}
ts_heal_pop <- wrap_plots(plots, ncol = 4) + 
  plot_annotation(title = 'Population hospital time series, across all locations', theme = theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)))

# CHECK ALL BY LOCATION
ts_pop_loc <- list()
for (i in location_levels){

  (ts <- wrap_plots(lapply(plots_loc, `[[`, i), ncol = 4) +
     plot_annotation(title = i, theme = theme(
       plot.title = element_text(hjust = 0.5, face = "bold", size = 16))))
  ts_pop_loc[[i]] <- ts

}
ts_pop_loc$`Peterborough City Hospital`
ts_pop_loc$`Addenbrooke's Hospital`
ts_pop_loc$`Hinchingbrooke Hospital`


#### SPATIAL CHARACTERISTICS ####
# FIG: Spatial distribution of weather variables (study_area) ####


provider_unique <- df %>% distinct(provider_name, .keep_all = TRUE)
provider_sf <- st_as_sf(provider_unique, coords = c("provider_lon", "provider_lat"), crs = 4326)
practice_unique <- df %>% distinct(practice_name, .keep_all = TRUE)  %>% filter(!is.na(practice_lon), !is.na(practice_lat))
practice_sf <- st_as_sf(practice_unique, coords = c("practice_lon", "practice_lat"), crs = 4326)


# Make sure CRS matches
country_sf <- st_transform(country_sf, st_crs(icb_sf))

# Mask = everything EXCEPT the ICB
icb_mask <- st_difference(
  st_union(country_sf),
  st_union(icb_sf)
)

(main_map <-
  tm_basemap("CartoDB.PositronNoLabels") +
  
  # Grey transparent mask
  tm_shape(icb_mask, bbox = st_bbox(st_buffer(icb_sf, dist = 15000))) +
  tm_polygons(
    col = "grey80",
    alpha = 0.7,
    border.col = NA
  ) +
  
  # Providers
  tm_shape(provider_sf) +
  tm_dots(col = "darkred", size = 0.5
  ) +
  tm_text(
    "provider_name",
    size = 0.7,
    col = "darkred",
    bg.color = "white",
    bg.alpha = 0.9,
    bg.frame = TRUE,
    bg.frame.col = "darkred",
    bg.frame.lwd = 0.3,
    just = "left",
    ymod = 1,
    remove.overlap = TRUE
  ) +
  
  tm_layout(
    main.title = "Cambridgeshire & Peterborough ICB Providers",
    legend.outside = TRUE,
  ))


inset_map <-
  tm_shape(country_sf) +
  tm_polygons(col = "white", border.col = "white") +
  
  tm_shape(icb_sf) +
  tm_polygons(col = "darkred", border.col = "black") +
  
  tm_layout(
    frame = TRUE,
    bg.color = "lightblue"
  )
main_map
print(
  inset_map,
  vp = viewport(
    x = 0.849, y = 0.664,
    width = 0.25, height = 0.25,
    just = c("right", "bottom")
  )
)


#### EXPORT ####

# Figures
figs <- list('1_study_area'=list(p=study_area, w=8, h=7*0.8),
             '2A_hist_tota_loc'=list(p=hist_tota_loc, w=10, h=7),
             '2B_cat_ward_loc'=list(p=cat_ward_loc, w=10, h=7),
             '2C_ts_heal_ward'=list(p=ts_heal_ward, w=7, h=5),
             '2D_ts_emergency'=list(p=ts_emergency, w=7, h=5),
             '3A_hist_weat'=list(p=hist_weat,  w=10, h=7), 
             '3B_ts_mx2t'=list(p=ts_mx2t, w=7, h=5),
             '3C_ts_tp'=list(p=ts_tp, w=7, h=5),
             '3D_ts_fg10'=list(p=ts_fg10, w=7, h=5)
)

for (f in names(figs)){
  ggsave(paste0(path_git, '2_DLNM/- Reporting/figures/descriptive/desc_', f, '.png'), figs[[f]]$p, width = figs[[f]]$w, height = figs[[f]]$h, dpi = 300, bg='white')
}

# TABLES
tabs <- c('base', 'weat_loc_tab', 'IMD_anova', 'IMD_interaction_tests', 'CLIM_anova', 'CLIM_interaction_tests')
for (t in tabs){
  write.csv(get(t), paste0(path_git, '2_DLNM/- Reporting/- tables/desc_', t, '.csv'), row.names = FALSE)
  
}
