### ------------------------------------------------------------------------ ### 
### HP Hospital cleaning - Controls                                          ###
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list = ls()); graphics.off(); gc()
setwd('C:/Users/claer14/Documents/GitHub/hisp-climate-health-forecast/')


#### DEPENDENCIES ####
source('1_Data/DATA - Dependencies.R')
setwd(path_data)

# Function to extract and save elevtaion raster according to bounding area
save_elev <- function(shp, loc_points, zoom = 7, buffer_km = 20){
  
  # Convert sf_points to sf sf_points
  cities <- st_as_sf(loc_points, coords = c("lon", "lat"), crs = 4326)
  cities <- cbind(cities, st_coordinates(cities))
  names(cities)[names(cities) %in% c("X","Y")] <- c("lon","lat")
  
  # Bounding box with buffer
  bbox <- st_as_sfc(st_bbox(shp)) |> st_buffer(buffer_km * 1000)
  
  # Elevation and terrain
  elev_raster <- get_elev_raster(locations = st_sf(bbox), z = zoom, clip = "bbox")
  
  # 4. Export
  writeRaster(elev_raster, "2_Active_Data/205_envi/ind_envi_s3_AWSOpenData_elevation_raster.tif", overwrite = TRUE)
  
}

#### MAIN CODE ####

hp <- st_read('2_Active_Data/200_admn/ind_admn_ad3_py_s1_GADM_pp_Himachal_Pradesh.shp')
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
loc_points_sf <- st_as_sf(loc_points, coords = c("lon", "lat"), crs = 4326)


# 1. Air pollution (LongPMInd via zendo.org) ####
# Source: https://essd.copernicus.org/articles/16/3565/2024/#section3
# https://zenodo.org/records/10073944

# Load nc files
nc_files <- list.files(path='1_Original_Data/115_LongPMInd/', pattern = '\\.nc$', full.names = TRUE)

# # Check time
# (r <- rast(nc_files[1]))
# time(r)
# names(r)
# (r_var <- r[[3]])
# 
# # Check CRS (same)
# st_crs(hp)
# crs(r)
# 
# # Plot
# r_crop <- crop(r_var, vect(hp))
# r_mask <- mask(r_crop, vect(hp))
# plot(r_mask)
# plot(hp$geometry, add = TRUE, border = 'blue')
# plot(loc_points_sf$geometry, add = TRUE, col = 'red', pch = 20)

# Extract data

out <- data.frame()
for (f in nc_files){
  
  # Load data as raster
  r <- rast(f)
  
  # Extract values at location points, across time
  vals <- terra::extract(r, vect(loc_points_sf)) %>%
    mutate(location = loc_points_sf$location) %>%
    tidyr::pivot_longer(
      cols = starts_with('pm25'),  # adjust if bands named differently
      names_to = 'band',
      values_to = 'pm25'
    ) %>%
    mutate(date = rep(time(r), times = nrow(loc_points_sf))) %>%
    dplyr::select(location, date, pm25)
  
  out <- rbind(out, vals)
  
}

# Convert location column to a factor with this order
pm25_res <- out %>%
  mutate(location = factor(location, levels=location_levels)) %>%
  arrange(location, date)

# 2. Elevation data (AWS Data via elevatr package) ####

# Source: https://github.com/USEPA/elevatr

# Mapzen terrain tiles (which supply the AWS source) contain 3DEP, SRTM, and 
# GMTED2010 content courtesy of the U.S. Geological Survey and ETOPO1 content 
# courtesy of U.S. National Oceanic and Atmospheric Administration. The Open 
# Topography API provide access to the SRTM and the ALOS World 3D datasets. 
# See https://opentopography.org/ for details.


# Load rasters
# save_elev(hp, loc_stats, ind)
elev_raster <- raster("2_Active_Data/204_misc/ind_envi_s3_AWSOpenData_elevation_raster.tif")

# Extract
loc_points$elev <- raster::extract(elev_raster, loc_points_sf)
elev_res <- dplyr::select(loc_points, location, elev)

# 3. Population data (WorldPop Counts via worldpop.org) ####

# Source: https://hub.worldpop.org/geodata/listing?id=136

# Load rasters (counts cover time period better than density, though they seem comparable enough to combine maybe?)
# raster_files <- list.files("1_Original_Data/117_WorldPop_Density", pattern = "\\.tif$", full.names = TRUE)
raster_files <- list.files("1_Original_Data/116_WorldPop_Counts.", pattern = "\\.tif$", full.names = TRUE)

# Loop over rasters: crop, mask, extract
out <- list()
for (f in raster_files) {
  # Load raster
  r <- raster(f)
  
  # Make sure shapefile CRS matches raster
  hp_r <- st_transform(hp, crs(r))
  
  # Crop + mask raster to HP
  r_hp <- crop(r, extent(hp_r))
  r_hp <- mask(r_hp, hp_r)
  
  # Transform points to raster CRS
  points_r <- st_transform(loc_points_sf, crs(r))
  
  # Extract raster values at points
  vals <- raster::extract(r_hp, points_r)
  
  # Extract year from filename
  year_val <- str_split(basename(f), "_")[[1]][3]  # 3rd element after splitting by "_"
  year_val <- str_remove(year_val, "\\.tif$")      # remove ".tif"
  
  # Store results in a data.frame
  df <- loc_points %>%
    mutate(year = as.numeric(year_val),
           population_density = vals)
  
  out[[f]] <- df
}

# Extract
(pop_res <- bind_rows(out))
pop_res <- dplyr::select(pop_res, location, year, population_density)

#### EXPORT ####

save(pm25_res, file = '2_Active_Data/204_misc/ind_adm4_LongPMInd_PM25_2000-2022.RData')
save(elev_res, file = '2_Active_Data/204_misc/ind_misc_cont_adm4_elev.RData')
save(pop_res, file = '2_Active_Data/204_misc/ind_misc_cont_adm4_pop_2000-2023.RData')

