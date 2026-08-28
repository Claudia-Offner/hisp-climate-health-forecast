### ------------------------------------------------------------------------ ### 
### HP Hospital cleaning - Controls                                          ###
### ------------------------------------------------------------------------ ### 

# Prep environment
rm(list = ls()); graphics.off(); gc()
setwd('C:/Users/claer14/Documents/GitHub/hisp-climate-health-forecast/')


#### DEPENDENCIES ####
source('1_Data/Data - Dependencies.R')
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
loc_points <- st_as_sf(
  data.frame(
    location = c('Chamba', 'Dharamshala', 'Solan',
                 'Hamirpur', 'Kullu', 'Kangra',
                 'Reckong Peo', 'Kyelang', 'Nahan',
                 'Bilaspur', 'Mandi', 'Shimla', 'Una'),
    fclass = c('District head', 'District head', 'District head', 
               'District head', 'District head', 'City',
               'District head', 'District head', 'District head',
               'District head', 'District head', 'District head', 'District head'),
    `Study Area` = c(1, 1, 1,
                     1, 1, 1, 
                     1, 0, 1, 
                     0, 1, 0, 1),
    lon = c(76.127017, 76.321482, 77.099506,
            76.522135, 77.107945, 76.274797,
            78.267181, 77.023559, 77.294546,
            76.761165, 76.928977, 77.170999, 76.276997),
    lat = c(32.554798, 32.216995, 30.906408,
            31.684233, 31.956550, 32.103220,
            31.539809, 32.571883, 30.564781,
            31.341244, 31.712559, 31.104221, 31.479129)
  ),
  coords = c("lon", "lat"),
  crs = 4326
)
loc_points_sf <- st_as_sf(loc_points, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(crs = 3857) |>
  st_buffer(dist = 5000) |>
  st_transform(crs = 4326)


# plot(loc_points)
# path <- paste0(path_data, '2_Active_Data/200_admn/ind_admn_gps_pt_s1_GADM_pp_study_sites.shp')
# st_write(loc_points, path, append=FALSE)


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
extract_pm25 <- function(f) {
  r <- rast(f)
  terra::extract(r, vect(loc_points_sf), fun = mean, na.rm = TRUE) %>%
    mutate(location = loc_points_sf$location) %>%
    pivot_longer(starts_with('pm25'), names_to = 'band', values_to = 'pm25') %>%
    mutate(date = rep(time(r), times = nrow(loc_points_sf))) %>%
    dplyr::select(location, date, pm25)
}

pm25_res <- map(nc_files, extract_pm25) %>%
  bind_rows() %>%
  mutate(location = factor(location, levels = location_levels)) %>%
  dplyr::arrange(location, date)

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
res(elev_raster) * 111000  # in metres

# Extract
loc_points$elev <- raster::extract(elev_raster, loc_points_sf, fun = mean, na.rm = TRUE)
elev_res <- dplyr::select(loc_points, location, elev) %>% st_drop_geometry()

# 3. Population data (WorldPop Counts via worldpop.org) ####

# Source: https://hub.worldpop.org/geodata/listing?id=136

# Load rasters (counts cover time period better than density, though they seem comparable enough to combine maybe?)
raster_files <- list.files("1_Original_Data/116_WorldPop_Counts.", pattern = "\\.tif$", full.names = TRUE)

extract_pop <- function(f) {
  r       <- raster(f)
  r_hp    <- mask(crop(r, extent(st_transform(hp, crs(r)))), st_transform(hp, crs(r)))
  vals    <- raster::extract(r_hp, st_transform(loc_points_sf, crs(r)), fun = sum, na.rm = TRUE)
  year_val <- str_remove(str_split(basename(f), "_")[[1]][3], "\\.tif$")
  loc_points %>% mutate(year = as.numeric(year_val), popu = vals)
}

pop_res <- map(raster_files, extract_pop) %>%
  bind_rows() %>%
  dplyr::select(location, year, popu) %>%
  st_drop_geometry()

pop_res$popu <- as.numeric(pop_res$popu)

# 4. Demographic & health data (WorldPop Counts via worldpop.org) ####

# Get all MEAN tif files across all folders
raster_folders <- list.files("1_Original_Data/121_WorldPop_Demo_Heal", full.names = TRUE) 
raster_files <- list.files(raster_folders, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
raster_files <- raster_files[grepl('MEAN', raster_files)]

# Extract values for each file, return named vector
extract_vals <- function(f) {
  print(f)
  r     <- raster(f)
  r_hp  <- crop(mask(r, st_transform(hp, crs(r))), extent(st_transform(hp, crs(r))))
  vals  <- raster::extract(r_hp, st_transform(loc_points_sf, crs(r)), fun = mean, na.rm = TRUE)
  measure <- str_remove(str_split(basename(f), "_")[[1]][2], "\\.tif$")
  setNames(as.numeric(vals), NULL) |> (\(v) tibble(!!measure := v))()
}

out <- map(raster_files, extract_vals)
demo_res <- bind_cols(loc_points, out) %>% st_drop_geometry()

# 5. Hospital information ####

districts <- st_read(paste0(path_data, "2_Active_Data/200_admn/ind_admn_ad2_py_s1_GADM_pp_india.shp")) %>%
  filter(NAME_1 == 'Himachal Pradesh')

pois <- st_read('1_Original_Data/122_OSM_geofabrik/northern-zone-260607-free.shp/gis_osm_pois_a_free_1.shp') %>%
  filter(fclass == 'hospital') %>%
  st_transform(st_crs(districts))

hosp <- districts %>%
  mutate(n_hospitals = lengths(st_intersects(., pois))) %>%
  group_by(NAME_2) %>%
  summarise(n_hospitals = sum(n_hospitals)) %>%
  st_drop_geometry() %>%
  add_row(NAME_2 = 'Dharamshala', n_hospitals = 5, .before = 3) %>%
  rename(location = NAME_2, hosp_n = n_hospitals) %>%
  mutate(
    location  = c("Bilaspur","Chamba","Dharamshala","Hamirpur","Kangra",
                  "Reckong Peo","Kullu","Kyelang","Mandi","Shimla","Nahan","Solan","Una"),
    hosp_teach = as.integer(location %in% c("Dharamshala", "Hamirpur"))
  )

demo_res <- demo_res %>%
  dplyr::select(-elev) %>%
  left_join(hosp, by = 'location') %>%
  left_join(elev_res, by = 'location') %>%
  rename_with(tolower)
demo_res$elev <- as.numeric(demo_res$elev)

#### EXPORT ####

save(pm25_res, file = '2_Active_Data/204_misc/ind_adm4_LongPMInd_PM25_2000-2022.RData')
save(pop_res, file = '2_Active_Data/204_misc/ind_misc_cont_adm4_pop_2000-2023.RData')
save(demo_res, file = '2_Active_Data/204_misc/ind_misc_cont_adm4_demo_elev.RData')

