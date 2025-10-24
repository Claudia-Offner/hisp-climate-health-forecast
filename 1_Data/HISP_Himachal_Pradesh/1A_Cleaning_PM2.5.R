### ------------------------------------------------------------------------ ### 
### HP Air pollution cleaning - LongPMInd                                    ###
### ------------------------------------------------------------------------ ### 

setwd(here())
source('1_Data/HISP_Himachal_Pradesh/HISP_Cleaning - Dependencies.R')

# Set data path
setwd(path_data)

####  FUNCTIONS  #####

#### MAIN CODE ####

# Load shape file
shp <- st_read('1_Original_Data/100_GADM/gadm41_IND_shp/gadm41_IND_3.shp')
hp <- shp[shp$NAME_1 == 'Himachal Pradesh', ]

# Get locations
locations <- data.frame(
  location = location_levels,
  lon = c(76.127017, 76.321482, 77.099506, 76.522135, 77.107945,
          78.267181, 77.023559, 77.294546, 76.761165, 76.928977,
          77.170999, 76.276997),
  lat = c(32.554798, 32.216995, 30.906408, 31.684233, 31.956550,
          31.539809, 32.571883, 30.564781, 31.341244, 31.712559,
          31.104221, 31.479129)
)
locations <- st_as_sf(locations, coords = c('lon', 'lat'), crs = 4326)

# Load nc files
nc_files <- list.files(path='1_Original_Data/115_LongPMInd/', pattern = '\\.nc$', full.names = TRUE)

# 1. Check spatial-temporal info ####
# Check time
(r <- rast(nc_files[1]))
time(r)
names(r)
(r_var <- r[[3]])

# Check CRS (same)
st_crs(hp)
crs(r)

# Plot
r_crop <- crop(r_var, vect(hp))
r_mask <- mask(r_crop, vect(hp))
plot(r_mask)
plot(hp$geometry, add = TRUE, border = 'blue')
plot(locations$geometry, add = TRUE, col = 'red', pch = 20)


# 2. Extract data #####

out <- data.frame()
for (f in nc_files){
  
  # Load data as raster
  r <- rast(f)
  
  # Extract values at location points, across time
  vals <- terra::extract(r, vect(locations)) %>%
    mutate(location = locations$location) %>%
    tidyr::pivot_longer(
      cols = starts_with('pm25'),  # adjust if bands named differently
      names_to = 'band',
      values_to = 'pm25'
    ) %>%
    mutate(date = rep(time(r), times = nrow(locations))) %>%
    select(location, date, pm25)
  
  out <- rbind(out, vals)
  
}

results <- out %>% arrange(location, date)

# Convert location column to a factor with this order
results <- out %>%
  mutate(location = factor(location, levels=location_levels)) %>%
  arrange(location, date)


#### EXPORT ####

save(results, file = '2_Active_Data/202_weat/ind_adm4_LongPMInd_PM25_2000-2022.RData')
