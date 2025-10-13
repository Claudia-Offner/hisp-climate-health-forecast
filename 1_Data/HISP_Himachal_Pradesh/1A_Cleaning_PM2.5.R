### ------------------------------------------------------------------------ ### 
### Air pollution cleaning - HISP Data                                       ###
### ------------------------------------------------------------------------ ### 

#### DEPENDENCIES ####

# Turn off scientific notation
options(scipen=999)

# Set data path
user <- 'claer14'
shp_path <- "C:/Users/claer14/OneDrive - University of Cambridge/III. Research/Data/1_Original_Data/100_GADM/gadm41_IND_shp/"
path <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/')
setwd(path)

# Packages
library(terra); library(sf); library(ggplot2); library(dplyr)

# Dependencies
location_order <- c("Chamba", "Dharamshala", "Solan", "Hamirpur", "Kullu",
                    "Reckong Peo", "Kyelang", "Nahan", "Bilaspur", "Mandi",
                    "Shimla", "Una")

#### Load data ####

# Get shape file
shp <- st_read(paste0(shp_path,'gadm41_IND_3.shp'))
hp <- shp[shp$NAME_1 == "Himachal Pradesh", ]

# Get locations
locations <- data.frame(
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
locations <- st_as_sf(locations, coords = c("lon", "lat"), crs = 4326)

# Load nc files
nc_files <- list.files(path='1_Original_Data/115_LongPMInd/', pattern = "\\.nc$", full.names = TRUE)

# Test file
r <- rast(nc_files[1])

#### Check temporal info ####
r
time(r)
names(r)

# Select time series like so...
(r_var <- r[[3]])

#### Check spatial info ####

# Check CRS (same)
st_crs(hp)
crs(r)

# Plot
r_crop <- crop(r_var, vect(hp))
r_mask <- mask(r_crop, vect(hp))
plot(r_mask)
plot(hp$geometry, add = TRUE, border = "blue")
plot(locations$geometry, add = TRUE, col = "red", pch = 20)

#### Extract data #####

out <- data.frame()
for (f in nc_files){
  
  # Load data as raster
  r <- rast(f)
  
  # Extract values at location points, across time
  vals <- terra::extract(r, vect(locations)) %>%
    mutate(location = locations$location) %>%
    tidyr::pivot_longer(
      cols = starts_with("pm25"),  # adjust if bands named differently
      names_to = "band",
      values_to = "pm25"
    ) %>%
    mutate(date = rep(time(r), times = nrow(locations))) %>%
    select(location, date, pm25)
  
  out <- rbind(out, vals)
  
}

results <- out %>% arrange(location, date)

# Convert location column to a factor with this order
results <- out %>%
  mutate(location = factor(location, levels = location_order)) %>%
  arrange(location, date)


#### EXPORT ####

# write.csv(results, file='2_Active_Data/202_weat/ind_adm4_LongPMInd_PM25_2000-2022.csv', row.names = FALSE)
save(results, file = '2_Active_Data/202_weat/ind_adm4_LongPMInd_PM25_2000-2022.RData')
