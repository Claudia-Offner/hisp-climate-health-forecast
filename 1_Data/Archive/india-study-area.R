### ------------------------------------------------------------------------ ### 
###  Study area map                                                          ###
### ------------------------------------------------------------------------ ### 

# Turn off scientific notation
options(scipen=999)

# Set data path
user <- 'claer14'
path1 <- paste0('C:/Users/', user, '/Documents/GitHub/HISP_TSGC/')
path2 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Data/2_Active_Data/')
path3 <- paste0('C:/Users/', user, '/OneDrive - University of Cambridge/III. Research/Chapter 0 - Research Planning/B. First Year Report')
setwd(path2)

####  PACKAGES & FUNCTIONS ####

library(ggplot2)
library(raster)
library(sf)
library(elevatr)
library(ggnewscale)
library(cowplot)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

####  MAIN CODE  #####


# 1. Load data ####

# Get admin data for india and filtered to Himachal Pradesh
shapefile_path <- "C:/Users/claer14/OneDrive - University of Cambridge/III. Research/Data/1_Original_Data/100_GADM/gadm41_IND_shp/"
ind <- st_read(paste0(shapefile_path,'gadm41_IND_1.shp'))
shp <- st_read(paste0(shapefile_path,'gadm41_IND_3.shp'))
hp <- shp[shp$NAME_1 == "Himachal Pradesh", ]



# Create a data frame with locations and coordinates
# locations <- data.frame(
#   location = c("Chamba", "Dharamshala", "Solan", "Hamirpur", "Kullu"),
#   lon = c(76.127017, 76.321482, 77.099506, 76.522135, 77.107945),
#   lat = c(32.554798, 32.216995, 30.906408, 31.684233, 31.956550)
# )
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
cities <- st_as_sf(locations, coords = c("lon", "lat"), crs = 4326)


# 2. Get elevation raster data ####

# Get bbox of state + 1 km buffer
bbox_geom <- st_as_sfc(st_bbox(hp))
bbox <- st_sf(st_buffer(bbox_geom, dist = 20000))

# Get elevation data (zoom level 9–11 is usually detailed enough)
elev_raster <- get_elev_raster(locations = bbox, z=7, clip = "bbox") #7=~1km res

# Inverse mask out Admin1 outside HP
# bbox_sf <- st_sf(geometry = bbox, crs = st_crs(hp))
bbox2 <- st_sf(st_buffer(bbox, dist = 10000))
outside_hp <- st_difference(bbox2, st_union(hp))

# Compute terrain derivatives
slope <- terrain(elev_raster, opt = 'slope')
aspect <- terrain(elev_raster, opt = 'aspect')
hillshade <- hillShade(slope, aspect, angle = 40, direction = 315)

# Convert for plotting
elev_df <- as.data.frame(elev_raster, xy = TRUE)
colnames(elev_df) <- c("x", "y", "elevation")
hill_df <- as.data.frame(hillshade, xy = TRUE)
colnames(hill_df) <- c("x", "y", "hill")

# 3. Plot ####

# Set colour scheme
arcgis_topo_colors <- c(
  "#5A7D24",  # dark green (low elevation)
  "#9DBF68",  # light green
  "#F3E99F",  # pale yellow
  "#D2B48C",  # tan/light brown
  "#8B7A55",  # brown (higher elevation)
  "#FFFFFF"   # white (highest peaks/snow)
)

main_map <- ggplot() +
  # Hillshade raster (full extent)
  geom_raster(data = hill_df, aes(x = x, y = y, fill = hill), show.legend = FALSE) +
  scale_fill_gradient(low = "grey20", high = "grey90", guide = "none") +
  
  new_scale_fill() +
  
  # Elevation raster (full extent)
  geom_raster(data = elev_df, aes(x = x, y = y, fill = elevation), alpha = 0.6) +
  scale_fill_gradientn(colours = arcgis_topo_colors, name = "Elevation (m)") +
  
  # Overlay semi-transparent grey mask outside HP
  geom_sf(data = outside_hp, fill = "grey50", alpha = 0.6) +  # white fill + black border

  # State, HP boundaries and cities
  geom_sf(data = ind, fill = NA, color = "grey50", size = 0.6) +
  geom_sf(data = hp, fill = NA, color = "black", size = 0.6) +
  geom_sf(data = cities, color = "darkred", size = 2) +
  geom_label(data = locations, aes(x = lon, y = lat, label = location),
    size = 3, color = "darkred", fill = "white", fontface = "bold",
    label.padding = unit(0.15, "lines"), label.r = unit(0.15, "lines"), 
    label.size = 0,   nudge_y = 0.1) +
  coord_sf(xlim = range(elev_df$x), ylim = range(elev_df$y), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal()+
  theme(
    plot.background = element_rect(fill = 'white', color = NA),
    panel.background = element_rect(fill = 'white', color = NA)
  )



# Create inset map highlighting HP
world <- ne_countries(scale = "medium", returnclass = "sf")
oceans <- ne_download(scale = "medium", type = "ocean", category = "physical", returnclass = "sf")
buffer <- 2
india_bbox <- st_bbox(ind)
xlim <- c(india_bbox$xmin - buffer, india_bbox$xmax + buffer)
ylim <- c(india_bbox$ymin - buffer, india_bbox$ymax + buffer)

inset_map <- ggplot() +
  geom_sf(data = oceans, fill = "lightblue") +
  geom_sf(data = world, fill = "white", color = "grey70") +
  geom_sf(data = ind, fill = "grey60", color = "grey60", size = 0.3) +
  geom_sf(data = hp, fill = "darkred", color = "darkred", size = 0.4) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +   # zoom in on India
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.6),
    plot.background = element_rect(fill = 'white', color = NA),
    panel.background = element_rect(fill = 'white', color = NA)
  )


final_plot <- main_map + 
  inset_element(
    inset_map,
    left = 0.75,  # fraction from left (0 = left, 1 = right)
    bottom = 0.75, # fraction from bottom (0 = bottom, 1 = top)
    right = 1.012,
    top = 1,
    align_to = 'panel'  # aligns to plot panel, not entire page
  )

final_plot

#### EXPORT RESULTS ####

# ggsave(paste0(path3, '/II. Figures/F2_studyarea.png'), final_plot, width = 10, height = 6, dpi = 300)
