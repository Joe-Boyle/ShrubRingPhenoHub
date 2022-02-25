####################################
### Tidy up GPS points for map   ###
### Sandra Angers-Blondin        ###
### 30-10-2020                   ###
####################################


# Packages ----------------------------------------------------------------

library(sf)
library(dplyr)
library(ggplot2)


# Read data ---------------------------------------------------------------

gps_raw <- read.csv("data/map_figure_data/JB_GPS_salix.csv", sep = ";")


# Convert to simple features ----------------------------------------------

# A points version
gps_pt = st_as_sf(gps, coords = c("lon", "lat"),
                 crs = 4326)

## Also create a quadrat polygon from bounding box
quadrat <- st_as_sfc(st_bbox(gps_pt))

# Check that everything looks right

ggplot() +
   geom_sf(data = quadrat, fill = "lightgrey") +
   geom_sf(data = gps_pt, shape = 16, size = 2, colour = "navy") +
   theme_bw()


# Export spatial features -------------------------------------------------

st_write(gps_pt, "data/GPS_points", layer = "transect_pts", driver = "ESRI shapefile")
st_write(quadrat, "data/GPS_points",layer = "quadrat", driver = "ESRI shapefile")
