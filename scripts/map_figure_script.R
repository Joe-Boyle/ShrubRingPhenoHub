# Quick script for Joe to generate map of North America with QHI location
# Jakob Assmann j.assmann@bio.au.dk 4 November 2020

# Dependencies
library(sf)
library(dplyr)
library(ggplot2)
library(cowplot)
library(rnaturalearth)
library(rnaturalearthdata)

# Define lambert conformic conic crs 
lambert_conformic <- 'PROJCS["North_America_Lambert_Conformal_Conic",
                            GEOGCS["GCS_North_American_1983",
                                   DATUM["North_American_Datum_1983",
                                         SPHEROID["GRS_1980",6378137,298.257222101]],
                                   PRIMEM["Greenwich",0],
                                   UNIT["Degree",0.017453292519943295]],
                            PROJECTION["Lambert_Conformal_Conic_2SP"],
                            PARAMETER["False_Easting",0],
                            PARAMETER["False_Northing",0],
                            PARAMETER["Central_Meridian",-140], 
                            PARAMETER["Standard_Parallel_1",20],
                            PARAMETER["Standard_Parallel_2",60],
                            PARAMETER["Latitude_Of_Origin",40],
                            UNIT["Meter",1],
                            AUTHORITY["EPSG","102009"]]'

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(lambert_conformic)
qhi_location <- data.frame(y = c(69.58), x = c(-139.05), label = "Qikiatruk") %>%
  st_as_sf(coords = c("x", "y"),crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  st_transform(lambert_conformic)
map_bounds <- data.frame(y = c(40, 70), x = c(-162.5, -97.5)) %>%
  st_as_sf(coords = c("x", "y"),crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  st_transform(lambert_conformic) %>%
  st_coordinates()
canada_map <- ggplot() + geom_sf(data = world, fill = "#ffffffFF", size = 0.5) + 
  geom_sf(data = qhi_location, colour = "white",
          fill = "#1e5c91FF", shape = 21, size = 8) +  
  coord_sf(xlim = c(min(map_bounds[,1]),
                    max(map_bounds[,1])),
           ylim = c(min(map_bounds[,2]),
                    max(map_bounds[,2])), expand = F) +
  theme(legend.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "mm"),
        axis.title = element_blank(),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white"), 
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(colour = "white", fill=NA, size=0.5)
  ) 
# theme(panel.background = element_rect(fill = "lightskyblue"))
canada_map
save_plot(canada_map, filename = "figures/map_figure/canada_map.png",
          base_aspect_ratio = 1, base_height = 5)  
