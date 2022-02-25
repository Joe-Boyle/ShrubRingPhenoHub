# Quick prep script to generate the QHI salix pheno transect coordiantes
# Jakob J Assmann 11 Feb 2022

library(tidyverse)
library(sf)

phenocoords <- data.frame(matrix(c(
"Phenocam1",	69.57560,	-138.90460,	110,
"Phenocam2",	69.57570,	-138.90525,	87,
"Phenocam3",	69.57781,	-138.91277,	83,
"Phenocam4",	69.57515,	-138.89453,	17,
"Phenocam5",	69.57651,	-138.86711,	67,
'Phenocam6',	69.57483,	-138.86305,	70,
"Joe Coords1", 69.5756,	-138.9046, NA,
"Herschel Manual Coords", 69.57571666666667, -138.90456666666667, NA), ncol = 4, byrow = T)) %>%
  setNames(c("Name", "Lat", "Long", "GPS height")) %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
  st_transform(crs = 32607)
write_sf(phenocoords, "data/map_figure_data/phenocoords/phenocoords.shp")
plot(phenocoords)
