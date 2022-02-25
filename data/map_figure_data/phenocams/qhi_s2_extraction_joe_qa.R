# Quick Quality Control script for Joe's QHI Phenocam extractions
# Jakob J. Assmann 8 March 2021 j.assmann@bio.au.dk

# Dependencies
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggrepel)

# Load phenocam extactions (from GEE script /user/imyersmith/ShrubHub/Phenology/QHI_S2_phenocams_JA)
  # A direct snapshot to the script can be found here: https://code.earthengine.google.com/2f006e66b09b02a2374a133288bd3480
pheno_data <- read_csv("data/map_figure_data/phenocams/S2QHIphenocam.csv")

# Prepare and plot data
plot_list <- pheno_data %>% 
  select(name, year, doi, CLP, CLP_20m, NDVI, NDVI_20m, kNDVI, kNDVI_20m, NDSI_20m) %>%
  mutate(CLP = CLP / 100,
         CLP_20m = CLP_20m / 100) %>%
  filter(CLP_20m <= 0.25) %>%
  pivot_longer(cols = 4:10, names_to = "Index") %>%
  group_by(name) %>%
  group_map(function(...){
    group_sub <- data.frame(...) 
    group_sub %>% ggplot(aes(x = doi, y = value, colour = Index)) +
      geom_line() +
      ggtitle(unique(group_sub$name)) +
      facet_wrap(vars(year)) +
      scale_y_continuous(limits = c(-0.6,1)) +
      theme_cowplot()
    })

# Prepare coordinates
phenocams <- pheno_data %>% 
  distinct(name, long, lat) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(unique(pheno_data$EPSG))

# Set bb for plotting
bb_for_plot <- st_bbox(st_buffer(phenocams, 5000))

# # Prepare coastline
# world <- ne_coastline(scale = "large", returnclass = "sf")  %>%
#   st_crop(phenocams %>% st_convex_hull() %>% st_buffer(500000) %>% st_transform(4326)) %>%
#   st_transform(unique(pheno_data$EPSG))

# Load detailed Yukon Coastline
yukon_coast <- read_sf("data/map_figure_data/phenocams/geometries/yukon_bounds_utmz7.shp") %>%
  st_crop(phenocams %>% st_convex_hull() %>% st_buffer(500000))

# Make map
phenocam_map <- ggplot() + geom_sf(data = yukon_coast, col = "grey") + 
  geom_sf(data = phenocams, aes(colour = name), size = 2) +
  coord_sf(xlim = c(bb_for_plot$xmin, bb_for_plot$xmax), 
           ylim = c(bb_for_plot$ymin, bb_for_plot$ymax), expand = F) +
  geom_label_repel(data = phenocams, aes(geometry = geometry, 
                                         label = name),
                   stat = "sf_coordinates",
                   min.segment.length = 0,
                   box.padding = 0.1,
                   nudge_x = c(1000,-2500,-2500,500, 2000,2000),
                   nudge_y = c(-2000, -500, 1000, 1500, 500, -1000)
             ) +
  theme_nothing() +
  theme(legend.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white")) 

# combine and make a blot
plot_grid(plotlist = plot_list,
          ncol = 2, 
          labels = letters[2:7]) %>%
  plot_grid(phenocam_map, .,
            labels = c("a",""),
            nrow = 2,
            rel_heights = c(1,3)) %>%
  save_plot("data/map_figure_data/phenocams/S2QHI_QA_plot.png", ., ncol = 2, nrow = 4)
