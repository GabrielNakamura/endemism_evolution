
library(ggplot2)
library(rcartocolor)
library(Herodotools)
library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth) # Modern replacement for rnaturalearthdata
library(rcartocolor)
# PE biogeographic reconstruction 

# coordinates
coords <- read.table(here::here("data", "coords.txt"), header = TRUE, sep = ";")

# bsm results
res_PE_one <- 
  readRDS(here::here("output", "PE_PD_insitu", "PD_PE_insitu_one_bsm.rds"))

# occurrence data and bioregion information
occ_birds <- read.table(here::here("data", "processed", "W_harvey.txt"), 
                        header = T)
bioregions_birds <- read.table(here::here("data", "matrixEco.txt"), 
                               header = T)
bioregions_birds2 <- data.frame(area = bioregions_birds[, 2])
rownames(bioregions_birds2) <- bioregions_birds$ID

# adding grid ID
res_PE_one1 <- data.frame(res_PE_one, grids = as.character(rownames(occ_birds)))


# plotting maps -----------------------------------------------------------

# creating limits to the maps

map_limits <- list(
  x = c(-170.2166, -13.21288),  # Longitude limits (xmin, xmax)
  y = c(-55.37714, 83.6236)     # Latitude limits (ymin, ymax)
)

# Spatializing the results

# Load data PEinsitu


# 2. Spatial Objects --------------------------------------------------

# sf coastline layer
coastline <- rnaturalearth::ne_coastline(scale = 50, returnclass = "sf")


# test 2 ------------------------------------------------------------------

# 3. Grid Matching & Data Preparation ----------------------------------------

# Initialize a blank raster for Americas
r_pe <- terra::rast(
  xmin = -170.2166, xmax = -13.21288, 
  ymin = -55.37714, ymax = 83.6236, 
  res = 1.0, 
  crs = "EPSG:4326",
  nlyr = 2                # Define 2 layers up front
)
names(r_pe) <- c("PEinsitu", "PE")

# Find the exact cell numbers for your coordinates
cells <- terra::cellFromXY(r_pe, coords[as.numeric(res_PE_one1$grids), ])

# Assign the PEinsitu values directly to those specific cells
values_to_assign <- res_PE_one1[, c("PEinsitu", "PE")]
r_pe[cells] <- values_to_assign


# Convert the raster grids into vector polygons so ggplot treats them natively
map_data_sf <- terra::as.polygons(r_pe, aggregate = FALSE) %>% 
  sf::st_as_sf()



# Map ---------------------------------------------------------------------

# 1. Calculate the spatial bounding box of ONLY your grid data
# This extracts the exact bounding polygon from the data you actually have
data_bbox <- sf::st_as_sfc(sf::st_bbox(map_data_sf))

# 2. Crop the coastline down to this data boundary
# This structurally deletes Africa and all distant small islands 
coastline_cropped <- sf::st_intersection(coastline, data_bbox)

# 3. Create dynamic map limits based tightly on your actual data
# This eliminates unnecessary white space
tight_bbox <- sf::st_bbox(map_data_sf)
map_limits_tight <- list(
  x = c(tight_bbox["xmin"], tight_bbox["xmax"]),
  y = c(tight_bbox["ymin"], tight_bbox["ymax"])
)

# Define the projection
equal_area_crs <- "+proj=cea +lon_0=-90 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# 4. Generate the clean figure
p_PE <- ggplot() +  
  # Plot the grid cells first
  geom_sf(data = map_data_sf, aes(fill = PE), color = NA) +  
  
  # Plot ONLY the cropped coastline on top
  geom_sf(data = coastline_cropped, fill = NA, color = "gray30", size = 0.3, inherit.aes = FALSE) +
  
  rcartocolor::scale_fill_carto_c(
    name = "Endemismo atual", 
    type = "quantitative", 
    palette = "SunsetDark",
    na.value = "transparent",
    guide = guide_colorbar(
      title.position = "top",   
      title.hjust = 0.5,        
      barwidth = unit(8, "cm"), 
      barheight = unit(0.3, "cm") 
    )
  ) +
  
  # Crop using the tight coordinates to drop blank spaces completely
  coord_sf(
    xlim = map_limits_tight$x, 
    ylim = map_limits_tight$y, 
    default_crs = sf::st_crs(4326), 
    crs = equal_area_crs,           
    expand = FALSE
  ) +
  
  labs(x = NULL, y = NULL) + 
  
  theme_minimal() + 
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "mm"),
    
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank(),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold")
  )

# Save the plot
ggsave(
  filename = here::here("Figures", "PE_atual.png"),
  plot = p_PE,
  width = 14, 
  height = 14, 
  units = "cm", 
  dpi = 500
)

p <- ggplot() +  
  # Plot the grid cells first
  geom_sf(data = map_data_sf, aes(fill = PEinsitu), color = NA) +  
  
  # Plot ONLY the cropped coastline on top
  geom_sf(data = coastline_cropped, fill = NA, color = "gray30", size = 0.3, inherit.aes = FALSE) +
  
  rcartocolor::scale_fill_carto_c(
    name = "Endemismo - manutenção de área", 
    type = "quantitative", 
    palette = "SunsetDark",
    na.value = "transparent",
    guide = guide_colorbar(
      title.position = "top",   
      title.hjust = 0.5,        
      barwidth = unit(8, "cm"), 
      barheight = unit(0.3, "cm") 
    )
  ) +
  
  # Crop using the tight coordinates to drop blank spaces completely
  coord_sf(
    xlim = map_limits_tight$x, 
    ylim = map_limits_tight$y, 
    default_crs = sf::st_crs(4326), 
    crs = equal_area_crs,           
    expand = FALSE
  ) +
  
  labs(x = NULL, y = NULL) + 
  
  theme_minimal() + 
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "mm"),
    
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank(),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold")
  )

ggsave(
  filename = here::here("Figures", "PEinsitu_bsm1.png"),
  plot = p,
  width = 14, 
  height = 14, 
  units = "cm", 
  dpi = 500
)
