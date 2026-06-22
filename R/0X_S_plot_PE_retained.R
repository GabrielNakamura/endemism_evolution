
library(dplyr)
library(sf)
library(ggplot2)
library(rcartocolor)

# PE biogeographic reconstruction 

# coordinates
coords <- read.table(here::here("data", "coords.txt"), header = TRUE, sep = ";")
coords_2 <- data.frame(coords, grids = rownames(coords))

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
res_PE_one1 <- data.frame(res_PE_one, grids = as.character(bioregions_birds$ID))
res_PE_one_2 <- 
  res_PE_one1 |> 
  left_join(coords_2, by = c(grids = "grids"))

# removing grids with lat and long NA
res_PE_one_3 <- 
  res_PE_one_2 |> 
  filter(!is.na(LAT))

AS_sf <- read_sf(
  here::here("data",  
             "shape_america2.shp"))

# transforming data frame into a spatial data frame
evo_metrics_sf <- st_as_sf(
  res_PE_one_3,
  coords = c("LON","LAT"), 
  crs = st_crs(AS_sf)
)

evo_metrics_AS <- st_filter(evo_metrics_sf, AS_sf)
site_xy_AS <- st_coordinates(evo_metrics_AS$geometry)

# just changing LON and LAT to x and y for simplicity and placing it before PE column
evo_metrics_AS_2 <- 
  evo_metrics_AS |> 
  mutate(
    x = site_xy_AS[,1],
    y = site_xy_AS[,2],
    .before = PE
  )

cont_cols <- rcartocolor::carto_pal(8, "Sunset")

map_PE <- 
  evo_metrics_AS_2 |> 
  ggplot() + 
  geom_tile(aes(x = x, y = y, fill = PEinsitu)) +
 # scale_fill_stepsn(
  #  breaks = seq(0, 17, 3),
  # limits = c(0, 17),
  #  name = "PE retained" ,
  # colours = cont_cols
  #) +
  labs(title = "Phylogenetic endemism") 

coastline <- rnaturalearth::ne_coastline(returnclass = "sf")
map_limits <- list(
  x = c(-95, -30),
  y = c(-55, 12)
)

map_PE <- 
  res_PE_one_3 |> 
  ggplot2::ggplot() + 
  ggplot2::geom_raster(ggplot2::aes(x = LON, y = LAT, fill = PEinsitu)) + 
  rcartocolor::scale_fill_carto_c(name = "PE", 
                                  type = "quantitative", 
                                  palette = "SunsetDark") +
  ggplot2::geom_sf(data = coastline) +
  ggplot2::coord_sf(xlim = map_limits$x, ylim = map_limits$y) +
  ggplot2::ggtitle("") +
  ggplot2::xlab("Longitude") +
  ggplot2::ylab("Latitude") +
  ggplot2::labs(fill = "PE") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "mm"),
    legend.text = element_text(size = 12), 
    axis.text = element_text(size = 7),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11)
  )


library(sf)
library(terra)
library(tidyterra)
library(ggplot2)

# 1. Convert sf object to a spatVector
sf_vector <- vect(evo_metrics_AS_2)

# 2. Extract your original data's spatial resolution
# (Finds the true spacing between coordinates along the X axis)
coords <- st_coordinates(evo_metrics_AS_2)
detected_res <- min(diff(sort(unique(coords[, "X"]))))

# 3. Create a clean grid template matching that exact resolution
template <- rast(sf_vector, res = detected_res)

# 4. Burn your data into the true raster grid structure
map_raster <- rasterize(sf_vector, template, field = "PEinsitu")

# 5. Plot instantly using tidyterra's specialized raster layer
ggplot() +
  geom_spatraster(data = map_raster, aes(fill = your_variable_column)) +
  scale_fill_viridis_c(na.value = "transparent") +
  theme_minimal()