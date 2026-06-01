library(rgdal)
library(raster)
library(rgeos)
library(dismo)
library(letsR)
source(here::here("R", "functions", "GridFilter.R"))

#######################################
#### SHAPE America
# Load shapefile
shape.america <- readOGR(here::here("data", "shape_america2.shp"))
shape.america
plot(shape.america)

#######################################
#### GRID
gridded <- GridFilter(shape.america, resol=1.0, prop=0.3)
saveRDS(object = gridded, file = here::here("data", "processed", "grid_tyranidae.rds"))
plot(gridded)
gridded
#str(gridded)
#writeOGR(gridded, dsn=getwd(), layer="final_shape", driver="ESRI Shapefile", overwrite_layer=T)

# coordenadas de cada celula
coords<-coordinates(gridded)
str(coords)
plot(coords)

# atribuir ID para cada celula
slot(gridded, "data") <- cbind("ID" = 1:length(gridded),
                            slot(gridded, "data"))
gridded$ID

#######################################
#### SHAPE Species
# load shapefiles
birds<-readOGR("shapes_tiranideos.shp")
birds
summary(birds)

#######################################
#### MATRIZ PA
# create presab from specified grid
pa.birds<-lets.presab.grid(shapes=birds,grid=gridded,sample.unit="ID")
names(pa.birds)
pa.birds$PAM # Presence/Absence Matrix
str(pa.birds$PAM) # sp

# Matriz W - Presence/Absence
write.table(pa.birds$PAM,"pa.birds.txt")

presab.birds<-pa.birds$PAM
fix(presab.birds)

#######################################
#### BIOMAS
# Load raster of Biomes
raster.biom.america<-raster("bioma_raster.asc") 
plot(raster.biom)
# crop raster to america
#raster.biom.america<-crop(raster.biom,shape.america)
#plot(raster.biom.america)

# extract raster values for each cell in the grid
biom_values<-extract(raster.biom.america,gridded,df=TRUE,fun=max)
str(biom_values) 
fix(biom_values)
levels(as.factor(biom_values$wwf_biomes_asc)) # 16 biomas
write.table(biom_values,file="a_biomes_values.txt")

