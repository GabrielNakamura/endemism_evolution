GridFilter<-function(shape, resol = 1, prop = 0){
grid <- raster(extent(shape))
res(grid) <- resol
proj4string(grid)<-proj4string(shape)
gridpolygon <- rasterToPolygons(grid)
drylandproj<-spTransform(shape, CRS("+proj=laea"))
gridpolproj<-spTransform(gridpolygon, CRS("+proj=laea"))
gridpolproj$layer <- c(1:length(gridpolproj$layer))
areagrid <- gArea(gridpolproj, byid=T)
dry.grid <- intersect(drylandproj, gridpolproj)
areadrygrid <- gArea(dry.grid, byid=T)
info <- cbind(dry.grid$layer, areagrid[dry.grid$layer], areadrygrid)
dry.grid$layer<-info[,3]/info[,2]
dry.grid <- spTransform(dry.grid, CRS(proj4string(shape)))
dry.grid.filtered <- dry.grid[dry.grid$layer >= prop,]}