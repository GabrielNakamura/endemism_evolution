#Spatializing the results

library(raster)
library(sp)
library(rnaturalearthdata)
costline <- coastline50
ages <- mean_age
nri<- ses.mpd_res_noNA
names(nri)<- names(mean_age)
rownames(coords)<- 1:5567
r <- raster(vals=NA, xmn = -170.2166 , xmx = -13.21288, ymn = -55.37714, ymx = 83.6236, resolution=1.0)

### NRI map ####
cell.r <- cellFromXY(r, coords[names(nri),])
values_cell <- rep(NA, ncell(r))
names(values_cell) <- 1:ncell(r)
nri.cells <- 1:ncell(r) %in% cell.r
values_cell[nri.cells] <- nri
r.nri <- raster::setValues(r, values = values_cell)

### AGE map ####
cell.r <- cellFromXY(r, coords[names(mean_age),])
values_cell <- rep(NA, ncell(r))
names(values_cell) <- 1:ncell(r)
age.cells <- 1:ncell(r) %in% cell.r
values_cell[age.cells] <- mean_age
r.age <- raster::setValues(r, values = values_cell)

### Figure ####
#exporting as tif image
tiff("Fig3_Ses_Age_06-07.tif", units = 'cm', width = 17.4, height = 12, res = 300)
par(mfrow = c(1,2), mar=c(6,4,4,5))
plot(r.nri, xlab = "Longitude", ylab = "Latitude")
plot(costline, add=T)
mtext("a",side = 3, line = 0.5, font = 2, adj = 0, cex = 1.5)
plot(r.age, xlab = "Longitude", ylab = "Latitude")
plot(costline, add=T)
mtext("b",side = 3, line = 0.5, font = 2, adj = 0, cex = 1.5)
dev.off()
