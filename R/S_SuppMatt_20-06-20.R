### Richness map for supplementary material ###

cell.r <- cellFromXY(r, coords[row.names(W),])


values_cell <- rep(NA, ncell(r))
names(values_cell) <- 1:ncell(r)
W.cells <- 1:ncell(r) %in% cell.r
values_cell[W.cells] <- rowSums(W)

r.W <- raster::setValues(r, values = values_cell)

windows(5,5)
plot(r.W, xlab = "Longitude", ylab = "Latitude")
plot(costline, add=T)

#### plot for supplementary material ####
nri_noNA<- nri[!is.na(nri[, 6]), "mpd.obs.z"]
coords_noNA<- coords[which(!is.na(nri[, 6]) == T),]
quartz()
plot(coords_noNA[, 2], nri_noNA, pch= 19, cex= 1, xlab= "Latitude", ylab= "NRI")
which(coords_noNA[,2] < 0)
