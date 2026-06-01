# Spatializing the results


# reading libraries and data ----------------------------------------------
ages_harvey <- readRDS(here::here("output", "agesResult_Harvey.rds"))
nri_harvey <- readRDS(here::here("output", "nriRes_harvey.rds"))

mean_age <- apply(ages_harvey, 1, function(x) mean(x[which(x != 0)])) #mean arrival age for each assemblage
ses_mpd_res <- nri_harvey[, 6]
names(ses_mpd_res) <- rownames(nri_harvey)
ses.mpd_res_noNA <- ifelse(!is.na(ses_mpd_res), ses_mpd_res, 0)
library(raster)
library(sp)
library(rnaturalearthdata)


# spatial objects ---------------------------------------------------------
costline <- coastline50
rownames(coords)<- 1:5567
r <- raster(vals=NA, xmn = -170.2166 , xmx = -13.21288, ymn = -55.37714, ymx = 83.6236, resolution=1.0)


# mapping NRI -------------------------------------------------------------

cell.r <- cellFromXY(r, coords[names(ses.mpd_res_noNA),])
values_cell <- rep(NA, ncell(r))
names(values_cell) <- 1:ncell(r)
nri.cells <- 1:ncell(r) %in% cell.r
values_cell[nri.cells] <- ses.mpd_res_noNA
r.nri <- raster::setValues(r, values = values_cell)


# mapping age -------------------------------------------------------------

cell.r <- cellFromXY(r, coords[names(mean_age),])
values_cell <- rep(NA, ncell(r))
names(values_cell) <- 1:ncell(r)
age.cells <- 1:ncell(r) %in% cell.r
values_cell[age.cells] <- mean_age
r.age <- raster::setValues(r, values = values_cell)


# generating and exporting figures ----------------------------------------
data_saving <- date()
data_edit <- unlist(strsplit(data_saving, split = " "))[c(1, 2, 3, 5)]
name_plot <- paste("Fig3_Ses_Age", 
                   data_edit[1],
                   data_edit[2], 
                   data_edit[3],
                   data_edit[4],
                   sep = "_")
tiff(filename = here::here("Figures", paste(name_plot, "tif", sep = ".")), units = 'cm', width = 25.4, height = 12, res = 500)
quartz()
par(mfrow = c(1,2), mar=c(6,4,4,5))
plot(r.nri, xlab = "Longitude", ylab = "Latitude")
plot(costline, add=T)
mtext("a",side = 3, line = 0.5, font = 2, adj = 0, cex = 1.5)
plot(r.age, xlab = "Longitude", ylab = "Latitude")
plot(costline, add=T)
mtext("b",side = 3, line = 0.5, font = 2, adj = 0, cex = 1.5)
dev.off()
