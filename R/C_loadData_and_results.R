# read libraries and functions --------------------------------------------
library(ape)
library(phytools)
library(SYNCSA)
library(picante)
library(geiger)
library(ade4)
library(phylobase)
library(rgdal)
library(raster)
library(rgeos)
library(dismo)
library(letsR)
source(here::here("R", "functions", "assembly_test_07-8-2019.R"))
source(here::here("R", "functions", "anova.1way.R"))
source(here::here("R", "functions", "GridFilter.R"))


# read data ---------------------------------------------------------------

tree_harvey <- ape::read.tree(here::here("data", "processed", "Tree_TF400Howard_Pruned.tre"))
W_harvey <- read.table(here::here("data", "processed", "W_harvey.txt"))
ancestral.area <- read.table(here::here("data", "processed", "Econodes_harvey.txt"), header= TRUE) #ancestral area data - from 
ancestral.area <- ancestral.area[, -1]
ancestral.area <- data.frame(ancestral.area)
biogeo<- read.table(here::here("data", "matrixEco.txt"), header= TRUE) #ecoregions of each point in the map
biogeo <- biogeo[, -1]
biogeo <- data.frame(biogeo)
Eco<- read.table(here::here("data", "matrixEco.txt"),header=TRUE) 
tree <- ape::read.tree(here::here("data", "processed", "Tree_TF400Howard_Pruned.tre")) #read tree
coords <- read.table("data/coords.txt", sep = ";") #coordinates for all points
temp_trop<- c(rep("temperate", length(1:2248)), rep("tropical", length(2249:nrow(W_harvey)))) #categorizing the coordinates

# read results NRI and Ages -----------------------------------------------

ages_harvey <- readRDS(here::here("output", "agesResult_Harvey.rds"))
nri<- readRDS(here::here("output", "nriRes_harvey.rds"))
