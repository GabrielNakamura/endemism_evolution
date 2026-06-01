
# reading data  -----------------------------------------------------------
tree_harvey <- ape::read.tree(here::here("data", "processed", "Tree_TF400Howard_Pruned.tre"))
W_harvey <- read.table(here::here("data", "processed", "W_harvey.txt"))
ancestral.area <- read.table(here::here("data", "processed", "Econodes_harvey.txt"), header= TRUE) #ancestral area data - from 
ancestral.area <- ancestral.area[, -1]
ancestral.area <- data.frame(ancestral.area)
biogeo<- read.table(here::here("data", "matrixEco.txt"), header= TRUE) #ecoregions of each point in the map
biogeo <- biogeo[, -1]
biogeo <- data.frame(biogeo)

# reading function --------------------------------------------------------
source(here::here("R", "functions", "DivB_metrics_18-08-20.R"))

# calculating age arrival -------------------------------------------------

res_age_arrival_Harvey <- diversification.assembly(W = W_harvey, tree = tree_harvey, ancestral.area = ancestral.area, biogeo = biogeo)

ages_harvey <- res_age_arrival_Harvey$age_arrival

saveRDS(ages_harvey, here::here("output", "Supp_agesResult_Harvey.rds"))

mean_age <- apply(ages_harvey, 1, function(x) mean(x[which(x != 0)])) #mean arrival age for each assemblage


# calculating NRI ---------------------------------------------------------

dis_harvey <- cophenetic(tree_harvey) #cophenetic distance matrix
org <- SYNCSA::organize.syncsa(comm = W_harvey, phylodist = dis_harvey) #organizing matrices
nri <- ses.mpd(org$community, org$phylodist, null.model = "taxa.labels", runs = 500) #nri calculation
ses_mpd_res <- nri[, 6]
ses.mpd_res_noNA <- ifelse(!is.na(ses_mpd_res), ses_mpd_res, 0)
saveRDS(nri, here::here("output", "nriRes_harvey.rds")) #saving nri results


# Anova age and NRI -------------------------------------------------------
coords_filter <- coords[match(names(mean_age), rownames(coords)), ]
anova_data_ses <- data.frame(mean_age= mean_age, NRI= -1*ses.mpd_res_noNA, local= temp_trop, coords = coords_filter)
res_anovaAge <- anova.1way(mean_age ~ temp_trop, data = anova_data_ses, nperm = 10000) #anova with age
res_anovaNRI <- anova.1way(NRI ~ temp_trop, data = anova_data_ses, nperm = 10000) #anova with NRI values

# saving processed data ---------------------------------------------------

write.csv(x = anova_data_ses, here::here("data", "processed", "anova_tableNRI_Age_Harvey.csv"))

