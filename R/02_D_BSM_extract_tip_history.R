

library(BioGeoBEARS)
library(Herodotools)
library(ape)
library(phytools) # helpful for extracting ancestral paths
library(progressr)
library(future)
library(furrr)

path_phyllip_file_birds <- here::here("data", "Phyllip_GeoAreas.txt")
path_tree <- here::here("data", "processed", "Tree_TF400Howard_Pruned.tre")
tree_pruned <- ape::read.tree(here::here("data",
                                         "processed",
                                         "Tree_TF400Howard_Pruned.tre"))

# load BioGeoBEARS output

res_DEC <- readRDS(here::here("output", "output_BioGeoBEARS", "model_DEC.rds"))

# Running biogeographic stochastic mapping

bsm_res_DEC <- 
  Herodotools::calc_bsm(BioGeoBEARS.data = res_DEC,
                        phyllip.file = path_phyllip_file_birds, 
                        tree.path = path_tree, 
                        max.maps = 100, 
                        n.maps.goal = 100,
                        seed = 42)

# saving bsm 
saveRDS(bsm_res_DEC, 
        here::here("output", "output_BioGeoBEARS", "bsm_res_DEC.rds"))


# Using functions from Herodotools to insert nodes and calculate some PE

# inserting anagenetic nodes -------------------------------------------------------
# here we produce a list with anagenetic nodes in the tree
insert_list <- 
  Herodotools::get_insert_df(
    bsm_res_DEC,
    phyllip.file = path_phyllip_file_birds,
    max.range.size = res_DEC$inputs$max_range_size
  )


# getting the ancestral range area for each node including anagenetic nodes
list_node_area <- 
  Herodotools::get_bsm_node_area(
    bsm = bsm_res_DEC, 
    BioGeoBEARS.data = res_DEC,
    phyllip.file = path_phyllip_file_birds,
    tree.path = path_tree,
    max.range.size = res_DEC$inputs$max_range_size
  )

# insert_nodes ----
# here the anagenetic nodes are inserted in the tree

bsm_tree <- insert_nodes(
  tree =  tree_pruned, 
  inserts = insert_list, 
  node_area = list_node_area)

# calculating PE with one tree

occ_birds <- read.table(here::here("data", "processed", "W_harvey.txt"), 
                        header = T)
bioregions_birds <- read.table(here::here("data", "matrixEco.txt"), 
                               header = T)
bioregions_birds2 <- data.frame(area = bioregions_birds[, 2])
rownames(bioregions_birds2) <- bioregions_birds$ID

# one bsm just for testing
PE_PD_insitu <- 
  Herodotools::calc_insitu_metrics(W = occ_birds, 
                                   tree = tree_pruned, 
                                   ancestral.area = bsm_tree[[42]]$node_area, 
                                   biogeo = bioregions_birds2)
# saving one output just for testing and downstream pipeline
saveRDS(PE_PD_insitu, 
        here::here("output", 
                   "PE_PD_insitu",
                   "PD_PE_insitu_one_bsm.rds")
        )

# subsetting trees
bsm_tree_sub <- bsm_tree[1:3]

# testing new function using only one realization of bsm
system.time(
  test_PD_PE <- 
    calc_PD_PE_insitu(W = occ_birds, 
                      tree = tree_pruned, 
                      ancestral.area = bsm_tree[[42]]$node_area, 
                      biogeo = bioregions_birds2)
)


# specifying the worker plan - multicore in this case
# number of cores
ncores <- parallel::detectCores() - 4 
plan(multisession, workers = ncores)
# Running for the 50 BSM - parallel computation 


# testing progress with future lapply
progressr::with_progress({
  p <- progressr::progressor(steps = length(bsm_tree_sub))
  group_perm <- 
      future.apply::future_lapply(bsm_tree_sub, function(x) {
        res_PE_PD <- 
          calc_PD_PE_insitu(W = occ_birds, 
                            tree = x$phylo, 
                            ancestral.area = x$node_area, 
                            biogeo = bioregions_birds2) 
        
        # update progress
        p(message = sprintf("max.nclust=%s perm=%s", x))
        return(res_PE_PD)
      }, future.seed=TRUE)
})

plan(sequential)

# using for
pb = txtProgressBar(min = 0, max = length(bsm_tree), initial = 0) 
res_all_bsm_PE <- vector(mode = "list", length = length(bsm_tree))
for(i in 1:length(bsm_tree)){
  res_all_bsm_PE[[i]] <- 
    calc_PD_PE_insitu(W = occ_birds, 
                      tree = tree_pruned, 
                      ancestral.area = bsm_tree[[i]]$node_area, 
                      biogeo = bioregions_birds2)
  setTxtProgressBar(pb,i)
}



