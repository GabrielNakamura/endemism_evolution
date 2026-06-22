

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

# Running for the 50 BSM - parallel computation 

# number of cores
ncores <- parallel::detectCores() - 1 

# subsetting trees
bsm_tree_sub <- bsm_tree[1:50]

# specifying the worker plan - multicore in this case
plan(multisession, workers = ncores)

# using progress package to monitor the progress of calculation for each bsm
with_progress({

  p <- progressor(steps = length(bsm_tree_sub))
  
  list_PE_PD_res <- 
    future_map(bsm_tree_sub, function(x){
      # phylogenetic tree
      tree <- x$phylo
      # ancestral area
      anc_area <- x$node_area
      # PE and PD calculation 
      res_PE_PD <- 
        Herodotools::calc_insitu_metrics(W = occ_birds, 
                                         tree = tree, 
                                         ancestral.area = anc_area, 
                                         biogeo = bioregions_birds2) 
      p() # progress bar
      return(res_PE_PD) # output
    }) # future worker
}) # progress counter

# end of multisession plan
plan(sequential)




