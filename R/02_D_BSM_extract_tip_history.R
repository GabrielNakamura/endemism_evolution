

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



BSMs_w_sourceAreas <- 
  simulate_source_areas_ana_clado(res_DEC, 
                                  bsm_res_DEC$RES_clado_events_tables, 
                                  bsm_res_DEC$RES_ana_events_tables, 
                                  names(tipranges@df))

# Simulated anagenetic and cladogenetic events 
ana_events_tables <- BSMs_w_sourceAreas$ana_events_tables
clado_events_tables <- BSMs_w_sourceAreas$clado_events_tables

# count the number of events
n_events <- 
  count_ana_clado_events(clado_events_tables, 
                         ana_events_tables, 
                         names(tipranges@df), 
                         names(tipranges@df))


all_counts_events <- n_events$summary_counts_BSMs 


# calculating event per tip -----------------------------------------------

# Assuming 'resDEC' contains your original ML fit (which holds the tree)
# and 'ana_events_tables' contains your BSM anagenetic results.

tip_labels <- tree_pruned$tip.label
num_sims <- length(ana_events_tables)

# Create a master dataframe to store the results
tip_events_summary <- 
  data.frame(
    tip = tip_labels,
    mean_dispersals = NA,
    mean_extirpations = NA,
    mean_vicariance = NA,
    mean_simpatry = NA,
    stringsAsFactors = FALSE
  )

# Loop through each tip to calculate its unique historical pathway
for(i in 1:length(tip_labels)) {
  # i = 1
  tip_node = i
  
  # Get all ancestral nodes leading from the root to this specific tip
  # nodepath() from 'ape' gives the sequence of nodes from root to tip
  root_node <- length(tree_pruned$tip.label) + 1
  path_nodes <- nodepath(tree_pruned, from = root_node, to = tip_node)
  
  # Accumulate event counts across all BSM simulation maps
  total_dispersals_for_tip = 0
  total_extirpations_for_tip = 0
  
  for(s in 1:num_sims) {
    # s = 10
    sim_table <- ana_events_tables[[s]]
    
    # In BioGeoBEARS ana_events_tables:
    # 'brlen' tells us branch lengths, and we can match branches using edge matrix indices
    # Or more directly, we can match rows where the 'headed_node' matches our path nodes
    
    # Filter the simulation table for branches that belong to this tip's lineage pathway
    # Note: The root node itself doesn't have an incoming branch, so we look at the tip's ancestral path nodes excluding the root
    lineage_branches <- sim_table[sim_table$node %in% path_nodes[-1], ]
    
    # Anagenetic events:
    # Sum the events along this pathway for this specific simulation map
    # 'd' and 'e' events are the standard naming for dispersal and extirpation respectively
    total_dispersals_for_tip <- 
      total_dispersals_for_tip + 
      sum(lineage_branches$event_type == "d", na.rm = TRUE) 
    total_extirpation_for_tip <- 
      total_dispersals_for_tip + 
      sum(lineage_branches$event_type == "e", na.rm = TRUE) 
    # Note: If checking pure anagenetic columns, verify column names via colnames(sim_table)
    # Frequently stored as specific event types depending on your BioGeoBEARS version:
    # e.g., sum(lineage_branches$anagenetic_events == "dispersal")
    
    # Cladogenetic events:
    # sum of the events along the species'pathway for a given specific bsm simulation
    # "subset (s)" and "vicariance (v)" are the names for the events in the standard DEC model
    # 
    total_dispersals_for_tip <- 
      total_dispersals_for_tip + 
      sum(lineage_branches$event_type == "d", na.rm = TRUE) 
    total_extirpation_for_tip <- 
      total_dispersals_for_tip + 
      sum(lineage_branches$event_type == "e", na.rm=TRUE) 
  }
  
  # Calculate the average across all stochastic maps
  tip_events_summary$mean_dispersals[i] <- 
    total_dispersals_for_tip / length(lineage_branches$event_type)
  tip_events_summary$mean_extirpations[i] <- 
    total_extirpations_for_tip / length(lineage_branches$event_type)
}

# View your final per-tip breakdown
print(tip_events_summary)

# Using functions from Herodotools to insert nodes and calculate some historical metrics

# assemblage age and in situ diversification ---------------------------------

# prepare the insertions -------------------------------------------------------
## get_insert_df ----

insert_list <- 
  Herodotools::get_insert_df(
    bsm_res_DEC,
    phyllip.file = path_phyllip_file_birds,
    max.range.size = res_DEC$inputs$max_range_size
  )


# getting the ancestral range area for each node 
list_node_area <- 
  Herodotools::get_bsm_node_area(
    bsm = bsm_res_DEC, 
    BioGeoBEARS.data = res_DEC,
    phyllip.file = path_phyllip_file_birds,
    tree.path = path_tree,
    max.range.size = res_DEC$inputs$max_range_size
  )

# insert_nodes ----

bsm_tree <- insert_nodes(
  tree =  phy_pruned, 
  inserts = insert_list, 
  node_area = list_node_area)

# calculating PE with one tree

occ_birds <- read.table(here::here("data", "processed", "W_harvey.txt"), 
                        header = T)
bioregions_birds <- read.table(here::here("data", "matrixEco.txt"), 
                               header = T)
bioregions_birds2 <- data.frame(area = bioregions_birds[, 2])
rownames(bioregions_birds2) <- bioregions_birds$ID
test_rds <- readRDS(here::here("data", "processed", "grid_tyranidae.rds"))

# one bsm just for testing
PE_PD_insitu <- 
  Herodotools::calc_insitu_metrics(W = occ_birds, 
                                   tree = phy_pruned, 
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




