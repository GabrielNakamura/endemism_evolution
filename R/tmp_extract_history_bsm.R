insert_list <- get_insert_df(
  bsm_res_DEC,
  phyllip.file = path_phyllip_file_birds,
  max.range.size = res_DEC$inputs$max_range_size
)


# getting the ancestral range area for each node 
list_node_area <- get_bsm_node_area(
  bsm = bsm_res, 
  BioGeoBEARS.data = resDECJ,
  phyllip.file = geog.path,
  tree.path = phy.path,
  max.range.size = resDECJ$inputs$max_range_size
)

# insert_nodes ----

bsm_tree <- insert_nodes(
  tree =  myrcia_tree, 
  inserts = insert_list, 
  node_area = list_node_area)


phangorn::Ancestors(x = phy_pruned, node = "Platyrinchus_flavigularis")
phy_pruned$tip.label
edg <- phy_pruned$edge[match(1:length(phy_pruned$tip.label), 
                                  phy_pruned$edge[, 2]), ]
dim(tip_edge)

edge2 <- phy_pruned$edge[match(tip_edge[, 1], phy_pruned$edge[, 2]), ]

edge3 <- phy_pruned$edge[match(edge2[, 1], phy_pruned$edge[, 2]), ]

edge4 <- phy_pruned$edge[match(edge3[, 1], phy_pruned$edge[, 2]), ]

lapply(1:phy_pruned$Nnode, function(x){
  
})

desc_mat <- matrix(NA,
                   nrow = length(phy_pruned$tip.label),
                   ncol = 901)

edg <- phy_pruned$edge[match(1:length(phy_pruned$tip.label), 
                             phy_pruned$edge[, 2]), ]
for(i in 1:900){
  # i = 2
  edg <- phy_pruned$edge[match(edg[, 1], 
                        phy_pruned$edge[, 2]), ]
  desc_mat[, i+1] <- edg[, 1]
}

desc_mat[, 1] <- phy_pruned$edge[match(1:length(phy_pruned$tip.label), 
                                       phy_pruned$edge[, 2]), ][, 1]

while(sum(edg, na.rm = TRUE) > 0){
  edg <- phy_pruned$edge[match(edg[, 1], 
                               phy_pruned$edge[, 2]), ]
  desc_mat[, ] <- edg[, 1]
  print(i)
}


1:phy_pruned$Nnode

phangorn::Ancestors(x = phy_pruned, node = phy_pruned$tip.label[384])

1:length(phy_pruned$tip.label)


### extracting biogeographic history of lineages

# using sampling 10 of bsm 
tmp_df_insert <- insert_list[[5]]

# using species 1 as an example
anc_list <- 
  lapply(1:length(phy_pruned$tip.label), function(x){
  c(x, phangorn::Ancestors(x = phy_pruned, node = phy_pruned$tip.label[x]))
})
# tmp_anc_sp1 <- phangorn::Ancestors(x = phy_pruned, node = phy_pruned$tip.label[1])
# tmp_anc_sp1_1 <- c(1, tmp_anc_sp1)

tbl_events_1 <- 
  lapply(anc_list, function(x){
  tmp_df_insert[unique(c(match(x, tmp_df_insert$child),
                         match(x, tmp_df_insert$parent)
  )), ]
})

tbl_events_2 <- 
  lapply(1:length(tbl_events_1), function(x){
  tbl_events_1[[x]] |> 
    filter(child %in% anc_list[[x]] & parent %in% anc_list[[x]])
})

names(tbl_events_2) <- phy_pruned$tip.label

tbl_one_bsm <- Map(cbind, tbl_events_2, species = names(tbl_events_2))

tbl_one_bsm <- do.call(rbind, tbl_one_bsm)

rownames(tbl_one_bsm) <- NULL

# all bsm 
lapply(1:length(insert_list), function(x){
  # getting all ancestors from each tip label
  tmp_anc_sp1 <- phangorn::Ancestors(x = phy_pruned, node = phy_pruned$tip.label[x])
  tmp_anc_sp1_1 <- c(x, tmp_anc_sp1)
  tmp_tbl_events_ex1 <- tmp_df_insert[unique(c(match(tmp_anc_sp1_1, tmp_df_insert$child),
                                               match(tmp_anc_sp1_1, tmp_df_insert$parent)
  )), ]
  tmp_tbl_events_ex2 <- 
    tmp_tbl_events_ex1 |> 
    filter(child %in% tmp_anc_sp1_1 & parent %in% tmp_anc_sp1_1)
  
  
})