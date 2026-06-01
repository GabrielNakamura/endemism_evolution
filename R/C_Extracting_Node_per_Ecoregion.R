# reading data and libraries ----------------------------------------------
load(here::here("output", "resBAYAREALIKE.RData"))


# extracting probabilites of states/range at each node --------------------
prob_state <- resBAYAREALIKE$ML_marginal_prob_each_state_at_branch_bottom_below_node

# teste <- resBAYAREALIKEj$ML_marginal_prob_each_state_at_branch_top_AT_node

# In this table:
# - columns are states/ranges
# - rows are nodes, in APE order (tips, then root, then internal)

#  You can see the node numbers in the same APE order with:
trtable = prt(tr, printflag=FALSE)

# Get your states list (assuming, say, 4-area analysis, with max. rangesize=4)
max_range_size = 3
areas = getareas_from_tipranges_object(tipranges)

# This is the list of states/ranges, where each state/range
# is a list of areas, counting from 0
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)

# Make the list of ranges
ranges_list = NULL
for (i in 1:length(states_list_0based)) {    
  if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
  {
    tmprange = "_"
  } else {
    tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
  }
  ranges_list = c(ranges_list, tmprange)
}

# Look at the ranges list
ranges_list

# Make the node numbers the row names
# Make the range_list the column names
range_probabilities = as.data.frame(prob_state)
row.names(range_probabilities) = trtable$node
names(range_probabilities) = ranges_list

# Write the table to a tab-delimited text file (for Excel etc.)
write.table(range_probabilities, file= here::here("data", "processed", 'range_probabilities.csv'), sep=';', dec=',', row.names=FALSE)

#############################################################################################################
#################################Fazer matrix eco Nodes######################################################
#############################################################################################################

# matrix Node x Biome -----------------------------------------------------

teste_max <- apply(range_probabilities,
                   MARGIN = 1, 
                   function(i) colnames(range_probabilities)[which(i == max(i))
                                                             ]
                   )

biome <- data.frame(rows_name = 1:768, biome = unlist(teste_max))


# saving table ------------------------------------------------------------
write.table(biome, file = here::here("data", "processed", 'EcoNodes.csv'), sep=';', dec=',', row.names=FALSE)

write.table(biome, file = here::here("data", "processed", "Econodes.txt"))
