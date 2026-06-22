library(dplyr)
W = occ_birds
tree = tree_pruned
ancestral.area = bsm_tree[[42]]$node_area
biogeo = bioregions_birds2
PD = TRUE
PE = TRUE

library(Herodotools)
W_toy <- matrix(c(0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0),
               nrow= 3,
               ncol= 5,
               dimnames=list(c("Comm 1", "Comm 2", "Comm 3"),
                             c(paste("s", 1:5, sep=""))))
data("toy_treeEx")
biogeo_toy <- data.frame(Ecoregion= c("A", "B", "C"))
ancestral_area_toy <- data.frame(state= c("ABC", "B", "C", "ABC"))
assemblage_phylo_metrics <- calc_insitu_metrics(W_toy,
                                                toy_treeEx,
                                                ancestral_area_toy, 
                                                biogeo_toy)
# edditing a bit occurrence matrix to create two assemblages in the same region
comm_add <- c(1, 1, 0, 0, 0)
W_toy2 <- rbind(W_toy, comm_add)
rownames(W_toy2)[4] <- "Comm 4"
biogeo_toy2 <- rbind(biogeo_toy, "C")
W = W_toy2
tree = toy_treeEx
ancestral.area = ancestral_area_toy
biogeo = biogeo_toy2

calc_PD_PE_insitu <- 
  function (W, tree, ancestral.area, biogeo, PD = TRUE, PE = TRUE) 
  {
    if (!is.matrix(W) == TRUE) {
      if (is.data.frame(W) == TRUE) {
        W <- as.matrix(W)
        rownames(W) <- 1:nrow(W)
      }
      else {
        stop("W must be a occurrence matrix with presences (1) or absences (0)")
      }
    }
    
    # reorganizing occurrence matrix to match phylogeny tip order
    W <- W[, match(tree$tip.label, colnames(W))]
    
    # calculating metrics without considering ancestral area 
    PDt <- picante::pd(samp = W, tree = tree)$PD
    PEt <- phyloregion::phylo_endemism(x = W, phy = tree, weighted = TRUE) # strict endemism sensu Rosauer
    # getting a list with each element representing a assemblage and each element in 
    # the asseblage representing a nodepath with insitu diversification 
    # if no insitu NA, $disp.anc.node represents the region from where the species came
    # in case of in situ diversification this element will be NA
    nodes.list <- get_nodes_info_core(W = W, tree = tree, ancestral.area = ancestral.area, 
                                      biogeo = biogeo)
    
    # keeping only species that represent insitu diversification path
    names_spp_noNull <- 
      lapply(
        lapply(
          lapply(nodes.list, function(x) {
            is.na(x$nodes_species)
          }), function(x) {
            which(x == FALSE)
          }), names)
    # list with each element as a assemblage and the insitu nodes for each species
    nodes_species_noNull <- vector(mode = "list", length = nrow(W))
    for (i in 1:length(names_spp_noNull)) {
      nodes_species_noNull[[i]] <- nodes.list[[i]]$nodes_species[names_spp_noNull[[i]]]
    }
    nodes_species_noNull_org <- nodes_species_noNull
    
    # a list with edge matrix for all species in each assemblage (elements of the list)
    list_matrix_nodes <- vector(mode = "list", length = length(nodes_species_noNull_org))
    for (i in 1:length(nodes_species_noNull)) {
      if (length(nodes_species_noNull[[i]]) == 0) {
        list_matrix_nodes[[i]] <- NA
      }
      else {
        names_spp <- names(nodes_species_noNull[[i]])
        list_nodes_org <- vector(mode = "list", length = length(names_spp))
        for (j in 1:length(names_spp)) {
          matrix_nodesSpp_nonull <- matrix(NA, nrow = round(length(unlist(nodes_species_noNull[[i]][j]))), 
                                           ncol = 2)
          nodes_org <- c(sort(nodes_species_noNull[[i]][j][[1]], 
                              decreasing = FALSE), which(tree$tip.label == 
                                                           names_spp[j]))
          for (k in 1:nrow(matrix_nodesSpp_nonull)) {
            matrix_nodesSpp_nonull[k, ] <- c(nodes_org[k], 
                                             nodes_org[k + 1])
          }
          list_nodes_org[[j]] <- matrix_nodesSpp_nonull
        }
        list_matrix_nodes[[i]] <- list_nodes_org
      }
    }
    
    # joining edges for each assemblage - this corresponds only to insitu edges
    list_matrix_edges <- vector(mode = "list", length = length(list_matrix_nodes))
    for (i in 1:length(list_matrix_nodes)) {
      if (any(is.na(list_matrix_nodes[[i]])) == TRUE) {
        list_matrix_edges[[i]] <- NA
      }
      else {
        list_matrix_edges[[i]] <- unique(do.call(rbind, list_matrix_nodes[[i]]))
      }
    }
    # naming the in situ edge list
    names(list_matrix_edges) <- rownames(W)
    list_df_edges <- 
      lapply(1:length(list_matrix_edges), function(x){
        if (any(is.na(list_matrix_edges[[x]])) == TRUE) {
          NA
        } else{
          data.frame(parent = list_matrix_edges[[x]][, 1],
                     child = list_matrix_edges[[x]][, 2],
                     site = rep(names(list_matrix_edges[x]), 
                                dim(list_matrix_edges[[x]])[1])
          )
        }
      })
    df_insitu_edges <- do.call(rbind, list_df_edges)
    
    # getting denominator for each site
    list_denominator_insitu <- 
      apply(df_insitu_edges, MARGIN = 1, function(x){
        which(as.numeric(x[c(1, 2)][1]) == df_insitu_edges[, 1] & 
                as.numeric(x[c(1, 2)][2]) ==  df_insitu_edges[, 2]) 
      })
    vec_denominator_insitu <- 
      unlist(lapply(list_denominator_insitu, 
                    function(x) length(x)))
    
    # data frame with all insitu edges and their incidence in grids 
    df_insitu_edges2 <- data.frame(df_insitu_edges,
                                   denominator = vec_denominator_insitu) 
    
    # getting branch lengths for in situ edges
    list_brLength_divLocal <- 
      vector(mode = "list", length = length(list_matrix_edges))
    for (i in 1:length(list_matrix_edges)) {
      # i = 500
      if (any(is.na(list_matrix_edges[[i]])) == TRUE) {
        list_brLength_divLocal[[i]] <- NA
      }
      else {
        pos_insitu_nodes <- 
          apply(tree$edge, MARGIN = 1, 
                function(x) {
                  which(x == apply(as.matrix(list_matrix_edges[[i]]), 
                                   MARGIN = 1, function(l) l))
                })
        list_brLength_divLocal[[i]] <- 
          tree$edge.length[which(unlist(lapply(pos_insitu_nodes, 
                                               function(x) length(x) > 1)) == TRUE)]
      }
      # print(i)
    }
    
    # adding insitu branch lenghts to the insitu edges data frame
    df_insitu_edges3 <- 
      data.frame(df_insitu_edges2, 
                 branch.lengths = unlist(list_brLength_divLocal))
    
    # calculating weighted phylogenetic endemism for each branch 
    df_insitu_edges4 <- 
      df_insitu_edges3 |> 
      dplyr::mutate(branch.lengths.weighted = branch.lengths/denominator)
    
    # calculating PD and PE insitu metrics for each assemblage/site
    # PD insitu
    PDinsitu <- unlist(lapply(list_brLength_divLocal, function(x) {
      sum(x)
    }))
    names(PDinsitu) <- rownames(W)
    
    # PE insitu
    df_insitu_edges5 <- 
      df_insitu_edges4 |> 
      group_by(site) |> 
      mutate(PE_insitu = sum(branch.lengths.weighted))
    
    # keeping only one PEinsitu value per site 
    df_insitu_sites <- 
      df_insitu_edges5 |> 
      distinct(site, .keep_all = T)
    
    # All PE 
    PEtotal <- data.frame(PEt = as.numeric(PEt), site = names(PEt))
    # PE insitu for sites with any insitu diversification
    PEinsitu_subset <- 
      df_insitu_sites |> 
      select("site", "PE_insitu")
    # joining PE total with PEinsitu to keep the size of the original occ table
    insitu_comm_metrics <- 
      PEtotal |> 
      left_join(PEinsitu_subset, by = c(site = "site"))
    # joining with PD and PDinsitu
    insitu_comm_metrics <- 
      data.frame(insitu_comm_metrics, PD = PDt, PDinsitu = PDinsitu)
    
    return(insitu_comm_metrics) # all metrics
  }