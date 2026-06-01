#' Diversification Based (DivB) Metrics for Community Phylogenetics
#' 
#' Compute in each site DivB metrics for diversification (Jetz and Frecklerton), 
#' diversity metrics (Phylogenenic diversity and phylogenetic endemism), 
#' age of arrival in the biogeographic area, and frequency of dispersal 
#' from other areas 
#' 
#' @param W community matrix (site x species)
#' @param tree phylogenetic tree
#' @param ancestral.area biogeographic area for ancestral. Must be a single 
#'   column data.frame in which rows are nodes in the tree. Each area must be a
#'   single letter that could be combined (eg. "A", "AB", "B" ).
#' @param biogeo sinlge collum data.frame with biogeographic area for each site. 
#' @param diversification character. Diversification index sensu Jetz et al. or 
#'   sensu Frecklerton et al. Defaut = c("jetz", "freck").
#' @param PD logical. calculate phylogenetic diversity and its diversification based version (PDlocal). Defaut = TRUE.
#' @param PE logical. calculate phylogenetic endemism and its diversification based version (Db-PE). Defaut = TRUE.
#' @param age.arrival logical. Calculate mean age of arrival of species in the 
#'   site's biogeographic area. Defaut = TRUE.
#' @param age.no.ancestor how to deal with age arrival for species in which the 
#'   late ancestor wasn't in the site's biogeographic area. Options: NA, 
#'   numeric(), or 'half.edge'. Defalt = NA.  See Details.
#' @param dispersal.from logical. Calculate frequecy of dipersal from each 
#'   biogeografic area of ancestors for each site 
#'   
#' 


DivB_metrics <-function(W,
                        tree,
                        ancestral.area,
                        biogeo, 
                        diversification = c("jetz", "freck"),
                        PD = TRUE,
                        PE = TRUE,
                        age.arrival = TRUE,
                        age.no.ancestor = NA,# 'half.edge' or numeric()
                        dispersal.from = TRUE, 
                        ED.type = "equal.splits"
){
  library(geiger)
  library(ape)
  library(ade4)
  library(phylobase)
  library(picante)
  library(sjmisc)
  
  s <- length(tree$tip.label)
  n <- tree$Nnode
  names_spComm <- colnames(W)
  
  
  # Total values for metrics ------------------------------------------------
  
  
  PDt<- picante::pd(samp = W, tree = tree, include.root = T)$PD
  PEt<- phyloregion::phylo_endemism(x = W, phy = tree, weighted = TRUE) #Phylogenetic endemism sensu Rosauer
  type<- ED.type #argument to compute equal splits metrics in ER calculation - internal function
  EDtotal<- evol.distinct(tree = tree, type = type)
  Jetz_total<- 1/EDtotal$w
  Freck_total<- numeric(length = length(tree$tip.label))
  names(Freck_total)<- tree$tip.label
  T_freck<- max(cophenetic(tree))/2
  for(l in 1:length(tree$tip.label)){
    nodes_Freck<- .get.nodes(tree, tree$tip.label[l]) 
    Freck_total[l]<- length(nodes_Freck)/T_freck
  }
  names(Jetz_total)<- tree$tip.label
  
  # species vs Node matrix --------------------------------------------------
  
  
  spxnode <- matrix(0,s,n)
  ages<-abs(node.depth.edgelength(phy = tree)
            -max(node.depth.edgelength(tree)))[-c(1:length(tree$tip.label))]
  ages<- data.frame(age= ages)
  rownames(ages)<- paste("N",(s+1):(s+(s-1)),sep="")
  tree_test1 <- suppressWarnings(phylo4(tree)) # precisa converter para esse formato para a funcao ancestors() funcionar
  for (i in 1:s){
    #i=1
    nos_ancestrais <- ancestors(tree_test1,i)
    nos_ancestrais <- nos_ancestrais - s  
    spxnode[i,nos_ancestrais] <- 1
  }
  rownames(spxnode)<-tree$tip.label #names for columns 
  colnames(spxnode)<- paste("N",(s+1):(s+(s-1)),sep="") #names for Nodes
  spxnode<- t(spxnode) #mesmo que Node
  
  # Ancestral species matrix [AS] -------------------------------------------
  
  
  AS<-spxnode #Ancestral State matrix 
  for(i in 1:nrow(spxnode)){
    pres<-which(spxnode[i,]==1)
    AS[i,pres]<-as.character(ancestral.area[i,1]) #matriz - a informacao em cada celula representa o estado do ancestral de cada especie
  }
  for(i in 1:nrow(AS)){
    zero<-which(AS[i,]==0)
    AS[i,zero]<-NA #NAs indicam nos que nao contem a especie 
  }
  
  
  # Decomposing local metrics -----------------------------------------------
  
  #matrix to receive the results of local Jetz metric
  matrix_XJetz<- matrix(0, 
                        nrow = nrow(W), 
                        ncol = ncol(W), 
                        dimnames = list(rownames(W), colnames(W))) 
  #matrix to receive the results of local Freckleton metric
  matrix_XFreck<- matrix(0, 
                         nrow = nrow(W), 
                         ncol = ncol(W), 
                         dimnames = list(rownames(W), colnames(W))) 
  #matrix to receive the age from which local diversification occurs in the area
  age_arrival<-  matrix(0, 
                        nrow = nrow(W), 
                        ncol = ncol(W), 
                        dimnames = list(rownames(W), colnames(W))) 
  
  #matrix to receive the area from which dispersal ocurred
  dispersal_from <-  matrix("-", 
                            nrow = nrow(W), 
                            ncol = ncol(W), 
                            dimnames = list(rownames(W), colnames(W))) 
  
  #matrix to receive the results of local PD metric
  PD_local<- matrix(NA, 
                    nrow= nrow(W), 
                    ncol= 1, 
                    dimnames= list(rownames(W), "PD_local"))
  
  
  # List with node path and dispersal node ----------------------------------
  
  nodes.list <- lapply(1:nrow(W), function(i){
    pres<- which(W[i,]==1)
    pres<- names_spComm[pres]
    nodes_species <- vector(mode= "list")
    disp.anc.node <- vector("numeric", length = length(pres))
    #j=1
    
    for(j in 1:length(pres)){
      nodes_sp<- AS[,pres[j]][!is.na(AS[,pres[j]])]
      nodes_sp<- nodes_sp[length(nodes_sp):1] #nodes for species j in community i
      
      nodes.T <- grep(biogeo[i,1], nodes_sp)
      nodes_all<- numeric(length = length(nodes_sp))
      nodes_all[nodes.T] <- 1
      
      if(all(nodes_all==1)){ #if all ancestors of species j are in the same ecoregion of local 1 this will be TRUE
        x<- names(nodes_sp[length(nodes_sp)]) # if TRUE, take basal node as the reference node for calculation of local diversification
        
        rec.anc.node <-as.numeric(substr(x, 2, nchar(x)))
        
        n.path <- nodepath(tree, 
                           rec.anc.node, 
                           which(tree$tip.label==pres[j]))
        
        nodes_species[[j]] <- sort(n.path[-length(n.path)])
        disp.anc.node[j] <- NA
        
        
      }else
      { # node for the ancenstor previous to dispersal
        out.situ.anc.pos <- which(nodes_all == 0)[1]
        x <- names(nodes_sp)[out.situ.anc.pos]
        disp.anc.node[j] <- as.numeric(substr(x, 2, nchar(x)))
        
        # node for the early ancestor in the same ecoregion
        x1 <- names(nodes_sp)[out.situ.anc.pos - 1]
        rec.anc.node <- ifelse(length(x1) == 1, 
                               yes = as.numeric(substr(x1, 2, nchar(x1))),
                               no = NA)
        if(!is.na(rec.anc.node)){
          n.path <- nodepath(tree, 
                             rec.anc.node, 
                             which(tree$tip.label==pres[j]))
          
          nodes_species[[j]] <- sort(n.path[-length(n.path)])
        }else{ nodes_species[[j]] <- NA}
        
      }
    }
    names(nodes_species) <- pres
    names(disp.anc.node) <- pres
    list(nodes_species = nodes_species,
         disp.anc.node = disp.anc.node)
  })
  
  
  
  # Jetz diversification index -----------------------------------------------
  
  
  if(any(diversification == "jetz")){
    
    ## Jetz local para cada espécie
    l.Jetz.local <- lapply(nodes.list , function(site){
      lapply(1:length(site$nodes_species), function(j){
        nodes_div <- site$nodes_species[[j]]
        sp <- names(site$nodes_species)[j]
        
        if(length(nodes_div)==1){
          internal.brlen_div <- tree$edge.length[which(tree$edge[, 
                                                                 2] %in% nodes_div)] #branch lenghts (in times) for internal branch lengts of the most ancient ancestral 
        } else{
          nodes_div <- nodes_div[-1] #internal nodes that form the path from most ancient ancestral that was presented at Ecoregion of local i to species j
          internal.brlen_div <- tree$edge.length[which(tree$edge[, 
                                                                 2] %in% nodes_div)] #ages for internal nodes of nodes_div object
        }
        if (length(internal.brlen_div) != 0) { #starts the calculation for ED measure from Redding et al.
          internal.brlen_div <- internal.brlen_div * switch(type, equal.splits = sort(rep(0.5, 
                                                                                          length(internal.brlen_div))^c(1:length(internal.brlen_div))), 
                                                            fair.proportion = {
                                                              for (j in 1:length(nodes_div)) {
                                                                sons <- .node.desc(tree, nodes_div[j])
                                                                n.descendents <- length(sons$tips)
                                                                if (j == 1) portion <- n.descendents else portion <- c(n.descendents, 
                                                                                                                       portion)
                                                              }
                                                              1/portion
                                                            })
          
        }
        
        if(any(is.na(nodes_div))){
          JetzLocal <- 0.00001
        }else{
          ED_div <- sum(internal.brlen_div, 
                        tree$edge.length[which.edge(tree, sp)]) #Local ED - modifyed ED considering only the edges of ancestros inside the biogeo of local i for species j
          
          EDtotal_spp <- EDtotal$w[which(EDtotal$Species==sp)]
          Jetz_total_spp<-Jetz_total[sp] #Diversificatio calculated according to Jetz for species j
          JetzLocal<- (Jetz_total_spp*(ED_div/EDtotal_spp)) #Modified Jetz to calculate only for local diversification
          
        }
        JetzLocal #Jetz local diversification
        
      })#lapply j
    })#lapply site
    
    # populate the matrix
    for(site in 1:length(l.Jetz.local)){
      for(sp in 1:length(l.Jetz.local[[site]])){
        pres<- which(W[site,]==1)
        pres<- names_spComm[pres]
        matrix_XJetz[site, pres[sp]] <- l.Jetz.local[[site]][[sp]]
      }
    }
    ## FIM Jetz local para cada espécie
    
    #### Summarised results for sites
    ##calculating harmonic mean for Jetz
    
    #sum of communities local diversification
    sum_localDiv<- apply(matrix_XJetz, 
                         1, 
                         function(x) sum(x[which(x!=0 & x!=0.00001)]))
    
    
    #objet to receive Jetz diversification values
    matrix_totalDiv_Jetz<- W 
    
    #substitui a matrix de ocorrencia pelos valores de div total do Jetz
    for(i in colnames(W)){ 
      #i=1
      matrix_totalDiv_Jetz[,i]<- ifelse(matrix_totalDiv_Jetz[,i]>=1,Jetz_total[i],0) #input Jetz total diversification values in occurence matrix
    }
    
    #sum of communities total diversification
    sum_totalDiv<- apply(matrix_totalDiv_Jetz, 
                         1, 
                         function(x) sum(x[which(x!=0)])) 
    
    numerator_values<- apply(matrix_totalDiv_Jetz, 1, function(x){
      sum(1/x[which(x!=0)])
    })
    denom_values<- apply(matrix_totalDiv_Jetz, 1, function(x){
      length(which(x!=0))
    })
    
    #harmonic mean for Jetz total diversification
    JetzTotalComm_harmonic <- denom_values/numerator_values 
    #harmonic mean for Jetz local diversification
    JetzLocalComm_harmonic <- (JetzTotalComm_harmonic*(sum_localDiv/sum_totalDiv)) 
    
  }#if(any(diversification == "jetz"))
  
  
  # Frecklerton diversification Index ----------------------------------------
  
  
  if(any(diversification == "freck")){
    ## Freckleton local
    # modifyed equation 4 from Freckleton et al (2008) considering only the nodes
    # that diversified in ecoregion of local i, plus one is only a correction of the 
    # previous step
    
    for(site in 1:length(nodes.list)){
      for(sp in 1:length(nodes.list[[site]]$nodes_species)){
        pres<- which(W[site,]==1)
        pres<- names_spComm[pres]
        nodes_div <- nodes.list[[site]]$nodes_species[[sp]]
        matrix_XFreck[site, pres[sp]] <- ifelse(is.na(nodes_div),
                                                0.00001,
                                                (length(nodes_div))/T_freck)[1]
      }
    }
    
    FreckLocalComm_mean <- sapply(1:nrow(matrix_XFreck), function(i){
      pres <- matrix_XFreck[i,][W[i,]==1]
      mean(pres, na.rm = T)
    })
    
    #objet to receive Jetz diversification values
    matrix_totalDiv_Freck <- W 
    
    #substitui a matrix de ocorrencia pelos valores de div total do Freck
    for(i in colnames(W)){ 
      #i=1
      matrix_totalDiv_Freck[,i]<- ifelse(matrix_totalDiv_Freck[,i]>=1,Freck_total[i],0) 
    }
    
    FreckTotalComm_mean <- sapply(1:nrow(matrix_totalDiv_Freck), function(i){
      pres <- matrix_totalDiv_Freck[i,][W[i,]==1]
      mean(pres, na.rm = T)
    })
    
  }#if(any(diversification == "freck"))
  
  
  # Age arrival -------------------------------------------------------------
  
  if(age.arrival){
    # site; species; fisrt node (root in the ecorregion)
    
    ## NA means the time of arrival is unknown
    for(site in 1:length(nodes.list)){
      for(sp in 1:length(nodes.list[[site]]$nodes_species)){
        pres<- which(W[site,]==1)
        pres<- names_spComm[pres] 
        node <- nodes.list[[site]]$nodes_species[[sp]][1]
        node.name  <- paste0("N", node)
        age_arrival[site,pres[sp]]<- ages[node.name,] #recebe a idade de chegada
      }
    }
    
    ## age.no.ancestor is an option to atribute an age time when diversification 
    ## wasn't local
    if(!is.na(age.no.ancestor)){
      if(is.numeric(age.no.ancestor)){
        age_arrival[is.na(age_arrival)] <- age.no.ancestor
      }
      
      if(age.no.ancestor == "half.edge"){
        for(site in 1:length(nodes.list)){
          for(sp in 1:length(nodes.list[[site]]$disp.anc.node)){
            
            pres<- which(W[site,]==1)
            pres<- names_spComm[pres] 
            if(is.na(age_arrival[site,pres[sp]])){
              node <- nodes.list[[site]]$disp.anc.node[sp]
              node.name  <- paste0("N", node)
              age.disp <- ages[node.name,]
              age.arri <- age_arrival[site,pres[sp]] #recebe a idade de chegada
              age.arri <- ifelse(is.na(age.arri), 0, age.arri)
              
              age_arrival[site,pres[sp]] <- mean(c(age.disp, age.arri))
            }
            
          }
        }
      }
      
    }
    
    mean_age_arrival <- sapply(1:nrow(age_arrival), function(i){
      pres <- age_arrival[i,][W[i,]==1]
      mean(pres, na.rm = T)
    })
    mean_age_arrival <- data.frame(mean_age_arrival = mean_age_arrival,
                                   row.names = row.names(W))
  }#if(age.arrival)
  
  
  # Dispersal From ----------------------------------------------------------
  
  if(dispersal.from){
    ### Age dispersal and the area from
    # site; species; fisrt node (root in the ecorregion)
    ## NA means the dispersal was previous to the root node
    for(site in 1:length(nodes.list)){
      for(sp in 1:length(nodes.list[[site]]$disp.anc.node)){
        pres<- which(W[site,]==1)
        pres<- names_spComm[pres] 
        node <- nodes.list[[site]]$disp.anc.node[sp]
        node.name  <- paste0("N", node)
        if(node.name == "NNA"){ dispersal_from[site,pres[sp]] <- NA
        }else{
          dispersal_from[site,pres[sp]] <- AS[node.name,pres[sp]]}
      }
    }
    
    l.freq.area <- lapply(1:nrow(dispersal_from), function(i){
      pres <- dispersal_from[i,][W[i,]==1]
      table(pres)/sum(table(pres))
    })
    
    
    areas <- unique(ancestral.area[,1])
    disp_from_frequency <- matrix(NA, nrow = nrow(W), ncol = length(areas))
    rownames(disp_from_frequency) <- row.names(W)
    colnames(disp_from_frequency) <- areas
    for(i in 1:nrow(disp_from_frequency)){
      temp.freq<- l.freq.area[[i]]
      temp.col <- colnames(disp_from_frequency) %in% names(temp.freq)
      
      disp_from_frequency[i,temp.col] <- temp.freq
    }
    
  }
  
  #### calculating dbPD and DbPE
  if(PD == TRUE | PE == TRUE){
    names_spp_noNull<- lapply(
      lapply(
        lapply(nodes.list, 
               function(x){
                 is.na(x$nodes_species)
               }
        ), 
        function(x){
          which(x == FALSE)
        }
      ),
      names)
    nodes_species_noNull<- vector(mode = "list", length= nrow(W))
    for(i in 1:length(names_spp_noNull)){
      nodes_species_noNull[[i]]<- nodes.list[[i]]$nodes_species[names_spp_noNull[[i]]]
    }
    
    ####organize node matrix
    nodes_species_noNull_org<- nodes_species_noNull
    list_matrix_nodes<- vector(mode = "list", length= length(nodes_species_noNull_org))
    for(i in 1:length(nodes_species_noNull)){
      #i= 3
      if(length(nodes_species_noNull[[i]]) == 0){
        list_matrix_nodes[[i]]<- NA
      } else{
        names_spp<- names(nodes_species_noNull[[i]])
        list_nodes_org<- vector(mode = "list", length = length(names_spp))
        for(j in 1:length(names_spp)){
          #j=1
          matrix_nodesSpp_nonull<- matrix(NA, nrow= round(length(unlist(nodes_species_noNull[[i]][j]))), ncol= 2)
          nodes_org<- c(sort(nodes_species_noNull[[i]][j][[1]], 
                             decreasing = FALSE), 
                        which(tree$tip.label == names_spp[j]))
          for(k in 1:nrow(matrix_nodesSpp_nonull)){
            matrix_nodesSpp_nonull[k, ]  <- c(nodes_org[k], nodes_org[k + 1])
          }
          list_nodes_org[[j]]<- matrix_nodesSpp_nonull
        }
        list_matrix_nodes[[i]]<- list_nodes_org
      }
    }
    ##########  
    
    #### extracting all internal insitu branches of the tree
    list_matrix_edges<- vector(mode = "list", length= length(list_matrix_nodes))
    for(i in 1:length(list_matrix_nodes)){
      #i= 1001
      if(any(is.na(list_matrix_nodes[[i]])) == TRUE){
        list_matrix_edges[[i]]<- NA
      } else{
        list_matrix_edges[[i]]<- unique(do.call(rbind, list_matrix_nodes[[i]]))
      }
    } 
    ####
    
    #####extracting internal insitu diversification branch length
    list_brLength_divLocal<- vector(mode = "list", length= length(list_matrix_edges))
    for(i in 1:length(list_matrix_edges)){
      if(any(is.na(list_matrix_edges[[i]])) == TRUE){
        list_brLength_divLocal[[i]]<- NA
      } else{
        pos_insitu_nodes<- apply(tree$edge, MARGIN = 1, function(x){
          which(x == apply(as.matrix(list_matrix_edges[[i]]), MARGIN = 1, function(l) l))
        })
        list_brLength_divLocal[[i]]<- tree$edge.length[which(unlist(lapply(pos_insitu_nodes, 
                                                                           function(x) length(x) > 1)
        ) == TRUE)
        ]
      }
      print(paste("run", i, sep= "_"))
    } 
    #####
    
    if(PD == TRUE){
      ######calculation of PDlocal based on internal branch lenghts with insitu diversification#####
      PDlocal<- unlist(lapply(list_brLength_divLocal, 
                              function(x){
                                sum(x)
                              }
      )
      )
      names(PDlocal)<- rownames(W) #PDlocal
    }
    
    if(PE == TRUE){
      ####### Db-PE calculation #####
      names(list_matrix_edges)<- paste("site", 1:nrow(W), sep= "_")
      all_comb<- gtools::permutations(length(list_matrix_edges), 2, 
                                      names(list_matrix_edges))
      list_occurence<- vector(mode= "list", length= length(list_matrix_edges))
      #list_occurence<- list()
      
      #calculating the number of branches co-occurence 
      for(i in 1:length(list_matrix_edges)){
        #i= 3000
        if(is.null(dim(list_matrix_edges[[i]])) == TRUE){
          list_occurence[[i]]<- NA
        } else {
          posit<- which(all_comb[, 1] == names(list_matrix_edges)[i])
          comb_occ<- matrix(NA, nrow = nrow(list_matrix_edges[[i]]), ncol = length(list_matrix_edges)-1)
          for(j in 1:length(posit)){
            if (is.na(list_matrix_edges[[all_comb[posit[j], 2]]]) == TRUE){
              comb_occ[,j]<- FALSE
            } else {
              comb_occ[,j]<- list_matrix_edges[[all_comb[posit[j], 1]]][, 1] %in% list_matrix_edges[[all_comb[posit[j], 2]]][, 1] & 
                list_matrix_edges[[all_comb[posit[j], 1]]][, 2] %in% list_matrix_edges[[all_comb[posit[j], 2]]][, 2]
            }
          }
          list_occurence[[i]]<- comb_occ 
        }
        print(paste("run", i, sep= "_"))
      }
      occurence_area<- lapply(list_occurence, function(x){
        if(is.null(dim(x)) == TRUE){
          NA
        } else{
          rowSums(x) + 1
        }
      }
      ) #numerator for PElocal calculation
      list_edge_brLen_insitu<- vector(mode = "list", length = length(occurence_area))
      for(i in 1:length(occurence_area)){
        list_edge_brLen_insitu[[i]]<- cbind( list_brLength_divLocal[[i]], occurence_area[[i]])
      }
      
      ####correcting denominator of PE for ancestral state occurrence
      node_ancestral<- seq((ncol(W) + 1), ncol(W) + (ncol(W) - 1), by= 1) #names of internal nodes
      ancestral_nodes<- cbind(ancestral.area, node_ancestral) #binding nodes name with ancestral state 
      
      list_edge_insitu<- lapply(list_matrix_edges, function(x){
        if(is.null(dim(x)) == TRUE){
          NA
        } else{
          cbind(x, 
                cbind(rep(NA, nrow(x)),
                      rep(NA,nrow(x))
                )
          ) 
        }
      }
      )
      
      ### ancestral state for each node that compose each community
      ancestral_comm_temp<- lapply(
        lapply(list_matrix_edges, 
               function(x){
                 if(is.null(dim(x)) == TRUE){
                   NA
                 } else {
                   apply(x, MARGIN = 1, function(l){
                     cbind(ancestral_nodes[match(l[1], ancestral_nodes[, 2]), 1], ancestral_nodes[match(l[2], ancestral_nodes[, 2]), 1])
                   }
                   )
                 }
               }
        ), 
        function(m) t(m)
      ) 
      
      ### binding nodes with their respective ancestral state
      list_matrix_edges_AS<- lapply(list_matrix_edges, function(x){
        if(is.na(x) == TRUE){
          NA
        } else{
          cbind(x, rep(NA, nrow(x)), rep(NA, nrow(x)))
        }
      }
      )
      for(i in 1:length(list_matrix_edges_AS)){
        if(is.null(dim(list_matrix_edges_AS[[i]])) == TRUE){
          list_matrix_edges_AS[[i]]<- NA
        } else{
          list_matrix_edges_AS[[i]][, c(3,4)]<- ancestral_comm_temp[[i]]
        }
      } 
      
      ### adding ancestor occurrence for branch lengths 
      for(i in 1:length(list_matrix_edges_AS)){
        if(is.null(dim(list_matrix_edges_AS[[i]])) == TRUE){
          list_matrix_edges_AS[[i]]<- NA
        } else{
          a<- strsplit(list_matrix_edges_AS[[i]][, 3], split= "")
          b<- strsplit(list_matrix_edges_AS[[i]][, 4], split= "")
          denom_tmp<- numeric()
          for(j in 1:length(a)){
            denom_tmp[j]<- length(which(!is.na(match(a[[j]], b[[j]])) == TRUE))
          }
          list_edge_brLen_insitu[[i]]<- cbind(list_edge_brLen_insitu[[i]], denom_tmp)
        }
      }
      ###
      
      ####Calculating Db-PE
      Db_PE <- lapply(
        lapply(list_edge_brLen_insitu, 
               function(x){
                 if(is.null(dim(x)) == TRUE){
                   NA
                 } else{
                   apply(x, MARGIN = 1, 
                         function(l){
                           if(is.na(l) == T){
                             NA
                           } else{
                             if(l[3] == 0){
                               l[1]/l[2]
                             } else{
                               ((l[1]/l[2])/l[3])
                             }
                           }
                         }
                   )
                 }
                 
               }
        ), 
        sum)
      Db_PEres<- matrix(unlist(Db_PE), nrow= length(list_matrix_edges_AS), ncol= 1, 
                        dimnames = list(names(list_matrix_edges_AS), "Db-PE"))
      ###
    }
    
    # Output for community metrics ---------------------------------------
    if(PD == TRUE & PE == TRUE){
      Db_comm_metrics<- data.frame(PE= PEt, Db_PE= Db_PEres, PD= PDt, PDlocal= PDlocal)
    } 
    if(PD == TRUE & PE == FALSE){
      Db_comm_metrics<- data.frame(PD= PDt, PDlocal= PDlocal)
    } 
    if(PE == TRUE & PD == FALSE){
      Db_comm_metrics<- data.frame(PE= PEt, Db_PE= Db_PEres)
    } 
    if(PE == FALSE & PD == FALSE){
      Db_comm_metrics <- NULL
    }
  }
  
  
  # Preparing the output ----------------------------------------------------
  
  
  div.opt <- c("jetz", "freck") %in% diversification
  
  if(div.opt[1]){
    mean.Jetz.Comm <- data.frame(JetzTotal = JetzTotalComm_harmonic, 
                                 JetzLocal = JetzLocalComm_harmonic,
                                 JetzLocalProp = JetzLocalComm_harmonic/
                                   JetzTotalComm_harmonic,
                                 row.names = row.names(W))
  }else{mean.Jetz.Comm <- NULL}
  
  if(div.opt[2]){
    mean.Freck.Comm <- data.frame(FreckTotal = FreckTotalComm_mean, 
                                  FreckLocal = FreckLocalComm_mean,
                                  FreckLocalProp = FreckLocalComm_mean/
                                    FreckTotalComm_mean,
                                  row.names = row.names(W))
  }else{mean.Freck.Comm <- NULL}
  
  if(!age.arrival){
    mean_age_arrival <- NULL
  }
  
  if(!dispersal.from){
    disp_from_frequency <- NULL
  }
  
  
  
  
  diversification <- list(mean.Jetz.Comm = mean.Jetz.Comm, 
                          mean.Freck.Comm = mean.Freck.Comm)
  
  age <- mean_age_arrival
  
  
  dispersal.freq <- disp_from_frequency
  
  
  return(list(diversification = diversification,
              age = age,
              dispersal.from = dispersal.freq,
              Db_comm_metrics = Db_comm_metrics))
  
  
}

