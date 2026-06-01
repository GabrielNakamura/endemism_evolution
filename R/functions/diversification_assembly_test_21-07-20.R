samp<- W_toy
tree<- toy_treeEx


function (samp, tree, include.root = TRUE) 
{
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths, cannot compute pd")
  }
  if (include.root) {
    if (!is.rooted(tree)) {
      stop("Rooted tree required to calculate PD with include.root=TRUE argument")
    }
    tree <- node.age(tree)
  }
  species <- colnames(samp)
  SR <- rowSums(ifelse(samp > 0, 1, 0))
  nlocations <- dim(samp)[1]
  nspecies <- dim(samp)[2]
  PDs <- rep(NA, nlocations)
  for (i in 1:nlocations) {
    #i= 1
    present <- species[samp[i, ] > 0]
    treeabsent <- tree$tip.label[which(!(tree$tip.label %in% 
                                           present))]
    if (length(present) == 0) {
      PDs[i] <- 0
    }
    else if (length(present) == 1) {
      if (!is.rooted(tree) || !include.root) {
        warning("Rooted tree and include.root=TRUE argument required to calculate PD of single-species communities. Single species community assigned PD value of NA.")
        PDs[i] <- NA
      }
      else {
        PDs[i] <- tree$ages[which(tree$edge[, 2] == which(tree$tip.label == 
                                                            present))]
      }
    }
    else if (length(treeabsent) == 0) {
      PDs[i] <- sum(tree$edge.length)
    }
    else {
      sub.tree <- drop.tip(tree, treeabsent)
      if (include.root) {
        if (!is.rooted(tree)) {
          stop("Rooted tree required to calculate PD with include.root=TRUE argument")
        }
        sub.tree.depth <- max(node.age(sub.tree)$ages)
        orig.tree.depth <- max(tree$ages[which(tree$edge[, 
                                                         2] %in% which(tree$tip.label %in% present))])
        PDs[i] <- sum(sub.tree$edge.length) + (orig.tree.depth - 
                                                 sub.tree.depth)
      }
      else {
        PDs[i] <- sum(sub.tree$edge.length)
      }
    }
  }
  PDout <- data.frame(PD = PDs, SR = SR)
  rownames(PDout) <- rownames(samp)
  return(PDout)
}



####migration assembly function with diversification####
#####input arguments#####
# W = composition matrix, in lines the sample units and columns species
# tree = phylogenetic tree from class phylo;
# ancestral.area = data.frame. Lines with node names from tree and one column containing the ancestral area of each node;
# biogeo = data.frame. Lines with names of sample units, must be the same names that in W and one column containing the Ecoregion that each sample unit belongs.

######outputs#######
# Ecoregion.per.node = a matrix object. lines containig the nodes and the colums the species. Each cell present the Ecorregion that each ancestral of species i was found;
# Diversification.Assembly_Jetz = numeric. Biome Diversification values calculated according to a modifyied version of Jetz's metric of diversification;
# JetzTotalComm_mean = numeric. A vector containing mean values of total diversification calculated as harmonic mean of Jetz diversification for each row of matrix W;
# JetzLocalComm_mean = numeric. A vector containing mean values of local diversification calculated as the portion of harmonic mean of total Jetz diversification for each row of matrix W;
# JetzLocalSpp_mean = numeric. A vector containing mean values of local species diversification calculated as the portion total diverfication of each species;
# Diversification.Assembly_Freck =  numeric. Biome Diversification values calculated according to a modifyied version of Freckleton's metric of diversification;
# Diversification_total = matrix. Diversification values calculate for species according to Jetz and Freckleton metrics.

W<- W_toy
tree<- toy_treeEx
ancestral.area<- ancestral_area_toy
biogeo<- biogeo_toy
diversification.assembly<- function(W, tree, ancestral.area, biogeo){
  library(geiger)
  library(ape)
  library(ade4)
  library(phylobase)
  library(picante)
  library(sjmisc)
  s <- length(tree$tip.label)
  n <- tree$Nnode
  names_spComm<- colnames(W)
  PDt<- picante::pd(samp = W, tree = tree, include.root = T)$PD
  EDtotal<- evol.distinct(tree = tree, type = "equal.splits")
  Jetz_total<- 1/EDtotal$w
  Freck_total<- numeric(length = length(tree$tip.label))
  names(Freck_total)<- tree$tip.label
  type<- "equal.splits" #argument to compute equal splits metrics in ER calculation - internal function
  T_freck<- max(cophenetic(tree))/2
  for(l in 1:length(tree$tip.label)){
    nodes_Freck<- .get.nodes(tree, tree$tip.label[l]) 
    Freck_total[l]<- length(nodes_Freck)/T_freck
  }
  names(Jetz_total)<- tree$tip.label
  spxnode <- matrix(0, s, n)
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
  AS<-spxnode #Ancestral State matrix 
  for(i in 1:nrow(spxnode)){
    pres<-which(spxnode[i,]==1)
    AS[i,pres]<-as.character(ancestral.area[i,1]) #matriz - a informacao em cada celula representa o estado do ancestral de cada especie
  }
  for(i in 1:nrow(AS)){
    zero<-which(AS[i,]==0)
    AS[i,zero]<-NA #NAs indicam nos que nao contem a especie 
  }
  matrix_XJetz<- matrix(0, nrow = nrow(W), ncol = ncol(W), dimnames = list(rownames(W), colnames(W))) #matrix to receive the results of local Jetz metric
  matrix_XFreck<- matrix(0, nrow = nrow(W), ncol = ncol(W), dimnames = list(rownames(W), colnames(W))) #matrix to receive the results of local Freckleton metric
  age_arrival<-  matrix(0, nrow = nrow(W), ncol = ncol(W), dimnames = list(rownames(W), colnames(W))) #matrix to reveive the age in which the ancestor arrive
  node_local<- vector(mode = "list") 
  PD_local<- matrix(NA, nrow= nrow(W), ncol= 1, dimnames= list(rownames(W), "PD_local"))
  #i= 1
  for(i in 1:nrow(W[1:2,])){
    pres<- which(W[i,]==1)
    pres<- names_spComm[pres]
    nodes_species<- vector(mode= "list")
    tip_edgeLength<- numeric()
    #j= 1
    for(j in 1:length(pres)){
      nodes_sp<- AS[,pres[j]][!is.na(AS[,pres[j]])]
      nodes_sp<- nodes_sp[length(nodes_sp):1] #nodes for species j in community i
      if(str_contains(nodes_sp[1],biogeo[i,1]) != TRUE)
      { #check if the most recent ancestor was in the same bioregion that the observed species
        matrix_XJetz[i,pres[j]]<- 0.00001 
        matrix_XFreck[i,pres[j]]<- 0.00001
        age_arrival[i,pres[j]]<- 0.00001
      } else{
        nodes_Freck<- .get.nodes(tree, pres[j]) #get the internal nodes from tip to root for species j
        nodes_Freck_internal<- nodes_Freck[1:(length(nodes_Freck))] #organize the internal nodes for species j
        Div_Freck<- length(nodes_Freck)/T_freck #modified equation to calculate local diversity based in Equation 4 Freckleton et al paper  Am.Nat 2008 
        nodes_all<- numeric(length = length(nodes_sp)) #test if all ancestors of species j are in the same ecoregion that local i
        for(m in 1:length(nodes_sp)){
          nodes_all[m]<- str_contains(nodes_sp[m], biogeo[i,1])
        }
        if(all(nodes_all==1)){ #if all ancestors of species j are in the same ecoregion of local 1 this will be TRUE
          x<- names(nodes_sp[length(nodes_sp)]) # if TRUE, take the most ancient ancestor (node 8 in this example) as the reference node for calculation of local diversification
        }else
        {
          if(length(nodes_sp)==1){
            x<- names(nodes_sp) # for the case in which the specie do not present internal branches
          }
          #k=1
          for(k in 1:length(nodes_sp)){ #find for the most ancient ancestor of species j that was present in the local i
            if(str_contains(nodes_sp[k],biogeo[i,1])!=TRUE){
              x<- names(nodes_sp)[k-1] # take the name of the node (ancestor) of species j that was present in the local i
              break
            } 
          } #extract the most ancient ancestral that was presented at Ecoregion of species j in local i
        }
        nodes_div<- sort(c(nodepath(tree, from = as.numeric(strsplit(x, split = "N")[[1]][2]), 
                                    to = which(tree$tip.label==pres[j])))[-length(nodepath(tree, from = as.numeric(strsplit(x, split = "N")[[1]][2]), to = which(tree$tip.label==pres[j])))],decreasing = TRUE) #organize the path of the most ancient ancestral that was presented at ecoregion of local i 
         
        nodes_species[[j]]<- c(nodepath(tree, from = as.numeric(strsplit(x, split = "N")[[1]][2]), 
                                        to = which(tree$tip.label==pres[j]))) #nodes for species j, this was done to compute PD local
        if(length(nodes_div) == 1){
          internal.brlen_div <- tree$edge.length[which(tree$edge[, 
                                                                 2] %in% nodes_div)] #branch lenghts (in times) for internal branch lengts of the most ancient ancestral 
        } else{
          nodes_div <- nodes_div[1:(length(nodes_div) - 1)] #internal nodes that form the path from most ancient ancestral that was presented at Ecoregion of local i to species j
          internal.brlen_div <- tree$edge.length[which(tree$edge[, 
                                                                 2] %in% unique(nodes_div))] #ages for internal nodes of nodes_div object
        }
        Div_Freck_local<- (length(nodes_div)+1)/T_freck #modifyed equation 4 from Freckleton et al (2008) considering only the nodes that diversified in ecoregion of local i, plus one is only a correction of the previous step
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
        tip_edgeLength[j]<- tree$edge.length[which.edge(tree, pres[j])]
        
        ED_div <- sum(internal.brlen_div, tree$edge.length[which.edge(tree, 
                                                                      pres[j])]) #Local ED - modifyed ED considering only the edges of ancestros inside the biogeo of local i for species j
        EDtotal_spp<- EDtotal$w[which(EDtotal$Species==pres[j])]
        JetzTot<- 1/EDtotal$w[which(EDtotal$Species==pres[j])] #Diversificatio calculated according to Jetz for species j
        JetzLocal<- (JetzTot*(ED_div/EDtotal_spp)) #Modified Jetz to calculate only for local diversification
        matrix_XJetz[i, pres[j]]<- JetzLocal #Jetz local diversification
        matrix_XFreck[i, pres[j]]<- Div_Freck_local #Freckleton local diversification
        age_arrival[i,pres[j]]<- ages[x,] #recebe a idade de chegada
      }
    }
    #age_arrival<- ifelse(age_arrival==1,NA,age_arrival)
    unique_internal_blength<- tree$edge.length[which(tree$edge[, 
                                                               2] %in% unique(unlist(nodes_species)))]
    
    
    nodes_species_noNull <- nodes_species[which(lapply(nodes_species, is.null) == FALSE)] #species with local PD
    spp_PDlocal_noNull <- vector(mode = "list", length = length(nodes_species_noNull)) #object to receive all values of branch lenghts of all species of community j
    recent_PD <- which(lapply(nodes_species, is.null) == TRUE) # species whith no local PD
    
    for(l in 1:length(nodes_species_noNull)){
      spp_PDlocal <- numeric(length = (length(unlist(nodes_species_noNull[l])) - 1))
      for(k in 1:(length(unlist(nodes_species_noNull[l])) - 1)){
        spp_PDlocal[k]<- tree$edge.length[which((tree$edge[, 1] == c(unlist(nodes_species_noNull[l])[k]) & 
                                                   tree$edge[, 2] == unlist(nodes_species_noNull[l])[k + 1]) == TRUE)]
      }
      spp_PDlocal_noNull[[l]] <- spp_PDlocal #PD for each species
    } #calculation of branch lengths for species with in situ diversification
    PD_local[i, ] <- sum(unlist(spp_PDlocal_noNull)) #local PD for local i
  }
  
  #calculating harmonic mean for Jetz
  sum_localDiv<- apply(matrix_XJetz, 1, function(x) sum(x[which(x!=0 & x!=0.00001)])) #sum of communities local diversification
  sum_localDivSpp<- apply(matrix_XJetz, 2, function(x) sum(x[which(x!=0 & x!=0.00001)])) #sum of species local diversification
  matrix_totalDiv_Jetz<- W #objet to receive Jetz diversification values
  for(i in colnames(W)){ #substitui a matrix de ocorrencia pelos valores de div total do Jetz
    #i=1
    matrix_totalDiv_Jetz[,i]<- ifelse(matrix_totalDiv_Jetz[,i]>=1,Jetz_total[i],0) #input Jetz total diversification values in occurence matrix
  }
  sum_totalDiv<- apply(matrix_totalDiv_Jetz, 1, function(x) sum(x[which(x!=0)])) #sum of communities total diversification
  sum_totalDivSpp<- apply(matrix_totalDiv_Jetz, 2, function(x) sum(x[which(x!=0)])) #sum of species total diversification
  prop_divLocal<- sum_localDiv/sum_totalDiv
  numerator_values<- apply(matrix_totalDiv_Jetz, 1, function(x){
    sum(1/x[which(x!=0)])
  })
  denom_values<- apply(matrix_totalDiv_Jetz, 1, function(x){
    length(which(x!=0))
  })
  JetzTotalComm_harmonic<- denom_values/numerator_values #harmonic mean for Jetz total diversification
  JetzLocalComm_harmonic<- (JetzTotalComm_harmonic*(sum_localDiv/sum_totalDiv)) #harmonic mean for Jetz local diversification
  JetzLocalSpp_harmonic<- (Jetz_total*(sum_localDivSpp/sum_totalDivSpp))
  Diversification_table<- data.frame(Jetz_total= Jetz_total, Freckleton_total= Freck_total) # original Jetz and Freckleton measures of diversification
  return(list(Ecoregion.per.node=AS,Diversification.Assembly_Jetz=matrix_XJetz, JetzTotalComm_mean= JetzTotalComm_harmonic, JetzLocalComm_mean= JetzLocalComm_harmonic, 
              JetzLocalSpp_mean= JetzLocalSpp_harmonic, Diversification.Assembly_Freck= matrix_XFreck, Diversification_total= Diversification_table, 
              PD_local= PD_local, PD_total= PDt, age_arrival= age_arrival))
}


#####testing migration assembly function######
test_migration<- FALSE
if(test_migration==TRUE){
  ####objetos necessarios######
  W<- matrix(c(0,0,0,0,1,
               0,0,0,0,1,
               0,0,1,0,0,
               0,0,0,1,0,
               0,0,1,1,0,
               1,1,1,1,0,
               1,1,0,1,0), nrow=5, ncol= 7, dimnames=list(c("A", "B", "C", "D", "E"),
                                                          paste("s", 1:7, sep="")))
  tree<- geiger::sim.bdtree(b = 1,d = 0,stop = "taxa",n = 7,extinct = FALSE)
  #tree<- makeNodeLabel(tree, method= "user", nodeList= list(N8AB=c("s"), N9AB= c("s7","s6","s5","s4","s3","s2"),
  #                                                          N10B= c("s6","s5","s4","s3","s2"),N11B=c("s4","s3","s2"),
  #                                                          N12A=c("s3","s4"), N13A=c("s6","s5")))
  ancestral.area<- data.frame(state= c("AB","AB", "B","B","A","A")) #example data frame derived from a reconstruction using BioGeoBEARS
  biogeo<- data.frame(Ecoregion= c("A", "B", "B", "A", "A")) # data frame containing the bioregions of each sample unit
  teste_migration<- diversification.assembly(W = W, tree = tree, ancestral.area = ancestral.area, biogeo = biogeo)
}
