####Supplementary material - Toy example########

toy_treeEx<- ape::read.tree(here::here("data", "toy_tree.new"))
toy_treeEx<- ape::makeNodeLabel(toy_treeEx, method= "user", nodeList= list(N6ABC= c("s"), 
                                                                           N7ABC= c("s4", "s5"),
                                                                           N8B= c("s1", "s2", "s3"),
                                                                           N9C= c("s1", "s2")
)
)
quartz()
plot(toy_treeEx,  show.node.label = T)
ape::axisPhylo()

### current occurrence area ####
W_toy<- matrix(c(0, 1, 1,
             0, 1, 1,
             0, 1, 0,
             1, 0, 0,
             1, 0, 0
             ), nrow= 3, ncol= 5, dimnames=list(c("Comm 1", "Comm 2", "Comm 3"),
                                                        c(paste("s", 1:5, sep=""))
                                                )
           )

biogeo_toy<- data.frame(Ecoregion= c("A", "B", "C"))
ancestral_area_toy<- data.frame(state= c("ABC", "B", "C", "ABC"))

####calculating age arrival with toy example#####
age_arrival_toy<- diversification.assembly(W = W_toy, tree = toy_treeEx, ancestral.area = ancestral_area_toy, biogeo = biogeo_toy)$age_arrival
apply(age_arrival_toy, MARGIN = 1, function(x) mean(x[which(x != 0)]))
abs(node.depth.edgelength(phy = toy_treeEx)
    -max(node.depth.edgelength(toy_treeEx)))[-c(1:length(toy_treeEx$tip.label))]
toy_treeEx$node.label
which(age_arrival_toy[1, ] != 0)
