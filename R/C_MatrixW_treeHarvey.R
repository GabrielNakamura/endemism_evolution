
# loading tree from Harvey ------------------------------------------------
trfn = np(paste("T400F_AOS_HowardMoore.tre", sep=""))
moref(trfn)
tr = ape::read.tree(here::here("data", trfn))

# data with species codes
spp_codes <- read.csv(here::here("data", "Species_name_map_uids.csv")) # species codes for phylogenetic tree
names_tipHoward <- spp_codes[match(tr$tip.label, spp_codes$tipnamecodes), 
                             "aos.howardmoore.species"]
tr$tip.label <- names_tipHoward
names_tipHoward_matchW <- unlist(lapply(strsplit(tr$tip.label, " "), function(x) paste(x[1], x[2], sep = "_")))
tr$tip.label <- names_tipHoward_matchW

moref(trfn)
tr = ape::read.tree(trfn)

# Editing matrix W for harvey´s tree --------------------------------------

W_edit <- W[, - match(c("Zimmerius_improbus", "Phylloscartes_flaviventris", "Phelpsia_inornatus"), colnames(W)
)
]
W_edit_sub <- W_edit[, match(c("Suiriri_islerorum", "Anairetes_agraphia", "Anairetes_agilis"), colnames(W_edit))]
colnames(W_edit_sub) <- c("Suiriri_affinis", "Uromyias_agraphia", "Uromyias_agilis") 
W_edit <- W_edit[, - match(c("Suiriri_islerorum", "Anairetes_agraphia", "Anairetes_agilis"), colnames(W_edit))]
W_edit <- cbind(W_edit, W_edit_sub)
Onychorhynchus_coronatus <- ifelse(rowSums(W_edit[, match(c("Onychorhynchus_coronatus", "Onychorhynchus_swainsoni", "Onychorhynchus_mexicanus", "Onychorhynchus_occidentalis"), colnames(W_edit))]) >=1, 1, 0)
W_edit <- W_edit[, - match(c("Onychorhynchus_coronatus", "Onychorhynchus_swainsoni", "Onychorhynchus_mexicanus", "Onychorhynchus_occidentalis"), colnames(W_edit))]
W_edit <- cbind(W_edit, Onychorhynchus_coronatus)
Xolmis_rubetra <- ifelse(rowSums(W_edit[, match(c("Xolmis_salinarum", "Xolmis_rubetra"), colnames(W_edit))]) >= 1, 1, 0)
W_edit <- W_edit[, - match(c("Xolmis_salinarum", "Xolmis_rubetra"), colnames(W_edit))]
W_edit <- cbind(W_edit, Xolmis_rubetra) 
colnames(W_edit)[which(is.na(match(colnames(W_edit), tr$tip.label)) == TRUE)] # checking
W_edit # matriz correta para analise com a arvore de Harvey


# saving results ----------------------------------------------------------

ct <- treedata(tr, geofile)
ct.tr<-ct$phy
write.tree(ct.tr,file="Tree_TF400Howard_Pruned.tre")
trfn = np(paste("Tree_TF400Howard_Pruned.tre", sep=""))
moref(trfn)
tr = ape::read.tree(here::here("data", trfn))
write.tree(tr,file="Tree_TF400Howard_tip_corrected.tre")
trfn = np(paste("Tree_TF400Howard_tip_corrected.tre", sep=""))
write.table(W_edit, here::here("data", "processed", "W_harvey.txt"))
