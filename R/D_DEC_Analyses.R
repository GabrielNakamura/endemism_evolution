###########################################################
###########################################################
###########################################################
# DEC ANALYSES
###########################################################
###########################################################
###########################################################
# # Install BioGeoBEARS from CRAN 0-cloud:
install.packages("BioGeoBEARS", dependencies=TRUE, repos="http://cran.rstudio.com")
require(BioGeoBEARS)
library(optimx)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(roxygen2)
library(BioGeoBEARS)
source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R") # (needed now that traits model added; source FIRST!)
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R") # added traits model
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
# slight speedup hopefully
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)



#######################################################
# SETUP: YOUR TREE FILE AND GEOGRAPHY FILE
#######################################################
# Phylogeny 
trfn = np(paste("T400F_AOS_HowardMoore.tre", sep=""))
moref(trfn)
tr = ape::read.tree(here::here("data", trfn))
spp_codes <- read.csv(here::here("data", "Species_name_map_uids.csv")) # species codes for phylogenetic tree
names_tipHoward <- spp_codes[match(tr$tip.label, spp_codes$tipnamecodes), 
                             "aos.howardmoore.species"]
tr$tip.label <- names_tipHoward
names_tipHoward_matchW <- unlist(lapply(strsplit(tr$tip.label, " "), function(x) paste(x[1], x[2], sep = "_")))
tr$tip.label <- names_tipHoward_matchW

# editing matrix W 
W_edit <- W[, - match(c("Zimmerius_improbus", "Phylloscartes_flaviventris", "Phelpsia_inornatus"), colnames(W)
                    )
            ]
W_edit_sub <- W_edit[, match(c("Suiriri_islerorum", "Anairetes_agraphia", "Anairetes_agilis"), colnames(W_edit))]
colnames(W_edit_sub) <- c("Suiriri_affinis", "Uromyias_agraphia", "Uromyias_agilis") 
W_edit <- W_edit[, - match(c("Suiriri_islerorum", "Anairetes_agraphia", "Anairetes_agilis"), colnames(W_edit))]
W_edit <- cbind(W_edit, W_edit_sub)
Onychorhynchus_coronatus <- rowSums(W_edit[, match(c("Onychorhynchus_coronatus", "Onychorhynchus_swainsoni", "Onychorhynchus_mexicanus", "Onychorhynchus_occidentalis"), colnames(W_edit))])
W_edit <- W_edit[, - match(c("Onychorhynchus_coronatus", "Onychorhynchus_swainsoni", "Onychorhynchus_mexicanus", "Onychorhynchus_occidentalis"), colnames(W_edit))]
W_edit <- cbind(W_edit, Onychorhynchus_coronatus)
Xolmis_rubetra <- rowSums(W_edit[, match(c("Xolmis_salinarum", "Xolmis_rubetra"), colnames(W_edit))])
W_edit <- W_edit[, - match(c("Xolmis_salinarum", "Xolmis_rubetra"), colnames(W_edit))]
W_edit <- cbind(W_edit, Xolmis_rubetra) 
colnames(W_edit)[which(is.na(match(colnames(W_edit), tr$tip.label)) == TRUE)] # checking
W_edit # matriz correta para analise com a arvore de Harvey

# Geography file - PHYLLIP FORMAT
colnames(presab.sigmodontinae.final)
presab.sigmodontinae.final[,1]
bior_values_factor

# extract bioregions values for presences of each species
t.2<-bior_values_factor[which(presab.sigmodontinae.final[,1]==1)]
t.2
unique(t.2)

MATRIXGEO<-matrix(NA,nrow=length(colnames(presab.sigmodontinae.final)),ncol=length(levels(bior_values_factor)),byrow=FALSE)
for(i in 1:length(presab.sigmodontinae.final[1,])){
  t.2<-bior_values_factor[which(presab.sigmodontinae.final[,i]==1)]
  MATRIXGEO[i,1]<-unique(t.2)[1]
  MATRIXGEO[i,2]<-unique(t.2)[2]
  MATRIXGEO[i,3]<-unique(t.2)[3]
  MATRIXGEO[i,4]<-unique(t.2)[4]
  MATRIXGEO[i,5]<-unique(t.2)[5]
  MATRIXGEO[i,6]<-unique(t.2)[6]
  MATRIXGEO[i,7]<-unique(t.2)[7]
  MATRIXGEO[i,8]<-unique(t.2)[8]
  MATRIXGEO[i,9]<-unique(t.2)[9]
  MATRIXGEO[i,10]<-unique(t.2)[10]
  MATRIXGEO[i,11]<-unique(t.2)[11]
}
rownames(MATRIXGEO)=colnames(presab.sigmodontinae.final)
fix(MATRIXGEO) # Os valores são a posição dentro da sequência de levels(bioregions): 1 2 4 8 9 10 11 12 13 14 15 
#levels(bioregions): 1 2 4 8 9 10 11 12 13 14 15
#sequence of file:   1 2 3 4 5 6  7  8  9  10 11
write.csv(MATRIXGEO,"MATRIXGEO.csv") # create phyllip matrix in excel
 

# Geography file - PHYLLIP FORMAT
geofile<-read_PHYLIP_data("MATRIXGEO.txt")
geogfn=geofile

# the location # the path to the file is the object to use, the the file itself
geogfn = np(paste("MATRIXGEO-3areas.txt", sep=""))
moref(geogfn)

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
max_range_size = 3
#max_tipsize=3
# To check the number of states for a given number of ranges, try:
numstates_from_numareas(numareas=11, maxareas=3, include_null_range=TRUE)


#######################################################
#######################################################
# DEC AND DEC+J ANALYSIS
#######################################################
#######################################################
#######################################################
# Run DEC
#######################################################
# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
# 1. Here, un-comment ONLY the files you want to use.
# 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
# 3. For example files see (a) extdata_dir, 
#  or (b) http://phylo.wikidot.com/biogeobears#files
#  and BioGeoBEARS Google Group posts for further hints)
#
# Uncomment files you wish to use in time-stratified analyses:
BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 4
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel" 
# is now default on Macs in R 3.0+, but apparently still 
# has to be typed on some Windows machines. Note: apparently 
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table


# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC model
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object


# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object,allow_huge_ranges = FALSE)

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
memory.limit(1000000)
memory.size(1000000)

runslow = TRUE
resfn = "Psychotria_DEC_M0_unconstrained_v1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}

str(tr) # 259 nodes
names(resDEC)
resDEC[4]

#######################################################
# PDF plots
#######################################################
pdffn = "Psychotria_DEC_vs_DEC+J_M0_unconstrained_v2.pdf"
pdf(pdffn, width=6, height=6)

#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on Psychotria M0_unconstrained"

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
x11()
plot(tr)
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.5, tipcex=0.5, statecex=0.5, splitcex=0.5, titlecex=0.5, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges,plotlegend=TRUE)
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.5, tipcex=0.5, statecex=0.5, splitcex=0.5, titlecex=0.5, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges,plotlegend=TRUE)

plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.35, tipcex=0.3, statecex=0.5, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges,plotlegend=TRUE)
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.3, tipcex=0.3, statecex=0.3, splitcex=0.3, titlecex=0.3, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges,plotlegend=TRUE)


################## OUTPUT ##################
names(res)
res[9] 
str(res[4]) #[1:519, 1:232]
res[12]
str(res[12]$outputs)

tr # 260 tips, 259 nodes


