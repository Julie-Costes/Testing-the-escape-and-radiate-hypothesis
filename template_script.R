#### FULL SCRIPT 


#### STEP 1 - UPDATE THE DATA ####


#### A - Update the insect tree ####

rm(list=ls())

library(phytools)
library(dplyr)

setwd("~/Documents/STAGE/Script")


## Import the most complete tree:

tree_nymphalidae <- read.nexus("~/Documents/STAGE/Data/nymphalini_data/intermediate_data/Chazot 2021 - complete tree.nex")


## Extract the clade:

taxa <- c("NW157_3_Symbrenthia_sinoides", "EW22_8_Polygonia_gracilis")
node <- getMRCA(tree_nymphalidae, taxa)
tree_B <- extract.clade(tree_nymphalidae, node)


## Modify the tree:

# Remove a duplicate species
new_tree_B <- drop.tip(tree_B, tip = "NW78_6_Nymphalis_vaualbum")

# Add a species
where <- which(new_tree_B$tip.label == "NW157_3_Symbrenthia_sinoides") # identify the insertion taxon
new_tree_B <- bind.tip(new_tree_B, tip.label = "__Symbrenthia_niphanda", where, edge.length = 3.2, position = 3.2)
new_tree_B = force.ultrametric(new_tree_B)


## Match species names:

# Display all species names
new_tree_B$tip.label 

# Targeted corrections 
new_tree_B$tip.label = gsub("jlr4_Vanessa_tameamea","jlr4__Vanessa_tameamea",new_tree_B$tip.label)
new_tree_B$tip.label = gsub("NW65_8_Polygonia_c_aureum","NW65_8_Polygonia_caureum",new_tree_B$tip.label)
new_tree_B$tip.label = gsub("NW70_3_Polygonia_c_album","NW70_3_Polygonia_calbum",new_tree_B$tip.label)
new_tree_B$tip.label = gsub("NW165_2_Polygonia_g_argenteum","NW165_2_Polygonia_gargenteum",new_tree_B$tip.label)
new_tree_B$tip.label = gsub("ZF_LY_001739_Polygonia_gongga","ZFLY_001739_Polygonia_gongga",new_tree_B$tip.label)
new_tree_B$tip.label = gsub("ZF_LY_001174_Mynes_talboti","ZFLY_001174_Mynes_talboti",new_tree_B$tip.label)

# General correction
new_tree_B$tip.label = paste0(sapply(strsplit(split = "_", new_tree_B$tip.label), "[[", 3), "_", sapply(strsplit(split = "_", new_tree_B$tip.label), "[[", 4)) # ne garder que Genre_espece


## Save the new tree:

# write.tree(new_tree_B, file = "~/M2/STAGE/Data/nymphalini_data/tree_chazot_nymphalini.tre")
test_tree <- read.tree("~/Documents/STAGE/Data/nymphalini_data/tree_chazot_nymphalini.tre")

## Max age:
max(node.depth.edgelength(new_tree_B))
# 28.01527



#### B - Update the network ####
#### to ensure species match

## Interactions from Braga (2020)

# Import the table:
network <- read.table("~/Documents/STAGE/Data/nymphalini_data/intermediate_data/network_Nymphaliny_plants.csv", sep=";", header=TRUE)
network[network == 1] <- 0
network[network == 2] <- 1
# network[network>0] <- 1

# Check if all the names in the table are present in the tree:
all(rownames(network) %in% new_tree_B$tip.label)

# Correct the species names:
print(rownames(network))
rownames(network) = gsub("Polygonia_c_aureum","Polygonia_caureum", rownames(network))
rownames(network) = gsub("Polygonia_c_album","Polygonia_calbum", rownames(network))
rownames(network) = gsub("Nymphalis_l_album","Nymphalis_vaualbum", rownames(network))

# Save the new network:
# write.table(network, file = "~/M2/STAGE/Data/nymphalini_data/network_nymphalini_plants.csv", sep = ";")

# Test if the table opens correctly:
test <- read.table("~/Documents/STAGE/Data/nymphalini_data/network_nymphalini_plants.csv", sep=";", header=TRUE)  
all(rownames(test) %in% new_tree_B$tip.label)


## All interactions 

# Import the table:
library(readxl)
network_all <- read_excel("~/Documents/STAGE/Data/nymphalini_data/intermediate_data/network_nymphalini_all.xlsx")
network_all <- as.data.frame(network_all)
row.names(network_all) <- network_all$row_names     # set the first column as row names
network_all <- network_all[, -1]                    # remove the now redundant first column
network_all[network_all != 0] <- 1

# Check if all the names in the table are present in the tree:
all(rownames(network_all) %in% new_tree_B$tip.label)


## Add interactions from HOSTS database:

full_dataset <- read.csv("~/Documents/STAGE/Data/1 - Research/HOSTS/resource.csv", header = TRUE, sep = ",", na.strings = c("", "NA"))
full_dataset$Insect.Genus_Species <- paste(full_dataset$Insect.Genus, full_dataset$Insect.Species, sep = "_")
full_dataset$Insect.Genus_Species[is.na(full_dataset$Insect.Species) | full_dataset$Insect.Species == "spp."] <- NA # to avoid names like "Colias_NA" and "Colias spp."

# Filter to keep only the Nymphalini
nymphalini_dataset <- subset(full_dataset, Insect.Family == "Nymphalidae")
genera_list <- unique(sapply(strsplit(new_tree_B$tip.label, "_"), `[`, 1))
nymphalini_dataset <- nymphalini_dataset[nymphalini_dataset$Insect.Genus %in% genera_list, ]

# Get unique insect species names and plant family names
insect_species <- unique(na.omit(nymphalini_dataset$Insect.Genus_Species))
plant_families <- unique(na.omit(nymphalini_dataset$Hostplant.Family))

# Create a zero matrix with insect species as rows and plant families as columns
network_HOSTS <- matrix(0, 
                             nrow = length(insect_species), 
                             ncol = length(plant_families), 
                             dimnames = list(insect_species, plant_families))

# Loop through the dataset to update the matrix based on interactions
for (i in 1:nrow(nymphalini_dataset)) {
  # Check if the species and family are not NA
  insect <- nymphalini_dataset$Insect.Genus_Species[i]
  hostplant <- nymphalini_dataset$Hostplant.Family[i]
  
  # If both insect species and plant family are present, set the corresponding matrix element to 1
  if (!is.na(insect) & !is.na(hostplant)) {
    network_HOSTS[insect, hostplant] <- 1
  }
}
network_HOSTS <- as.data.frame(network_HOSTS)

# Remove the rows and columns where there is no interaction
network_HOSTS <- network_HOSTS[,which(colSums(network_HOSTS)>0)]
network_HOSTS <- network_HOSTS[which(rowSums(network_HOSTS)>0),]

# Comparison 
setdiff(colnames(network_all), colnames(network_HOSTS))
setdiff(colnames(network_HOSTS), colnames(network_all))
setdiff(rownames(network_all), rownames(network_HOSTS))
setdiff(rownames(network_HOSTS), rownames(network_all))

# Rename plant families
colnames(network_HOSTS)[colnames(network_HOSTS) == "Leguminosae (P)"] <- "Fabaceae"
colnames(network_HOSTS)[colnames(network_HOSTS) == "Compositae"] <- "Asteraceae"

# Rename insect species
rownames(network_HOSTS)[rownames(network_HOSTS) == "Polygonia_c-album"] <- "Polygonia_calbum"
rownames(network_HOSTS)[rownames(network_HOSTS) == "Polygonia_c-aureum"] <- "Polygonia_caureum"
rownames(network_HOSTS)[rownames(network_HOSTS) == "Nymphalis_io"] <- "Aglais_io"
rownames(network_HOSTS)[rownames(network_HOSTS) == "Antanartia_abyssinica"] <- "Vanessa_abyssinica"
rownames(network_HOSTS)[rownames(network_HOSTS) == "Antanartia_dimorphica"] <- "Vanessa_dimorphica"
rownames(network_HOSTS)[rownames(network_HOSTS) == "Antanartia_hippomene"] <- "Vanessa_hippomene"

# Add species present in the tree
setdiff(rownames(network_HOSTS), new_tree_B$tip.label)
setdiff(new_tree_B$tip.label, rownames(network_HOSTS))
rownames(network_HOSTS)[rownames(network_HOSTS) == "Aglais_kaschmirensis"] <- "Aglais_caschmirensis"
rownames(network_HOSTS)[rownames(network_HOSTS) == "Hypanartia_kefersteinii"] <- "Hypanartia_kefersteini"

new_rows <- c(intersect(setdiff(rownames(network_HOSTS), rownames(network_all)), new_tree_B$tip.label))
zeros <- as.data.frame(matrix(0,
                              nrow = length(new_rows),
                              ncol = ncol(network_all),
                              dimnames = list(new_rows, colnames(network_all))))
network_all <- rbind(network_all, zeros)

# Subset comparison
common_cols <- intersect(colnames(network_all), colnames(network_HOSTS))
common_rows <- intersect(rownames(network_all), rownames(network_HOSTS))
modif_count <- 0
for (row in common_rows) {
  for (col in common_cols) {
    if (network_HOSTS[row, col] == 1 && network_all[row, col] != 1) {
      network_all[row, col] <- 1
      print(paste("Modif", row, col))
      modif_count <- modif_count + 1
    }
  }
}
print(paste("Modif count:", modif_count))


## Add interactions from LepTraits database:

LepTraits_full <- read.csv("~/Documents/STAGE/Data/1 - Research/LepTraits/consensus/consensus.csv", header = TRUE, sep = ",", na.strings = c("", "NA"))
LepTraits_full$Species <- gsub(" ", "_", LepTraits_full$Species)
LepTraits_full <- subset(LepTraits_full, select = c(Family, Genus, Species, NumberOfHostplantFamilies, SoleHostplantFamily, PrimaryHostplantFamily, SecondaryHostplantFamily, EqualHostplantFamily, NumberOfHostplantAccounts))

# Filter to keep only the Nymphalini
nymphalini_dataset <- subset(LepTraits_full, Genus %in% genera_list)

# Get unique insect species names and plant family names
insect_species <- unique(na.omit(nymphalini_dataset$Species))
plant_families <- sort(unique(na.omit(c(nymphalini_dataset$SoleHostplantFamily, nymphalini_dataset$PrimaryHostplantFamily, nymphalini_dataset$SecondaryHostplantFamily))))

# Create a zero matrix with insect species as rows and plant families as columns
network_LepTraits <- matrix(0, 
                             nrow = length(insect_species), 
                             ncol = length(plant_families), 
                             dimnames = list(insect_species, plant_families))

# Create a list of plant families by species in nymphalini_dataset
hostplant_families <- nymphalini_dataset %>%
  dplyr::select(Species, SoleHostplantFamily, PrimaryHostplantFamily, 
                SecondaryHostplantFamily) %>%
  tidyr::pivot_longer(cols = -Species, values_to = "Family") %>%
  dplyr::filter(!is.na(Family)) %>%
  dplyr::distinct(Species, Family)

# Update network_LepTraits by replacing 0 with 1 if the family is a host for the species
for (i in 1:nrow(hostplant_families)) {
  species <- hostplant_families$Species[i]
  family <- hostplant_families$Family[i]
  
  if (species %in% rownames(network_LepTraits) && family %in% colnames(network_LepTraits)) {
    network_LepTraits[species, family] <- 1
  }
}
network_LepTraits <- as.data.frame(network_LepTraits)

# Remove the rows and columns where there is no interaction
network_LepTraits <- network_LepTraits[,which(colSums(network_LepTraits)>0)]
network_LepTraits <- network_LepTraits[which(rowSums(network_LepTraits)>0),]

# Comparison 
setdiff(colnames(network_all), colnames(network_LepTraits))
setdiff(colnames(network_LepTraits), colnames(network_all))
setdiff(rownames(network_all), rownames(network_LepTraits))
setdiff(rownames(network_LepTraits), rownames(network_all))

# Rename insect species
rownames(network_LepTraits)[rownames(network_LepTraits) == "Polygonia_c-aureum"] <- "Polygonia_caureum"
rownames(network_LepTraits)[rownames(network_LepTraits) == "Nymphalis_l-album"] <- "Nymphalis_vaualbum"
rownames(network_LepTraits)[rownames(network_LepTraits) == "Symbrenthia_geoffroyii"] <- "Mynes_geoffroyi"

# Add species present in the tree
setdiff(rownames(network_LepTraits), new_tree_B$tip.label)
setdiff(new_tree_B$tip.label, rownames(network_LepTraits))
rownames(network_LepTraits)[rownames(network_LepTraits) == "Symbrenthia_websteri"] <- "Mynes_websteri"
rownames(network_LepTraits)[rownames(network_LepTraits) == "Symbrenthia_woodfordi"] <- "Mynes_woodfordi"

new_rows <- c(intersect(setdiff(rownames(network_LepTraits), rownames(network_all)), new_tree_B$tip.label))
zeros <- as.data.frame(matrix(0,
                              nrow = length(new_rows),
                              ncol = ncol(network_all),
                              dimnames = list(new_rows, colnames(network_all))))
network_all <- rbind(network_all, zeros)

# Subset comparison
common_cols <- intersect(colnames(network_all), colnames(network_LepTraits))
common_rows <- intersect(rownames(network_all), rownames(network_LepTraits))
modif_count <- 0
for (row in common_rows) {
  for (col in common_cols) {
    if (network_LepTraits[row, col] == 1 && network_all[row, col] != 1) {
      network_all[row, col] <- 1
      print(paste("Modif", row, col))
      modif_count <- modif_count + 1
    }
  }
}
print(paste("Modif count:", modif_count))


# Save the new network:
# write.table(network_all, file = "~/Documents/STAGE/Data/nymphalini_data/network_nymphalini_all.csv", sep = ";")

# Test if the table opens correctly:
test <- read.table("~/Documents/STAGE/Data/nymphalini_data/network_nymphalini_all.csv", sep=";", header=TRUE)  
all(rownames(test) %in% new_tree_B$tip.label)



#### C - Update the plant tree ####

# Import the angiosperm tree (Zuntini 2024):
tree_angiosperm <- read.tree("~/Documents/STAGE/Data/nymphalini_data/intermediate_data/3_calibrated_family_level_young_tree.tre")

# Check for duplicated names:
anyDuplicated(tree_angiosperm$tip.label)

# Check that the tree is binary:
is.binary(tree_angiosperm)

# Check that the tree is ultrametric:
is.ultrametric(tree_angiosperm)
tree_angiosperm <- force.ultrametric(tree_angiosperm) 

# Keep only the plant families that are in the network:
tips_to_keep <- c(colnames(network_all))
new_tree_angiosperm <- drop.tip(tree_angiosperm, setdiff(tree_angiosperm$tip.label, tips_to_keep))
plot(new_tree_angiosperm)
all(colnames(network_all) %in% new_tree_angiosperm$tip.label)

# Save the new tree:
# write.tree(new_tree_angiosperm, file = "~/Documents/STAGE/Data/nymphalini_data/tree_angiosperm_families.tre")


# scp -r /mnt/c/Users/julie/Documents/M2/STAGE/Data/nymphalini_data/tree* jcostes@genobioinfo.toulouse.inrae.fr:/home/jcostes/work/nymphalini/nymphalini_data/
# scp -r /mnt/c/Users/julie/Documents/M2/STAGE/Data/nymphalini_data/network_nymphalini_plants.csv jcostes@genobioinfo.toulouse.inrae.fr:/home/jcostes/work/nymphalini/nymphalini_data/
# scp -r /Users/bperez/Documents/STAGE/Data/nymphalini_data/network_nymphalini_all.csv jcostes@genobioinfo.toulouse.inrae.fr:/home/jcostes/work/nymphalini/nymphalini_data/



#### STEP 2 - BDD ####

#### julia


# need to install them once only (using SLURM)
# import Pkg; Pkg.add("ArgParse")
# import Pkg; Pkg.add("Tapestree")
# import Pkg; Pkg.add("Plots")


using DelimitedFiles, Distributed, ArgParse, Base
using Tapestree
using Plots


cd("/home/jcostes/work/nymphalini/BDD/")

names = ["nymphalini_test_1_sf_0.83", "nymphalini_test_2_sf_0.76", "nymphalini_test_3_sf_0.70"]

sampling_fractions = [0.83, 0.76, 0.70]


output_dir = "/home/jcostes/work/nymphalini/BDD/"


#### Open tree file  --> has to be a newick file
tree = read_newick("../nymphalini_data/tree_chazot_nymphalini.tre")



#### Run BDD with constant extinction

for (sampling_fraction, name) in zip(sampling_fractions, names)

path_results = string(output_dir, "/", name)

r, tv = insane_gbmce(tree,
                     nburn    = 500_000,
                     niter    = 6_250_000,
                     nthin    = 25_000,
                     nflush   = 25_000,
                     ofile    = path_results,
                     α_prior  = (0.0, 10.0),       # Normal distribution on drift (trend)
                     σλ_prior = (0.1, 0.01),       # Inverse gamma on diffusion (stochasticity)
                     μ_prior  = (1.0, 1.0),        # Gamma prior on extinction
                     δt       = 1e-3,
                     survival = true,              # condition on survival
                     mxthf    = 0.2,               # this is for better mixing, keep like this
                     tρ       = Dict("" => sampling_fraction)) # fixed sampling fraction



### save augmented trees as newick

treesDA = iread(string(path_results, ".txt"))

label_trees = [sT_label(ti, tree) for ti in treesDA]

write_newick(label_trees, path_results)



### estimating posterior average rates along the tree

tv0 = remove_unsampled(treesDA)

# we can then estimate the average tree
treeDA = imean(tv0)

write_nexus(treeDA, tree, string(path_results, "_mean_rates"))



### plot


ENV["TMPDIR"] = "/tmp"


gr()
plot(ltt(treeDA), linewidth = 2.0)
file_name = string("plot_ltt_", name, ".pdf")
savefig(joinpath(output_dir, file_name))

gr()
plot_julia = plot(treeDA, b)
file_name = string("plot_tree_speciation_", name, ".pdf")
savefig(joinpath(output_dir, file_name))

gr()
plot(treeDA, b, type = :rates)
file_name = string("plot_rates_speciation_", name, ".pdf")
savefig(joinpath(output_dir, file_name))

gr()
plot(treeDA, b, 0.1)
file_name = string("plot_rates_speciation_quantile_", name, ".pdf")
savefig(joinpath(output_dir, file_name))

end



#### STEP 3 - ELEFANT ####


#### A - Run ELEFANT ####


rm(list=ls())


#_________________________________________________________________________________________________

path <- "/home/jcostes/work/nymphalini/ELEFANT/"

names <- c("nymphalini_test_1_sf_0.83", "nymphalini_test_2_sf_0.76", "nymphalini_test_3_sf_0.70") 

nb_cores <- 10 # must match cpus (file .sh)

#_________________________________________________________________________________________________


setwd(path)

library(phytools, lib.loc="~/work/R")
library(mvMORPH, lib.loc="~/work/R")
library(dplyr, lib.loc="~/work/R")
library(tidyr, lib.loc="~/work/R")
library(igraph, lib.loc="~/work/R")
library(reshape2, lib.loc="~/work/R")
library(ggplot2, lib.loc="~/work/R")
library(parallel, lib.loc="~/work/R")
library(ggpubr, lib.loc="~/work/R")
library(RPANDA, lib.loc="~/work/R")
library(pryr, lib.loc="~/work/R")
options(dplyr.summarise.inform = FALSE)

dir.create(file.path(getwd(), "figures/"), showWarnings = FALSE)

source("/home/jcostes/work/function_network_diversification.R") 

path_results <- paste0(getwd(),"/")

args <- commandArgs()
print(args)

nb_recon <- 250 # 250
stochastic_mapping <- TRUE

plot=FALSE
seed=3

global_metrics=FALSE



for (name in names) {
  
  
  
  #### Load data
  
  
  #____________________________________________________________________________________________
  
  network <- read.table("../nymphalini_data/network_nymphalini_plants.csv", sep=";", header=TRUE)  
  
  tree_B <- read.tree("../nymphalini_data/tree_chazot_nymphalini.tre")    # butterfly tree
  tree_B <- force.ultrametric(tree_B, method = "extend")
  tree_A <- read.tree("../nymphalini_data/tree_angiosperm_families.tre")  # plant tree
  tree_A <- force.ultrametric(tree_A, method = "extend")
  
  list_ages <- c(seq(0, 28, 1)) 
  
  #____________________________________________________________________________________________
  
  
  treesDA_B <- read.tree(paste0("../BDD/", name, ".trees")) # augmented trees from BDD
  
  dim(network) 
  Ntip(tree_B)
  Ntip(tree_A)
  
  tree_B_full <- tree_B # the most complete tree possible
  
  # remove species not present in the network 
  tree_A <- drop.tip(tree_A, tip=tree_A$tip.label[!tree_A$tip.label %in% colnames(network)])
  tree_B <- drop.tip(tree_B, tip=tree_B$tip.label[!tree_B$tip.label %in% rownames(network)])
  
  # organize the interaction matrix to match the trees
  network <- as.matrix(network)[tree_B$tip.label, tree_A$tip.label]
  
  all(tree_A$tip.label %in% colnames(network))
  all(colnames(network) %in% tree_A$tip.label)
  
  network <- network[,which(colSums(network)>0)]
  network <- network[which(rowSums(network)>0),]
  
  tree_A <- drop.tip(tree_A, tip=tree_A$tip.label[!tree_A$tip.label %in% colnames(network)])
  tree_B <- drop.tip(tree_B, tip=tree_B$tip.label[!tree_B$tip.label %in% rownames(network)])
  
  
  dim(network) 
  Ntip(tree_B) 
  Ntip(tree_A) 
  
  
  # write final trees
  write.tree(tree_B, file = file.path(path_results, "tree_B_nymphalini_run_ELEFANT.tre"))
  write.tree(tree_A, file = file.path(path_results, "tree_A_plants_run_ELEFANT.tre"))
  
  
  
  #### Get the traits at present
  
  # singular value decomposition 
  s <- svd(network) 
  L <- s$u
  D <- s$d 
  R <- s$v
  
  # percentage variance explained by each singular value 
  explained_variance <- D^2/sum(D^2)*100
  round(explained_variance)
  
  d <- min(which(cumsum(explained_variance)>=90)) # at least 90% of the present-day variance
  d <- max(c(d,3)) # at least 3 traits
  
  
  E <- matrix(0, nrow = d,  ncol = d)
  diag(E) <- D[1:d]
  l <- L[,1:d]
  r <- t(R[,1:d])
  V <- l%*%E%*%r
  
  
  rownames(V) <- rownames(l) <- rownames(network)
  colnames(V) <- colnames(r) <- colnames(network)
  plot(network, V)
  
  traits_extant_B <- l[tree_B$tip.label,]
  traits_extant_A <- t(r)[tree_A$tip.label,]
  
  
  # Thresholding approach maximizing Youden's J statistic 
  thresh_proba <- threshold_proba_interaction(network, V, tol=0.001)
  
  
  results_summary <- c()
  nb_A <- ncol(network)
  nb_B <- nrow(network)
  
  
  
  #### Perform cross validation
  
  perc_cv_A=0
  perc_cv_B=0.1 # percentage removed
  res_cv <- cross_validation(tree_extant_A=tree_A, tree_extant_B=tree_B, network=network, threshold=thresh_proba,
                             d=d, perc_cv_A, perc_cv_B, nb_cv=100, seed=1, echo=FALSE, obligate_A = FALSE, obligate_B = TRUE)
  print(res_cv)
  
  
  
  results_summary <- rbind(results_summary, c("cv_true_positive", 0, nb_A, nb_B,  NA, res_cv[1]))
  results_summary <- rbind(results_summary, c("cv_false_positive",  0, nb_A, nb_B, NA, res_cv[2]))
  results_summary <- rbind(results_summary, c("cv_f1score",  0, nb_A, nb_B, NA, res_cv[3]))
  
  
  tree_extant_A <- tree_A
  tree_extant_B <- tree_B
  
  
  
  #### Load augmented trees
  
  for (i in 1:nb_recon){
    tree_DA <- treesDA_B[[i]]
    if (!is.ultrametric(tree_DA)) {
      f = file() # silence
      sink(file = f)
      tree_DA <- force.ultrametric(tree_DA, method = "extend")
      sink()
      close(f)
    }
    tree_DA <- root(tree_DA, node=which(node.depth.edgelength(tree_DA)==0))
    tree_DA$tip.label[grep("^t",tree_DA$tip.label)] <- paste0("DA_B", 1:length(grep("^t",tree_DA$tip.label)))
    tree_DA$tip.label <- gsub("DA_Bt","DA_B",tree_DA$tip.label)
    treesDA_B[[i]] <- tree_DA
  }
  
  
  
  #### Measure global metrics and phylogenetic signal
  
  if (global_metrics){
    original_network_metrics <- get_network_metrics(network)
    print(original_network_metrics)
    mantel_test <- phylosignal_network(network, tree_A = tree_A, tree_B=tree_B, method = "Jaccard_binary")
    mantel_test_unifrac <- phylosignal_network(network, tree_A = tree_A, tree_B=tree_B, method = "UniFrac_unweighted")
    
    results_summary <- rbind(results_summary, c("d",  0, nb_A, nb_B,  NA, d))
    results_summary <- rbind(results_summary, c("thresh_proba",  0, nb_A, nb_B,  NA, thresh_proba))
    results_summary <- rbind(results_summary, c("nb_A",  "0_observed", nb_A, nb_B, original_network_metrics[1], original_network_metrics[1]))
    results_summary <- rbind(results_summary, c("nb_B",  "0_observed", nb_A, nb_B,  original_network_metrics[2], original_network_metrics[2]))
    results_summary <- rbind(results_summary, c("connectance",  "0_observed", nb_A, nb_B, original_network_metrics[3], original_network_metrics[3]))
    results_summary <- rbind(results_summary, c("nestedness",  "0_observed", nb_A, nb_B, original_network_metrics[4], original_network_metrics[4]))
    results_summary <- rbind(results_summary, c("modularity",  "0_observed", nb_A, nb_B, original_network_metrics[5], original_network_metrics[5]))
    results_summary <- rbind(results_summary, c("phylosig_cladeA",  "0_observed", nb_A, nb_B,  mantel_test[c(3,4)]))
    results_summary <- rbind(results_summary, c("phylosig_cladeB",  "0_observed", nb_A, nb_B,  mantel_test[c(6,7)]))
    results_summary <- rbind(results_summary, c("phylosig_unifrac_cladeA",  "0_observed", nb_A, nb_B, mantel_test_unifrac[c(3,4)]))
    results_summary <- rbind(results_summary, c("phylosig_unifrac_cladeB",  "0_observed", nb_A, nb_B, mantel_test_unifrac[c(6,7)]))
  }
  
  
  
  ##### Run ELEFANT
  
  data_augmentation_A=FALSE
  data_augmentation_B=TRUE
  treesDA_A <- NULL
  
  
  obligate_A <- FALSE
  obligate_B <- TRUE
  evolution_A <- TRUE
  
  
  # List extant species in the network
  list_extant_A <- tree_extant_A$tip.label
  list_extant_B <- tree_extant_B$tip.label
  
  
  # List extant species used for the diversification analyses (may not be in the network)
  list_full_extant_A <- NULL
  list_full_extant_B <- tree_B_full$tip.label
  
  
  null_model=FALSE
  save_DA=TRUE
  
  
  list_network_metrics <- run_ELEFANT(name=name, nb_recon, list_ages=list_ages, tree_extant_A=tree_extant_A, tree_extant_B=tree_extant_B, 
                                      list_extant_A = list_extant_A, list_extant_B = list_extant_B,
                                      list_full_extant_A = list_full_extant_A, list_full_extant_B = list_full_extant_B,
                                      traits_extant_A=traits_extant_A, traits_extant_B=traits_extant_B, 
                                      d=d, E=E, thresh_proba=thresh_proba, stochastic_mapping=stochastic_mapping, 
                                      obligate_A=obligate_A, obligate_B=obligate_B, evolution_A=evolution_A, 
                                      data_augmentation_A=data_augmentation_A, data_augmentation_B=data_augmentation_B, 
                                      treesDA_A=treesDA_A, treesDA_B=treesDA_B, 
                                      global_metrics=global_metrics, null_model=null_model, save_DA=save_DA,
                                      seed=3, nb_cores=nb_cores, path_results=path_results)
  
  
  # Temporary saved
  save(list = ls(), file=paste0(path_results, "results_", name, ".Rdata"))
  write.table(results_summary, paste0(path_results, "results_",name, ".csv"), sep=";", quote=FALSE, row.names = FALSE)
  # load(paste0(path_results, "results_", name, ".Rdata"))
  
  # Process outputs
  
  
  
  #### Look at metrics
  
  
  if (global_metrics){
    for (age in list_ages){
      results_summary <- rbind(results_summary, c("nb_A", age, nb_A, nb_B, NA, mean(list_network_metrics$nb_A[which(list_network_metrics$age==age)], na.rm=TRUE)))
      results_summary <- rbind(results_summary, c("nb_B", age, nb_A, nb_B, NA, mean(list_network_metrics$nb_B[which(list_network_metrics$age==age)], na.rm=TRUE)))
      results_summary <- rbind(results_summary, c("connectance", age, nb_A, nb_B, NA, mean(list_network_metrics$connectance[which(list_network_metrics$age==age)], na.rm=TRUE)))
      results_summary <- rbind(results_summary, c("nestedness", age, nb_A, nb_B, NA, mean(list_network_metrics$nestedness[which(list_network_metrics$age==age)], na.rm=TRUE)))
      results_summary <- rbind(results_summary, c("modularity", age, nb_A, nb_B, NA, mean(list_network_metrics$modularity[which(list_network_metrics$age==age)], na.rm=TRUE)))
      
      results_summary <- rbind(results_summary, c("mantel_cor_A", age, nb_A, nb_B, NA, mean(list_network_metrics$mantel_cor_A[which(list_network_metrics$age==age)], na.rm=TRUE)))
      results_summary <- rbind(results_summary, c("pvalue_A", age, nb_A, nb_B, NA, mean(list_network_metrics$pvalue_A[which(list_network_metrics$age==age)], na.rm=TRUE)))
      results_summary <- rbind(results_summary, c("mantel_cor_B", age, nb_A, nb_B, NA, mean(list_network_metrics$mantel_cor_B[which(list_network_metrics$age==age)], na.rm=TRUE)))
      results_summary <- rbind(results_summary, c("pvalue_B", age, nb_A, nb_B, NA, mean(list_network_metrics$pvalue_B[which(list_network_metrics$age==age)], na.rm=TRUE)))
      
      results_summary <- rbind(results_summary, c("rate_coext_A", age, nb_A, nb_B, NA, mean(list_network_metrics$rate_coext_A[which(list_network_metrics$age==age)], na.rm=TRUE)))
      results_summary <- rbind(results_summary, c("rate_coext_B", age, nb_A, nb_B, NA, mean(list_network_metrics$rate_coext_B[which(list_network_metrics$age==age)], na.rm=TRUE)))
      
      
      ### 95% confidence interval
      results_summary <- rbind(results_summary, c("nb_A_CI", age, nb_A, nb_B, quantile(list_network_metrics$nb_A[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$nb_A[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
      results_summary <- rbind(results_summary, c("nb_B_CI", age, nb_A, nb_B, quantile(list_network_metrics$nb_B[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$nb_B[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
      results_summary <- rbind(results_summary, c("connectance_CI", age, nb_A, nb_B, quantile(list_network_metrics$connectance[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$connectance[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
      results_summary <- rbind(results_summary, c("nestedness_CI", age, nb_A, nb_B, quantile(list_network_metrics$nestedness[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$nestedness[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
      results_summary <- rbind(results_summary, c("modularity_CI", age, nb_A, nb_B, quantile(list_network_metrics$modularity[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$modularity[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
      
      results_summary <- rbind(results_summary, c("mantel_cor_A_CI", age, nb_A, nb_B, quantile(list_network_metrics$mantel_cor_A[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$mantel_cor_A[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
      results_summary <- rbind(results_summary, c("mantel_cor_B_CI", age, nb_A, nb_B, quantile(list_network_metrics$mantel_cor_B[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$mantel_cor_B[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
      
      results_summary <- rbind(results_summary, c("rate_coext_A_CI", age, nb_A, nb_B, quantile(list_network_metrics$rate_coext_A[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$rate_coext_A[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
      results_summary <- rbind(results_summary, c("rate_coext_B_CI", age, nb_A, nb_B, quantile(list_network_metrics$rate_coext_B[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$rate_coext_B[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
      
      
      #### Null model
      
      if (null_model){
        results_summary <- rbind(results_summary, c("nb_A_NM", age, nb_A, nb_B, NA, mean(list_network_metrics$nb_A_NM[which(list_network_metrics$age==age)], na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("nb_B_NM", age, nb_A, nb_B, NA, mean(list_network_metrics$nb_B_NM[which(list_network_metrics$age==age)], na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("connectance_NM", age, nb_A, nb_B, NA, mean(list_network_metrics$connectance_NM[which(list_network_metrics$age==age)], na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("nestedness_NM", age, nb_A, nb_B, NA, mean(list_network_metrics$nestedness_NM[which(list_network_metrics$age==age)], na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("modularity_NM", age, nb_A, nb_B, NA, mean(list_network_metrics$modularity_NM[which(list_network_metrics$age==age)], na.rm=TRUE)))
        
        results_summary <- rbind(results_summary, c("mantel_cor_A_NM", age, nb_A, nb_B, NA, mean(list_network_metrics$mantel_cor_A_NM[which(list_network_metrics$age==age)], na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("pvalue_A_NM", age, nb_A, nb_B, NA, mean(list_network_metrics$pvalue_A_NM[which(list_network_metrics$age==age)], na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("mantel_cor_B_NM", age, nb_A, nb_B, NA, mean(list_network_metrics$mantel_cor_B_NM[which(list_network_metrics$age==age)], na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("pvalue_B_NM", age, nb_A, nb_B, NA, mean(list_network_metrics$pvalue_B_NM[which(list_network_metrics$age==age)], na.rm=TRUE)))
        
        results_summary <- rbind(results_summary, c("rate_coext_A_NM", age, nb_A, nb_B, NA, mean(list_network_metrics$rate_coext_A_NM[which(list_network_metrics$age==age)], na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("rate_coext_B_NM", age, nb_A, nb_B, NA, mean(list_network_metrics$rate_coext_B_NM[which(list_network_metrics$age==age)], na.rm=TRUE)))
        
        
        ### 95% confidence interval
        results_summary <- rbind(results_summary, c("nb_A_CI_NM", age, nb_A, nb_B, quantile(list_network_metrics$nb_A_NM[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$nb_A_NM[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("nb_B_CI_NM", age, nb_A, nb_B, quantile(list_network_metrics$nb_B_NM[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$nb_B_NM[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("connectance_CI_NM", age, nb_A, nb_B, quantile(list_network_metrics$connectance_NM[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$connectance_NM[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("nestedness_CI_NM", age, nb_A, nb_B, quantile(list_network_metrics$nestedness_NM[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$nestedness_NM[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("modularity_CI_NM", age, nb_A, nb_B, quantile(list_network_metrics$modularity_NM[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$modularity_NM[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
        
        results_summary <- rbind(results_summary, c("mantel_cor_A_CI_NM", age, nb_A, nb_B, quantile(list_network_metrics$mantel_cor_A_NM[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$mantel_cor_A_NM[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("mantel_cor_B_CI_NM", age, nb_A, nb_B, quantile(list_network_metrics$mantel_cor_B_NM[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$mantel_cor_B_NM[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
        
        results_summary <- rbind(results_summary, c("rate_coext_A_CI_NM", age, nb_A, nb_B, quantile(list_network_metrics$rate_coext_A_NM[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$rate_coext_A_NM[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
        results_summary <- rbind(results_summary, c("rate_coext_B_CI_NM", age, nb_A, nb_B, quantile(list_network_metrics$rate_coext_B_NM[which(list_network_metrics$age==age)], 0.025, na.rm=TRUE), quantile(list_network_metrics$rate_coext_B_NM[which(list_network_metrics$age==age)], 0.975, na.rm=TRUE)))
      }
      
    }
    
    save(list = ls(), file=paste0(path_results, "results_", name, ".Rdata"))
    write.table(results_summary, paste0(path_results, "results_",name, ".csv"), sep=";", quote=FALSE, row.names = FALSE)
  }
  
  
  
  #### Plot ancestral networks (following Braga)
  
  
  
  #### Consider an interaction as true is present in >XX% (keep % reconstructions to plot confidence)
  
  for (past in list_ages){
    
    print(past)
    
    for (threshold_prop_inter in c(0.10, 0.25, 0.50, thresh_proba)){
      
      print(threshold_prop_inter)
      
      list_interactions_past <- read.table(file = paste0(path_results, "list_interactions_past_",name,"_",past,".csv"), sep=";", header=TRUE)
      
      # Convert back with one row per replicate 
      list_interactions_past <- list_interactions_past %>%
        separate_rows(rep, sep = "-")
      
      colnames(list_interactions_past) <- c("species_B", "species_A", "rep")
      list_interactions_past$species_B <- as.character(list_interactions_past$species_B)
      list_interactions_past$species_A <- as.character(list_interactions_past$species_A)
      list_interactions_past$rep <- as.character(list_interactions_past$rep)
      
      list_interactions_past_extant <- list_interactions_past
      if (length(grep("DA_A", list_interactions_past_extant$species_A))>0) {
        list_interactions_past_extant <- list_interactions_past_extant[-grep("DA_A", list_interactions_past_extant$species_A),]
        list_interactions_past$species_A[grep("DA_A", list_interactions_past$species_A)] <- paste0(list_interactions_past$species_A[grep("DA_A", list_interactions_past$species_A)], "_", list_interactions_past$rep[grep("DA_A", list_interactions_past$species_A)])
      }
      if (length(grep("DA_B", list_interactions_past_extant$species_B))>0) {
        list_interactions_past_extant <- list_interactions_past_extant[-grep("DA_B", list_interactions_past_extant$species_B),]
        list_interactions_past$species_B[grep("DA_B", list_interactions_past$species_B)] <- paste0(list_interactions_past$species_B[grep("DA_B", list_interactions_past$species_B)], "_", list_interactions_past$rep[grep("DA_B", list_interactions_past$species_B)])
      }
      
      
      list_interactions_past_extant <- as.data.frame(list_interactions_past_extant %>%
                                                       group_by(species_B, species_A) %>%
                                                       summarize(count = n()))
      list_interactions_past_extant$species_B <- as.character(list_interactions_past_extant$species_B)
      list_interactions_past_extant$species_A <- as.character(list_interactions_past_extant$species_A)
      list_interactions_past_extant$count <- as.numeric(as.character(list_interactions_past_extant$count))
      
      list_interactions_past_extant$count <- list_interactions_past_extant$count/nb_recon
      
      # only keep interactions > threshold_prop_inter %
      list_interactions_past_extant <- list_interactions_past_extant[which(list_interactions_past_extant$count>=threshold_prop_inter),]
      
      
      plot_augmeted <- FALSE
      
      ### Add augmented partners 
      if (plot_augmeted==TRUE) {
        list_interactions_past_augmented <- add_augmented_partners(list_interactions_past_extant, list_interactions_past, threshold_prop_inter)
      } else {
        list_interactions_past_augmented <- list_interactions_past_extant
      }
      
      
      # Plot augmented graph:
      
      if (nrow(list_interactions_past_augmented)>0){
        
        new_tips_B <- paste0("b", 1:length(list_full_extant_B))
        names(new_tips_B) <- sort(list_full_extant_B)
        new_tips_A <- paste0("a", 1:length(list_full_extant_A))
        names(new_tips_A) <- sort(list_full_extant_A)
        
        g <- graph_from_edgelist(as.matrix(list_interactions_past_augmented[,1:2]), directed = FALSE)
        edge.attributes(g)$weight <- as.numeric(list_interactions_past_augmented[,3])
        
        V(g)$type <- rep("Clade A", length(V(g)))
        V(g)$type[which(names(V(g)) %in% c(tree_B$tip.label, new_tips_B))] <- "Clade B"
        # for (tip in c(tree_B$tip.label)) { # wrong ? 
        #   V(g)$type[grep(paste0(tip,"_"), names(V(g)))] <- "Clade B"
        #   V(g)$type[grep(paste0("_",tip), names(V(g)))] <- "Clade B"
        # }
        for (tip in c(tree_B$tip.label, new_tips_B)) { 
          V(g)$type[grep(paste0(tip,"-"), names(V(g)))] <- "Clade B"
          V(g)$type[grep(paste0("-",tip), names(V(g)))] <- "Clade B"
        }
        V(g)$type[grep("DA_B", names(V(g)))] <- "Clade B"
        V(g)$type_col <- V(g)$type
        #V(g)$type_col[names(V(g))=="A_1"] <- "A_1"
        V(g)$augmented <- "FALSE"
        V(g)$augmented[grep("DA_A|augmented_B", names(V(g)))] <- "TRUE"
        
        # define color and shape mappings.
        shape <- c("square", "circle")
        size <- c(5, 2.5)
        col <- c( "#229954", "#6e2c00")
        names(shape) <- c("Clade A", "Clade B")
        names(col) <- c("Clade A", "Clade B")
        names(size) <- c("FALSE", "TRUE")
        
        color_simul=FALSE
        color_edge <- rep("grey",nrow(list_interactions_past_augmented))
        if (color_simul==TRUE){
          list_interactions_past_augmented$pair <- paste0(list_interactions_past_augmented$species_B," - ", list_interactions_past_augmented$species_A)
          color_edge[list_interactions_past_augmented$pair %in% simulate_interactions$pair] <- "#212f3d"
        }
        
        pdf(paste0(path_results, "/figures/ancestral_network_",name, "_thresh_", round(threshold_prop_inter, 2), "_t",past,".pdf"), width=7, height = 7)
        set.seed(1)
        plot(g, vertex.color = col[V(g)$type_col], vertex.shape = shape[V(g)$type], vertex.label=NA, vertex.size=size[V(g)$augmented],
             edge.width=5*E(g)$weight, edge.color=color_edge, layout= layout_with_fr(g, weights = 5*E(g)$weight,niter = 5000 ) )
        plot(1,1,col="white")
        legend("bottomleft", pch=19,names(col),col=col,cex=0.8)
        dev.off()
      } else {
        pdf(paste0(path_results, "/figures/ancestral_network_",name, "_thresh_", round(threshold_prop_inter, 2), "_t",past,".pdf"), width=7, height = 7)
        dev.off()
      }
    }
    
  }
  
  
  save(list = ls(), file=paste0(path_results, "results_", name, ".Rdata"))
  write.table(results_summary, paste0(path_results, "results_",name, ".csv"), sep=";", quote=FALSE, row.names = FALSE)
  
}



#### B - Link with diversification ####


rm(list=ls())


#_______________________________________________________________________________

path <- "/home/jcostes/work/nymphalini/ELEFANT/"

path_BDD <- "/home/jcostes/work/nymphalini/BDD/"

nb_cores <- 10

#_______________________________________________________________________________


setwd(path)

library(phytools, lib.loc="~/work/R")
library(mvMORPH, lib.loc="~/work/R")
library(dplyr, lib.loc="~/work/R")
library(tidyr, lib.loc="~/work/R")
library(igraph, lib.loc="~/work/R")
library(reshape2, lib.loc="~/work/R")
library(ggplot2, lib.loc="~/work/R")
library(parallel, lib.loc="~/work/R")
library(ggpubr, lib.loc="~/work/R")
library(RPANDA, lib.loc="~/work/R")
library(pryr, lib.loc="~/work/R")
options(dplyr.summarise.inform = FALSE)

path_results <- path

nb_digits <- 3 # nb significant digits

args <- commandArgs()
print(args)


#_________________________________________________________________________________________________________

names_tests <- c("nymphalini_test_1_sf_0.83", "nymphalini_test_2_sf_0.76", "nymphalini_test_3_sf_0.70")
# constant extinction ; global sampling fraction 

# open the tree used for summarizing the data
tree_B <- read.tree("tree_B_nymphalini_run_ELEFANT.tre") # tree generated by ELEFANT  

list_ages <- c(seq(0, 28, 1))

#_________________________________________________________________________________________________________


for (name_original in names_tests) {
  
  consensus_tree_name <- paste0(name_original, "_mean_rates.nex")
  # consensus tree = complete tree with speciation rates information (tree generated by BDD)
  
  
  list_names <- c(name_original, paste0(name_original, "_random_", 1:50), 
                  paste0(name_original, "_random_nb_partners_", 1:50))
  
  list_names <- c(name_original,
                  sample(paste0(name_original, "_random_nb_partners_", 1:50)),
                  sample(paste0(name_original, "_random_", 1:50)))
  
  
  stats <- TRUE
  stats <- FALSE
  
  echo=FALSE
  long=TRUE
  
  
  list_alpha <- unique(c(seq(0.01, 0.1, 0.01), seq(0.2, 2, 0.1), seq(2, 10, 1)))
  
  
  for (name in list_names){
    
    # verification of necessary files
    
    if (file.exists(paste0(path_results, "results_", name, ".csv"))){
      
      if (!file.exists(paste0("list_rates_ancestral_interactions_",name,".csv"))){
        
        # display the name being processed
        print(name) 
        
        file_name <- paste0(path_BDD, consensus_tree_name)
        
        # loading the tree
        tree_rates <- read.nexus(file_name) 
        
        # extracting tree metadata
        read.nexus.metadata <- function(file){
          
          ph <- read.nexus(file)
          X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
          tab	<- X[grep("tree 1 =", X)]
          tab <- gsub("tree 1 = ", "", tab)
          tab <- unlist(strsplit(tab, "\\["))[-1]
          tab <- gsub("&|;|\\]", "", tab)
          tab <- gsub(":.+$", "", tab)
          tab <- lapply(tab, function(x) unlist(strsplit(x, "="))	)
          tab	<- lapply(tab, function(x){
            x	<- gsub('{','',x,fixed=1)
            x	<- gsub(',}','',x,fixed=1)
            x	<- gsub(',dt','',x,fixed=1)
            x	<- gsub(',fdt','',x,fixed=1)
            x	<- gsub(',da','',x,fixed=1)
            x <- x[x!="sr"]
            return(x)
          })
          
          tab <- data.frame(matrix(unlist(tab), ncol = 4, byrow = TRUE), stringsAsFactors = FALSE)
          
          tree <- X[grep("tree 1", X)]
          tree <- gsub("tree 1 = ", "", tree)
          tree <- gsub("\\[.*?\\]", "", tree)
          tree <- substr(tree, 1, nchar(tree)-1)
          
          tmp <- unlist(strsplit(tree, ":"))
          tmp <- sapply(seq_along(tmp), function(i) paste0(tmp[i],"_xxx",i))
          tmp	<- paste(paste(tmp, collapse=':'),';',sep='')
          tmp	<- read.tree(text=tmp)
          
          tab$name <- NA
          tab$min_age <- NA
          tab$name[as.numeric(sapply(strsplit(split="_xxx", tmp$tip.label), "[[", 2))] <- sapply(strsplit(split="_xxx", tmp$tip.label), "[[", 1)
          tab$min_age[as.numeric(sapply(strsplit(split="_xxx", tmp$tip.label), "[[", 2))] <- 0
          tmp$tip.label <- sapply(strsplit(split="_xx", tmp$tip.label), "[[", 1)
          
          for (node in tmp$node.label){
            node_number <- as.numeric(gsub("_xxx", "", node))
            tab$name[node_number] <- paste(sort(extract.clade(tmp, node=node)$tip.label), collapse ="-")
            tab$min_age[node_number] <- max(node.depth.edgelength(extract.clade(tmp, node=node)))
          }
          
          colnames(tab) <- c("list_rates", "dt", "fdt", "da", "name", 'min_age')
          
          tab$dt <- as.numeric(as.character(tab$dt))
          tab$fdt <- as.numeric(as.character(tab$fdt))
          tab$min_age <- as.numeric(as.character(tab$min_age))
          
          return(tab[,c("name", "list_rates", "dt", "fdt", "min_age")])
        }
        
        table_rates <- read.nexus.metadata(file_name)
        
        
        # Age max each node
        table_rates$max_age <- 0
        
        for (i in 1:nrow(table_rates)){
          species <- unlist(strsplit(split="-", table_rates$name[i]))
          if (length(species)<Ntip(tree_rates)){
            if (length(species)>1) {table_rates$max_age[i] <- node.depth.edgelength(tree_rates)[tree_rates$edge[which(tree_rates$edge[,2]==getMRCA(tree_rates, species)), 1]]
            }else{ table_rates$max_age[i] <- node.depth.edgelength(tree_rates)[tree_rates$edge[which(tree_rates$edge[,2]==which(tree_rates$tip.label==species)), 1]] }
          }
        }
        table_rates$max_age <-  max(node.depth.edgelength(tree_rates)) - table_rates$max_age 
        
        
        # Update tree names
        new_tip_names <- read.table(paste0(path_results, "list_interactions_past_storage_name_",name,".csv"), sep=";", header = TRUE)
        new_tips_B <- new_tip_names$storage_name[grep("b", new_tip_names$storage_name)]
        names(new_tips_B) <- new_tip_names$species_name[grep("b", new_tip_names$storage_name)]
        
        
        ### Removed species not in the network from table_rates$name 
        additional_species <- names(new_tips_B)[!names(new_tips_B) %in% tree_B$tip.label]
        for (i in 1:nrow(table_rates)){
          list_species <- unlist(strsplit(split="-", table_rates$name[i]))
          if (length(which(list_species %in% additional_species))>0){
            if (!all(list_species %in% additional_species)){
              table_rates$name[i] <- paste0(list_species[!list_species %in% additional_species], collapse = "-")
            }
          }
        }
        
        
        ##### Link between diversification and species interactions
        
        extract_inter_rates <- function(name, age, new_tips_B, 
                                        table_rates, list_alpha, 
                                        echo=FALSE, long=TRUE){
          
          print(age)
          
          list_inter_store <- read.table(paste0(path_results, "list_interactions_past_",name,"_",age,".csv"), sep=";", header=TRUE)
          
          # Update the name
          for (species in sort(new_tips_B, decreasing = TRUE)) {
            list_inter_store$species_B <- gsub(species, names(new_tips_B)[which(new_tips_B==species)],  list_inter_store$species_B, fixed=TRUE)
          }
          list_lineages <- names(table(list_inter_store$species_B))
          
          # Remove DA lineages (as we cannot compare them)
          if (length(grep("DA_B", list_lineages))>0) list_lineages <- list_lineages[-grep("DA_B", list_lineages)]
          
          if (!all(list_lineages %in% table_rates$name)){
            print(paste0("Problem with wrong lineages in ", name))
          } # should be TRUE 
          
          # Convert back with one row per replicate 
          list_inter <- list_inter_store %>%
            separate_rows(rep, sep = "-")
          rm(list_inter_store)
          
          ### Extract information
          
          nb_recon <- length(unique(list_inter$rep))
          problem <- FALSE
          if (nb_recon<250) {  # REMOVE
            print(paste0("Problem with less than 250 rep in ", name))
            problem <- TRUE
          }
          
          # Create table
          if (long) {
            table <- matrix(nrow = length(list_lineages), ncol=4+length(list_alpha), 0)
          } else {
            table <- matrix(nrow = length(list_lineages), ncol=4, 0)
          }
          rownames(table) <- list_lineages
          
          for (species in list_lineages){
            
            if (echo) print(species)
            
            # Mean number of partners
            mean_nb_partners <- length(which(list_inter$species_B==species))/nb_recon
            
            # Log speciation rate
            rate_species <- table_rates[which(table_rates$name==species),,drop=FALSE] # can correspond to several rows if species are in the tree but not in the network
            rate_species <- rate_species[intersect(which(rate_species$min_age<age+1e-10),which(rate_species$max_age>age-1e-10)),,drop=FALSE]
            list_rates <- as.numeric(unlist(strsplit(split=",", rate_species$list_rates)))
            time_steps <-  c(rate_species$max_age   - (0:(length(list_rates)-2))*rate_species$dt, rate_species$min_age)
            log_rate <- list_rates[which.min(abs(age-time_steps))[1]]
            
            # Ehrlich and Raven hypothesis 
            
            if (long){
              nb_inter_partners <- integer(0)
              for (recon in 1:nb_recon) {
                list_partners <- list_inter$species_A[list_inter$species_B == species & list_inter$rep == recon]
                partner_counts <- table(list_inter$species_A[list_inter$rep == recon])
                nb_inter_partners <- c(nb_inter_partners, partner_counts[list_partners])
              }
            }
            
            
            list_info <- c(age, species, signif(mean_nb_partners, nb_digits), signif(log_rate,nb_digits))
            
            if (long){
              list_info <- c(list_info, sapply(list_alpha, function(alpha) signif(sum((nb_inter_partners)^(-alpha))/length(nb_inter_partners), nb_digits)) )
            }
            
            table[species,] <- list_info
          }
          
          table <- data.frame(table)
          rownames(table) <- NULL
          colnames(table) <- c("age", "species", "mean_nb_partners", "log_rate", paste0("sum_exp_alpha_", list_alpha))
          table$mean_nb_partners <- as.numeric(as.character(table$mean_nb_partners))
          table$log_rate <- as.numeric(as.character(table$log_rate))
          for (var in grep("sum_exp", colnames(table))){table[,var] <- as.numeric(as.character(table[,var]))}
          
          if (!problem) {
            return(table)
          } else {
            table <- rep(NA, length(table))
          }
          
        }
        
        list_all_table <- do.call(rbind, mclapply(list_ages, extract_inter_rates, mc.preschedule=F, mc.cores = nb_cores,
                                                  name=name, new_tips_B=new_tips_B, table_rates=table_rates, list_alpha=list_alpha))
        
        
        if (!any(is.na(list_all_table[, 1]))){ # don't save if NA (meaning that same run have not been saved)
          write.table(list_all_table, paste0("list_rates_ancestral_interactions_",name,".csv"), sep=";", row.names=FALSE, quote=FALSE)
        } else {
          print(paste0("Problem with NA in ", name))
        }
      }
    }
  }
  
}



#### STEP 4 - ANALYZE THE RESULTS ####

#### directly on the computer


# scp jcostes@genobioinfo.toulouse.inrae.fr:/home/jcostes/work/nymphalini/ELEFANT/list_rates_ancestral_interactions_*csv /mnt/c/Users/julie/Documents/M2/STAGE/Results/Nymphalini/ELEFANT/
# scp jcostes@genobioinfo.toulouse.inrae.fr:/home/jcostes/work/nymphalini/BDD/plot* /mnt/c/Users/julie/Documents/M2/STAGE/Results/Nymphalini/BDD/


rm(list=ls())

library(phytools)
library(mvMORPH)
library(dplyr)
library(igraph)
library(reshape2)
library(ggplot2)
library(parallel)
library(ggpubr)
library(RPANDA)
options(dplyr.summarise.inform = FALSE)
library(nlme)
library(tidyr)


# setwd("~/Documents/STAGE/Results/Nymphalini/ELEFANT/")
setwd("~/Documents/STAGE/Results/Nymphalini/ELEFANT/all_interactions/")


source("~/Documents/STAGE/Script/function_network_diversification.R")


dir.create(file.path(getwd(), "/plot_ER_score_div/"), showWarnings = FALSE)


list_names_original <- list_names <- c("nymphalini_test_1_sf_0.83", "nymphalini_test_2_sf_0.76", "nymphalini_test_3_sf_0.70")



### Prepare the covariance matrix

name <- list_names[1] # just need one name to prepare the tree

list_all_table <-  read.table(paste0("list_rates_ancestral_interactions_",name,".csv"), sep=";", head=TRUE)

colnames(list_all_table)


tree <- read.tree("~/Documents/STAGE/Data/nymphalini_data/tree_chazot_nymphalini.tre") # the most complete tree possible

# network <- read.table("~/Documents/STAGE/Data/nymphalini_data/network_nymphalini_plants.csv", sep=";", header=TRUE)
network <- read.table("~/Documents/STAGE/Data/nymphalini_data/network_nymphalini_all.csv", sep=";", header=TRUE)
network <- network[rownames(network) %in% tree$tip.label,]


### Make tree with all the tips corresponding to time steps
list_species_network <- rownames(network)
list_ages <- sort(unique(list_all_table$age), decreasing = TRUE)
list_ages <- list_ages[list_ages>0]
list_ages <- list_ages[list_ages!=28]
print(list_ages)
# list_ages <- list_ages[list_ages %% 2 ==0]
# list_ages <- list_ages[list_ages %% 5 ==0]

tree_fit_link <- add_past_tips(tree, list_ages, list_species_network)


#### Compute -- Individual plots  ####

list_all_models <- c()

for (name in list_names){
  
  if
  (file.exists(paste0("list_rates_ancestral_interactions_",name,".csv"))){
    
    print(name)
    
    list_all_table <- read.table(paste0("list_rates_ancestral_interactions_",name,".csv"), sep=";", head=TRUE)
    list_all_table <- subset(list_all_table, age != 28)
    
    list_all_table$species_age <- paste0(list_all_table$age, "-", list_all_table$species)
    
    # Add 0 = constant rates
    list_all_table$sum_exp_alpha_0 <- 1
    
    # Fit link BDD and ELEFANT model (with PGLS)
    list_aic <- fit_link_BDD_ELEFANT(list_all_table, tree_fit_link)
    
    alpha_best <- as.numeric(list_aic$alpha[which.min(list_aic$AIC)])
    best_model <- (list_aic$model[which.min(list_aic$AIC)])
    AIC_best <- min(list_aic$AIC)
    
    
    ### Look at best E&R model (excluding constant)
    list_aic_ER <- list_aic[grep("sum_exp", list_aic$name),]
    list_aic_ER <- list_aic_ER[list_aic_ER$name!="sum_exp_alpha_0",]
    alpha_best_ER <- as.numeric(list_aic_ER$alpha[which.min(list_aic_ER$AIC)])
    list_all_table$sum_exp_best <- list_all_table[,grep("sum_exp_alpha", colnames(list_all_table))[which.min(list_aic_ER$AIC)]]
    
    AIC_constant <- list_aic$AIC[list_aic$name=="sum_exp_alpha_0"]
    AIC_ER <- list_aic$AIC[which(list_aic$model=="E&R")][which.min(list_aic$AIC[which(list_aic$model=="E&R")])]
    AIC_temporal <- list_aic$AIC[which(list_aic$model=="temporal")][which.min(list_aic$AIC[which(list_aic$model=="temporal")])]
    AIC_nbpartners <- list_aic$AIC[which(list_aic$model=="nbpartners")][which.min(list_aic$AIC[which(list_aic$model=="nbpartners")])]
    c(AIC_constant, AIC_ER, AIC_temporal, AIC_nbpartners)
    
    # AIC weights ###
    
    aic_values <- c(AIC_constant, AIC_temporal, AIC_nbpartners, AIC_ER) # list of AIC values
    
    AIC_min <- min(aic_values)  # find the minimal AIC
    delta_AIC <- aic_values - AIC_min  # calculate the AIC differences (ΔAIC) from the best AIC
    
    weights <- exp(-delta_AIC / 2)  # calculate the probability of each model
    sum_weights <- sum(weights)  # sum of the weights
    AIC_weights <- weights / sum_weights  # normalization of the weights
    
    # assign the AIC weights to the variables
    AICweight_constant <- AIC_weights[1]
    AICweight_temporal <- AIC_weights[2]
    AICweight_nbpartners <- AIC_weights[3]
    AICweight_ER <- AIC_weights[4]
    
    ###
    
    R2_constant <- list_aic$R2[list_aic$name=="sum_exp_alpha_0"]
    R2_ER <- list_aic$R2[which(list_aic$model=="E&R")][which.min(list_aic$AIC[which(list_aic$model=="E&R")])]
    R2_temporal <- list_aic$R2[which(list_aic$model=="temporal")][which.min(list_aic$AIC[which(list_aic$model=="temporal")])]
    R2_nbpartners <- list_aic$R2[which(list_aic$model=="nbpartners")][which.min(list_aic$AIC[which(list_aic$model=="nbpartners")])]
    
    lambda0_constant <- list_aic$lambda0[list_aic$name=="sum_exp_alpha_0"]
    lambda0_ER <- list_aic$lambda0[which(list_aic$model=="E&R")][which.min(list_aic$AIC[which(list_aic$model=="E&R")])]
    lambda0_temporal <- list_aic$lambda0[which(list_aic$model=="temporal")][which.min(list_aic$AIC[which(list_aic$model=="temporal")])]
    lambda0_nbpartners <- list_aic$lambda0[which(list_aic$model=="nbpartners")][which.min(list_aic$AIC[which(list_aic$model=="nbpartners")])]
    
    #if ((length(grep("_1$", name))==1)|(length(grep("_2$", name))==1)|(length(grep("_3$", name))==1)|(length(grep("_4$", name))==1)){
      pdf(paste0("plot_ER_score_div/link_speciation_rates_interactions_",name,".pdf"), width = 5, height=3.5)
      p <- ggplot(list_all_table, aes(x=sum_exp_best, y=exp(log_rate), col=age))+geom_point()+
        ylab("Speciation rate (/Myr)") + xlab("Ehrlich & Raven score")+theme_bw()+
        ggtitle( bquote(lambda[0] == .(round(lambda0_ER, 2)) ~ ", " ~ alpha == .(round(alpha_best_ER,2)) ) )
      if (alpha_best_ER>0) p <- p + geom_abline(slope = lambda0_ER, intercept = 0, col="orange", lwd=1.25)
      if (alpha_best_ER==0) p <- p + geom_hline(yintercept = lambda0_ER, col="orange", lwd=1.25)
      print(p)
      dev.off()
    #}
    
    # plot: tree with E&R score
    tree_df <- fortify(tree)
    tree_df <- tree_df %>%
      mutate(species = ifelse(isTip, label, NA_character_))
    internal_nodes <- tree_df %>%
      filter(!isTip) %>%
      arrange(desc(x)) %>%
      pull(node)
    for (node_id in internal_nodes) {
      children <- tree_df %>% filter(parent == node_id)
      children_names <- children$species
      children_names <- ifelse(is.na(children_names), children$label, children_names)
      children_names <- children_names[!is.na(children_names)]
      if (length(children_names) > 0) {
        new_name <- paste(sort(children_names), collapse = "-")
        tree_df$species[tree_df$node == node_id] <- new_name
      }
    }
    tree_df <- tree_df %>%
      left_join(list_all_table %>% select(species, sum_exp_best), by = "species")
    plot_tree_ER_score <- ggtree(tree) %<+% tree_df + 
      geom_tree(aes(color = sum_exp_best), size = 1) +  # colorier selon sum_exp_best
      scale_color_viridis_c(option = "plasma") +        # palette couleur continue
      theme_tree2()                                     # thème classique pour arbres
    print(plot_tree_ER_score)
    ggsave("~/Documents/STAGE/Results/Nymphalini/ELEFANT/all_interactions/fig/plot_tree_ER_score.pdf", plot_tree_ER_score)
    
    list_all_models <- rbind(list_all_models, c(name, best_model, AIC_best, alpha_best, 
                                                AIC_constant, AIC_temporal, AIC_nbpartners, AIC_ER, 
                                                AICweight_constant, AICweight_temporal, AICweight_nbpartners, AICweight_ER,
                                                alpha_best_ER, R2_constant, R2_temporal, R2_nbpartners, R2_ER,
                                                lambda0_constant, lambda0_temporal, lambda0_nbpartners, lambda0_ER))
  }
}


list_all_models <- data.frame(list_all_models)
colnames(list_all_models) <- c("name", "best_model", "AIC_best", "alpha_best",
                               "AIC_constant", "AIC_temporal", "AIC_nbpartners", "AIC_ER", 
                               "AICweight_constant", "AICweight_temporal", "AICweight_nbpartners", "AICweight_ER",
                               "alpha_best_ER", "R2_constant", "R2_temporal", "R2_nbpartners", "R2_ER",
                               "lambda0_constant", "lambda0_temporal", "lambda0_nbpartners", "lambda0_ER")

table(list_all_models$best_model)
list_all_models$AIC_best <- as.numeric(as.character(list_all_models$AIC_best))
list_all_models$AIC_constant <- as.numeric(as.character(list_all_models$AIC_constant))
list_all_models$AIC_temporal <- as.numeric(as.character(list_all_models$AIC_temporal))
list_all_models$AIC_nbpartners <- as.numeric(as.character(list_all_models$AIC_nbpartners))
list_all_models$AIC_ER <- as.numeric(as.character(list_all_models$AIC_ER))
list_all_models$alpha_best <- as.numeric(as.character(list_all_models$alpha_best))
list_all_models$alpha_best_ER <- as.numeric(as.character(list_all_models$alpha_best_ER))
list_all_models$AICweight_constant <- as.numeric(as.character(list_all_models$AICweight_constant))
list_all_models$AICweight_temporal <- as.numeric(as.character(list_all_models$AICweight_temporal))
list_all_models$AICweight_nbpartners <- as.numeric(as.character(list_all_models$AICweight_nbpartners))
list_all_models$AICweight_ER <- as.numeric(as.character(list_all_models$AICweight_ER))

list_all_models$type <- "original"
list_all_models$type[grep("random",list_all_models$name)] <- "random"
list_all_models$type[grep("nb_partners",list_all_models$name)] <- "random_nbpartners"


# Colors 
colors_model <- c("constant" = "blue", 
                  "temporal" = "red", 
                  "nb partners" = "green", 
                  "E&R" = "purple")


pdf("results_model_selection_nymphalini.pdf")
ggplot(list_all_models, aes(x=type, fill=best_model))+
  geom_bar(position = "stack") + # or use "stack" for a stacked bar chart
  labs(title = "",x = "",y = "Count") +
  theme_bw() +
  scale_fill_manual(values=colors_model, name="Model") +
  facet_wrap(~ name, scales = "free")
dev.off()


pdf("results_R2_nymphalini.pdf")
list_R2 <- data.frame(name=c(list_all_models$name,list_all_models$name,list_all_models$name),
             R2=as.numeric(as.character(c(list_all_models$R2_temporal,
                                          list_all_models$R2_ER, list_all_models$R2_nbpartners))),
             type=rep(c("temporal", "E&R", "nb partners"),
                      nrow(list_all_models)))
ggplot(list_R2, aes(x=type, y=R2))+
  geom_boxplot() +
  labs(title = "",x = "",y = "R2") +
  theme_bw() 
dev.off()


write.table(list_all_models, "results_models_all_nymphalini.csv",
            sep=";", row.names=FALSE, quote=FALSE)


# Plot AIC weights

# Convert AIC weights to long format
list_all_models_long <- pivot_longer(list_all_models, 
                                     cols = c(AICweight_constant, AICweight_temporal, 
                                              AICweight_nbpartners, AICweight_ER),
                                     names_to = "model", values_to = "AIC_weight")
# Rename the models for clarity
list_all_models_long$model <- factor(list_all_models_long$model, 
                                     levels = c("AICweight_constant", "AICweight_temporal", 
                                                "AICweight_nbpartners", "AICweight_ER"),
                                     labels = c("constant", "temporal", "nb partners", "E&R"))
# Plot
pdf("results_AIC_weights_nymphalini.pdf")
p <- ggplot(list_all_models_long, aes(x=type, y=AIC_weight, fill=model)) +
  geom_bar(stat="identity", position="stack") +
  labs(title="AIC Weights", x="", y="AIC Weight") +
  theme_bw() +
  scale_fill_manual(values=colors_model, name="Model") +
  facet_wrap(~ name, scales = "free")
print(p)
dev.off()



#### STOP ####
#### randomisations plus tard éventuellement 


#### restrict on E&R hypothesis 

list_all_models$lambda_0_original <- NA
list_all_models$AIC_original <- NA
list_all_models$alpha_original <- NA
list_all_models$name_original <- NA
for (name in list_names_original){
  list_all_models$lambda_0_original[grep(paste0(name,"_random_"), list_all_models$name)] <- list_all_models$lambda_0_ER[which(list_all_models$name==name)]
  list_all_models$AIC_original[grep(paste0(name,"_random_"), list_all_models$name)] <- list_all_models$AIC_ER[which(list_all_models$name==name)]
  list_all_models$alpha_original[grep(paste0(name,"_random_"), list_all_models$name)] <- list_all_models$alpha_best_ER[which(list_all_models$name==name)]
  list_all_models$name_original[grep(paste0(name,"_random_"), list_all_models$name)] <- name
}
list_all_random <- list_all_models[grep("random", list_all_models$name),]
list_all_random$delta_AIC <- list_all_random$AIC_ER - list_all_random$AIC_original
list_all_random$type <- "random"
list_all_random$type[grep("nb_partners", list_all_random$name)] <- "random_nb_partners"

list_all_models[-grep("random", list_all_models$name),]

list_all_random$new_name <- "constant extinction (sf 1)"
list_all_random$new_name[grep("test_2", list_all_random$name)] <- "constant extinction (sf 2)"
list_all_random$new_name[grep(paste0("test_3"), list_all_random$name)] <- "constant turnover (sf 1)"
list_all_random$new_name[grep(paste0("test_4"), list_all_random$name)] <- "constant turnover (sf 2)"


library(dplyr)
# Calculating mean, lower, and upper CI
list_all_random_CI <- list_all_random %>%
  group_by(new_name, type) %>%
  summarize(lower = quantile(delta_AIC, probs = 0.025),
            upper = quantile(delta_AIC, probs = 0.975),
            median = median(delta_AIC), # Optional: For central tendency
            .groups = "drop")

pdf("delta_AIC_random_model_nointercept.pdf", width=6, height=4)

ggplot(list_all_random, aes(x=new_name, y=delta_AIC))+geom_boxplot(fill="orange")+
  theme_bw()+theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ylab("Diff. AIC (random - model)")+xlab("")+geom_hline(yintercept = 0, lwd=0.75)+
  geom_hline(yintercept = 2, lwd=0.5)

ggplot(list_all_random, aes(x=new_name, y=delta_AIC))+geom_violin(fill="orange")+
  theme_bw()+theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ylab("Diff. AIC (random - model)")+xlab("")+geom_hline(yintercept = 0, lwd=0.75)+
  geom_hline(yintercept = 2, lwd=0.5)

# ggplot(list_all_random, aes(x=new_name, y=delta_AIC, fill=type))+geom_boxplot()+
#   theme_bw()+theme(axis.text.x = element_text(angle = 60, hjust = 1))+
#   ylab("Diff. AIC (random - model)")+xlab("")+geom_hline(yintercept = 0, lwd=0.75)+
#   geom_hline(yintercept = 2, lwd=0.5)
ggplot(list_all_random, aes(x=new_name, y=delta_AIC, fill=type))+geom_violin()+
  theme_bw()+theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ylab("Diff. AIC (random - model)")+xlab("")+geom_hline(yintercept = 0, lwd=0.75)+
  geom_hline(yintercept = 2, lwd=0.5)

ggplot(list_all_random_CI, aes(x = new_name, y = median, color = type)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylab("Diff. AIC (random - model)") +
  xlab("") +
  geom_hline(yintercept = 0, lwd = 0.75) +
  geom_hline(yintercept = 2, lwd = 0.5)

dev.off()

pdf("est_lambda0_random_model_nointercept.pdf", width=6, height=4)
ggplot(list_all_random, aes(x=lambda_0_original, y=lambda_0_ER, col=new_name))+geom_point()+geom_abline(slope=1, intercept = 0)+
  theme_bw()+ylab("lambda0 (random)")+xlab("lambda0 (model)")

ggplot(list_all_random, aes(x=lambda_0_original, y=lambda_0_ER, col=new_name, shape=type))+geom_point()+geom_abline(slope=1, intercept = 0)+
  theme_bw()+ylab("lambda0 (random)")+xlab("lambda0 (model)")
dev.off()


pdf("est_alpha_random_model_nointercept.pdf", width=6, height=4)
p <- ggplot(list_all_random, aes(x=as.factor(alpha_original), y=alpha_best_ER, col=new_name))+geom_boxplot()+geom_abline(slope=1, intercept = 0)+
  theme_bw()+ylab("alpha (random)")+xlab("alpha (model)")
ind <- 1
for (alpha in sort(unique(list_all_random$alpha_original))){
  list_all_random
  p <- p + annotate("point", x = ind, y = alpha, col="black", shape=18, size=5)
  ind <- ind+1
}
print(p)
dev.off()






#### B - Compute per age slice  ####
#### plus tard éventuellement


list_all_models <- c()
#time_slice <- c(0-5, 6-10, 11-15, 16-20, 21-25, 26-30, 31-35, 35-38)
time_slice <- c("0-10", "11-20", "21-30", "31-38")

for (name in list_names){
  
  if (file.exists(paste0("list_rates_ancestral_interactions_",name,".csv"))){
    
    print(name)
    
    for (slice in time_slice){
      
      list_all_table <- read.table(paste0("list_rates_ancestral_interactions_",name,".csv"), sep=";", head=TRUE)
      
      list_all_table$species_age <- paste0(list_all_table$age, "-", list_all_table$species)
      
      lower=as.numeric(strsplit(split="-", slice)[[1]][1])
      upper=as.numeric(strsplit(split="-", slice)[[1]][2])
      
      list_all_table <- list_all_table[which(list_all_table$age>=lower),]
      list_all_table <- list_all_table[which(list_all_table$age<=upper),]
      
      inv_L <- inv_L_all[list_all_table$species_age, list_all_table$species_age]
      
      # Add 0 = constant rates
      list_all_table$sum_exp_alpha_0 <- 1
      
      # Fit link BDD and ELEFANT model :
      list_aic <- fit_link_BDD_ELEFANT(list_all_table, inv_L)
      
      alpha_best <- as.numeric(list_aic$alpha[which.min(list_aic$AIC_cholesky)])
      best_model <- (list_aic$model[which.min(list_aic$AIC_cholesky)])
      lambda_0_best <- list_aic$lambda_0_cholesky[which.min(list_aic$AIC_cholesky)]
      AIC_best <- min(list_aic$AIC_cholesky)
      
      
      ### Look at best E&R model (including constant)
      list_aic_ER <- list_aic[grep("sum_exp", list_aic$name),]
      alpha_best_ER <- as.numeric(list_aic_ER$alpha[which.min(list_aic_ER$AIC_cholesky)])
      lambda_0_ER <- as.numeric(list_aic_ER$lambda_0_cholesky[which.min(list_aic_ER$AIC_cholesky)])
      list_all_table$sum_exp_best <- list_all_table[,grep("sum_exp_alpha", colnames(list_all_table))[which.min(list_aic_ER$AIC_cholesky)]]
      
      AIC_constant <- list_aic$AIC_cholesky[list_aic$name=="sum_exp_alpha_0"]
      AIC_ER <- list_aic$AIC_cholesky[which.min(list_aic$AIC_cholesky[which(list_aic$model=="E&R")])]
      AIC_temporal <- list_aic$AIC_cholesky[which(list_aic$model=="temporal")][which.min(list_aic$AIC_cholesky[which(list_aic$model=="temporal")])]
      AIC_nbpartners <- list_aic$AIC_cholesky[which(list_aic$model=="nbpartners")][which.min(list_aic$AIC_cholesky[which(list_aic$model=="nbpartners")])]
      
      list_all_models <- rbind(list_all_models, c(name, slice, best_model, AIC_best, lambda_0_best, alpha_best,
                                                  AIC_constant, AIC_temporal, AIC_nbpartners, AIC_ER, lambda_0_ER, alpha_best_ER))
      
    }
  }
}

list_all_models <- data.frame(list_all_models)
colnames(list_all_models) <- c("name", "slice", "best_model", "AIC_best", "lambda_0_best", "alpha_best",
                               "AIC_constant", "AIC_temporal", "AIC_nbpartners", "AIC_ER", "lambda_0_ER", "alpha_best_ER")
table(list_all_models$slice, list_all_models$best_model)
list_all_models$lambda_0_best <- as.numeric(as.character(list_all_models$lambda_0_best))
list_all_models$lambda_0_ER <- as.numeric(as.character(list_all_models$lambda_0_ER))
list_all_models$AIC_best <- as.numeric(as.character(list_all_models$AIC_best))
list_all_models$AIC_constant <- as.numeric(as.character(list_all_models$AIC_constant))
list_all_models$AIC_temporal <- as.numeric(as.character(list_all_models$AIC_temporal))
list_all_models$AIC_nbpartners <- as.numeric(as.character(list_all_models$AIC_nbpartners))
list_all_models$AIC_ER <- as.numeric(as.character(list_all_models$AIC_ER))
list_all_models$alpha_best <- as.numeric(as.character(list_all_models$alpha_best))
list_all_models$alpha_best_ER <- as.numeric(as.character(list_all_models$alpha_best_ER))


write.table(list_all_models, "results_models_all_nymphalini_slice.csv", sep=";", row.names=FALSE, quote=FALSE)


list_all_models$lambda_0_original <- NA
list_all_models$AIC_original <- NA
list_all_models$alpha_original <- NA
list_all_models$name_original <- NA
for (name in list_names_original){
  list_all_models$lambda_0_original[grep(paste0(name,"_random_"), list_all_models$name)] <- list_all_models$lambda_0_ER[which(list_all_models$name==name)]
  list_all_models$AIC_original[grep(paste0(name,"_random_"), list_all_models$name)] <- list_all_models$AIC_ER[which(list_all_models$name==name)]
  list_all_models$alpha_original[grep(paste0(name,"_random_"), list_all_models$name)] <- list_all_models$alpha_best_ER[which(list_all_models$name==name)]
  list_all_models$name_original[grep(paste0(name,"_random_"), list_all_models$name)] <- name
}
list_all_random <- list_all_models[grep("random", list_all_models$name),]
list_all_random$delta_AIC <- list_all_random$AIC_ER - list_all_random$AIC_original
list_all_random$type <- "random"
list_all_random$type[grep("nb_partners", list_all_random$name)] <- "random_nb_partners"


list_all_original <- list_all_models[-grep("random", list_all_models$name),]
list_all_original


list_all_random$new_name <- "constant extinction (sf 1)"
list_all_random$new_name[grep("test_2", list_all_random$name)] <- "constant extinction (sf 2)"
list_all_random$new_name[grep(paste0("test_3"), list_all_random$name)] <- "constant turnover (sf 1)"
list_all_random$new_name[grep(paste0("test_4"), list_all_random$name)] <- "constant turnover (sf 2)"


library(dplyr)
# Calculating mean, lower, and upper CI
list_all_random_CI <- list_all_random %>%
  group_by(new_name, type, slice) %>%
  summarize(lower = quantile(delta_AIC, probs = 0.025),
            upper = quantile(delta_AIC, probs = 0.975),
            median = median(delta_AIC), # Optional: For central tendency
            .groups = "drop")

pdf("delta_AIC_random_model_nointercept_slice.pdf", width=7.5, height=4)
ggplot(list_all_random, aes(x=new_name, y=delta_AIC))+geom_boxplot(fill="orange")+
  theme_bw()+theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ylab("Diff. AIC (random - model)")+xlab("")+geom_hline(yintercept = 0, lwd=0.75)+
  geom_hline(yintercept = 2, lwd=0.5) + facet_grid(~slice)

ggplot(list_all_random, aes(x=new_name, y=delta_AIC, fill=type))+geom_boxplot()+
  theme_bw()+theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  ylab("Diff. AIC (random - model)")+xlab("")+geom_hline(yintercept = 0, lwd=0.75)+
  geom_hline(yintercept = 2, lwd=0.5) + facet_grid(~slice)

ggplot(list_all_random_CI, aes(x = new_name, y = median, color = type)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylab("Diff. AIC (random - model)") +
  xlab("") +
  geom_hline(yintercept = 0, lwd = 0.75) +
  geom_hline(yintercept = 2, lwd = 0.5)+ facet_grid(~slice)

dev.off()


pdf("est_lambda0_random_model_nointercept_slice.pdf", width=6, height=4)
ggplot(list_all_random, aes(x=lambda_0_original, y=lambda_0_ER, col=new_name))+geom_point()+geom_abline(slope=1, intercept = 0)+
  theme_bw()+ylab("lambda0 (random)")+xlab("lambda0 (model)")

ggplot(list_all_random, aes(x=lambda_0_original, y=lambda_0_ER, col=new_name, shape=type))+geom_point()+geom_abline(slope=1, intercept = 0)+
  theme_bw()+ylab("lambda0 (random)")+xlab("lambda0 (model)")

ggplot(list_all_random, aes(x=lambda_0_original, y=lambda_0_ER, col=slice, shape=type))+geom_point()+geom_abline(slope=1, intercept = 0)+
  theme_bw()+ylab("lambda0 (random)")+xlab("lambda0 (model)")
dev.off()


pdf("est_alpha_random_model_nointercept_slice.pdf", width=6, height=4)
p <- ggplot(list_all_random, aes(x=as.factor(alpha_original), y=alpha, col=new_name))+geom_boxplot()+ #geom_abline(slope=1, intercept = 0)+
  theme_bw()+ylab("alpha (random)")+xlab("alpha (model)")
ind <- 1
for (alpha in sort(unique(list_all_random$alpha_original))){
  list_all_random
  p <- p + annotate("point", x = ind, y = alpha, col="black", shape=18, size=5)
  ind <- ind+1
}
print(p)

p <- ggplot(list_all_random, aes(x=as.factor(alpha_original), y=alpha, col=slice))+geom_boxplot()+ #geom_abline(slope=1, intercept = 0)+
  theme_bw()+ylab("alpha (random)")+xlab("alpha (model)")
ind <- 1
for (alpha in sort(unique(list_all_random$alpha_original))){
  list_all_random
  p <- p + annotate("point", x = ind, y = alpha, col="black", shape=18, size=5)
  ind <- ind+1
}
print(p)
dev.off()


