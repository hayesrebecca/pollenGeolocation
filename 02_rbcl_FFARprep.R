## from this
## tutorial:https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121

## tutorial also have some good pcoa plots and heatmaps, which arent
## on here
#dir.bombus <- '/Volumes/bombus/Dropbox (University of Oregon)'
rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path

setwd(local.path)
setwd('sunflower_saved')



if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
##BiocManager::install("TreeSummarizedExperiment")
library(TreeSummarizedExperiment)
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)

library(tidyr)
library(dplyr)
library(bipartite)
library(phyloseq)
library(TreeTools)
library(devtools)
library(ape)
library(picante)

source("../pnw_survey/data_prep/sequence_prep/src/misc.R")

# Run again with merged runs!
## reading artifacts
# qza.rbcl.path  <- "SI_pipeline/merged/RBCL/final"

# For now, just run 2023 until we have the rest of the seq data
qza.rbcl.path  <- "ffar_pipeline/merged/RBCL/final"


## function from: https://stackoverflow.com/questions/38570074/phylogenetics-in-r-collapsing-descendant-tips-of-an-internal-node

drop_dupes <- function(tree, thres=1e-5){
  tips <- which(tree$edge[,2] %in% 1:Ntip(tree))
  toDrop <- tree$edge.length[tips] < thres
  drop.tip(tree,tree$tip.label[toDrop])
}


## ***********************************************************************
## RBCL
## ***********************************************************************
# ## R0
# weightedUFrbclqzaR0 <- read_qza(file.path(qza.rbcl.path,
#                                           'core_metricsRBCL/weighted_unifrac_distance_matrix.qza'))
# 
# unweightedUFrbclqzaR0 <- read_qza(file.path(qza.rbcl.path,
#                                             'core_metricsRBCL/unweighted_unifrac_distance_matrix.qza'))
# # 
# # ## access data inside artifacts
# # 
# # # RBCL
# # # R0
# wphylo.dist.rbclR0 <- weightedUFrbclqzaR0$data
# phylo.dist.rbclR0 <- unweightedUFrbclqzaR0$data
# 
# # 
# # #Note, when taxonomy is imported, a single string is returned along
# # #with a confidence score.  For many analysis we will want to break up
# # #this string and for that purpose the parse_taxonomy() function is
# # #provided:
# # 
# # # R0
# taxonomyRBCLR0 <-
#   read_qza("pipeline/R2023/RBCL/final/core_metricsRBCL/rarefied_table.qza")
# taxonomyRBCLR0 <- taxonomyRBCLR0$data
# 
# # 
# # #Phyloseq tutorials here: https://joey711.github.io/phyloseq/
# # 
# # # ****CHECK PATHS BELOW ONCE FINISH RBCL PIPELINE****
# # 
# # ## NEED TO CHECK whether this actually has all of the samples and tips
# # ## TODO: I think I need to rerun the pipeline steps to make sure the correct qzas are used
# # 
# # #RBCL 2023 phylogeny
# physeqRBCLR0 <- qza_to_phyloseq(
#   features="pipeline/R2023/RBCL/final/core_metricsRBCL/rarefied_table.qza",
#   tree="pipeline/R2023/RBCL/rooted-treeRBCL2023.qza",
#   "pipeline/R2023/RBCL/taxonomyRBCL.qza",
#   metadata = "pipeline/R2023/maps/pnw2023mapRBCL.txt"
# )
# physeqRBCLR0
# plot(physeqRBCLR0@phy_tree, show.tip.label = FALSE)
# # 

ffar <- read.csv("../../Documents/1ROTATION kern lab/prepBeeRBCL/FFARspec.csv")



# ## ***********************************************************************
# 


## read in rep-seqs to get DNA seq and feature id. right now tip labels are feature ids,
## we want to switch to actual DNA seqs

library(Biostrings)

# Replace with your actual file path
fasta_file <- "pipeline/R2023/RBCL/rep-seqs-RBCL.fasta"

# Read the fasta file
fasta_seqs <- readDNAStringSet(fasta_file)

# Convert to data.frame
fasta_df <- data.frame(
  FeatureID = names(fasta_seqs),
  Sequence = as.character(fasta_seqs),
  stringsAsFactors = FALSE
)


# ## add 16s/rbcl to make the columns easy to find when they are added
# ## to the specimen data
# 
#feature.2.tax.rbcl$Taxon  <- gsub(" ", "_", feature.2.tax.rbcl$Taxon)
fasta_df$Sequence <- paste("RBCL", fasta_df$Sequence,
                                  sep=':')
## R0
tree.rbclR0 <- phy_tree(physeqRBCLR0, errorIfNULL=TRUE)

# First, filter fasta_df to only those in the tree
fasta_df_filtered <- fasta_df[fasta_df$FeatureID %in% tree.rbclR0$tip.label, ]

# Then reassign the tip labels using match on the filtered set
tree.rbclR0$tip.label <- fasta_df_filtered$Sequence[
  match(tree.rbclR0$tip.label, fasta_df_filtered$FeatureID)
]


#fasta_df <- fasta_df[fasta_df$FeatureID %in% tree.rbclR0$tip.label,]

# ## ***********************************************************************
# ## rbcl networks
# ## ***********************************************************************
# 
#R0
indiv.comm.rbclR0 <-
    bipartite::empty(catchDups(makeComm(taxonomyRBCLR0,
                                        fasta_df_filtered,
                                        feature.col="FeatureID")))
                                        # changed to "Feature.ID" to match 'makeComm' function
colnames(indiv.comm.rbclR0) <- fasta_df_filtered$Sequence
indiv.comm.rbclR0 <- indiv.comm.rbclR0/rowSums(indiv.comm.rbclR0)

bees.rbcl <- c(rownames(indiv.comm.rbclR0))
# 
comms <- list(indiv.comm.rbclR0)
# 

species.rbcl <- unique(unlist(sapply(comms, colnames)))
# 
merged.comm.rbcl <- plyr::rbind.fill(lapply(comms, as.data.frame))
# 
# ## check with number of columns against number of unique species
dim(merged.comm.rbcl)
length(species.rbcl)
# 
rownames(merged.comm.rbcl) <- bees.rbcl
# 
megaRBCLdata <- read_qza("pipeline/R2023/RBCL/rooted-treeRBCL2023.qza")
tree.rbcl <- megaRBCLdata$data
# 
tree.rbcl <- tip_glom(tree.rbcl, h=0.1)
# 
tree.rbcl <- drop_dupes(tree.rbcl, thres=1e-5)
# 
tree.rbcl$tip.label <- fasta_df_filtered$Sequence[
  match(tree.rbcl$tip.label, fasta_df_filtered$FeatureID)
]

## ***********************************************************************
## Make mega dataset
## ***********************************************************************

## spec already includes parasite data
load('../pnw_survey/data/spec_net_fire.Rdata')

indiv.comm.rbcl <- as.data.frame(merged.comm.rbcl)
pollen <- colnames(indiv.comm.rbcl)
indiv.comm.rbcl$SampleID <- rownames(indiv.comm.rbcl)
indiv.comm.rbcl$SampleID <- sub("^", "NCASI-S", indiv.comm.rbcl$SampleID)


spec.net <-cbind(spec.net, indiv.comm.rbcl[, pollen][match(spec.net$SampleID,
                           indiv.comm.rbcl$SampleID),])


## check for duplicates
#any(duplicated(spec.net$UniqueID))

save(spec.net, file= "../pollenGeolocation/data/NCASIpollen.Rdata")

write.csv(spec.net, file= "../pollenGeolocation/data/NCASIpollen.csv",
          row.names=FALSE)



