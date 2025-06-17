## from this
## tutorial:https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121
## tutorial also have some good pcoa plots and heatmaps, which arent
## on here

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


## reading artifacts
#qza.16s.path  <- "ffar_pipeline/merged/16s/final/"
qza.rbcl.path  <- "ffar_pipeline/merged/RBCL/final/"

# 16s
#weightedUF16sqzaR0 <- read_qza(file.path(qza.16s.path,
#                             'core_metrics16sR0/weighted_unifrac_distance_matrix.qza'))
#weightedUF16sqzaR1 <- read_qza(file.path(qza.16s.path,
#                             'core_metrics16sR1/weighted_unifrac_distance_matrix.qza'))
#weightedUF16sqzaR2 <- read_qza(file.path(qza.16s.path,
#                              'core_metrics16sR2/weighted_unifrac_distance_matrix.qza'))
#weightedUF16sqzaR3 <- read_qza(file.path(qza.16s.path,
#                              'core_metrics16sR3/weighted_unifrac_distance_matrix.qza'))
#weightedUF16sqzaR4 <- read_qza(file.path(qza.16s.path,
#                              'core_metrics16sR4/weighted_unifrac_distance_matrix.qza'))

# unweightedUF16sqzaR0 <- read_qza(file.path(qza.16s.path,
#                             'core_metrics16sR0/unweighted_unifrac_distance_matrix.qza'))
# unweightedUF16sqzaR1 <- read_qza(file.path(qza.16s.path,
#                             'core_metrics16sR1/unweighted_unifrac_distance_matrix.qza'))
# unweightedUF16sqzaR2 <- read_qza(file.path(qza.16s.path,
#                             'core_metrics16sR2/unweighted_unifrac_distance_matrix.qza'))
# unweightedUF16sqzaR3 <- read_qza(file.path(qza.16s.path,
#                             'core_metrics16sR3/unweighted_unifrac_distance_matrix.qza'))
# unweightedUF16sqzaR4 <- read_qza(file.path(qza.16s.path,
#                             'core_metrics16sR4/unweighted_unifrac_distance_matrix.qza'))
# 
# # 16s
# wphylo.dist.16sR0 <- weightedUF16sqzaR0$data
# phylo.dist.16sR0 <-unweightedUF16sqzaR0$data
# 
# wphylo.dist.16sR1 <- weightedUF16sqzaR1$data
# phylo.dist.16sR1 <-unweightedUF16sqzaR1$data
# 
# wphylo.dist.16sR2 <- weightedUF16sqzaR2$data
# phylo.dist.16sR2 <-unweightedUF16sqzaR2$data
# 
# wphylo.dist.16sR3 <- weightedUF16sqzaR3$data
# phylo.dist.16sR3 <-unweightedUF16sqzaR3$data
# 
# wphylo.dist.16sR4 <- weightedUF16sqzaR4$data
# phylo.dist.16sR4 <-unweightedUF16sqzaR4$data
# 
# ## Note, when taxonomy is imported, a single string is returned along
# ## with a confidence score.  For many analysis we will want to break
# ## up this string and for that purpose the parse_taxonomy() function
# ## is provided:
# 
# tax.path <- "ffar_pipeline/merged/16s/final/"
# taxonomy16sR0 <- read_qza(
#     file.path(tax.path, "core_metrics16sR0/rarefied_table.qza"))
# taxonomy16sR0 <- taxonomy16sR0$data
# 
# taxonomy16sR1 <- read_qza(
#     file.path(tax.path, "core_metrics16sR1/rarefied_table.qza"))
# taxonomy16sR1 <- taxonomy16sR1$data
# 
# taxonomy16sR2 <- read_qza(
#     file.path(tax.path, "core_metrics16sR2/rarefied_table.qza"))
# taxonomy16sR2 <- taxonomy16sR2$data
# 
# taxonomy16sR3 <- read_qza(
#     file.path(tax.path, "core_metrics16sR3/rarefied_table.qza"))
# taxonomy16sR3 <- taxonomy16sR3$data
# 
# taxonomy16sR4 <- read_qza(
#     file.path(tax.path, "core_metrics16sR4/rarefied_table.qza"))
# taxonomy16sR4 <- taxonomy16sR4$data
# 
# ## 16s R0 phylogeny
# physeq16sR0 <- qza_to_phyloseq(
#     features=
#         file.path(tax.path, "core_metrics16sR0/rarefied_table.qza"),
#     tree="ffar_pipeline/merged/16s/rooted-tree16s.qza",
#     "ffar_pipeline/merged/16s/taxonomy16s.qza",
#     metadata = "ffar_pipeline/merged/16s/maps/ffar2018map16s.txt"
# )
# physeq16sR0
# ##plot(physeq16sR0@phy_tree)
# 
# ## plot(physeq16sR0@phy_tree, show.tip.label = FALSE)
# 
# ##16s R1 phylogeny
# physeq16sR1 <- qza_to_phyloseq(
#     features=
#         file.path(tax.path, "core_metrics16sR1/rarefied_table.qza"),
#     tree="ffar_pipeline/merged/16s/rooted-tree16s.qza",
#     "ffar_pipeline/merged/16s/taxonomy16s.qza",
#     metadata = "ffar_pipeline/merged/16s/maps/ffar2019_R1map16s.txt"
# )
# physeq16sR1
# ##plot(physeq16sR1@phy_tree)
# 
# ## plot(physeq16sR1@phy_tree, show.tip.label = FALSE)
# 
# 
# ##16s R2 phylogeny
# physeq16sR2 <- qza_to_phyloseq(
#     features=
#         file.path(tax.path, "core_metrics16sR2/rarefied_table.qza"),
#     tree="ffar_pipeline/merged/16s/rooted-tree16s.qza",
#     "ffar_pipeline/merged/16s/taxonomy16s.qza",
#     metadata = "ffar_pipeline/merged/16s/maps/ffar2019_R2map16s.txt"
# )
# physeq16sR2
# ##plot(physeq16sR2@phy_tree)
# 
# ## plot(physeq16sR2@phy_tree, show.tip.label = FALSE)
# 
# 
# ##16s R3 phylogeny
# physeq16sR3 <- qza_to_phyloseq(
#     features=
#         file.path(tax.path, "core_metrics16sR3/rarefied_table.qza"),
#     tree="ffar_pipeline/merged/16s/rooted-tree16s.qza",
#     "ffar_pipeline/merged/16s/taxonomy16s.qza",
#     metadata = "ffar_pipeline/merged/16s/maps/ffar2019_R3map16s.txt"
# )
# physeq16sR3
# ## plot(physeq16sR3@phy_tree)
# 
# ## plot(physeq16sR3@phy_tree, show.tip.label = FALSE)
# 
# ##16s R4 phylogeny
# physeq16sR4 <- qza_to_phyloseq(
#     features=
#         file.path(tax.path, "core_metrics16sR4/rarefied_table.qza"),
#     tree="ffar_pipeline/merged/16s/rooted-tree16s.qza",
#     "ffar_pipeline/merged/16s/taxonomy16s.qza",
#     metadata = "ffar_pipeline/merged/16s/maps/ffar2019_R4map16s.txt"
# )
# physeq16sR4
# 
# 
# ## plot(physeq16sR4@phy_tree)
# ## plot(physeq16sR4@phy_tree, show.tip.label = FALSE)
# 
# feature.2.tax.16s <-
#     read.table("ffar_pipeline/merged/16s/taxonomy16s.txt", sep="\t",
#                header=TRUE)
# 
# feature.2.tax.16s$Taxon <- paste("16s", feature.2.tax.16s$Taxon, sep=':')
# 
# ## convert to a phylo class which is more useful downstream
# tree.16sR0 <- phy_tree(physeq16sR0, errorIfNULL=TRUE)
# tree.16sR1 <- phy_tree(physeq16sR1, errorIfNULL=TRUE)
# tree.16sR2 <- phy_tree(physeq16sR2, errorIfNULL=TRUE)
# tree.16sR3 <- phy_tree(physeq16sR3, errorIfNULL=TRUE)
# tree.16sR4 <- phy_tree(physeq16sR4, errorIfNULL=TRUE)
# 
# ## match the tip labs to the table with feature ID and Taxon
# tree.16sR0$tip.label  <-  feature.2.tax.16s$Taxon[match(tree.16sR0$tip.label,
#                                            feature.2.tax.16s$Feature.ID)]
# tree.16sR1$tip.label  <-  feature.2.tax.16s$Taxon[match(tree.16sR1$tip.label,
#                                                         feature.2.tax.16s$Feature.ID)]
# tree.16sR1$tip.label  <-  feature.2.tax.16s$Taxon[match(tree.16sR2$tip.label,
#                                                         feature.2.tax.16s$Feature.ID)]
# tree.16sR3$tip.label  <-  feature.2.tax.16s$Taxon[match(tree.16sR3$tip.label,
#                                                         feature.2.tax.16s$Feature.ID)]
# tree.16sR4$tip.label  <-  feature.2.tax.16s$Taxon[match(tree.16sR4$tip.label,
#                                                         feature.2.tax.16s$Feature.ID)]
# 
# ## ***********************************************************************
# ## 16s networks
# ## ***********************************************************************
# 
# indiv.comm.16sR0 <-
#     bipartite::empty(catchDups(makeComm(taxonomy16sR0,
#                                         feature.2.tax.16s)))
# indiv.comm.16sR0 <- indiv.comm.16sR0/rowSums(indiv.comm.16sR0)
# 
# indiv.comm.16sR1 <-
#     bipartite::empty(catchDups(makeComm(taxonomy16sR1,
#                                         feature.2.tax.16s)))
# indiv.comm.16sR1 <- indiv.comm.16sR1/rowSums(indiv.comm.16sR1)
# 
# indiv.comm.16sR2 <-
#     bipartite::empty(catchDups(makeComm(taxonomy16sR2,
#                                         feature.2.tax.16s)))
# indiv.comm.16sR2 <- indiv.comm.16sR2/rowSums(indiv.comm.16sR2)
# 
# indiv.comm.16sR3 <-
#     bipartite::empty(catchDups(makeComm(taxonomy16sR3,
#                                         feature.2.tax.16s)))
# indiv.comm.16sR3 <- indiv.comm.16sR3/rowSums(indiv.comm.16sR3)
# 
# indiv.comm.16sR4 <-
#     bipartite::empty(catchDups(makeComm(taxonomy16sR4,
#                                         feature.2.tax.16s)))
# indiv.comm.16sR4 <- indiv.comm.16sR4/rowSums(indiv.comm.16sR4)
# 
# 
# bees.16s <- c(rownames(indiv.comm.16sR0),
#               paste0("SF",  rownames(indiv.comm.16sR1)),
#               paste0("SF", rownames(indiv.comm.16sR2)),
#                      paste0("SF", rownames(indiv.comm.16sR3)),
#               paste0("SF", rownames(indiv.comm.16sR4)))
# 
# comms <- list(indiv.comm.16sR0, indiv.comm.16sR1,
#               indiv.comm.16sR2, indiv.comm.16sR3,
#               indiv.comm.16sR4)
# 
# species.16s <- unique(unlist(sapply(comms, colnames)))
# 
# merged.comm.16s <- plyr::rbind.fill(lapply(comms, as.data.frame))
# 
# ## check with number of columns against number of unique species
# dim(merged.comm.16s)
# length(species.16s)
# 
# rownames(merged.comm.16s) <- bees.16s

## ***********************************************************************
## working with a merged 16s tree
## ***********************************************************************
#upload our mega 16s phylogenetic tree

# mega16sdata <- read_qza("ffar_pipeline/merged/16s/rooted-tree16s.qza")
# ## mega16sdata <- read_qza("../../sunflower_saved/ffar_pipeline/merged/16s/rooted-tree16s.qza")
# tree.16s <- mega16sdata$data

drop_dupes <- function(tree,thres=1e-5){
  tips <- which(tree$edge[,2] %in% 1:Ntip(tree))
  toDrop <- tree$edge.length[tips] < thres
  drop.tip(tree,tree$tip.label[toDrop])
}

## tree.16s <- tip_glom(tree.16s, h=0.05)

#tree.16s <- drop_dupes(tree.16s, thres=1e-5)


# make distance matrix

#tree.16s$tip.label  <-  feature.2.tax.16s$Taxon[match(tree.16s$tip.label,
#                                           feature.2.tax.16s$Feature.ID)]

## ***********************************************************************
## RBCL
## ***********************************************************************
weightedUFrbclqzaR0 <- read_qza(file.path(qza.rbcl.path,
                  'core_metricsRBCLR0/weighted_unifrac_distance_matrix.qza'))
weightedUFrbclqzaR1 <- read_qza(file.path(qza.rbcl.path,
                  'core_metricsRBCLR1/weighted_unifrac_distance_matrix.qza'))
weightedUFrbclqzaR2 <- read_qza(file.path(qza.rbcl.path,
                  'core_metricsRBCLR2/weighted_unifrac_distance_matrix.qza'))
weightedUFrbclqzaR3 <- read_qza(file.path(qza.rbcl.path,
                  'core_metricsRBCLR3/weighted_unifrac_distance_matrix.qza'))
weightedUFrbclqzaR4 <- read_qza(file.path(qza.rbcl.path,
                  'core_metricsRBCLR4/weighted_unifrac_distance_matrix.qza'))

unweightedUFrbclqzaR0 <- read_qza(file.path(qza.rbcl.path,
                 'core_metricsRBCLR0/unweighted_unifrac_distance_matrix.qza'))
unweightedUFrbclqzaR1 <- read_qza(file.path(qza.rbcl.path,
                 'core_metricsRBCLR1/unweighted_unifrac_distance_matrix.qza'))
unweightedUFrbclqzaR2 <- read_qza(file.path(qza.rbcl.path,
                 'core_metricsRBCLR2/unweighted_unifrac_distance_matrix.qza'))
unweightedUFrbclqzaR3 <- read_qza(file.path(qza.rbcl.path,
                  'core_metricsRBCLR3/unweighted_unifrac_distance_matrix.qza'))
unweightedUFrbclqzaR4 <- read_qza(file.path(qza.rbcl.path,
                  'core_metricsRBCLR4/unweighted_unifrac_distance_matrix.qza'))

## access data inside artifacts

# RBCL
wphylo.dist.rbclR0 <- weightedUFrbclqzaR0$data
phylo.dist.rbclR0 <- unweightedUFrbclqzaR0$data

wphylo.dist.rbclR1 <- weightedUFrbclqzaR1$data
phylo.dist.rbclR1 <- unweightedUFrbclqzaR1$data

wphylo.dist.rbclR2 <- weightedUFrbclqzaR2$data
phylo.dist.rbclR2 <- unweightedUFrbclqzaR2$data

wphylo.dist.rbclR3 <- weightedUFrbclqzaR3$data
phylo.dist.rbclR3 <- unweightedUFrbclqzaR3$data

wphylo.dist.rbclR4 <- weightedUFrbclqzaR4$data
phylo.dist.rbclR4 <- unweightedUFrbclqzaR4$data

#Note, when taxonomy is imported, a single string is returned along with a confidence score.
#For many analysis we will want to break up this string and for that
                                        #purpose the parse_taxonomy() function is provided:

taxonomyRBCLR0 <- read_qza("ffar_pipeline/merged/RBCL/final/core_metricsRBCLR0/rarefied_table.qza")
taxonomyRBCLR0 <- taxonomyRBCLR0$data

taxonomyRBCLR1 <- read_qza("ffar_pipeline/merged/RBCL/final/core_metricsRBCLR1/rarefied_table.qza")
taxonomyRBCLR1 <- taxonomyRBCLR1$data

taxonomyRBCLR2 <- read_qza("ffar_pipeline/merged/RBCL/final/core_metricsRBCLR2/rarefied_table.qza")
taxonomyRBCLR2 <- taxonomyRBCLR2$data

taxonomyRBCLR3 <- read_qza("ffar_pipeline/merged/RBCL/final/core_metricsRBCLR3/rarefied_table.qza")
taxonomyRBCLR3 <- taxonomyRBCLR3$data

taxonomyRBCLR4 <- read_qza("ffar_pipeline/merged/RBCL/final/core_metricsRBCLR4/rarefied_table.qza")
taxonomyRBCLR4 <- taxonomyRBCLR4$data


#Phyloseq tutorials here: https://joey711.github.io/phyloseq/

# ****CHECK PATHS BELOW ONCE FINISH RBCL PIPELINE****

#RBCL R0 phylogeny
physeqRBCLR0 <- qza_to_phyloseq(
  features="ffar_pipeline/merged/RBCL/final/core_metricsRBCLR0/rarefied_table.qza",
  tree="ffar_pipeline/merged/RBCL/rooted-treeRBCL.qza",
  "ffar_pipeline/merged/RBCL/taxonomyRBCL.qza",
  metadata = "ffar_pipeline/merged/RBCL/maps/ffar2018mapRBCL.txt"
)
physeqRBCLR0

## plot(physeqRBCLR0@phy_tree, show.tip.label = FALSE)

#RBCL R1 phylogeny
physeqRBCLR1 <- qza_to_phyloseq(
  features="ffar_pipeline/merged/RBCL/final/core_metricsRBCLR1/rarefied_table.qza",
  tree="ffar_pipeline/merged/RBCL/rooted-treeRBCL.qza",
  "ffar_pipeline/merged/RBCL/taxonomyRBCL.qza",
  metadata = "ffar_pipeline/merged/RBCL/maps/ffar2019_R1mapRBCL.txt"
)
physeqRBCLR1

## plot(physeqRBCLR1@phy_tree, show.tip.label = FALSE)

#RBCL R2 phylogeny
physeqRBCLR2 <- qza_to_phyloseq(
  features="ffar_pipeline/merged/RBCL/final/core_metricsRBCLR2/rarefied_table.qza",
  tree="ffar_pipeline/merged/RBCL/rooted-treeRBCL.qza",
  "ffar_pipeline/merged/RBCL/taxonomyRBCL.qza",
  metadata = "ffar_pipeline/merged/RBCL/maps/ffar2019_R2mapRBCL.txt"
)
physeqRBCLR2

## plot(physeqRBCLR2@phy_tree, show.tip.label = FALSE)

#RBCL R3 phylogeny
physeqRBCLR3 <- qza_to_phyloseq(
  features="ffar_pipeline/merged/RBCL/final/core_metricsRBCLR3/rarefied_table.qza",
  tree="ffar_pipeline/merged/RBCL/rooted-treeRBCL.qza",
  "ffar_pipeline/merged/RBCL/taxonomyRBCL.qza",
  metadata = "ffar_pipeline/merged/RBCL/maps/ffar2019_R3mapRBCL.txt"
)
physeqRBCLR3

## plot(physeqRBCLR3@phy_tree, show.tip.label = FALSE)

#RBCL R4 phylogeny
physeqRBCLR4 <- qza_to_phyloseq(
  features="ffar_pipeline/merged/RBCL/final/core_metricsRBCLR4/rarefied_table.qza",
  tree="ffar_pipeline/merged/RBCL/rooted-treeRBCL.qza",
  "ffar_pipeline/merged/RBCL/taxonomyRBCL.qza",
  metadata = "ffar_pipeline/merged/RBCL/maps/ffar2019_R4mapRBCL.txt"
)
physeqRBCLR4

## plot(physeqRBCLR4@phy_tree, show.tip.label = FALSE)



## ***********************************************************************

## read in rep-seqs to get DNA seq and feature id. right now tip labels are feature ids,
## we want to switch to actual DNA seqs

library(Biostrings)

# Replace with your actual file path
fasta_file <- "ffar_pipeline/merged/RBCL/rep-seqs-RBCL.fasta"

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



# feature.2.tax.rbcl <-
#     read.table("ffar_pipeline/merged/RBCL/taxonomyRBCL.txt", sep="\t",
#                header=TRUE)
# 
# ## ## R knows to ignore the row in the txt with a # (see defaults read.table)
# 
# ## add 16s/rbcl to make the columns easy to find when they are added
# ## to the specimen data
# 
# feature.2.tax.rbcl$Taxon  <- gsub(" ", "_", feature.2.tax.rbcl$Taxon)
# feature.2.tax.rbcl$Taxon <- paste("RBCL", feature.2.tax.rbcl$Taxon,
#                                   sep=':')

tree.rbclR0 <- phy_tree(physeqRBCLR0, errorIfNULL=TRUE)
tree.rbclR1 <- phy_tree(physeqRBCLR1, errorIfNULL=TRUE)
tree.rbclR2 <- phy_tree(physeqRBCLR2, errorIfNULL=TRUE)
tree.rbclR3 <- phy_tree(physeqRBCLR3, errorIfNULL=TRUE)
tree.rbclR4 <- phy_tree(physeqRBCLR4, errorIfNULL=TRUE)

## 
# First, filter fasta_df to only those in the tree
all_tips <- c(tree.rbclR0$tip.label,
              tree.rbclR1$tip.label,
              tree.rbclR2$tip.label,
              tree.rbclR3$tip.label,
              tree.rbclR4$tip.label)

fasta_df_filtered <- fasta_df[fasta_df$FeatureID %in% all_tips, ]
fasta_df_filteredR0 <- fasta_df[fasta_df$FeatureID %in% tree.rbclR0$tip.label, ]
fasta_df_filteredR1 <- fasta_df[fasta_df$FeatureID %in% tree.rbclR1$tip.label, ]
fasta_df_filteredR2 <- fasta_df[fasta_df$FeatureID %in% tree.rbclR2$tip.label, ]
fasta_df_filteredR3 <- fasta_df[fasta_df$FeatureID %in% tree.rbclR3$tip.label, ]
fasta_df_filteredR4 <- fasta_df[fasta_df$FeatureID %in% tree.rbclR4$tip.label, ]

## Then reassign the tip labels using match on the filtered set
## match the tip labs to the table with feature ID and Taxon
tree.rbclR0$tip.label <- fasta_df_filtered$Sequence[
  match(tree.rbclR0$tip.label, fasta_df_filtered$FeatureID)
]
tree.rbclR1$tip.label <- fasta_df_filtered$Sequence[
  match(tree.rbclR1$tip.label, fasta_df_filtered$FeatureID)
]
tree.rbclR2$tip.label <- fasta_df_filtered$Sequence[
  match(tree.rbclR2$tip.label, fasta_df_filtered$FeatureID)
]
tree.rbclR3$tip.label <- fasta_df_filtered$Sequence[
  match(tree.rbclR3$tip.label, fasta_df_filtered$FeatureID)
]
tree.rbclR4$tip.label <- fasta_df_filtered$Sequence[
  match(tree.rbclR4$tip.label, fasta_df_filtered$FeatureID)
]


## ***********************************************************************
## rbcl networks
## ***********************************************************************

#R0
indiv.comm.rbclR0 <-
  bipartite::empty(catchDups(makeComm(taxonomyRBCLR0,
                                      fasta_df_filteredR0,
                                      feature.col="FeatureID")))
# changed to "Feature.ID" to match 'makeComm' function
colnames(indiv.comm.rbclR0) <- fasta_df_filteredR0$Sequence
indiv.comm.rbclR0 <- indiv.comm.rbclR0/rowSums(indiv.comm.rbclR0)

#R1
indiv.comm.rbclR1 <-
  bipartite::empty(catchDups(makeComm(taxonomyRBCLR1,
                                      fasta_df_filteredR1,
                                      feature.col="FeatureID")))
# changed to "Feature.ID" to match 'makeComm' function
colnames(indiv.comm.rbclR1) <- fasta_df_filteredR1$Sequence
indiv.comm.rbclR1 <- indiv.comm.rbclR1/rowSums(indiv.comm.rbclR1)

#R2
indiv.comm.rbclR2 <-
  bipartite::empty(catchDups(makeComm(taxonomyRBCLR2,
                                      fasta_df_filteredR2,
                                      feature.col="FeatureID")))
# changed to "Feature.ID" to match 'makeComm' function
colnames(indiv.comm.rbclR2) <- fasta_df_filteredR2$Sequence
indiv.comm.rbclR2 <- indiv.comm.rbclR2/rowSums(indiv.comm.rbclR2)

#R3
indiv.comm.rbclR3 <-
  bipartite::empty(catchDups(makeComm(taxonomyRBCLR3,
                                      fasta_df_filteredR3,
                                      feature.col="FeatureID")))
# changed to "Feature.ID" to match 'makeComm' function
colnames(indiv.comm.rbclR3) <- fasta_df_filteredR3$Sequence
indiv.comm.rbclR3 <- indiv.comm.rbclR3/rowSums(indiv.comm.rbclR3)

#R4
indiv.comm.rbclR4 <-
  bipartite::empty(catchDups(makeComm(taxonomyRBCLR4,
                                      fasta_df_filteredR4,
                                      feature.col="FeatureID")))
# changed to "Feature.ID" to match 'makeComm' function
colnames(indiv.comm.rbclR4) <- fasta_df_filteredR4$Sequence
indiv.comm.rbclR4 <- indiv.comm.rbclR4/rowSums(indiv.comm.rbclR4)
###

bees.rbcl <- c(rownames(indiv.comm.rbclR0),
              paste0("SF",  rownames(indiv.comm.rbclR1)),
              paste0("SF", rownames(indiv.comm.rbclR2)),
                     paste0("SF", rownames(indiv.comm.rbclR3)),
              paste0("SF", rownames(indiv.comm.rbclR4)))

comms <- list(indiv.comm.rbclR0, indiv.comm.rbclR1,
              indiv.comm.rbclR2, indiv.comm.rbclR3,
              indiv.comm.rbclR4)

species.rbcl <- unique(unlist(sapply(comms, colnames)))

merged.comm.rbcl <- plyr::rbind.fill(lapply(comms, as.data.frame))

## check with number of columns against number of unique species
dim(merged.comm.rbcl)
length(species.rbcl)

rownames(merged.comm.rbcl) <- bees.rbcl


megaRBCLdata <- read_qza("ffar_pipeline/merged/RBCL/rooted-treeRBCL.qza")
tree.rbcl <- megaRBCLdata$data

## tree.rbcl <- tip_glom(tree.rbcl, h=0.1)

## tree.rbcl <- drop_dupes(tree.rbcl, thres=1e-5)

tree.rbcl$tip.label <- fasta_df_filtered$Sequence[
  match(tree.rbcl$tip.label, fasta_df_filtered$FeatureID)
]

## ***********************************************************************
## Make mega dataset
## ***********************************************************************
## spec already includes parasite data
ffar <- read.csv("../../Documents/1ROTATION kern lab/prepBeeRBCL/FFARspec.csv") %>%
  select(-starts_with("X16s"), -starts_with("RBCL"))


indiv.comm.rbcl <- as.data.frame(merged.comm.rbcl)
pollen <- colnames(indiv.comm.rbcl)
indiv.comm.rbcl$UniqueID <- rownames(indiv.comm.rbcl)

# indiv.comm.16s <- as.data.frame(merged.comm.16s)
# bact <- colnames(indiv.comm.16s)
# indiv.comm.16s$UniqueID <- rownames(indiv.comm.16s)
# 
# spec <-cbind(spec, indiv.comm.16s[, bact][match(spec$UniqueID,
#                            indiv.comm.16s$UniqueID),])

ffar.rbcl <-cbind(ffar, indiv.comm.rbcl[, pollen][match(ffar$UniqueID,
                                                        indiv.comm.rbcl$UniqueID),])

save(ffar.rbcl, file= "../pollenGeolocation/data/FFARpollen.Rdata")

write.csv(ffar.rbcl, file= "../pollenGeolocation/data/FFARpollen.csv",
          row.names=FALSE)

