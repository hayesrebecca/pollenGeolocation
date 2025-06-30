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
setwd('skyIslands_saved')

## pollen_type = 'raw' if you want raw dna sequences as features
## pollen_type = 'taxonomy' if you want the classified dna seqs

pollen_type = 'taxonomy'

#if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
##BiocManager::install("TreeSummarizedExperiment")
library(TreeSummarizedExperiment)
# if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
# devtools::install_github("jbisanz/qiime2R")
library(qiime2R)

library(tidyr)
library(dplyr)
library(bipartite)
library(phyloseq)
library(TreeTools)
#library(devtools)
library(ape)
library(picante)

source("../skyIslands/dataPrep/src/misc.R")

# 2018 Samples!!
## reading artifacts
#qza.16s.path  <- "SI_pipeline/merged/16s/final"
qza.rbcl.path  <- "SI_pipeline/merged/RBCL/final"


## ***********************************************************************
## RBCL
## ***********************************************************************
## R0
weightedUFrbclqzaR0 <- read_qza(file.path(qza.rbcl.path,
                                          'core_metricsRBCL-2025/weighted_unifrac_distance_matrix.qza'))

unweightedUFrbclqzaR0 <- read_qza(file.path(qza.rbcl.path,
                                            'core_metricsRBCL-2025/unweighted_unifrac_distance_matrix.qza'))
# 
# ## access data inside artifacts
# 
# # RBCL
# # R0
wphylo.dist.rbclR0 <- weightedUFrbclqzaR0$data
phylo.dist.rbclR0 <- unweightedUFrbclqzaR0$data
# 
# 
# #Note, when taxonomy is imported, a single string is returned along
# #with a confidence score.  For many analysis we will want to break up
# #this string and for that purpose the parse_taxonomy() function is
# #provided:
# 
# # R0
taxonomyRBCLR0 <-
  read_qza("SI_pipeline/merged/RBCL/final/core_metricsRBCL-2025/rarefied_table.qza")
taxonomyRBCLR0 <- taxonomyRBCLR0$data
# 
# 

# 
# 
# #Phyloseq tutorials here: https://joey711.github.io/phyloseq/
# 
# # ****CHECK PATHS BELOW ONCE FINISH RBCL PIPELINE****
# 
# ## NEED TO CHECK whether this actually has all of the samples and tips
# ## TODO: I think I need to rerun the pipeline steps to make sure the correct qzas are used
# 
#RBCL 2023 phylogeny
physeqRBCLR0 <- qza_to_phyloseq(
  features="SI_pipeline/merged/RBCL/final/core_metricsRBCL-2025/rarefied_table.qza",
  tree="SI_pipeline/merged/RBCL/rooted-treeRBCL.qza",
  "SI_pipeline/merged/RBCL/taxonomyRBCL-2025.qza",
  metadata = "SI_pipeline/merged/RBCL/maps/skyMap2018-2021RBCL.txt"
)
physeqRBCLR0
plot(physeqRBCLR0@phy_tree, show.tip.label = FALSE)
# 

# ## ***********************************************************************
## read in rep-seqs to get DNA seq and feature id. right now tip labels are feature ids,
## we want to switch to actual DNA seqs

library(Biostrings)

if(pollen_type=="raw"){
  # Replace with your actual file path
  fasta_file <- "SI_pipeline/merged/RBCL/rep-seqs-RBCL.fasta"
  
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
  # tree.rbclR1 <- phy_tree(physeqRBCLR1, errorIfNULL=TRUE)
  # tree.rbclR2 <- phy_tree(physeqRBCLR2, errorIfNULL=TRUE)
  # tree.rbclR3 <- phy_tree(physeqRBCLR3, errorIfNULL=TRUE)
  # tree.rbclR4 <- phy_tree(physeqRBCLR4, errorIfNULL=TRUE)
  
  ## 
  # First, filter fasta_df to only those in the tree
  # all_tips <- c(tree.rbclR0$tip.label,
  #               tree.rbclR1$tip.label,
  #               tree.rbclR2$tip.label,
  #               tree.rbclR3$tip.label,
  #               tree.rbclR4$tip.label)
  
  #fasta_df_filtered <- fasta_df[fasta_df$FeatureID %in% all_tips, ]
  fasta_df_filteredR0 <- fasta_df[fasta_df$FeatureID %in% tree.rbclR0$tip.label, ]
  # fasta_df_filteredR1 <- fasta_df[fasta_df$FeatureID %in% tree.rbclR1$tip.label, ]
  # fasta_df_filteredR2 <- fasta_df[fasta_df$FeatureID %in% tree.rbclR2$tip.label, ]
  # fasta_df_filteredR3 <- fasta_df[fasta_df$FeatureID %in% tree.rbclR3$tip.label, ]
  # fasta_df_filteredR4 <- fasta_df[fasta_df$FeatureID %in% tree.rbclR4$tip.label, ]
  
  ## Then reassign the tip labels using match on the filtered set
  ## match the tip labs to the table with feature ID and Taxon
  tree.rbclR0$tip.label <- fasta_df_filteredR0$Sequence[
    match(tree.rbclR0$tip.label, fasta_df_filteredR0$FeatureID)
  ]
  # tree.rbclR1$tip.label <- fasta_df_filtered$Sequence[
  #   match(tree.rbclR1$tip.label, fasta_df_filtered$FeatureID)
  # ]
  # tree.rbclR2$tip.label <- fasta_df_filtered$Sequence[
  #   match(tree.rbclR2$tip.label, fasta_df_filtered$FeatureID)
  # ]
  # tree.rbclR3$tip.label <- fasta_df_filtered$Sequence[
  #   match(tree.rbclR3$tip.label, fasta_df_filtered$FeatureID)
  # ]
  # tree.rbclR4$tip.label <- fasta_df_filtered$Sequence[
  #   match(tree.rbclR4$tip.label, fasta_df_filtered$FeatureID)
  # ]
  
  
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
  
  # #R1
  # indiv.comm.rbclR1 <-
  #   bipartite::empty(catchDups(makeComm(taxonomyRBCLR1,
  #                                       fasta_df_filteredR1,
  #                                       feature.col="FeatureID")))
  # # changed to "Feature.ID" to match 'makeComm' function
  # colnames(indiv.comm.rbclR1) <- fasta_df_filteredR1$Sequence
  # indiv.comm.rbclR1 <- indiv.comm.rbclR1/rowSums(indiv.comm.rbclR1)
  # 
  # #R2
  # indiv.comm.rbclR2 <-
  #   bipartite::empty(catchDups(makeComm(taxonomyRBCLR2,
  #                                       fasta_df_filteredR2,
  #                                       feature.col="FeatureID")))
  # # changed to "Feature.ID" to match 'makeComm' function
  # colnames(indiv.comm.rbclR2) <- fasta_df_filteredR2$Sequence
  # indiv.comm.rbclR2 <- indiv.comm.rbclR2/rowSums(indiv.comm.rbclR2)
  # 
  # #R3
  # indiv.comm.rbclR3 <-
  #   bipartite::empty(catchDups(makeComm(taxonomyRBCLR3,
  #                                       fasta_df_filteredR3,
  #                                       feature.col="FeatureID")))
  # # changed to "Feature.ID" to match 'makeComm' function
  # colnames(indiv.comm.rbclR3) <- fasta_df_filteredR3$Sequence
  # indiv.comm.rbclR3 <- indiv.comm.rbclR3/rowSums(indiv.comm.rbclR3)
  # 
  # #R4
  # indiv.comm.rbclR4 <-
  #   bipartite::empty(catchDups(makeComm(taxonomyRBCLR4,
  #                                       fasta_df_filteredR4,
  #                                       feature.col="FeatureID")))
  # # changed to "Feature.ID" to match 'makeComm' function
  # colnames(indiv.comm.rbclR4) <- fasta_df_filteredR4$Sequence
  # indiv.comm.rbclR4 <- indiv.comm.rbclR4/rowSums(indiv.comm.rbclR4)
  # ###
  
  bees.rbcl <- c(rownames(indiv.comm.rbclR0))
  
  comms <- list(indiv.comm.rbclR0)
  
  species.rbcl <- unique(unlist(sapply(comms, colnames)))
  
  merged.comm.rbcl <- plyr::rbind.fill(lapply(comms, as.data.frame))
  
  ## check with number of columns against number of unique species
  dim(merged.comm.rbcl)
  length(species.rbcl)
  
  rownames(merged.comm.rbcl) <- bees.rbcl
  
  
  megaRBCLdata <- read_qza("SI_pipeline/merged/RBCL/rooted-treeRBCL.qza")
  tree.rbcl <- megaRBCLdata$data
  
  tree.rbcl <- tip_glom(tree.rbcl, h=0.1)
  
  tree.rbcl <- drop_dupes(tree.rbcl, thres=1e-5)
  
  tree.rbcl$tip.label <- fasta_df_filteredR0$Sequence[
    match(tree.rbcl$tip.label, fasta_df_filteredR0$FeatureID)
  ]
  
  ## ***********************************************************************
  ## Make mega dataset
  ## ***********************************************************************
  
  si <- read.csv("../skyIslands/data/spec_RBCL_16s.csv") %>%
    select(-starts_with("X16s"), -starts_with("RBCL"))
  
  
  indiv.comm.rbcl <- as.data.frame(merged.comm.rbcl)
  pollen <- colnames(indiv.comm.rbcl)
  indiv.comm.rbcl$UniqueID <- rownames(indiv.comm.rbcl)
  
  
  si.rbcl <-cbind(si, indiv.comm.rbcl[, pollen][match(si$UniqueID,
                                                          indiv.comm.rbcl$UniqueID),])
  
  save(si.rbcl, file= "../pollenGeolocation/data/SIpollen_raw.Rdata")
  
  write.csv(si.rbcl, file= "../pollenGeolocation/data/SIpollen_raw.csv",
            row.names=FALSE)
} else {
  # 
  feature.2.tax.rbcl <-
      read.table("SI_pipeline/merged/RBCL/taxonomyRBCL-2025.txt", sep="\t",
                 header=TRUE)
  # 
  # ## ## R knows to ignore the row in the txt with a # (see defaults read.table)
  # 
  # ## add 16s/rbcl to make the columns easy to find when they are added
  # ## to the specimen data
  # 
  feature.2.tax.rbcl$Taxon  <- gsub(" ", "_", feature.2.tax.rbcl$Taxon)
  feature.2.tax.rbcl$Taxon <- paste("RBCL", feature.2.tax.rbcl$Taxon,
                                    sep=':')
  ## R0
  
  tree.rbclR0 <- phy_tree(physeqRBCLR0, errorIfNULL=TRUE)

  # ## match the tip labs to the table with feature ID and Taxon
  ## R0
  tree.rbclR0$tip.label  <-  feature.2.tax.rbcl$Taxon[
                             match(tree.rbclR0$tip.label,
                             feature.2.tax.rbcl$Feature.ID)]
  
  # ## ***********************************************************************
  # ## rbcl networks
  # ## ***********************************************************************
  # 
  #R0
  indiv.comm.rbclR0 <-
      bipartite::empty(catchDups(makeComm(taxonomyRBCLR0,
                                          feature.2.tax.rbcl,
                                          feature.col="Feature.ID")))
                                          # changed to "Feature.ID" to match 'makeComm' function

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
  megaRBCLdata <- read_qza("SI_pipeline/merged/RBCL/rooted-treeRBCL.qza")
  tree.rbcl <- megaRBCLdata$data
  # 
  # tree.rbcl <- tip_glom(tree.rbcl, h=0.1)
  # 
  # ## tree.rbcl <- drop_dupes(tree.rbcl, thres=1e-5)
  # 
  tree.rbcl$tip.label  <-  feature.2.tax.rbcl$Taxon[match(tree.rbcl$tip.label,
                                              feature.2.tax.rbcl$Feature.ID)]
  
  
  ## ***********************************************************************
  ## Make mega dataset
  ## ***********************************************************************
  
  ## spec already includes parasite data
  load('../skyIslands/data/spec_net_fdiv.Rdata')
  
  indiv.comm.rbcl <- as.data.frame(merged.comm.rbcl)
  pollen <- colnames(indiv.comm.rbcl)
  indiv.comm.rbcl$UniqueID <- rownames(indiv.comm.rbcl)
  
  
  si.pollen <-cbind(spec.net, indiv.comm.rbcl[, pollen][match(spec.net$UniqueID,
                             indiv.comm.rbcl$UniqueID),])
  
  si.pollen$`RBCL:BACTERIA` <- NULL
  
  save(si.pollen, file= "../pollenGeolocation/data/SIpollen_tax.Rdata")
  
  write.csv(si.pollen, file= "../pollenGeolocation/data/SIpollen_tax.csv",
            row.names=FALSE)
  
}
