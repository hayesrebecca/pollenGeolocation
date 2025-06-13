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
setwd('pnw_survey_saved')



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
qza.rbcl.path  <- "pipeline/R2023/RBCL/final"


# 16s
#weightedUF16sqzaR0 <- read_qza(file.path(qza.16s.path,
#                      'core_metrics16s/weighted_unifrac_distance_matrix.qza'))

#unweightedUF16sqzaR0 <- read_qza(file.path(qza.16s.path,
#                      'core_metrics16s/unweighted_unifrac_distance_matrix.qza'))

# 16s
#wphylo.dist.16sR0 <- weightedUF16sqzaR0$data
#phylo.dist.16sR0 <-unweightedUF16sqzaR0$data

## Note, when taxonomy is imported, a single string is returned along
## with a confidence score.  For many analysis we will want to break
## up this string and for that purpose the parse_taxonomy() function
## is provided:

#taxonomy16sR0 <- read_qza(
#    file.path(qza.16s.path,
#              "core_metrics16s/rarefied_table.qza"))
#
#taxonomy16sR0 <- taxonomy16sR0$data

## ## 16s R0 phylogeny

#physeq16sR0 <- qza_to_phyloseq(
#    features=
#        file.path(qza.16s.path, "core_metrics16s/rarefied_table.qza"),
#    tree="pipeline/R2023/rooted-tree16s.qza",
#    "pipeline/R2023/taxonomy16s.qza",
#    metadata = "pipeline/R2023/maps/pnw2023map16s.txt"
#)



#physeq16sR0
#plot(physeq16sR0@phy_tree, show.tip.label = FALSE)

# To get taxonomy.tsv, convert taxonomy16s.qza to a zip file in terminal
# cp taxonomy16s.qza taxonomy16s.zip
# unzip taxonomy16s.zip
# then open the resulting folder and cd into data then copy metadata.tsv to taxonomy16s.tsv in terminal
# cp metadata.tsv ../../taxonomy16s.tsv

#feature.2.tax.16s <-
#    read.table("pipeline/R2023/taxonomy16s.tsv", sep="\t",
#               header=TRUE)

#feature.2.tax.16s$Taxon <- paste("16s", feature.2.tax.16s$Taxon, sep=':')

## convert to a phylo class which is more useful downstream
#tree.16sR0 <- phy_tree(physeq16sR0, errorIfNULL=TRUE)

## match the tip labs to the table with feature ID and Taxon
#tree.16sR0$tip.label  <-  feature.2.tax.16s$Taxon[match(tree.16sR0$tip.label,
#                                           feature.2.tax.16s$Feature.ID)]
#

## 10-24-2023 Rebecca is dropping all sequences that are only resolved to the first level D_0__Bacteria

# Identify tips with labels exactly matching '16s:D_0__Bacteria'
#matching_tips <- grep('^16s:d__Bacteria$', tree.16sR0$tip.label)

# Drop the matching tips
#tree.16sR0 <- drop.tip(tree.16sR0, matching_tips)

# # Drop the tips that are NA
#tree.16sR0 <- drop.tip(tree.16sR0, tree.16sR0$tip.label[is.na(tree.16sR0$tip.label)])

#drop unassigned tips
#unassigned_tips <- grep('^16s:Unassigned$', tree.16sR0$tip.label)

# Drop the tips that are NA
#tree.16sR0 <- drop.tip(tree.16sR0, tree.16sR0$tip.label[is.na(tree.16sR0$tip.label)])

#plot(tree.16sR0, show.tip.label = FALSE)
#save(physeq16sR0, tree.16sR0, file= "../pnw_survey/data/physeq16s.Rdata")
## ***********************************************************************
## 16s networks
## ***********************************************************************

#finalASVtable <- read.csv(file.path(qza.16s.path, "asvTable2023_final.csv"), header=TRUE)

#rownames(finalASVtable) <- finalASVtable[,1]

#finalASVtable[,1] <- NULL

## drop the barcode cols
#finalASVtable <- finalASVtable %>%
#  dplyr::select(starts_with("d__"))

#fixing naming inconsistency
#colnames(finalASVtable) <- paste("16s", colnames(finalASVtable), sep=':')

#colnames(finalASVtable) <- gsub("\\.__", "", colnames(finalASVtable))

#colnames(finalASVtable) <- gsub("\\.", "; ", colnames(finalASVtable))

#indiv.comm.16sR0 <-
#    bipartite::empty(catchDups(makeComm(taxonomy16sR0,
#                                        feature.2.tax.16s)))
## saving this out so we can make the tree visualizations
#save(indiv.comm.16sR0, file="../pnw_survey/data/indiv.comm16sR0.Rdata")

#finalASVtable <- finalASVtable/rowSums(finalASVtable)

#finalASVtable <- as.matrix(finalASVtable)

#bees.16s <- c(rownames(finalASVtable))

#comms <- list(finalASVtable)

#species.16s <- unique(unlist(sapply(comms, colnames)))

#merged.comm.16s <- plyr::rbind.fill(lapply(comms, as.data.frame))

## check with number of columns against number of unique species
#dim(merged.comm.16s)
#length(species.16s)

#rownames(merged.comm.16s) <- bees.16s


#length(tree.16sR0$tip.label[tree.16sR0$tip.label %in% colnames(merged.comm.16s)])
#labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
#print(length(labels_mismatch)) #337


# # fixing g__IS-44 issue
# tree.16sR0$tip.label[grep("16s:d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Nitrosomonadaceae;", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("16s:d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Nitrosomonadaceae;", colnames(merged.comm.16s))]
# 
# colnames(merged.comm.16s) <- gsub("g__IS; 44;", "g__IS-44;", colnames(merged.comm.16s))
# colnames(merged.comm.16s)[grep("g__IS-44;", colnames(merged.comm.16s))]
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #336
# 
# # fixing g__Burkholderia-Caballeronia-Paraburkholderia. issue
# tree.16sR0$tip.label[grep("g__Burkholderia-Caballeronia-Paraburkholderia", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("g__Burkholderia; Caballeronia", colnames(merged.comm.16s))]
# 
# colnames(merged.comm.16s) <- gsub("g__Burkholderia; Caballeronia; Paraburkholderia", "g__Burkholderia-Caballeronia-Paraburkholderia", colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #308
# 
# # fixing f__TRA3 issue
# tree.16sR0$tip.label[grep("f__TRA3", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("f__TRA3", colnames(merged.comm.16s))]
# 
# 
# colnames(merged.comm.16s) <- gsub("f__TRA3; 20; g__TRA3; 20", "f__TRA3-20; g__TRA3-20", colnames(merged.comm.16s))
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #304
# 
# # fixing f__SC issue
# tree.16sR0$tip.label[grep("f__SC", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("f__SC", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("f__SC; I; 84; g__SC; I; 84", "f__SC-I-84; g__SC-I-84", colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #297
# 
# # fixing s__Burkholderia_sp. issue
# tree.16sR0$tip.label[grep("s__Burkholderia_sp.", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("s__Burkholderia_sp.", colnames(merged.comm.16s))]
# 
# colnames(merged.comm.16s) <- gsub("s__Burkholderia_sp; ",
#                                 "s__Burkholderia_sp.",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #296
# 
# # fixing s__Lampropedia_sp. issue
# tree.16sR0$tip.label[grep("s__Lampropedia_sp.", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("s__Lampropedia_sp.", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("s__Lampropedia_sp; ",
#                                 "s__Lampropedia_sp.",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #295
# 
# # fixing s__Bdellovibrio_sp. issue
# tree.16sR0$tip.label[grep("s__Bdellovibrio_sp.", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("s__Bdellovibrio_sp; ", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("s__Bdellovibrio_sp; ",
#                                 "s__Bdellovibrio_sp.",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #294
# 
# # fixing s__Chlamydia_sp. issue
# tree.16sR0$tip.label[grep("s__Chlamydia_sp.", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("s__Chlamydia_sp; ", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("s__Chlamydia_sp; ",
#                                 "s__Chlamydia_sp.",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #293
# 
# # fixing s__Neochlamydia_sp. issue
# tree.16sR0$tip.label[grep("s__Neochlamydia_sp.", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("s__Neochlamydia_sp.", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("s__Neochlamydia_sp; ",
#                                 "s__Neochlamydia_sp.",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #288
# 
# # fixing terminal _sp. issue
# tree.16sR0$tip.label[grep("_sp.", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("_sp.", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("_sp; ",
#                                 "_sp.",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #258
# 
# # fixing g__Amb issue
# tree.16sR0$tip.label[grep("g__Amb", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("g__Amb", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__Amb; 16S; 1323",
#                                 "__Amb-16S-1323",
#                                 colnames(merged.comm.16s))
# 
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #257
# 
# # fixing g__Allorhizobium issue
# tree.16sR0$tip.label[grep("g__Allorhizobium", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("g__Allorhizobium", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("g__Allorhizobium; Neorhizobium; Pararhizobium; Rhizobium",
#                                 "g__Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #203
# 
# # fixing g__Methylobacterium issue
# tree.16sR0$tip.label[grep("g__Methylobacterium", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("g__Methylobacterium", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__Methylobacterium; Methylorubrum",
#                                 "__Methylobacterium-Methylorubrum",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #180
# 
# # fixing g__1174 issue
# tree.16sR0$tip.label[grep("g__1174", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("g__1174", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__1174; 901; 12",
#                                 "__1174-901-12",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #163
# 
# # fixing g__Hafnia issue
# tree.16sR0$tip.label[grep("g__Hafnia", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("g__Hafnia", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__Hafnia; Obesumbacterium",
#                                 "__Hafnia-Obesumbacterium",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #158
# 
# # fixing g__wb1 issue
# tree.16sR0$tip.label[grep("g__wb1", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("g__wb1", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__wb1; P19",
#                                 "__wb1-P19",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #157
# 
# # fixing coprostanoligenes issue
# tree.16sR0$tip.label[grep("coprostanoligenes", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("coprostanoligenes", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("f__; Eubacterium; _coprostanoligenes_group; g__; Eubacterium; _coprostanoligenes_group",
#                                 "f__[Eubacterium]_coprostanoligenes_group; g__[Eubacterium]_coprostanoligenes_group",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #156
# 
# # fixing __Clostridia issue
# tree.16sR0$tip.label[grep("__Clostridia_UCG", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__Clostridia_UCG", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__Clostridia_UCG; 014",
#                                 "__Clostridia_UCG-014",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #155
# 
# # fixing __Peptostreptococcales issue
# tree.16sR0$tip.label[grep("__Peptostreptococcales", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__Peptostreptococcales", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__Peptostreptococcales; Tissierellales",
#                                 "__Peptostreptococcales-Tissierellales",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #134
# 
# # fixing _coagulans issue
# tree.16sR0$tip.label[grep("_coagulans", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("_coagulans", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("s__; Bacteroides; _coagulans",
#                                 "s__[Bacteroides]_coagulans",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #133
# 
# # fixing g__Christensenellaceae issue
# tree.16sR0$tip.label[grep("g__Christensenellaceae", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("g__Christensenellaceae", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("g__Christensenellaceae_R; 7_group",
#                                 "g__Christensenellaceae_R-7_group",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #132
# 
# # fixing ic1294 issue
# tree.16sR0$tip.label[grep("ic1294", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("ic1294", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__; Ruminococcus; _torques_group",
#                                 "__[Ruminococcus]_torques_group",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #128
# 
# # fixing __KF issue
# tree.16sR0$tip.label[grep("__KF", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__KF", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("f__KF; JG30; B3; g__KF; JG30; B3",
#                                 "f__KF-JG30-B3; g__KF-JG30-B3",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #127
# 
# # fixing __Gitt issue
# tree.16sR0$tip.label[grep("__Gitt", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__Gitt", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("c__Gitt; GS; 136; o__Gitt; GS; 136; f__Gitt; GS; 136; g__Gitt; GS; 136",
#                                 "c__Gitt-GS-136; o__Gitt-GS-136; f__Gitt-GS-136; g__Gitt-GS-136",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #126
# 
# # fixing __Gitt issue
# tree.16sR0$tip.label[grep("__Gitt", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__Gitt", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("c__Gitt; GS; 136; o__Gitt; GS; 136; f__Gitt; GS; 136; g__Gitt; GS; 136",
#                                 "c__Gitt-GS-136; o__Gitt-GS-136; f__Gitt-GS-136; g__Gitt-GS-136",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #126
# 
# # fixing __CL500 issue
# tree.16sR0$tip.label[grep("__CL500", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__CL500", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("g__CL500; 29_marine_group",
#                                 "g__CL500-29_marine_group",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #121
# 
# # fixing __0319 issue
# tree.16sR0$tip.label[grep("__0319", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__0319", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__0319; 7L14",
#                                 "__0319-7L14",
#                                 colnames(merged.comm.16s))
# colnames(merged.comm.16s) <- gsub("__0319; 6G20",
#                                 "__0319-6G20",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #106
# 
# # fixing __ADurb issue
# tree.16sR0$tip.label[grep("__ADurb", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__ADurb", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("g__ADurb; Bin063; 1",
#                                 "g__ADurb.Bin063-1",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #105
# 
# # fixing __Elev issue
# tree.16sR0$tip.label[grep("__Elev", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__Elev", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__Elev; 16S; 1166",
#                                 "__Elev-16S-1166",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #104
# 
# # fixing __RCP2 issue
# tree.16sR0$tip.label[grep("__RCP2", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__RCP2", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("p__RCP2; 54; c__RCP2; 54; o__RCP2; 54; f__RCP2; 54; g__RCP2; 54",
#                                 "p__RCP2-54; c__RCP2-54; o__RCP2-54; f__RCP2-54; g__RCP2-54",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #103
# 
# # fixing __mle1 issue
# tree.16sR0$tip.label[grep("__mle1", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__mle1", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__mle1; 27; f__mle1; 27; g__mle1; 27",
#                                 "__mle1-27; f__mle1-27; g__mle1-27",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #99
# 
# # fixing __P3OB issue
# tree.16sR0$tip.label[grep("__P3OB", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__P3OB", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("g__P3OB; 42",
#                                 "g__P3OB-42",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #91
# 
# # fixing __YC issue
# tree.16sR0$tip.label[grep("__YC", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__YC", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("g__YC; ZSS; LKJ147",
#                                 "g__YC-ZSS-LKJ147",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #90
# 
# # fixing __NS11 issue
# tree.16sR0$tip.label[grep("__NS11", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__NS11", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__NS11; 12_marine_group",
#                                 "__NS11-12_marine_group",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #86
# 
# # fixing __LiUU issue
# tree.16sR0$tip.label[grep("__LiUU", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__LiUU", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__LiUU; 11; 161",
#                                 "__LiUU-11-161",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #85
# 
# # fixing __env issue
# tree.16sR0$tip.label[grep("__env", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__env", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__env; OPS_17",
#                                 "__env.OPS_17",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #76
# 
# # fixing __KD3 issue
# tree.16sR0$tip.label[grep("__KD3", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__KD3", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__KD3; 93",
#                                 "__KD3-93",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #74
# 
# # fixing __37 issue
# tree.16sR0$tip.label[grep("__37", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__37", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__37; 13",
#                                 "__37-13",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #71
# 
# # fixing VC2 issue
# tree.16sR0$tip.label[grep("VC2", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("VC2", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__Bacteroidetes_VC2; 1_Bac22; f__Bacteroidetes_VC2; 1_Bac22; g__Bacteroidetes_VC2; 1_Bac22",
#                                 "__Bacteroidetes_VC2.1_Bac22; f__Bacteroidetes_VC2.1_Bac22; g__Bacteroidetes_VC2.1_Bac22",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #68
# 
# # fixing __cf issue
# tree.16sR0$tip.label[grep("__cf", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("__cf", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__cf; _Chryseobacterium",
#                                 "__cf._Chryseobacterium",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #67
# 
# # fixing Absconditabacteriales issue
# tree.16sR0$tip.label[grep("Absconditabacteriales", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("Absconditabacteriales", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("__Absconditabacteriales_; SR1",
#                                 "__Absconditabacteriales_(SR1)",
#                                 colnames(merged.comm.16s))
# colnames(merged.comm.16s) <- gsub("; ;",
#                                 ";",
#                                 colnames(merged.comm.16s))
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #66
# 
# # fixing Subgroup_1) issue
# tree.16sR0$tip.label[grep("Subgroup_1", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("Subgroup_1", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("f__Acidobacteriaceae_; Subgroup_1;",
#                                 "f__Acidobacteriaceae_(Subgroup_1);",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #23
# 
# # fixing Subgroup_1) issue
# tree.16sR0$tip.label[grep("f__Acidobacteriaceae", tree.16sR0$tip.label)]
# colnames(merged.comm.16s)[grep("f__Acidobacteriaceae", colnames(merged.comm.16s))]
# colnames(merged.comm.16s) <- gsub("; $",
#                                 "",
#                                 colnames(merged.comm.16s))
# 
# # check labels again
# labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(merged.comm.16s))]
# print(length(labels_mismatch)) #0 YAY!!
# 
# colnames(merged.comm.16s) <- gsub("; ;",
#                                 ";",
#                                 colnames(merged.comm.16s))
# 
# save(merged.comm.16s, file= "../pnw_survey/data/presAbsTable.Rdata")
# 


## ***********************************************************************
## working with a merged 16s tree
## ***********************************************************************

# upload our mega 16s phylogenetic tree

#mega16sdata <- read_qza("pipeline/R2023/rooted-tree16s.qza")

#tree.16s <- mega16sdata$data

## function from: https://stackoverflow.com/questions/38570074/phylogenetics-in-r-collapsing-descendant-tips-of-an-internal-node

drop_dupes <- function(tree, thres=1e-5){
  tips <- which(tree$edge[,2] %in% 1:Ntip(tree))
  toDrop <- tree$edge.length[tips] < thres
  drop.tip(tree,tree$tip.label[toDrop])
}

## tree.16s <- tip_glom(tree.16s, h=0.05)

## what is a good cutoff?
#tree.16s <- drop_dupes(tree.16s, thres=1e-3)

# make distance matrix

#tree.16s$tip.label  <-  feature.2.tax.16s$Taxon[match(tree.16s$tip.label,
#                                           feature.2.tax.16s$Feature.ID)]

#tree.16s <- drop.tip(tree.16s, which(duplicated(tree.16s$tip.label)))


## check for duplicates
#any(duplicated(spec.net$SpecimenID))

## ***********************************************************************
## RBCL
## ***********************************************************************
## R0
weightedUFrbclqzaR0 <- read_qza(file.path(qza.rbcl.path,
                                          'core_metricsRBCL/weighted_unifrac_distance_matrix.qza'))

unweightedUFrbclqzaR0 <- read_qza(file.path(qza.rbcl.path,
                                            'core_metricsRBCL/unweighted_unifrac_distance_matrix.qza'))
# 
# ## access data inside artifacts
# 
# # RBCL
# # R0
wphylo.dist.rbclR0 <- weightedUFrbclqzaR0$data
phylo.dist.rbclR0 <- unweightedUFrbclqzaR0$data

# 
# #Note, when taxonomy is imported, a single string is returned along
# #with a confidence score.  For many analysis we will want to break up
# #this string and for that purpose the parse_taxonomy() function is
# #provided:
# 
# # R0
taxonomyRBCLR0 <-
  read_qza("pipeline/R2023/RBCL/final/core_metricsRBCL/rarefied_table.qza")
taxonomyRBCLR0 <- taxonomyRBCLR0$data

# 
# #Phyloseq tutorials here: https://joey711.github.io/phyloseq/
# 
# # ****CHECK PATHS BELOW ONCE FINISH RBCL PIPELINE****
# 
# ## NEED TO CHECK whether this actually has all of the samples and tips
# ## TODO: I think I need to rerun the pipeline steps to make sure the correct qzas are used
# 
# #RBCL 2023 phylogeny
physeqRBCLR0 <- qza_to_phyloseq(
  features="pipeline/R2023/RBCL/final/core_metricsRBCL/rarefied_table.qza",
  tree="pipeline/R2023/RBCL/rooted-treeRBCL2023.qza",
  "pipeline/R2023/RBCL/taxonomyRBCL.qza",
  metadata = "pipeline/R2023/maps/pnw2023mapRBCL.txt"
)
physeqRBCLR0
plot(physeqRBCLR0@phy_tree, show.tip.label = FALSE)
# 

# ## ***********************************************************************
# 
#feature.2.tax.rbcl <-
#    read.table("pipeline/R2023/RBCL/taxonomyRBCL.txt", sep="\t",
#               header=TRUE)
# 
# ## ## R knows to ignore the row in the txt with a # (see defaults read.table)
# 

## read in rep-seqs to get DNA seq and feature id. right now tip labels are feature ids,
## we want to switch to actual DNA seqs

#feature.2.seq.rbcl <- read_qza("pipeline/R2023/dada2-RBCL/rep-seqs-dada2-RBCL.qza")$data


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



#16s
#indiv.comm.16s <- as.data.frame(merged.comm.16s)
#bact <- colnames(indiv.comm.16s)
#indiv.comm.16s$SampleID <- rownames(indiv.comm.16s)
#indiv.comm.16s$SampleID <- sub("^", "NCASI-S", indiv.comm.16s$SampleID)
  

## Seems to be working up to here

#spec.net <-cbind(spec.net, indiv.comm.16s[, bact][match(spec.net$SampleID,
#                           indiv.comm.16s$SampleID),])

spec.net <-cbind(spec.net, indiv.comm.rbcl[, pollen][match(spec.net$SampleID,
                           indiv.comm.rbcl$SampleID),])


## check for duplicates
#any(duplicated(spec.net$UniqueID))

save(spec.net, file= "../pollenGeolocation/data/NCASIpollen.Rdata")

write.csv(spec.net, file= "../pollenGeolocation/data/NCASIpollen.csv",
          row.names=FALSE)



