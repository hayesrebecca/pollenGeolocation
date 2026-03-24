#!/usr/bin/env bash 

# Pipeline Code for Processing 16s and RBCL reads from fastq files. 

# 1: #1: Download docker on your computer, allows you to run qiime and
# other software without downloading them.

# On whatever computer will be running the tasks, innitiate a docker
# container for qiime1

#only need to do this docker pull step once per user
docker pull mbari/qiime1

#Helpful notes:
# 1. if running on the shared lab computer, you may get errors if your zipped sequence results are not synched 
#     to your dropbox account. Double clicking the file will usually make it appear. If you look at the filesize and 
# it is 0, this is likely because the file is not synched to the local.

# 2. If you get an error about memory or space issues, check that the other users do not have their dropbox files
#     on the local computer, you can fix this by logging in and making all dropbox files online only. This should
#     clear up enough space to fix this error.

## m1 mac specific flag: --platform linux/amd64 ???

#docker run -itv /Volumes/bombus/Dropbox\ \(University\ of\ Oregon\)/skyIslands_saved/SI_pipeline:/mnt/SI_pipeline mbari/qiime1
docker run -itv /Volumes/bombus/rhayes/Dropbox\ \(University\ of\ Oregon\)/pnw_survey_saved/pipeline:/mnt/pipeline mbari/qiime1
source activate qiime1 

# 2 for 2018 cd into the correct folder within the container, verify that the files you need are there

cd ../../mnt/pipeline/
ls

#3: If your reads are in *fasta.gz format, unzip them with qiime, then
#rename them "forward.fastq" and "reverse.fastq". Quinn designed the
#barcodes in reverse so flowcell1 is the reverse and flowcell2 is the
#forward

## Step 3 only need to be run once ever.

#cd R2024/raw
#gunzip *.gz

#mv GC3F-RH-11142---9462_S1_L001_R1_001.fastq rawreverse.fastq
#mv GC3F-RH-11142---9462_S1_L001_R2_001.fastq rawforward.fastq


#4a: parse the barcodes in the files, putting our data into a format
#qiime2 will be able to use. -- NOTE this step takes ~40 minutes to run!
cd R2024
extract_barcodes.py -f raw/rawforward.fastq -r raw/rawreverse.fastq  -c barcode_paired_end --bc1_len 8 --bc2_len 8 -o parsed_barcodes



#5a: re-zip the output files
cd parsed_barcodes
gzip *.fastq

mv reads1.fastq.gz forward.fastq.gz
mv reads2.fastq.gz reverse.fastq.gz


#6: Exit Qiime1 and use docker to open the environment for Qiime 2. 
exit
 
 
docker run -itv /Volumes/bombus/rhayes/Dropbox\ \(University\ of\ Oregon\)/pnw_survey_saved/pipeline:/mnt/pipeline qiime2/core
source activate qiime2

#6: Test that the container for Qiime 2 is properly associated, then
#make sure you are in the root directory using ls and/or pwd, then set working
#directory to the mounted volume

#7: Import your parser barcodes into an object you can demultiplex in
#Qiime 2. Make sure working directory is set correctly.

cd ../../
cd mnt/pipeline/

## Run 1 2023
# cd R2023/
# qiime tools import --type EMPPairedEndSequences --input-path parsed_barcodes/ --output-path seqs.qza

## Run 2 2024
cd R2024/
qiime tools import --type EMPPairedEndSequences --input-path parsed_barcodes/ --output-path seqs.qza


#8: Examine your mapping file, which is metadata you create associated
#with the project. Use qiime to make it into a qzv object.

#If you are examining multiple amplicon types, pick a map associated
#with one to start with (e.g. 16s)

# ## 2023 run 1
# cd ../../R2023/
# qiime metadata tabulate --m-input-file maps/pnw2023map16s.txt --o-visualization pnw2023map16s.qzv
# qiime tools view pnw2023map16s.qzv
# 
# 
# exit

## 2024 run 1
#cd ../../R2024/
qiime metadata tabulate --m-input-file maps/pnw2024map16s.txt --o-visualization pnw2024map16s.qzv
qiime tools view pnw2024map16s.qzv


exit


#9: Demultiplex 16s reads first. Only works in version Qiime2 2019.1
#note: this step takes ~2 hours!

docker run -itv /Volumes/bombus/rhayes/Dropbox\ \(University\ of\ Oregon\)/pnw_survey_saved/pipeline:/mnt/pipeline qiime2/core:2019.1

cd ../../mnt/pipeline

# Run 1
# cd R2023/
# qiime demux emp-paired --i-seqs seqs.qza --m-barcodes-file maps/pnw2023map16s.txt --m-barcodes-column barcodesequence --o-per-sample-sequences demux16s.qza 


cd R2024/



qiime demux emp-paired --i-seqs seqs.qza --m-barcodes-file maps/pnw2024map16s.txt --m-barcodes-column barcodesequence --o-per-sample-sequences demux16s.qza 


## visualize the results
qiime demux summarize --i-data demux16s.qza --o-visualization demux16s.qzv
qiime tools view demux16s.qzv

## click the tab Interactive Quality from the main visualization tab

## when interpretating the quality boxes, you can use the bottom of
## the black box as a conservative measure for the phred score (not the
## whiskers and not the middle of the box). Quality score around 30 is
## acceptable, lower not so good. Select the number before the
## quality degings to decline. 

## watch out for mnt

## The truncation length will vary between each run! Make sure to adjust the numbers pasted below. 
## R2023 16s: f = 173, r = 215
## R2024 16s: f = 206, r = 252

#note this step takes hours!

# ## Run 1 2023
# cd R2023/
# qiime dada2 denoise-paired  \
# --i-demultiplexed-seqs demux16s.qza  \
# --p-trunc-len-f 173  \
# --p-trunc-len-r 215  \
# --p-trim-left-f 0  \
# --p-n-threads 2  \
# --output-dir dada2-16s  \
#  --o-representative-sequences dada2-16s/rep-seqs-dada2-16s.qza  \
#  --o-table dada2-16s/table16s.qza
# 
# 
# qiime feature-table tabulate-seqs --i-data dada2-16s/rep-seqs-dada2-16s.qza --o-visualization dada2-16s/rep-seqs-dada2-16s.qzv
# qiime feature-table summarize --i-table dada2-16s/table16s.qza --o-visualization dada2-16s/table16s.qzv
# 
# 
# exit

## Run 2 2024

qiime dada2 denoise-paired  \
--i-demultiplexed-seqs demux16s.qza  \
--p-trunc-len-f 206  \
--p-trunc-len-r 252  \
--p-trim-left-f 0  \
--p-n-threads 2  \
--output-dir dada2-16s  \
 --o-representative-sequences dada2-16s/rep-seqs-dada2-16s.qza  \
 --o-table dada2-16s/table16s.qza


qiime feature-table tabulate-seqs --i-data dada2-16s/rep-seqs-dada2-16s.qza --o-visualization dada2-16s/rep-seqs-dada2-16s.qzv
qiime feature-table summarize --i-table dada2-16s/table16s.qza --o-visualization dada2-16s/table16s.qzv


exit

## repeat the steps for RBCL

#8: Examine your mapping file, which is metadata you create associated
#with the project. Use qiime to make it into a qzv object.


## 2023

docker run -itv /Volumes/bombus/rhayes/Dropbox\ \(University\ of\ Oregon\)/pnw_survey_saved/pipeline:/mnt/pipeline qiime2/core:2019.1

cd ../../mnt/pipeline
cd R2023
qiime metadata tabulate --m-input-file maps/pnw2023mapRBCL.txt --o-visualization pnw2023mapRBCL.qzv
qiime tools view pnw2023mapRBCL.qzv


exit

#9: Demultiplex RBCL reads. Only works in version Qiime2 2019.1
#note: this step takes ~2 hours!

docker run -itv /Volumes/bombus/rhayes/Dropbox\ \(University\ of\ Oregon\)/pnw_survey_saved/pipeline:/mnt/pipeline qiime2/core:2019.1

cd ../../mnt/pipeline

# Run 1
cd R2023/
qiime demux emp-paired --i-seqs seqs.qza --m-barcodes-file maps/pnw2023mapRBCL.txt --m-barcodes-column barcodesequence --o-per-sample-sequences demuxRBCL.qza 


## visualize the results
qiime demux summarize --i-data demuxRBCL.qza --o-visualization demuxRBCL.qzv
qiime tools view demuxRBCL.qzv

## click the tab Interactive Quality from the main visualization tab

## when interpretating the quality boxes, you can use the bottom of
## the black box as a conservative measure for the phred score (not the
## whiskers and not the middle of the box). Quality score around 30 is
## acceptable, lower not so good. Select the number before the
## quality degings to decline. 

## watch out for mnt

## The truncation length will vary between each run! Make sure to adjust the numbers pasted below. 
## R2023 RBCL: f = 109, r = 111

#note this step takes hours!

## Run 1 2023
# cd R2023/
# qiime dada2 denoise-paired  \
# --i-demultiplexed-seqs demuxRBCL.qza  \
# --p-trunc-len-f 109  \
# --p-trunc-len-r 111  \
# --p-trim-left-f 0  \
# --p-n-threads 2  \
# --output-dir dada2-RBCL  \
#  --o-representative-sequences dada2-RBCL/rep-seqs-dada2-RBCL.qza  \
#  --o-table dada2-RBCL/tableRBCL.qza

## no pollen for run 2024


# qiime feature-table tabulate-seqs --i-data dada2-RBCL/rep-seqs-dada2-RBCL.qza --o-visualization dada2-RBCL/rep-seqs-dada2-RBCL.qzv
# qiime feature-table summarize --i-table dada2-RBCL/tableRBCL.qza --o-visualization dada2-RBCL/tableRBCL.qzv


exit

# check outputs to make sure you didn't lose too many samples. 
# We found being more conservative and doing shorter truncations gives you the same number of sequences, 
# but sorted into fewer features, likely cause trimmed off poor quality reads

## *****************************************************************************
##       MERGE files from runs ONCE THEY EXIST
## *****************************************************************************
# time to merge the files from your different runs.
# NOTE: much of this comes from https://john-quensen.com/tutorials/merging-dada2-results-in-qiime2/ 

# cd back to your main folder (in this case pipeline) that has the separate run folders
# make a couple new directories

## cd back into before you do into R2023
cd ../../ 
mkdir merged
cd merged
mkdir 16s
mkdir RBCL
cd ../

#### 1. 16s ###
# 1a. first merge the table files. do this from your main SI_pipeline folder
qiime feature-table merge \
      --i-tables R2023/dada2-16s/table16s.qza \
      --i-tables R2024/dada2-16s/table16s.qza \
      --o-merged-table merged/16s/table16s.qza

# 1b: next merge the rep-seqs
qiime feature-table merge-seqs \
      --i-data R2023/dada2-16s/rep-seqs-dada2-16s.qza \
      --i-data R2024/dada2-16s/rep-seqs-dada2-16s.qza \
      --o-merged-data merged/16s/rep-seqs-16s.qza
