## ---- include = FALSE---------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)


## ----setup, results = "hide"--------------------------------------------------------------------------------------------------------------------
# Attach the packages you will need for the analysis.
library(blaseRtools)
library(blaseRdata)
library(GenomicRanges)
library(tidyverse)


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Read in the data

vignette_CXCL8_ape <- bb_parseape(system.file("extdata/hg38_CXCL8.ape", package = "blaseRdata"))



## -----------------------------------------------------------------------------------------------------------------------------------------------
# Show the Ape Data 

vignette_CXCL8_ape



## -----------------------------------------------------------------------------------------------------------------------------------------------

vignette_CXCL1_ape <- bb_hg38_ape("CXCL1")

vignette_CXCL1_ape



## -----------------------------------------------------------------------------------------------------------------------------------------------
# Get genomic sequence and extend 100 bp left and right.
vignette_CXCL1_ape <- bb_hg38_ape("CXCL1", extend_left = 100, extend_right = 100)

vignette_CXCL1_ape



## -----------------------------------------------------------------------------------------------------------------------------------------------
# Get genomic sequence and extend 100 bp left and right.
# Now add a new custom feature based on original coordinates:  chr4 73869293-73871408
vignette_CXCL1_ape <- bb_hg38_ape("CXCL1", 
                                  extend_left = 100, 
                                  extend_right = 100, 
                                  additional_granges = GenomicRanges::makeGRangesFromDataFrame(data.frame(
                                    seqname = "chr4", 
                                    start = 73869293,
                                    end = 73871408,
                                    strand = "+",
                                    gene_name = "CXCL1",
                                    type = "custom_feature",
                                    label = "custom_feature_1",
                                    fwdcolor = "red",
                                    revcolor = "blue"
                                  ), keep.extra.columns = T))

vignette_CXCL1_ape



## ----eval = FALSE-------------------------------------------------------------------------------------------------------------------------------
## # Save as a genbank/Ape file
## Ape.save(vignette_CXCL1_ape, out = "/path/to/file/filename.ape")
## 
## # Save as fasta
## Ape.fasta(vignette_CXCL1_ape, out = "/path/to/file/filename.fa")
## 


## -----------------------------------------------------------------------------------------------------------------------------------------------
# get the sequence
Ape.DNA(vignette_CXCL1_ape)

# get the features
Ape.granges(vignette_CXCL1_ape)



## -----------------------------------------------------------------------------------------------------------------------------------------------
# define the new feature set
old_features <- Ape.granges(vignette_CXCL1_ape)
new_features <- old_features[mcols(old_features)$type == "gene"]

new_vignette_CXCL1_ape <- Ape.setFeatures(vignette_CXCL1_ape, gr = new_features)
new_vignette_CXCL1_ape


## ----eval=FALSE---------------------------------------------------------------------------------------------------------------------------------
## Ape.fimo(vignette_CXCL1_ape, fimo_feature = "CXCL1_gene")

