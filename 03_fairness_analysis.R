---------------------------------------------------------------------
  # FAIRNESS VIA CONFORMAL UNCERTAINTY QUANTIFICATION FOR TRANSCRIPTOMICS
  #
  # this mini-project is described in more detail here:
  # https://www.overleaf.com/read/nxpvdqkmsgnp#b69020
  #
  # second script to extract the Oncotype Dx and Mammaprint scores from
  # the TCGA transcriptomic profiles
  # ---------------------------------------------------------------------
#
# 01_TGx_SCORES-CALCULATION
#
### SET UP -----------------------------------------------------------

## TODO
# finish this to make sure it can run both independently and with the 
# other scripts

# rm(lis=ls())
# setwd("")
# read.

# DESeq2 
# not available for this version of R
#  "R version 4.3.1 (2023-06-16 ucrt)"
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# load("DESeq2")



# install.packages("ggplot2")
# install.packages("patchwork")
library(ggplot2)
library(patchwork)
library(dplyr)

# load("data/TCGA data/TCGA-BRCA_LB.Rd")
## 60660 x 1095 matrix of unstranded count data of RNAseq
## rows are genes, cols are patients
# 
## load("TCGA data/gene_annos.Rd")

# load the output from the previous script
# 02_conformal_ODX.R


### FAIRNESS ANALYSIS --------------------------------------------------
mean(lengths(prediction_sets))

mean_by_subgroup <- aggregate(lengths(prediction_sets) ~ race, data = clin_TGx, FUN = mean)
#                               race lengths(prediction_sets)
# 1 american indian or alaska native                 2.000000
# 2                            asian                 2.803279
# 3        black or african american                 2.873626
# 4                     not reported                 2.894737
# 5                            white                 2.922975
# right there's essentially no difference?
# run again after proper data split

