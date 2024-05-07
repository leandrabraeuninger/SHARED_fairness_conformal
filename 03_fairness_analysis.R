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
# install.packages("dplyr")
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
# 2.818681


mean_by_subgroup <- data.frame(aggregate(lengths(prediction_sets) ~ race, data = test, FUN = mean))
# race lengths(prediction_sets)
# 1                     asian                 2.521739
# 2 black or african american                 2.783333
# 3              not reported                 2.696970
# 4                     white                 2.870968

# larger average mean for white subpopulation
# that sounds like the uncertainty would be higher for the people of that subpopulation
colnames(mean_by_subgroup)[2] <- "mean_length"


# !!
# I'm not saying anything about whether or not they are more or less right

ggplot(mean_by_subgroup, aes(x = race, y = mean_length, fill = race)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean Lengths of Prediction Sets by Race",
       x = "Race",
       y = "Mean Length",
       fill = "Race") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# seems very minimally informative
# the uncertainty is quite high. with 3 risk bins available, and almost always 
# including all 3 risk bins (mean is above 2.5 in all subgroups), that means the uncertainty is really high
# my bins are a bit useless probably

range_by_subgroup <- data.frame(aggregate(lengths(prediction_sets) ~ race, data = test, FUN = range))
# all subgroups have predictions in all risk bins

# count(lengths(test$prediction_sets))
table(lengths(test$prediction_sets))

# Count occurrences of unique lengths of prediction sets by race
occurrences_by_race <- test %>%
  group_by(race, lengths_prediction_sets = lengths(prediction_sets)) %>%
  count() %>%
  ungroup() %>%
  rename(occurrences = n)

ggplot(occurrences_by_race, aes(x = lengths_prediction_sets, y = occurrences, fill = race)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Occurrences of Unique Lengths of Prediction Sets by Race",
       x = "Length of Prediction Sets",
       y = "Occurrences",
       fill = "Race") +
  theme_minimal()