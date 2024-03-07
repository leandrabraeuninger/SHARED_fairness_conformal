# ---------------------------------------------------------------------
# FAIRNESS VIA CONFORMAL UNCERTAINTY QUANTIFICATION FOR TRANSCRIPTOMICS
#
# this mini-project is described in more detail here:
# https://www.overleaf.com/read/nxpvdqkmsgnp#b69020
#
# third script to use conformal uncertainty quantification on the rough
# approximation of Oncotype DX scores calculated in the previous script
# ---------------------------------------------------------------------
#
# 02_conformal_ODX_SECOND
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
# install.packages("tcltk2")

library(ggplot2)
library(patchwork)
library(dplyr)
library(survival)
library(survminer)
library(tcltk2)


# load("TCGA data/TCGA-BRCA_LB.Rd")
## 60660 x 1095 matrix of unstranded count data of RNAseq
## rows are genes, cols are patients
# 
## load("TCGA data/gene_annos.Rd")


#
# clin_TGx$days_to_death <- as.numeric(clin_TGx$days_to_death)
clin_TGx$event_time <- as.numeric(clin_TGx$event_time)
clin_TGx$event <- as.numeric(clin_TGx$event)

# TODO create a new smaller dataframe that contains only the coloumns that we need
# data <- as_tibble(cbind(clin_TGx$case_id, clin_TGx$case_submitter_id, clin_TGx$race, clin_TGx$ODX,
#               clin_TGx$event, clin_TGx$event_time))



### conformal analysis -----------------------------------------------------------
# Chris B's recidivism script does all of this in a loop in order to do it with 30 different samples from the data
# TODO do that

# we split the data into train, calibration and test sets

# so for the split, Chris B. did a 3rd split
# TODO should I remove the right-censored patients?
# because right now I simply treated them as low-risk

# for now we are going to do the same and use 364 for train, test and calib each


## TRAIN
# with the train set we calculate survival probability based on the ODX scores
# HOW?

# survfit(Surv(time = event_time, event =  event) ~ race, clin_TGx)
# wait this works:
summary(coxph(Surv(time = event_time, event =  event) ~ ODX, clin_TGx))

# Call:
#   coxph(formula = Surv(time = event_time, event = event) ~ clin_TGx$ODX, 
#         data = clin_TGx)
# 
# n= 1092, number of events= 151 
# 
# coef exp(coef) se(coef)    z Pr(>|z|)
# clin_TGx$ODX 0.5189    1.6802   0.5523 0.94    0.347
# 
# exp(coef) exp(-coef) lower .95 upper .95
# clin_TGx$ODX      1.68     0.5952    0.5692      4.96
# 
# Concordance= 0.573  (se = 0.029 )
# Likelihood ratio test= 0.83  on 1 df,   p=0.4
# Wald test            = 0.88  on 1 df,   p=0.3
# Score (logrank) test = 0.88  on 1 df,   p=0.3

# Survfit is Kaplan Meier
# only takes categorical values
fit <- survfit(Surv(time = event_time, event =  event) ~ ODX, clin_TGx)
# this def doesn't do what I want, there's no point in stratifying by ODX
plot(fit)

# TODO play with ideal cutpoints
cutpoint_1 <- surv_cutpoint(data = clin_TGx, time = "event_time", event = "event", variables = "ODX")
# creates one cutpoint, I want two
# Chris B. in his recidivism example simply found the second cutpoint by filtering 
# for event (here death) rows with events after 1y and repeating the surv_cutpoint()
# he created bins at <1 year, 1-2 years, and >2 years
# I am doing similarly (see test data bit) with the following bins:
# <1 year, 1-5 years, and >5 years
temp <- clin_TGx[clin_TGx$event == 1, ] # filtering for death events to get the cutpoints
temp$event[temp$event_time > 365] <- 0 # setting the event to 0 for the deaths that occured after the first year 
cutpoint_2 <- surv_cutpoint(data = temp, time="event_time", event = "event", variables ="ODX")
# extracting the cutpoints and putting them togethe:
cutpoints <- unlist(c(cutpoint_1$cutpoint[1],cutpoint_2$cutpoint[1]))
# these are the cutpoints for the ODX scores to be sorted into risk bins
# TODO check that this is entirely correct, because we haven't used the 5y cutoff here, 
# just two somewhat random cutoff points
# TODO check that this makes sense! is a higher ODX score really a higher risk?

# with this we have 3 bins, let's fill them:
clin_TGx <- clin_TGx %>%
  mutate(pred_risk_bin = ifelse(ODX <= cutpoints[1], 1, # high risk 
                                ifelse(ODX <= cutpoints[2], 2, # medium risk
                                       3))) # low risk
hist(clin_TGx$pred_risk_bin)
# > sum(clin_TGx$pred_risk_bin == 3)
# [1] 112
# > sum(clin_TGx$pred_risk_bin == 2)
# [1] 630
# > sum(clin_TGx$pred_risk_bin == 1)
# [1] 350

## TEST
# for the test set we calculate the true survival probability according to the data

# no for now we go with risk bins
# somewhat arbitrarily we define survival risk categories as follows:
# 1 - high risk: patient died within 1 year of diagnosis
# 2 - medium risk: patient died 1-5 years after diagnosis OR had a positive (alive) follow up less than 1 year after diagnosis
# 3 - low risk: patient died more than 5 years after diagnosis OR had a positive (alive) follow up more than 1 year after diagnosis

clin_TGx <- clin_TGx %>%
  mutate(true_risk_bin = ifelse(event_time <= 365 & event == 1 , 1, # high risk
                                ifelse(event_time <= 5*365 & event == 1, 2, # medium risk (dead)
                                       ifelse(event_time <= 365 & event == 0, 2, # medium risk (alive)
                                              ifelse(event_time > 5*365 & event == 1, 3, # low risk (dead)
                                                     ifelse(event_time > 365 & event == 0, 3, 4)))))) # low risk (alive) and others: 4

# just to check that the bins are assigned correctly, replacing them with 0:
# clin_TGx <- clin_TGx %>%
#   mutate(true_risk_bin = 0)
hist(clin_TGx$true_risk_bin)
# > sum(clin_TGx$true_risk_bin == 3)
# [1] 827
# > sum(clin_TGx$true_risk_bin == 2)
# [1] 243
# > sum(clin_TGx$true_risk_bin == 1)
# [1] 22

## COMMENT
# of course this bin allocation could be done differently, especially the patients with an alive follow up
# could be classified into higher risk bins 


## CALIB 
# for the calibration we calculate a conformal score that represents the 
# inaccuracy / difference between the TRUE survival probability and the one based on ODX

# for now the ODX predicted bin is correct for only 236 instances
# TODO play with the cutoff values to make them more accurate

# provide a conformal score for the calibration set:
# if the bin is correctly assigned, the conformal score is 0 
# if not, assign the score with the "smallest positive distance between the ODX score and the 
# prediction value that would have assigned the appropriate risk category" (taken from Chris B's recidivism task)

clin_TGx <- clin_TGx %>%
  mutate(conf_score = ifelse(pred_risk_bin == true_risk_bin, 0, # if correct bin assigned the conformal score is 0, otherwise
                             ifelse(pred_risk_bin == 1 & true_risk_bin == 2, abs(ODX - cutpoints[1]), 
                                    ifelse(pred_risk_bin == 1 & true_risk_bin == 3, abs(ODX - cutpoints[2]),
                                           ifelse(pred_risk_bin == 2 & true_risk_bin == 1, abs(ODX - cutpoints[1]),
                                                  ifelse(pred_risk_bin == 2 & true_risk_bin == 3, abs(ODX - cutpoints[2]),
                                                         ifelse(pred_risk_bin == 3 & true_risk_bin == 1, abs(ODX - cutpoints[1]),
                                                                ifelse(pred_risk_bin == 3 & true_risk_bin == 2, abs(ODX - cutpoints[2]),
                                                                       NA # else NA
                                                  ))))))))
# IT WORKS!
# ok what do I do with this now?
