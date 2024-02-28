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
# 02_conformal_ODX
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
library(survival)
library(survminer)

# load("TCGA data/TCGA-BRCA_LB.Rd")
## 60660 x 1095 matrix of unstranded count data of RNAseq
## rows are genes, cols are patients
# 
## load("TCGA data/gene_annos.Rd")



### OLD - toy recurrence --------------------------------------------------

## I thought for a while that I will use recurrence instead of survival
## and have ODX directly predict recurrence, so tried to create some toy data
## because we don't have recurrence data 
## this is no longer needed however, because we are looking at survival
# 
# 
# # the approximate scores calculated are continuous and normalised
# # discretise them into bins to allow for a simple conformal prediction test 
# 
# # ODX_cats <- c("Bin 1", "Bin 2", "Bin 3", "Bin 4", "Bin 5", "Bin 6", "Bin 7", "Bin 8", "Bin 9", "Bin 10")
# clin_TGx$ODX_cat <- cut(clin_TGx$ODX, breaks = 10, labels = FALSE)
# 
# # create toy data based on vital status:
# # recurrence probability, if they are dead we will assume recurrence prob = 1
# clin_TGx <- clin_TGx %>%
#   mutate(vital_status_num = ifelse(vital_status == "Alive", 1, ifelse(vital_status == "Dead", 0, NA)))
# 
# # creating a toy probability that is taken from a uniform distribution if the patient is alive
# # and is 0 if the person is dead
# # could use different distributions, this is just a toy example
# set.seed(123)
# toy_recurrence_prob <- (runif(nrow(clin_TGx), min = -1, max = 1)) * clin_TGx$vital_status_num
# 
# clin_TGx <- cbind(clin_TGx, toy_recurrence_prob)
# 
# # let's only look at the chemo therapy treated patients, as this is what ODX aims at
# chemo_patients <- clin_TGx[grepl("Pharma", clin_TGx$treatment_type), ]
# # 543 remaining


### conformal uncertainty ----------------------------------------------------------
# conformal bins

# survival analysis:
# event: death
# event_time: either days_to_death for those who died or days_to_last_follow_up
# t_0 is the day of the diagnosis (days_to_diagnosis is 0 for everyone)

# we are then going to look at survival bins
# i) patient has died within >5 years of diagnosis (event_time < 5*365)
# ii) patient has died within 5-10 years of diagnosis (5*365 =< event_time < 10*365)
# iii) patient has NOT died within 10 years of diagnosis (event_time =< 10*365)


# > clin_TGx |>
#   + as_tibble() |>
#   + mutate(days_to_death = as.numeric(days_to_death)) |>
#   + summarise(longest_life = max(days_to_death, na.rm = T))
# # A tibble: 1 × 1
# longest_life
# <dbl>
#   1         7455

# change the data to be numerical, then try the binning again
# bin as intended and ignore the censoring for now.
# change both the event time coloumns into numericals

clin_TGx$days_to_death <- as.numeric(clin_TGx$days_to_death)



# keep in mind that only 151 people have died, so most people are alive after 3 years
# if days_days_to_last_follow_up
# censoring is when their last follow up was under 3 years, i.e. less than the entire time
clin_TGx <- clin_TGx %>%
  mutate(surv_bin = ifelse(is.na(days_to_death), 3, # all of the patients which have no recorded death event in the dataset # careful: censoring!
                          ifelse(days_to_death < 5*365, 1, # i) patient has died within >5 years of diagnosis
                           ifelse(days_to_death < 10*365, 2, # ii) patient has died within 1.5-3 years of diagnosis
                                  3)))) # iii) patient has died later than 3 years of diagnosis
                                  # if time_to_last_follow_up is under 3 years, we cannot be sure that they are still alive -> censoring
                                  # they could either be in the truly alive bin (category 3) or in one of the others
                                  # how do we deal with censoring?
                                  # for now we just WRONGLY treat everyone alive as alive at 3 years -> category 3

hist(clin_TGx$surv_bin)

## SURVIVAL ANALYSIS -----------------------------------------------------------
# for status we are using the survival bins
# nah this doesn't work, let's try with 1 for event (death) and 0 for no event recorded (survival or censoring)
clin_TGx <- clin_TGx %>%
  mutate(event = ifelse(!is.na(clin_TGx$days_to_death), 1, 0))
                        # ifelse(!is.na(clin_TGx$days_to_last_follow_up), 0, NA)))

# for times we are using event time, which is either days_to_death for those who died or days_to_last_follow_up
# for ease of use we are merging the event_times into one coloumn
clin_TGx <- clin_TGx %>%
  mutate(event_time = coalesce(clin_TGx$days_to_death, clin_TGx$days_to_last_follow_up))

# clin_TGx <- clin_TGx[, -ncol(clin_TGx)]

clin_TGx$event_time <- as.numeric(clin_TGx$event_time)

summary(coxph(Surv(time = event_time, event =  event) ~ ODX + strata(race), clin_TGx))

# Call:
#   coxph(formula = Surv(time = event_time, event = event) ~ ODX + 
#           strata(race), data = clin_TGx)
# 
# n= 1092, number of events= 151 
# 
# coef exp(coef) se(coef)    z Pr(>|z|)
# ODX 0.4397    1.5522   0.5634 0.78    0.435
# 
# exp(coef) exp(-coef) lower .95 upper .95
# ODX     1.552     0.6443    0.5145     4.683
# 
# Concordance= 0.546  (se = 0.033 )
# Likelihood ratio test= 0.58  on 1 df,   p=0.4
# Wald test            = 0.61  on 1 df,   p=0.4
# Score (logrank) test = 0.61  on 1 df,   p=0.4

## it seems ODX is not significant for survival?

fit <- survfit(Surv(time = event_time, event =  event) ~ ODX + strata(race), clin_TGx)
# ggsurvplot(fit)
# plot(fit)

###################
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  #conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in years",   # customize X axis label
  xlim=c(0,2),
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  # risk.table = "abs_pct",  # absolute number and percentage at risk.
  # risk.table.y.text.col = T,# colour risk table text annotations.
  # risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  # ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  # legend.labs = 
  #   c("low","med","high"),    # change legend labels.
  palette = 
    c("pink","#E7B800", "#2E9FDF") # custom color palettes.
)
