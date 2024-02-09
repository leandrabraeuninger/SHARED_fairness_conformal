# ---------------------------------------------------------------------
# FAIRNESS VIA CONFORMAL UNCERTAINTY QUANTIFICATION FOR TRANSCRIPTOMICS
#
# this mini-project is described in more detail here:
# https://www.overleaf.com/read/nxpvdqkmsgnp#b69020
#
# this is the first script, which is a bunch of notes so far and will
# be developed and then split into sub-scripts
# ---------------------------------------------------------------------
#
# 00_DATA-EXPLORATION
#
### SET UP -----------------------------------------------------------

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

# load("TCGA data/TCGA-BRCA_LB.Rd")
## 60660 x 1095 matrix of unstranded count data of RNAseq
## rows are genes, cols are patients
# 
## load("TCGA data/gene_annos.Rd")





## put these sections into different scripts


### DATA EXPLORATION --------------------------------------------------

colnames(clin)
# 157 clinical and sociodemographic variables

unique(clin$ethnicity)
# very limited: 
# [1] "not hispanic or latino" "not reported"           "hispanic or latino"     "'--"    
unique(clin$race)
# [1] "white"                            "black or african american"        "not reported"                    
# [4] "asian"                            "american indian or alaska native" "'--"  

unique(clin$country_of_residence_at_enrollment)

range(as.numeric(clin$age_at_index), na.rm = TRUE)
unique(clin$occupation_duration_years)

unique(gene_annos$gene_type)


# interesting cols:
# race, vital_status, age_at_diagosis, age_at_index, occupation_duration_years, 
# gender, initial_disease_status, last_known_disease_status, tumor_stage, 
# treatment_type, treatment_arm


table(clin$vital_status)
table(clin$race)

## let's create some visualisations to explore the data
## pies

# # Function to create a pie chart for a specific column in the dataset
# create_pie_chart <- function(data, column_name) {
#   # Count occurrences for each category
#   counts <- table(data$column_name)
# 
#   # Create a data frame from the counts
#   counts_df <- as.data.frame(counts)
# 
#   # # Rename columns for better aesthetics
#   # colnames(counts_df) <- c("Category", "Count")
# 
#   # ggplot(counts_df) +
#   #   geom_bar(stat = "identity", width = 1, color = "white") +
#   #   coord_polar("y") +
#   #   ggtitle(paste("Pie Chart for", column_name)) +
#   #   theme_void()
# 
#   ggplot(counts_df, aes(x = "", y = Freq)) +
#     geom_bar(stat = "identity", width = 1, color = "white") +
#     coord_polar("y") +
#     ggtitle(paste("Pie Chart for", column_name)) +
#     geom_text(aes(label = Var1), position = position_stack(vjust = 0.5)) +
#     theme_void()
# }

# Loop through columns (excluding the first two column which are ids)
# for (col in colnames(clin)[-2]) {
#   pie_chart <- create_pie_chart(clin, col)
#   print(pie_chart)
# }



test <- ggplot(as.data.frame(table(clin$vital_status)), aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  ggtitle(paste("Pie Chart for TEST")) +
  geom_text(aes(label = Var1), position = position_stack(vjust = 0.5)) +
  theme_void()
# Create an empty patchwork layout
combined_plots <- test



# # for all the vars:
# for (col in colnames(clin)) {
#   counts <- as.data.frame(table(clin[[col]]))
#   # print(counts)
#   plot <- ggplot(counts, aes(x = "", y = Freq, fill = Var1)) +
#     geom_bar(stat = "identity", width = 1, color = "white") +
#     coord_polar("y") +
#     ggtitle(paste("Pie Chart for", col)) +
#     geom_text(aes(label = Var1), position = position_stack(vjust = 0.5)) +
#     theme_void()
#     combined_plots <- combined_plots + plot
# }

# # for the interesting vars:

interesting_cols <- c( "race", "vital_status", "ethnicity",
                       # "age_at_diagnosis", "age_at_index",
                       "occupation_duration_years", "gender", "initial_disease_status", 
                       "chemo_concurrent_to_radiation", "classification_of_tumor",
                       "year_of_birth",
                       "last_known_disease_status", "tumor_stage", "treatment_type", "treatment_arm")
# int_vars_plots <- interesting_cols

for (col in interesting_cols) {
  counts <- as.data.frame(table(clin[[col]]))
  # print(counts)
  plot <- ggplot(counts, aes(x = "", y = Freq, fill = Var1)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y") +
    ggtitle(paste("Pie Chart for", col)) +
    geom_text(aes(label = paste(Var1, Freq)), position = position_stack(vjust = 0.5)) +
    theme_void()

  # int_vars_plots <- int_vars_plots + plot
  print(plot)
  ggsave(paste0(Sys.Date(), "_TCGA_data_",  col, ".png"), plot = plot, width = 40, height = 20)
}

# print(int_vars_plots) # doesn't work, must combine with the method used in the review

# the age variables are not great
# year_of_birth seems too equally distributed -> what's the reason for this? how were the patients selected for this dataset?
# treatmen_arm has no info


# ggplot(as.data.frame(table(clin$vital_status)), aes(x = "", y = Freq, fill = Var1)) +
#   geom_bar(stat = "identity", width = 1, color = "white") +
#   coord_polar("y") +
#   ggtitle(paste("Pie Chart for vital status")) +
#   geom_text(aes(label = Var1), position = position_stack(vjust = 0.5)) +
#   theme_void()



### EXTRACTING THE Oncotype DX scores -------------------------------------

### different models for survival prediction (just use all of the inputs?)
## create survival at 5 years after diagnosis and 10 years after diagnosis labels

