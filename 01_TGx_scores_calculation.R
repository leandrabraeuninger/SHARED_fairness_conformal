# ---------------------------------------------------------------------
# FAIRNESS VIA CONFORMAL UNCERTAINTY QUANTIFICATION FOR TRANSCRIPTOMICS
#
# script to extract the Oncotype Dx and Mammaprint scores from
# the TCGA transcriptomic profiles
# ---------------------------------------------------------------------
#
# 01_TGx_SCORES-CALCULATION
#
### SET UP -----------------------------------------------------------

# rm(lis=ls())
# setwd("")
# read.

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
# load("data/TCGA data/gene_annos.Rd")



### SCORE CALCULATION --------------------------------------------------

## ONCYTYPE DX ---------------------------------------------------------
## taken from Chris Banerji's script "initial_notes.r" lines 175-199
## slightly adapted
## I'm unsure if the coefficients are available in the initial paper
#####################
prolif <- c("MKI67","AURKA","BIRC5","CCNB1","MYBL2")
invasion <- c("MMP11","CTSV")
HER2 <- c("GRB7","ERBB2")
oest <- c("ESR1","PGR","BCL2","SCUBE2")
g <- "GSTM1"
b <- "BAG1"
cd <- "CD68"
og <- c(prolif,invasion,HER2,oest,g,b,cd)
refs <- c("GAPDH","RPLP0","GUSB","ACTB","TFRC")
annos <- gene_annos

onco_un <- 0.47*apply(data[match(HER2,annos$gene_name),],2,sum)-
            0.34*apply(data[match(oest,annos$gene_name),],2,sum)+
            1.04*apply(data[match(prolif,annos$gene_name),],2,sum)+
            0.1*apply(data[match(invasion,annos$gene_name),],2,sum)+
            0.05*data[match(cd,annos$gene_name),]-
            0.08*data[match(g,annos$gene_name),]-
            0.07*data[match(b,annos$gene_name),]
## RT-qPCR based , so we normalise
norms <- apply(data[match(refs,annos$gene_name),],2,sum)

hist(onco_un)
hist(norms)
hist(onco_un/norms)
onco_n <- onco_un/norms
## TODO check this, bc the ODX paper mentions that the sccores are between 0-100
## here the normalised scores lie between -0.7402087 and 0.8885420
#########################

ODX <- onco_n #  copying this to a different name for clarity

## should be able to match these scores to the clinical data by ids
# checking that they are in the same order:
all.equal(rownames(clin) , names(ODX)) # they are

clin_TGx <- cbind(clin, ODX) # add the ODX scores to the full table


# the table does not use NA for empty values, so let's replace that
clin_TGx[clin_TGx == "'--"] <- NA

# finding the empty cols (this assigns true or false to each colname)
empty_columns <- colSums(is.na(clin_TGx)) == nrow(clin_TGx)
# listing only the empty ones
empty_columns <- names(clin_TGx)[empty_columns] 
# removing them
clin_TGx <- clin_TGx[, colSums(is.na(clin_TGx)) != nrow(clin_TGx)] # subsetting all rows and only the non-empty cols
# 36 vars left


## Mammaprint (MP) ---------------------------------------------------------
## do the same for the Mammaprint score

## choose the genes
## based on this paper I extracted the 70 gene names used for MP
# https://www.futuremedicine.com/doi/10.2217/fon-2018-0221?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed#_i2
# also through the original MP paper, but it didn't provide a full list easily extractable
# https://www.nature.com/articles/415530a#Sec1

# I'm still unsure about the scoring algorithm itself and whether other gene expression data is used to inform the scores
# the full info might be proprietary
# taken the gene list from the figure description in Figure 2
# this lists some genes multiple times due to the content of the paper
# it doesn't mean anything for us
MP_genes <- c("FLT1", "HRASLS", "STK32B", "RASSF7", "DCK", "MELK", 
              "EXT1", "GNAZ", "EBF4", "MTDH", "PITRM1", "QSCN6L1", 
              "TGFB3", "FLT1", "HRASLS", "STK32B", "RASSF7", "DCK", 
              "MELK", "EXT1", "GNAZ", "EBF4", "MTDH", "PITRM1", 
              "QSCN6L1", "COL4A2", "GPR180", "MMP9", "GPR126", "RTN4RL1", 
              "DIAPH3", "CDC42BPA", "PALM2", "CCNE2", "ECT2", "CENPA", 
              "LIN9", "KNTC2", "MCM6", "NUSAP1", "ORC6L", "TSPYL5", 
              "RUNDC1", "PRC1", "RFC4", "RECQL5", "CDCA7", "DIL", 
              "ALDH4A1", "AYTL2", "OXCT1", "PEC1", "GMPS", "GSTM3", 
              "SLC2A3", "FLT1", "FGF18", "COL4A2", "GPR180", "EGLN1", 
              "MMP9", "BBC3", "EGLN1", "FLT1", "HRASLS", "STK32B", 
              "RASSF7", "DCK", "MELK", "EXT1", "GNAZ", "EBF4", 
              "MTDH", "PITRM1", "QSCN6L1", "LOC100288906", "C9orf30", 
              "ZNF533", "C16orf61", "SERF1A", "C20orf46", "LOC730018", "LOC100131053", 
              "AA555029_RC", "LGP2", "NMU", "UCHL5", "JHDM1D", "AP2B1", 
              "MS4A7", "RAB6B"
)
length(MP_genes)
length(unique(MP_genes)) # but why are there 8 missing for the full 70 genes?

# only keeping each gene once
MP_genes <- unique(MP_genes)

# trying to match the MP genes to all the gene expression profiles I have from TCGA
# something like this
unweightedMPgenes <- data[match(MP_genes,annos$gene_name), ]

# some genes are not found in the TCGA set
# could be because of different naming conventions
missing <- MP_genes[is.na(rownames(unweightedMPgenes))] # get input from Chris Banerji on whether these are actually missing in the gene_annos or just slightly differently named

## define the algorithm (normalised if needed)
## add it to the the main dataframe (clin_TGx)
clin_TGx <- cbind(clin, MP)


## the MP part was abandoned for the continuing analysis
