# Install all the packages
#install.packages("imputeLCMD")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("pcaMethods",force = TRUE)
#BiocManager::install("impute", force = TRUE)

setwd("/Users/francescocaperna/Desktop/Github_Tesi/Machine-Learning-approaches-to-predict-Drug-Sensitivity-in-tumoral-cell-lines") #set dir

library(imputeLCMD)
library(dplyr)
library(tidyverse)
library(readr)
library(biomaRt)
library(data.table)
source('Preprocessing/Gene_expression_Sign/compute_zscore.R')
#########################
#                       #
#### TRANSCRIPTOMICS ####
#                       #
#########################

transcriptomics <- fread(paste0('pazienti/GSE199455_AML_Fedratinib_preTreatment_TPM.txt'))



# perform log transformation and replace low values with zeros
transcriptomic_clean_log <- transcriptomics %>%
  mutate(across(where(is.numeric), ~ log(. + 1)))

# replace all values ≤ 0.1 with 0 (only for numeric columns)
transcriptomic_clean_log <- transcriptomic_clean_log %>%
mutate(across(where(is.numeric), ~ ifelse(. <= 0.1, 0, .)))


#connect to Ensembl and specify the human dataset
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# retrieve protein-coding genes with their external gene names
coding_genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                      filters = 'biotype',
                      values = 'protein_coding',
                      mart = ensembl)


# keep only protein-coding genes in the dataset
transcriptomic_clean_log <- transcriptomic_clean_log %>%
  filter(Gene %in% coding_genes$external_gene_name)


dt <- as.data.table(transcriptomic_clean_log)
transcriptomic_final_agg <- dt[!is.na(Gene), lapply(.SD, mean, na.rm = TRUE), by = Gene]


# remove the first row (if necessary, check why it is being removed)
transcriptomic_final_agg <- transcriptomic_final_agg[-1,]

# convert the dataset to a matrix for further processing
transcriptomic_clean_log_matrix <- as.matrix(transcriptomic_final_agg[, -c(1)])

# compute the z-score normalization by column using the mean as a metric
transcriptomic_clean_zscore <- compute_zscore(omic_matrix = transcriptomic_clean_log_matrix, by = "column", metric = "mean")
transcriptomic_clean_zscore <- as.data.frame(transcriptomic_clean_zscore)

transcriptomic_clean_zscore$gene.name <- transcriptomic_final_agg$Gene
transcriptomic_clean_zscore <- transcriptomic_clean_zscore[, c(ncol(transcriptomic_clean_zscore), 1:(ncol(transcriptomic_clean_zscore)-1))]

write_csv(transcriptomic_clean_zscore, 'pazienti/transcriptomic_zscore.csv')

