# Install all the packages
#install.packages("imputeLCMD")

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

setwd("/Users/francescocaperna/Desktop/Github_Tesi/Machine-Learning-approaches-to-predict-Drug-Sensitivity-in-tumoral-cell-lines/Preprocessing/Gene_expression_Sign") # set the right dir
#BiocManager::install("pcaMethods",force = TRUE)
#BiocManager::install("impute", force = TRUE)
library(imputeLCMD)
library(dplyr)
library(tidyverse)
library(readr)
library(data.table)
library(biomaRt)
source('compute_zscore.R')
#########################
#                       #
#### TRANSCRIPTOMICS ####
#                       #
#########################


transcriptomics <- fread(paste0('CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt'))
transcriptomics_num <- transcriptomics %>%
  mutate(across(3:ncol(transcriptomics), as.numeric))  # convert all numeric columns starting from the third column

transcriptomics_clean <- transcriptomics_num %>%
  filter(rowSums(dplyr::select(., -c(gene_id, transcript_id)) == 0) < (ncol(.) - 3)) # remove rows that contain only zeros in expression values

transcriptomics_clean$transcript_id <- gsub("\\..*", "", transcriptomics_clean$transcript_id) # remove transcript version from 'transcript_id' to facilitate mapping
transcriptomics_clean$gene_id <- gsub("\\..*", "", transcriptomics_clean$gene_id)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 110)

# use cleaned 'transcript_id' values to retrieve gene names from Ensembl
mapping <- getBM(
  attributes = c('ensembl_transcript_id', 'external_gene_name'),
  filters = 'ensembl_transcript_id',
  values = transcriptomics_clean$transcript_id,
  mart = ensembl
)


# merge the retrieved mapping with the original dataframe using 'transcript_id'
transcriptomics_clean <- transcriptomics_clean %>%
  left_join(mapping, by = c("transcript_id" = "ensembl_transcript_id"))

# rename the column to store the gene name
names(transcriptomics_clean)[names(transcriptomics_clean) == 'external_gene_name'] <- 'gene.name'

#move 'gene.name' column to the first position
transcriptomics_clean <- transcriptomics_clean %>%
  dplyr::relocate(gene.name)

# perform log transformation and replace low values with zeros
transcriptomic_clean_log <- transcriptomics_clean %>%
  mutate(across(where(is.numeric), ~ log(. + 1)))


#remove unnecessary columns ('gene_id', 'transcript_id') after filtering
transcriptomic_filtered <- transcriptomic_clean_log %>%
  dplyr::select(-gene_id, -transcript_id)


# retrieve protein-coding genes with their external gene names
coding_genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                      filters = 'biotype',
                      values = 'protein_coding',
                      mart = ensembl)


# keep only protein-coding genes in the dataset
transcriptomic_clean_cod <- transcriptomic_filtered %>%
  filter(gene.name %in% coding_genes$external_gene_name)


dt <- as.data.table(transcriptomic_clean_cod)
transcriptomic_final_agg <- dt[!is.na(gene.name), lapply(.SD, mean, na.rm = TRUE), by = gene.name]


# remove the first row (if necessary, check why it is being removed)
transcriptomic_final_agg <- transcriptomic_final_agg[-1,]

# convert the dataset to a matrix for further processing
transcriptomic_clean_log_matrix <- as.matrix(transcriptomic_final_agg[, -c(1)])

# compute the z-score normalization by column using the mean as a metric
transcriptomic_clean_zscore <- compute_zscore(omic_matrix = transcriptomic_clean_log_matrix, by = "column", metric = "mean")
transcriptomic_clean_zscore <- as.data.frame(transcriptomic_clean_zscore)

transcriptomic_clean_zscore$gene.name <- transcriptomic_final_agg$gene.name
transcriptomic_clean_zscore <- transcriptomic_clean_zscore[, c(ncol(transcriptomic_clean_zscore), 1:(ncol(transcriptomic_clean_zscore)-1))]

write_csv(transcriptomic_clean_zscore, 'transcriptomic_zscore.csv')

