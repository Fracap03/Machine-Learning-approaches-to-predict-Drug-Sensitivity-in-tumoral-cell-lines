# Install all the packages
#install.packages("imputeLCMD")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("pcaMethods",force = TRUE)
#BiocManager::install("impute", force = TRUE)
library(imputeLCMD)
library(dplyr)
library(tidyverse)
library(readr)
library(org.Hs.eg.db)
library("AnnotationDbi")


setwd("/Users/francescocaperna/Desktop/Tesi biologia/Francesco/")
source('file_puliti_ele/compute_zscore.R')
#########################
#                       #
#### TRANSCRIPTOMICS ####
#                       #
#########################

#transcriptomics <- read_delim('CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt', delim = '\t') 
library(data.table)
#transcriptomics <- fread(paste0('file_iniziali/CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt'))
transcriptomics <- fread(paste0('pazienti_2/GSE81259_log2rpkmTable.txt'))
#transcriptomics_num <- transcriptomics %>%
  #mutate(across(3:ncol(transcriptomics), as.numeric))  # convert all numeric columns starting from the third column
valid_keys <- na.omit(as.character(transcriptomics$geneid))

transcriptomics$Gene <- mapIds(org.Hs.eg.db,
                           keys = valid_keys,
                           column = "SYMBOL",
                           keytype = "ENTREZID",
                           multiVals = "first")
transcriptomics[, geneid := NULL]







#transcriptomics_clean <- transcriptomics_num %>%
  #filter(rowSums(dplyr::select(., -c(gene_id, transcript_id)) == 0) < (ncol(.) - 3)) # remove rows that contain only zeros in expression values


#transcriptomics_clean$transcript_id <- gsub("\\..*", "", transcriptomics_clean$transcript_id) # remove transcript version from 'transcript_id' to facilitate mapping

#transcriptomics_clean$gene_id <- gsub("\\..*", "", transcriptomics_clean$gene_id)

library(biomaRt)

# Prova con il mirror degli Stati Uniti Est
#ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")

# Prova con il mirror dell'Asia
#ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "asia")

#Prova co funzione precedente 
#ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 110)

# use cleaned 'transcript_id' values to retrieve gene names from Ensembl
#mapping <- getBM(
  #attributes = c('ensembl_transcript_id', 'external_gene_name'),
  #filters = 'ensembl_transcript_id',
  #values = transcriptomics_clean$transcript_id,
  #mart = ensembl
#)
# Salva l'oggetto ensembl in un file RDS
#saveRDS(ensembl, "ensembl_object.rds")
# Salva il mapping in un file CSV
#write.csv(mapping, "mapping_ensembl.csv", row.names = FALSE)


# merge the retrieved mapping with the original dataframe using 'transcript_id'
#transcriptomics_clean <- transcriptomics_clean %>%
  #left_join(mapping, by = c("transcript_id" = "ensembl_transcript_id"))

# rename the column to store the gene name
#names(transcriptomics_clean)[names(transcriptomics_clean) == 'external_gene_name'] <- 'gene.name'

#move 'gene.name' column to the first position
#transcriptomics_clean <- transcriptomics_clean %>%
  #dplyr::relocate(gene.name)

# perform log transformation and replace low values with zeros
#transcriptomic_clean_log <- transcriptomics %>%
  #mutate(across(where(is.numeric), ~ log(. + 1)))

# replace all values ≤ 0.1 with 0 (only for numeric columns)
transcriptomic_clean_log <- transcriptomics %>%
mutate(across(where(is.numeric), ~ ifelse(. <= 0.1, 0, .)))

# define a threshold for filtering
#threshold <- 1e-10

# filter out rows where more than 80% of values are below the threshold
#transcriptomic_clean_log_filtered <- transcriptomic_clean_log %>%
#mutate(percent_zero = rowMeans(across(where(is.numeric), ~ . <= threshold), na.rm = TRUE) * 100) %>%
#filter(percent_zero < 80) %>%
#dplyr::select(-percent_zero)  # Remove the auxiliary column

#remove unnecessary columns ('gene_id', 'transcript_id') after filtering
#transcriptomic_filtered <- transcriptomic_clean_log %>%
  #dplyr::select(-gene_id, -transcript_id)

#connect to Ensembl and specify the human dataset
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# retrieve protein-coding genes with their external gene names
coding_genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                      filters = 'biotype',
                      values = 'protein_coding',
                      mart = ensembl)

# Salva coding_genes in un file CSV
#write.csv(coding_genes, "coding_genes.csv", row.names = FALSE)

# display the first few results
#head(coding_genes)

# keep only protein-coding genes in the dataset
transcriptomic_clean_log <- transcriptomic_clean_log %>%
  filter(Gene %in% coding_genes$external_gene_name)

# aggregate data by 'gene.name' (if necessary) using the mean of expression values
#transcriptomic_final_agg <- transcriptomic_clean_cod %>%
#group_by(gene.name) %>%
#summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = 'drop') %>%
#filter(!is.na(gene.name))  # Remove NA gene names

dt <- as.data.table(transcriptomic_clean_log)

# Rimuove righe con NA in gene.name e aggrega calcolando la media
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

write_csv(transcriptomic_clean_zscore, 'pazienti_2/transcriptomic_zscore.csv')

