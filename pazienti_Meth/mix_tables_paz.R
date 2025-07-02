setwd("/Users/francescocaperna/Desktop/Tesi biologia/Francesco/")

## Install dependencies
library(readr)
library(readxl)

source("matrice_finale/Drug_screening_cleaned.R")
source("matrice_finale/Combining function.R")

## Read files
#transcriptomic_clean_zscore <- read_csv("file_puliti_ele/transcriptomic_clean_zscore.csv")
#proteomics_zscore <- read_csv("proteomic.csv")

library(readxl)

# Percorso del file
drug_screening <- read_excel("pazienti_2/Responder.xlsx")
merged_results_trans <- read.csv("merged_results_trans.csv", header = TRUE, stringsAsFactors = FALSE)



k <- 3000
for (i in seq(1, nrow(drug_screening) , by = k)) {
  start <- i
  end <- min(i + k - 1, nrow(drug_screening))  
  #debug(combine_datasets)
  #debug(combine_datasets)
  drug_proteom <- combine_datasets(drug_screening, merged_results, 
                                   'Sample', 'gene_name', 
                                   batch_start = start, batch_end = end)
  
  filename <- sprintf("pazienti_2/trans_zscore%d_%d.csv", start, end)
  
  #gzfile_conn <- gzfile(filename, "wt")  # Apertura connessione compressa
  write.csv(drug_proteom, filename, row.names = TRUE)  
  #close(gzfile_conn)  # Chiudere la connessione per evitare problemi
}

