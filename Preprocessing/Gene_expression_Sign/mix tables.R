setwd("/Users/francescocaperna/Desktop/Github_Tesi/Machine-Learning-approaches-to-predict-Drug-Sensitivity-in-tumoral-cell-lines/Preprocessing/Gene_expression_Sign") #set the right dir

## Install dependencies
library(readr)
library(readxl)

source("Drug_screening_cleaned.R")
source("Combining function.R")

## Read files
transcriptomic_clean_zscore <- read_csv("transcriptomic_zscore.csv")
merged_results<- read_csv("merged_results_trans.csv")


k <- 3000
for (i in seq(1, nrow(drug_screening) , by = k)) {
  start <- i
  end <- min(i + k - 1, nrow(drug_screening))  
  drug_proteom <- combine_datasets(drug_screening, merged_results, 
                                   'cell.name', 'gene_name', 
                                   batch_start = start, batch_end = end)
  
  filename <- sprintf("transcrittomica_finale/trans_zscore%d_%d.csv", start, end)

  write.csv(drug_proteom, filename, row.names = TRUE)  

}

