## Install dependencies
library(readr)
library(readxl)

setwd("/Users/francescocaperna/Desktop/Github_Tesi/Machine-Learning-approaches-to-predict-Drug-Sensitivity-in-tumoral-cell-lines") #set dir
## Read files
drug_screening <- read.delim("pazienti/GSE199455_AML_Fedratinib_preTreatment_annotations.txt", header = TRUE, sep = "\t")
merged_results<- read_csv("pazienti/merged_results_trans.csv")

setwd("/Users/francescocaperna/Desktop/Github_Tesi/Machine-Learning-approaches-to-predict-Drug-Sensitivity-in-tumoral-cell-lines/Preprocessing/Gene_expression_Sign") #set dir
source("Combining function.R")
setwd("/Users/francescocaperna/Desktop/Github_Tesi/Machine-Learning-approaches-to-predict-Drug-Sensitivity-in-tumoral-cell-lines") 


k <- 3000
for (i in seq(1, nrow(drug_screening) , by = k)) {
  start <- i
  end <- min(i + k - 1, nrow(drug_screening))  
  drug_proteom <- combine_datasets(drug_screening, merged_results, 
                                   'Sample', 'gene_name', 
                                   batch_start = start, batch_end = end)
  filename <- sprintf("pazienti/trans_zscore%d_%d.csv", start, end)
  write.csv(drug_proteom, filename, row.names = TRUE)  
}

