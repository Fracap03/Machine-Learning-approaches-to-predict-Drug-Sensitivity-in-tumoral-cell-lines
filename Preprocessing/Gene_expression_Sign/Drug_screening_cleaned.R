setwd("/Users/francescocaperna/Desktop/Github_Tesi/Machine-Learning-approaches-to-predict-Drug-Sensitivity-in-tumoral-cell-lines/Preprocessing/Gene_expression_Sign") # set the right directory

## Install dependencies
library(readxl)

## Read files
drug_screening <- read_excel("Pharmacogenomic_1001_long_GLOBAL.xlsx")

# Remove Nan from Drug screening IC50
dataset_filtrato_IC50 <- drug_screening[!is.na(drug_screening$IC50),]
dataset_filtrato_IC50 <- dataset_filtrato_IC50[dataset_filtrato_IC50$IC50!='NA',]

# Remove Nan from Drug screening Sensitivity
dataset_filtrato_sens <- drug_screening[!is.na(drug_screening['Sensitivity']),]


