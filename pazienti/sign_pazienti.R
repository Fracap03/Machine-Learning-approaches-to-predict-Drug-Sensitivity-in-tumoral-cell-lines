setwd("/Users/francescocaperna/Desktop/Github_Tesi/Machine-Learning-approaches-to-predict-Drug-Sensitivity-in-tumoral-cell-lines") #set dir

#devtools::install_github('https://github.com/SaccoPerfettoLab/SignalingProfiler/')
library(dplyr)
library(tidyr)
## Import dependencies
library(readr)
library(readxl)
library(SignalingProfiler)
source('Preprocessing/Gene_expression_Sign/footprint_functions_Francesco.R')

#Reading file
transcriptomic <- read_csv('pazienti/transcriptomic_zscore.csv')

gene_name <- transcriptomic[['gene.name']]

#   gene_name difference logpval significant
# 1    A4GALT -0.4963227      NA        <NA>
# 2      AAAS -0.8501827 2.39387           +
# 3      AACS -0.6911841      NA        <NA>
# 4     AAED1  1.3261344      NA        <NA>
# 5     AAGAB  0.5938689      NA        <NA>
# 6     AAMDC  0.5814601      NA        <NA>
# Initialize an empty data frame to store the merged results
merged_results <- data.frame()
lengths <- c()

for (i in 2:ncol(transcriptomic)) {
  difference <- transcriptomic[[i]]
  dataset <- data.frame(
    gene_name = gene_name,
    difference = difference,
    logpval = NA,
    significant = ifelse(abs(difference) > 1.96, "+", NA)  
  )
  tf_activity_foot <- run_footprint_based_analysis_Francesco(omic_data = dataset, 
                                                             analysis = 'tfea', 
                                                             organism = 'human', 
                                                             reg_minsize = 10, 
                                                             exp_sign = FALSE,
                                                             collectri = TRUE,
                                                             hypergeom_corr = TRUE,
                                                             GO_annotation = TRUE, 
                                                             correct_proteomics = FALSE,
                                                             significant = FALSE
  )
  
  cell <- colnames(transcriptomic)[i]  # Get the column name dynamically
  print(cell)
  nes_dataset <- data.frame(
    gene_name = tf_activity_foot$gene_name  # or `gene_name` if available directly
  )
  
  # Assign the NES values to a column with the name stored in `cell`
  nes_dataset[[cell]] <- tf_activity_foot$weightedNES  # Correct way to dynamically assign column name
  
  # Perform a full join with the existing merged_results
  if (nrow(merged_results) == 0) {
    merged_results <- nes_dataset  # Initialize if empty
  } else {
    merged_results <- full_join(merged_results, nes_dataset, by = "gene_name")  # Full join
  }
  lengths <- c(lengths, nrow(merged_results)) 
  print(lengths)
}



write.csv(merged_results, file = "pazienti/merged_results_trans.csv", row.names = FALSE)




