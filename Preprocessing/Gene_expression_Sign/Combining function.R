combine_datasets <- function(dataset1, dataset2, row_dataset1, row_dataset2, batch_start = 1, batch_end = nrow(dataset1)) {
  dataset3 <- dataset1[batch_start:batch_end, ]
  gene_names <- dataset2[[row_dataset2]]
  for (gene in gene_names) {
    if (!gene %in% colnames(dataset3)) {
      dataset3[[gene]] <- NA 
    }
  }
  
  for (i in batch_start:batch_end) {

    cell_name <- dataset1[[i,row_dataset1]]
    if (cell_name %in% colnames(dataset2)) {
      cell_data <- dataset2[[cell_name]]
      row_index <- which(dataset3[[row_dataset1]] == cell_name) 
      dataset3[row_index, gene_names] <- as.list(as.numeric(cell_data))
      
    }
  }
  na_percentage <- rowMeans(is.na(dataset3[, 14:ncol(dataset3)]))
  dataset3_cleaned <- dataset3[na_percentage < 1, ]
  return(dataset3_cleaned)
}





