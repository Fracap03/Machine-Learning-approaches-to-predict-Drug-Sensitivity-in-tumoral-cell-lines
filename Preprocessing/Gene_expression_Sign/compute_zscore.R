calculate_zscore_col <- function(mat, metric = 'median') {
  # Verifica che metric sia 'median' o 'mean'
  if (!metric %in% c("median", "mean")) {
    stop("Il parametro metric deve essere 'median' o 'mean'")
  }
  
  scaled <- matrix(nrow = nrow(mat), ncol = ncol(mat))
  
  for (i in seq_len(ncol(mat))) {
    column <- mat[, i]
    
    if (metric == "median") {
      centered <- column - median(column, na.rm = TRUE)
    } else {
      centered <- column - mean(column, na.rm = TRUE)
    }
    
    scaled[, i] <- centered / sd(centered, na.rm = TRUE)
  }
  
  colnames(scaled) <- colnames(mat)
  rownames(scaled) <- rownames(mat)
  
  return(scaled)
}

# Function to calculate z-score for each ROW of a matrix

calculate_zscore_row <- function(mat, metric = 'median') {
  centering_func <- ifelse(metric == "median", median, mean)
  centered <- t(apply(mat, 1, function(x) x - centering_func(x, na.rm = TRUE)))
  # centered <- t(apply(mat, 1, function(x) x - median(x, na.rm = TRUE)))
  scaled <- t(apply(centered, 1, function(x) x / sd(x, na.rm = TRUE)))
  return(scaled)

}

#' compute_zscore
#'
#' @param omic_matrix omic data matrix
#' @param by string, 'row' or 'column' according to the orientation of the z-score computation
#' @param metric string, 'median' (default) or 'mean' according to the metric for centering data
#'
#' @return omic data matrix with zscore
#' @export
#'
#' @examples
compute_zscore <- function(omic_matrix, by, metric = 'median'){

  if(by == 'row'){
    if(metric %in% c('median', 'mean')){
      omic_z_mat <- calculate_zscore_row(omic_matrix, metric)
    }else{
      stop('Please provide a valid metric: mean or median')
    }

  }else if( by == 'column'){
    if(metric %in% c('median', 'mean')){
      omic_z_mat <- calculate_zscore_col(omic_matrix, metric)
    }else{
      stop('Please provide a valid metric: mean or median')
    }
  }else{
    stop('Please provide a valid type of analysis: column or row')
  }

  return(omic_z_mat)
}


