
#' prepare_comet_data
#'
#' Prepare comet assay data for one sample.
#' 
#' @param dataset input dataframe storing comet kinetic data, requires columns: Sample and F_t values at different time points
#' @param sample_name input individual sample name
#'
#' @return a formatted dataframe with time, F_t, and sample name
#' @export

prepare_comet_data <- function(dataset, sample_name) {
  # Identify columns that start with 'c_'
  columns_to_use <- grep("^c_", names(dataset), value = TRUE)

  # Extract time vector from column names
  time <- as.numeric(gsub("c_", "", columns_to_use))

  # Subset the dataset
  F_t <- dataset[dataset$Sample == sample_name, columns_to_use]

  # Prepare the data frame
  data <- data.frame(
    time = rep(time, nrow(F_t)),
    F_t = as.vector(t(F_t)),
    Sample_Name = rep(sample_name, length(time) * nrow(F_t))
  )

  return(data)
}

#' RE_to_zscore
#' 
#' Log transform and z-score normalize FM-HCR data.
#' 
#' @param dataset input dataframe storing raw reporter expression data - make sure the data is appropriately ordered
#' @param var_inverse input pathway names with inverse relationship between raw reporter expression and pathway activity
#'
#' @return list with standardized data and mean/sd
#' @export

RE_to_zscore <- function(dataset, var_inverse = NULL) {
  dataset_log <- log(dataset)
  if (!is.null(var_inverse)) {
    dataset_log[, var_inverse] <- -dataset_log[, var_inverse]
  }
  mean_sd_priorscaling <- data.frame(
    mean = apply(dataset_log, 2, mean, na.rm = TRUE),
    sd = apply(dataset_log, 2, sd, na.rm = TRUE)
  )
  dataset_zscore <- scale(dataset_log)
  return(list(dataset_zscore = dataset_zscore, mean_sd_priorscaling = mean_sd_priorscaling))
}