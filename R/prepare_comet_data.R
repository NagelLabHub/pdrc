# 1. data prep for single individual dataset

#' prepare_comet_data
#'
#' @param dataset input dataframe storing comet kinetic data, requires columns: Sample and F_t values at different time points
#' @param sample_name input individual sample name
#'
#' @return a formatted dataframe with time, F_t, and sample name
#' @export
#'
#' @examples
#' in_data <- data.frame(Sample = c("001","002","003"),
#' c_0 = c(40,50,60), c_15=c(25,27,30), c_30=c(20,21,22), c_60=c(15,17,18), c_120=c(12,11,10))
#' out_data <- prepare_comet_data(in_data, "001")

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
