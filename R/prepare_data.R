# 2.Data prep - separate for replicates input vs averaged input

#' prepare_data
#'
#' @param dataset input dataframe name
#' @param sample_name input individual sample name
#'
#' @return a formatted dataframe with time, F_t, and sample name
#' @export
#'
#' @examples
#' in_data <- data.frame(Sample = c("001","002","003"),
#' c_0 = c(40,50,60), c_15=c(25,27,30), c_30=c(20,21,22), c_60=c(15,17,18), c_120=c(12,11,10))
#' out_data <- prepare_data(in_data, "001")
prepare_data <- function(dataset, sample_name) {
  time <- c(0, 15, 30, 60, 120)

  columns_to_use <-c("c_0", "c_15", "c_30", "c_60", "c_120")

  F_t <- dataset[dataset$Sample == sample_name, columns_to_use]

  data <- data.frame(
    time = rep(time, nrow(F_t)),
    F_t = as.vector(t(F_t)),
    Sample_Name = rep(sample_name, length(time) * nrow(F_t))
  )

  return(data)
}
