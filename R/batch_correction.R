# 1. batch correction

#' batch_correction
#'
#' @import sva
#'
#' @param dataset input dataframe storing standardized data for batch correction - make sure the rows are appropriately ordered
#' @param pheno input dataframe storing phenotype data for batch correction - make sure the rows are appropriately ordered
#' @param covariate input covariate names in pheno whose variations need to be retained
#' @param batch input batch column name in pheno
#'
#' @return a dataframe with batch corrected data after ComBat
#' @export
#'
#' @examples
#' dat_og <- data.frame(UG = c(0.05,0.08,0.10), NHEJ = c(23,16,19), HR = c(2,4,6))
#' dat_std <- RE_to_zscore(dat_og, c("UG"))
#' pheno <- data.frame(Sample = c("001","002","003"), age = c(34,36,57), sex = c("M","F","M"), batch = c("A","B","A"))
#' dataset_batch <- batch_correction(dat_std, pheno, batch="batch")
#'

batch_correction <- function(dataset, pheno, covariate=NULL, batch) {

  # Check if the batch column exists in the pheno dataframe
  if (!(batch %in% names(pheno))) {
    stop("The batch column does not exist in the phenotype data.")
  }

  # Extract the batch column from the pheno dataframe
  batch_col <- pheno[[batch]]

  # Convert input dataset to matrix
  t_dataset <- t(dataset)

  # Perform batch correction using ComBat
  if (is.null(covariate)) {
    t_corrected_data <- ComBat(dat = as.matrix(t_dataset), batch = batch_col, par.prior = TRUE, prior.plots = FALSE)
  } else {
    formula <- paste("~", paste(covariate, collapse = " + "))
    model <- model.matrix(as.formula(formula), data = pheno)
    t_corrected_data <- ComBat(dat = as.matrix(t_dataset), batch = batch_col, mod = model, par.prior = TRUE, prior.plots = FALSE)
  }

  corrected_data = as.data.frame(t(t_corrected_data))

  return(corrected_data)
}
