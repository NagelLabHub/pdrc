
#' batch_evaluation
#' 
#' Evaluate batch effects across different batches
#'
#' @param dataset input dataframe storing data for batch effect evaluation - make sure the rows are appropriately ordered
#' @param pheno input dataframe storing phenotype data - make sure the rows are appropriately ordered
#' @param batch_col column name in pheno dataframe that stores batch information
#'
#' @return a list of different evaluation metrics for batch effect: if input data is on original scale, use cv; if standardized, use sd
#' @export
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @import ggthemes

batch_evaluation <- function(dataset, pheno, batch_col) {

  # check if the batch column is in pheno dataframe
  if (!(batch_col %in% colnames(pheno))) {
    stop("The batch column must be in pheno dataframe")
  }

  data <- cbind(dataset, pheno) # join dataset and pheno for ease of evaluation
  variables <- colnames(dataset)

  data_long <- pivot_longer(data, cols = variables, names_to = "variable", values_to = "value")
  plot_batch <- ggplot(data_long, aes(x = as.factor(.data[[batch_col]]), y = value)) +
    geom_boxplot() + theme_few() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~ variable, nco=1, scales = "free_y") + xlab("Batch") + ylab("")

  # calculate mean and sd of each batch
  each_batch_stats <- data_long %>%
    group_by(.data[[batch_col]], variable) %>%
    summarize(mean = mean(value), sd = sd(value), cv = sd(value) / mean(value) )

  # calculate the sd of batch mean across all batches
  across_batch_stats <- each_batch_stats %>%
    group_by(variable) %>%
    summarize(sd_across_batch = sd(mean), cv_across_batch = sd(mean) / mean(mean))

  eval_list = list(plot_batch = plot_batch, each_batch_stats = each_batch_stats, across_batch_stats = across_batch_stats)

  return(eval_list)
}

#' batch_correction
#'
#' Batch correction using ComBat
#' 
#' @importFrom sva ComBat
#'
#' @param dataset input dataframe storing standardized data for batch correction - make sure the rows are appropriately ordered
#' @param pheno input dataframe storing phenotype data for batch correction - make sure the rows are appropriately ordered
#' @param covariate input covariate names in pheno whose variations need to be retained
#' @param batch input batch column name in pheno
#'
#' @return a dataframe with batch corrected data after ComBat
#' @export

batch_correction <- function(dataset, pheno, covariate=NULL, batch) {

  if (!(batch %in% names(pheno))) {
    stop("The batch column does not exist in the phenotype data.")
  }

  batch_col <- pheno[[batch]]

  t_dataset <- t(dataset)

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
