# 1. batch effect evaluation

#' batch_evaluation
#'
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @import ggthemes
#'
#' @param dataset input dataframe storing data for batch effect evaluation - make sure the rows are appropriately ordered
#' @param pheno input dataframe storing phenotype data - make sure the rows are appropriately ordered
#' @param batch_col column name in pheno dataframe that stores batch information
#'
#' @return a list of different evaluation metrics for batch effect: if input data is on original scale, use cv; if standardized, use sd
#' @export
#'
#' @examples
#' dat_og <- data.frame(UG = c(0.05,0.08,0.10), NHEJ = c(23,16,19), HR = c(2,4,6))
#' list_std <- RE_to_zscore(dat_og, c("UG"))
#' dat_std <- list_std$dataset_zscore
#' pheno <- data.frame(Sample = c("001","002","003"),
#' age = c(34,36,57), sex = c("M","F","M"), batch = c("A","B","A"))
#' list = batch_evaluation(dataset, pheno, "batch")
#'

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

