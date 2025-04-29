
#' loop_comet_data
#'
#' Loop through each sample in the dataset and fit the biphasic model.
#' 
#' @param dataset_name full dataset name
#' @param method Character; choose between "nls", "bayes", or "adaptive" (default = "adaptive").
#'
#' @return a large list of result lists ordered by sample name
#' @export
#' @importFrom dplyr group_by summarise count
#' @importFrom purrr map
loop_comet_data <- function(dataset_name, method = "adaptive") {
  results_list <- list()

  for (sample_name in unique(dataset_name$Sample)) {
    data <- prepare_comet_data(dataset_name, sample_name)

    if (method == "nls") {
      initial_params <- c(F = data$F_t[1], S = 0, k_f = 0.04, k_s = 0.01)
      fit_results <- fit_biphasic_nls(data, initial_params)
    } else if (method == "bayes") {
      fit_results <- fit_biphasic_bayes(data)
    } else if (method == "adaptive") {
      fit_results <- fit_adaptive_spike_model(data)
    } else {
      stop("Invalid method: must be 'nls', 'bayes', or 'adaptive'.")
    }

    results_list[[sample_name]] <- fit_results
  }

  # If method = "adaptive", summarize model selection
  if (method == "adaptive") {
    model_names <- purrr::map_chr(results_list, "best_model_name")
    model_summary <- dplyr::count(tibble::tibble(Model = model_names), Model) %>%
      dplyr::mutate(Percent = n / sum(n) * 100)

    message("Summary of selected models across samples:")
    print(model_summary)
  }

  return(results_list)
}

#' comet_list_to_df
#'
#' Summarize list of comet fitting results into a single data frame.
#'
#' @param results_list A list generated from `loop_comet_data()`.
#' @param method Character; method used to fit: "nls", "bayes", or "adaptive".
#'
#' @return A data frame summarizing fitted parameters per sample.
#' @export
#' @importFrom purrr map_dfr
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr mutate
comet_list_to_df <- function(results_list, method = "adaptive") {
  rows <- purrr::map_dfr(names(results_list), function(sample_name) {
    res <- results_list[[sample_name]]

    if (method == "nls" || method == "bayes") {
      # Direct extraction for nls and bayes
      tibble::tibble(
        Sample = sample_name,
        F_amp = res$params["F"],
        S_amp = res$params["S"],
        k_f = res$params["k_f"],
        k_s = res$params["k_s"],
        Half_life_overall = res$half_life$overall,
        Half_life_fast = res$half_life$fast_phase,
        Half_life_slow = res$half_life$slow_phase
      )
    } else if (method == "adaptive") {
      # Use summarize_and_plot() to extract full summary for adaptive models
      summary_result <- summarize_and_plot(res)  # uses res$data internally

      summary_wide <- summary_result$summary_table %>%
        tidyr::pivot_wider(names_from = Parameter, values_from = Median) %>%
        dplyr::mutate(Sample = sample_name)

      summary_wide
    } else {
      stop("Invalid method: must be 'nls', 'bayes', or 'adaptive'.")
    }
  })

  return(rows)
}

#' extract_comet_plots
#'
#' Extract and optionally save plots from comet fitting results by first running summarize_and_plot.
#'
#' @param results_list A list of fitting results from `loop_comet_data()`.
#' @param pdf_file Optional; file name to save the plots as a single PDF.
#' @param ncol Integer; number of columns for plot layout (default = 2).
#' @param width Numeric; width of PDF (inches) if saving (default = 12).
#' @param height_per_row Numeric; height per row of plots (default = 4).
#'
#' @return A list of ggplot2 plot objects.
#' @export
#' @importFrom purrr map compact
#' @importFrom patchwork wrap_plots
#' @importFrom grDevices pdf dev.off
extract_comet_plots <- function(results_list, pdf_file = NULL, ncol = 2, width = 12, height_per_row = 4) {
  plot_list <- purrr::map(results_list, ~ {
    if (is.null(.x)) {
      NULL
    } else {
      summary <- summarize_and_plot(.x)
      summary$plot
    }
  }) %>%
    purrr::compact()

  num_plots <- length(plot_list)
  nrow <- ceiling(num_plots / ncol)
  pdf_height <- nrow * height_per_row

  if (!is.null(pdf_file)) {
    grDevices::pdf(pdf_file, width = width, height = pdf_height)
    p <- patchwork::wrap_plots(plot_list, ncol = ncol)
    print(p)  # <-- IMPORTANT
    grDevices::dev.off()
    message(sprintf("Saved plots to %s (%d plots, %d rows, height %d inches)", pdf_file, num_plots, nrow, pdf_height))
  }

  return(plot_list)
}

