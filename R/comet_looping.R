
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
#' Summarise list of comet-fit results into a tidy data-frame
#'
#' @param results_list Output from `loop_comet_data()`
#' @param method       "nls", "bayes", or "adaptive"
#' @return             Tibble: one row per sample, one block of
#'                     <estimate / low / high> columns per parameter
#' @export
#' @importFrom dplyr mutate relocate all_of
#' @importFrom purrr imap_dfr
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_replace str_replace_all
comet_list_to_df <- function(results_list, method = "adaptive") {

  purrr::imap_dfr(results_list, function(res, sample_name) {

    ## ------------------------------------------------------------- ##
    ##  NLS / Bayes – keep original output, add empty CI placeholders
    ## ------------------------------------------------------------- ##
    if (method %in% c("nls", "bayes")) {

      out <- tibble::tibble(
        Sample          = sample_name,
        F_amp           = res$params["F" ],
        F_amp_low       = NA_real_,
        F_amp_high      = NA_real_,
        S_amp           = res$params["S" ],
        S_amp_low       = NA_real_,
        S_amp_high      = NA_real_,
        k_f             = res$params["k_f"],
        k_f_low         = NA_real_,
        k_f_high        = NA_real_,
        k_s             = res$params["k_s"],
        k_s_low         = NA_real_,
        k_s_high        = NA_real_,
        Half_life_overall = res$half_life$overall,
        Half_life_fast    = res$half_life$fast_phase,
        Half_life_slow    = res$half_life$slow_phase
      )

      return(out)
    }

    ## ------------------------------------------------------------- ##
    ##  Adaptive – pull everything from `summarize_and_plot()`
    ## ------------------------------------------------------------- ##
    if (method == "adaptive") {

      sum_tbl <- summarize_and_plot(res)$summary_table |>
        dplyr::select(Parameter, Median, CI_lower, CI_upper) |>
        ## make syntactically safe names (t1/2_fast -> t1_2_fast etc.)
        dplyr::mutate(Parameter = stringr::str_replace_all(Parameter, "[/ ]", "_"))

      wide <- sum_tbl |>
        tidyr::pivot_wider(
          names_from  = Parameter,
          values_from = c(Median, CI_lower, CI_upper),
          names_glue  = "{Parameter}_{.value}"
        ) |>
        ## nicer suffixes
        dplyr::rename_with(~stringr::str_replace(., "_Median$",  ""),
                           dplyr::ends_with("_Median")) |>
        dplyr::rename_with(~stringr::str_replace(., "_CI_lower$", "_low"),
                           dplyr::ends_with("_CI_lower")) |>
        dplyr::rename_with(~stringr::str_replace(., "_CI_upper$", "_high"),
                           dplyr::ends_with("_CI_upper")) |>
        dplyr::mutate(Sample = sample_name) |>
        dplyr::relocate(Sample, .before = tidyselect::everything())

      ## Optional: order columns so every “estimate/low/high” block stays together
      param_order <- sum_tbl$Parameter
      wanted_cols <- unlist(
        purrr::map(param_order, \(p) c(p, paste0(p, "_low"), paste0(p, "_high")))
      )
      wide <- wide |>
        dplyr::relocate(dplyr::all_of(intersect(wanted_cols, names(wide))),
                        .after = Sample)

      return(wide)
    }

    stop("`method` must be one of 'nls', 'bayes', or 'adaptive'.")
  })
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

#' extract_model_scores
#'
#' Extract model_scores tables from all samples in an adaptive fitting result list.
#'
#' @param results_list A named list of results (from loop_comet_data(..., method = "adaptive")).
#'
#' @return A data frame combining model_scores for all samples, with `Sample` as the first column.
#' @export
#' @importFrom purrr map_dfr
#' @importFrom dplyr select relocate mutate
extract_model_scores <- function(results_list) {
  purrr::map_dfr(names(results_list), function(sample) {
    res <- results_list[[sample]]
    if (!is.null(res$model_scores)) {
      res$model_scores %>%
        dplyr::mutate(Sample = sample) %>%
        dplyr::relocate(Sample, .before = Model)  # Move Sample to the first column
    } else {
      NULL
    }
  })
}

