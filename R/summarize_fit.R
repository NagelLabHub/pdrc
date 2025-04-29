#' summarize_and_plot
#'
#' Summarize fitted model parameters and half-lives, and generate a plot with annotations.
#'
#' @param res A list output from `fit_adaptive_spike_model()` (or similar structure).
#' @param data Optional; data used for fitting (if needed).
#'
#' @return A list containing the summary data and the plot object.
#' @export
#' 
#' @import brms
#' @importFrom ggthemes theme_few
#' @importFrom dplyr group_by summarise
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_line labs annotate
#' @importFrom purrr map_chr compact

summarize_and_plot <- function(res, data = NULL) {
  if (is.null(data) && !is.null(res$data)) {
    data <- res$data
  }
  if (is.null(data)) stop("Data must be provided or stored inside 'res' object.")

  fit <- res$best_model
  best_model_name <- res$best_model_name
  looic_value <- res$model_scores %>%
    dplyr::filter(Model == best_model_name) %>%
    dplyr::pull(LOOIC)

  samples <- brms::as_draws_df(fit)

  # Helper to extract parameters
  extract_or_na <- function(samples, name) {
    varname <- paste0("b_", name, "_Intercept")
    if (varname %in% names(samples)) {
      return(samples[[varname]])
    } else {
      return(rep(NA, nrow(samples)))
    }
  }

  # Extract amplitudes and rates
  F_amp <- extract_or_na(samples, "F")
  S_amp <- extract_or_na(samples, "S")
  A_amp <- extract_or_na(samples, "A")

  k1 <- extract_or_na(samples, "k1")
  k2 <- extract_or_na(samples, "k2")
  k3 <- extract_or_na(samples, "k3")

  # Half-life calculation
  half_life <- function(k) log(2) / k

  t1_2_k1 <- half_life(k1)
  t1_2_k2 <- half_life(k2)
  t1_2_k3 <- half_life(k3)

  # Overall decay half-life
  overall_decay_half_life <- if (all(is.na(S_amp))) {
    t1_2_k2
  } else {
    total_F_S <- F_amp + S_amp
    weighted_k <- (F_amp * k2 + S_amp * k3) / total_F_S
    log(2) / weighted_k
  }

  accumulation_included <- !all(is.na(A_amp))

  # Peak calculations (only if accumulation included)
  if (accumulation_included) {
    peak_times <- sapply(1:nrow(samples), function(i) {
      optimize(function(t) {
        -(F_amp[i] * exp(-k2[i] * t) +
            ifelse(!is.na(S_amp[i]), S_amp[i] * exp(-k3[i] * t), 0) +
            A_amp[i] * (1 - exp(-k1[i] * t)))
      }, c(0, max(data$time)))$minimum
    })

    peak_values <- sapply(1:nrow(samples), function(i) {
      F_amp[i] * exp(-k2[i] * peak_times[i]) +
        ifelse(!is.na(S_amp[i]), S_amp[i] * exp(-k3[i] * peak_times[i]), 0) +
        A_amp[i] * (1 - exp(-k1[i] * peak_times[i]))
    })
  }

  # Build summary table
  parameters <- list(
    if (accumulation_included) list(name = "A_amp", values = A_amp),
    list(name = "F_amp", values = F_amp),
    if (!all(is.na(S_amp))) list(name = "S_amp", values = S_amp),
    if (accumulation_included) list(name = "k1_accum", values = k1),
    list(name = "k2_fast", values = k2),
    if (!all(is.na(S_amp))) list(name = "k3_slow", values = k3),
    if (accumulation_included) list(name = "t1/2_accum", values = t1_2_k1),
    list(name = "t1/2_fast", values = t1_2_k2),
    if (!all(is.na(S_amp))) list(name = "t1/2_slow", values = t1_2_k3),
    list(name = "Overall_decay_t1/2", values = overall_decay_half_life),
    if (accumulation_included) list(name = "Peak_time", values = peak_times),
    if (accumulation_included) list(name = "Peak_value", values = peak_values)
  ) %>% purrr::compact()

  summary_table <- tibble::tibble(
    Parameter = purrr::map_chr(parameters, "name"),
    Median = purrr::map_chr(parameters, ~ format(median(.x$values, na.rm = TRUE), scientific = FALSE, digits = 3)),
    CI_lower = purrr::map_chr(parameters, ~ format(quantile(.x$values, 0.025, na.rm = TRUE), scientific = FALSE, digits = 3)),
    CI_upper = purrr::map_chr(parameters, ~ format(quantile(.x$values, 0.975, na.rm = TRUE), scientific = FALSE, digits = 3))
  )

  # Predictions for plotting
  new_data <- tibble(time = seq(min(data$time), max(data$time), length.out = 1000))
  pred <- fitted(fit, newdata = new_data, re_formula = NA)
  new_data$predicted <- pred[, 1]
  new_data$lower <- pred[, 3]
  new_data$upper <- pred[, 4]

  summary_data <- data %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(avg = mean(F_t), sem = sd(F_t) / sqrt(n()))

  sample_name <- unique(data$Sample_Name)
  sample_line <- paste0("Sample: ", sample_name, "\n")
  looic_line <- sprintf("LOOIC: %.2f\n", looic_value)

  # Dynamic depending on best model
  if (best_model_name == "fast") {
    fast_line <- sprintf("Decay t½: %.1f min\n", median(t1_2_k2, na.rm = TRUE))
    accum_line <- ""
    slow_line <- ""
    overall_line <- ""
  } else if (best_model_name == "fast_accum") {
    accum_line <- sprintf("Accum t½: %.1f min\n", median(t1_2_k1, na.rm = TRUE))
    fast_line <- sprintf("Decay t½: %.1f min\n", median(t1_2_k2, na.rm = TRUE))
    slow_line <- ""
    overall_line <- ""
  } else if (best_model_name %in% c("fast_slow", "fast_slow_accum")) {
    accum_line <- if (accumulation_included) sprintf("Accum t½: %.1f min\n", median(t1_2_k1, na.rm = TRUE)) else ""
    fast_line <- sprintf("Fast t½: %.1f min\n", median(t1_2_k2, na.rm = TRUE))
    slow_line <- if (!all(is.na(S_amp))) sprintf("Slow t½: %.1f min\n", median(t1_2_k3, na.rm = TRUE)) else ""
    overall_line <- sprintf("Overall decay t½: %.1f min", median(overall_decay_half_life, na.rm = TRUE))
  } else {
  # fallback just in case
    fast_line <- sprintf("Fast t½: %.1f min\n", median(t1_2_k2, na.rm = TRUE))
    accum_line <- ""
    slow_line <- ""
    overall_line <- sprintf("Overall decay t½: %.1f min", median(overall_decay_half_life, na.rm = TRUE))
  }

  if (accumulation_included && best_model_name != "fast") {
    peak_time_med <- median(peak_times, na.rm = TRUE)
    peak_value_med <- median(peak_values, na.rm = TRUE)
    peak_line <- paste0(
      "Peak Time: ", format(peak_time_med, scientific = FALSE, digits = 3), " min\n",
      "Peak Value: ", format(peak_value_med, scientific = FALSE, digits = 3), "\n"
    )
  } else {
    peak_line <- ""
  }

  # Assemble annotation text
  annotation_text <- paste0(
    sample_line,
    looic_line,
    peak_line,
    accum_line,
    fast_line,
    slow_line,
    overall_line
  )

  # Create plot
  plot <- ggplot2::ggplot(summary_data, aes(time, avg)) +
    ggplot2::geom_point(color = "blue") +
    ggplot2::geom_errorbar(aes(ymin = avg - sem, ymax = avg + sem), width = 1, color = "blue") +
    ggplot2::geom_line(data = new_data, aes(time, predicted), linetype = "dashed") +
    ggthemes::theme_few() +
    ggplot2::labs(x = "Repair time (min)", y = "% DNA in tail") +
    ggplot2::annotate("text", x = Inf, y = Inf, label = annotation_text,
             hjust = 1.1, vjust = 1.2, size = 4)

  return(list(
    summary_table = summary_table,
    plot = plot
  ))
}
