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

  # ───────────────────────── 0. checks ─────────────────────────
  if (is.null(data) && !is.null(res$data)) data <- res$data
  if (is.null(data)) stop("Data must be provided or embedded in `res`.")

  fit         <- res$best_model
  best_model  <- res$best_model_name
  looic_value <- dplyr::filter(res$model_scores, Model == best_model)$LOOIC

  draws <- brms::as_draws_df(fit)

  pull_d <- function(df, nm) {
    v <- paste0("b_", nm, "_Intercept")
    if (v %in% names(df)) df[[v]] else rep(NA_real_, nrow(df))
  }

  # ────────────────── 1. raw posterior vectors ─────────────────
  F_amp <- pull_d(draws, "F")
  S_amp <- pull_d(draws, "S")
  A_amp <- pull_d(draws, "A")
  k1    <- pull_d(draws, "k1")
  k2    <- pull_d(draws, "k2")   # fast decay
  k3    <- pull_d(draws, "k3")   # slow decay (may be NA)

  accumulation  <- !all(is.na(A_amp))
  slow_included <- !all(is.na(S_amp))

  # ───────────── 2. enforce k2 > k3 for two-phase models ───────
  if (slow_included) {
    swap_idx <- which(!is.na(k3) & k3 > k2)
    if (length(swap_idx)) {
      temp <- k2[swap_idx];  k2[swap_idx] <- k3[swap_idx];  k3[swap_idx] <- temp
      temp <- F_amp[swap_idx]; F_amp[swap_idx] <- S_amp[swap_idx]; S_amp[swap_idx] <- temp
    }
  }

  hl   <- function(k) log(2) / k
  t12_1 <- hl(k1)
  t12_2 <- hl(k2)
  t12_3 <- hl(k3)

  # Overall half-life (weighted rate) ‒ now safe after swap
  overall_t12 <- if (slow_included) {
    w_k <- (F_amp * k2 + S_amp * k3) / (F_amp + S_amp)
    hl(w_k)
  } else t12_2

  # ───────────── 3. peak stats (only if accumulation) ──────────
  if (accumulation) {
    peak_times  <- vapply(seq_along(k2), function(i) {
      optimise(function(t) {
        -((F_amp[i] + A_amp[i]*(1-exp(-k1[i]*t))) * exp(-k2[i]*t) +
            if (!is.na(S_amp[i])) S_amp[i]*exp(-k3[i]*t) else 0)
    }, c(0, max(data$time)))$minimum }, numeric(1))

    peak_values <- vapply(seq_along(k2), function(i) {
      (F_amp[i] + A_amp[i]*(1-exp(-k1[i]*peak_times[i]))) * exp(-k2[i]*peak_times[i]) +
        if (!is.na(S_amp[i])) S_amp[i]*exp(-k3[i]*peak_times[i]) else 0
    }, numeric(1))
  }

  # ───────────── 4. assemble summary table ─────────────────────
  blocks <- list(
    list(name = "F_amp", values = F_amp),
    if (slow_included) list(name = "S_amp", values = S_amp),
    if (accumulation)  list(name = "A_amp", values = A_amp),

    list(name = "k2_fast", values = k2),
    if (slow_included) list(name = "k3_slow", values = k3),
    if (accumulation)  list(name = "k1_accum", values = k1),

    list(name = "t1/2_fast", values = t12_2),
    if (slow_included) list(name = "t1/2_slow", values = t12_3),
    if (accumulation)  list(name = "t1/2_accum", values = t12_1),

    list(name = "Overall_decay_t1/2", values = overall_t12),
    if (accumulation)  list(name = "Peak_time",  values = peak_times),
    if (accumulation)  list(name = "Peak_value", values = peak_values)
  ) |> purrr::compact()

  summary_table <- tibble::tibble(
    Parameter = purrr::map_chr(blocks, "name"),
    Median    = purrr::map_chr(blocks, ~ format(median(.x$values, na.rm = TRUE),
                                                digits = 3, scientific = FALSE)),
    CI_lower  = purrr::map_chr(blocks, ~ format(quantile(.x$values, 0.025, na.rm = TRUE),
                                                digits = 3, scientific = FALSE)),
    CI_upper  = purrr::map_chr(blocks, ~ format(quantile(.x$values, 0.975, na.rm = TRUE),
                                                digits = 3, scientific = FALSE))
  )

  # ───────────── 5. prediction grid & plot data ────────────────
  grid <- tibble::tibble(time = seq(min(data$time), max(data$time), length.out = 1000))
  fit_pred <- fitted(fit, newdata = grid, re_formula = NA)
  grid$predicted <- fit_pred[, 1]

  summary_data <- data |> dplyr::group_by(time) |>
    dplyr::summarise(avg = mean(F_t),
                     sem = stats::sd(F_t)/sqrt(dplyr::n()),
                     .groups = "drop")

  # ───────────── 6. dynamic annotation block ───────────────────
  ann <- switch(best_model,
    "fast" = sprintf("Decay t½: %.1f min\n", median(t12_2, na.rm = TRUE)),
    "fast_accum" = paste0(
        sprintf("Accum t½: %.1f min\n", median(t12_1, na.rm = TRUE)),
        sprintf("Decay t½: %.1f min\n",  median(t12_2, na.rm = TRUE))),
    "fast_slow" = paste0(
        sprintf("Fast t½: %.1f min\n", median(t12_2, na.rm = TRUE)),
        sprintf("Slow t½: %.1f min\n", median(t12_3, na.rm = TRUE)),
        sprintf("Overall decay t½: %.1f min\n", median(overall_t12, na.rm = TRUE))),
    "fast_slow_accum" = paste0(
        sprintf("Accum t½: %.1f min\n", median(t12_1, na.rm = TRUE)),
        sprintf("Fast t½: %.1f min\n",  median(t12_2, na.rm = TRUE)),
        sprintf("Slow t½: %.1f min\n",  median(t12_3, na.rm = TRUE)),
        sprintf("Overall decay t½: %.1f min\n", median(overall_t12, na.rm = TRUE)))
  )

  peak_block <- if (accumulation && best_model != "fast") {
    sprintf("Peak Time: %.2f min\nPeak Value: %.1f\n",
            median(peak_times,  na.rm = TRUE),
            median(peak_values, na.rm = TRUE))
  } else ""

  annotation <- paste0(
    "Sample: ", unique(data$Sample_Name), "\n",
    sprintf("LOOIC: %.2f\n", looic_value),
    peak_block,
    ann
  )

  # ───────────── 7. plot ───────────────────────────────────────
  plt <- ggplot2::ggplot(summary_data, ggplot2::aes(time, avg)) +
    ggplot2::geom_point(colour = "blue") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = avg-sem, ymax = avg+sem),
                           colour = "blue", width = 1) +
    ggplot2::geom_line(data = grid,
                       ggplot2::aes(time, predicted),
                       linetype = "dashed") +
    ggthemes::theme_few() +
    ggplot2::labs(x = "Repair time (min)", y = "% DNA in tail") +
    ggplot2::annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.2,
                      label = annotation, size = 4)

  list(summary_table = summary_table,
       plot          = plt)
}

