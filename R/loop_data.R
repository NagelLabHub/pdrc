# 6. loop through all samples and apply the model fitting in step 5.

#' loop_data
#'
#' @param dataset_name full dataset name
#' @param method choose between "nls" and "brms" for model fitting algrithm
#'
#' @return a large list of result lists ordered by sample name
#' @export
#'
#' @examples
#' in_data <- data.frame(Sample = c("001","002","003"),
#' c_0 = c(40,50,60), c_15=c(25,27,30), c_30=c(20,21,22), c_60=c(15,17,18), c_120=c(12,11,10))
#' sample_list <- loop_data(in_data, "nls")
loop_data <- function(dataset_name, method) {

  results_list <- list()

  for (sample_name in unique(dataset_name$Sample)) {
    data <- prepare_data(dataset_name, sample_name)

    if (method == "nls") {
      initial_params_list <- list(
        c(F = data$F_t[1], S = 0, k_f = 0.04, k_s = 0.02),
        optimize_parameters(data, c(F = data$F_t[1], S = 0, k_f = 0.04, k_s = 0.02)),
        c(F = data$F_t[1], S = 0, k_f = 0.4, k_s = 0.04),
        c(F = 0.5 * data$F_t[1], S = 0.5 * data$F_t[1], k_f = 0.4, k_s = 0.04),
        c(F = data$F_t[1], S = 0, k_f = 0.1, k_s = 0.01),
        c(F = 0.5 * data$F_t[1], S = 0.5 * data$F_t[1], k_f = 0.1, k_s = 0.01),
        c(F = data$F_t[1], S = 0, k_f = 0.04, k_s = 0)
      )

      fitted <- FALSE
      for (initial_params in initial_params_list) {
        tryCatch({
          fit_results <- biphasic_fit_nls(data, initial_params)

          results_list[[sample_name]] <- fit_results

          fitted <- TRUE
          break  # Exit the loop if successful fit
        }, error = function(e) {
          # Print error message and try next set of initial_params
          cat(paste("Error processing sample:", sample_name, "- Trying next set of initial parameters.\n"))
        })

        if (fitted) break  # Exit the loop if fitted successfully
      }

      if (!fitted) {
        # If none of the initial_params combinations worked, add an empty row
        cat(paste("Error processing sample:", sample_name, "- None of the initial parameters combinations worked.\n"))
      }
    } else if (method == "brms") {
      # Suppress output of biphasic_fit_bayesian()
      suppressMessages({
        fit_results <- biphasic_fit_bayesian(data)
      })

      # Store the results
      results_list[[sample_name]] <- fit_results
    } else {
      stop("Invalid method argument. Method must be either 'nls' or 'brms'.")
    }
  }

  # Return the results_list
  return(results_list)
}
