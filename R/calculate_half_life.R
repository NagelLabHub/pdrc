# 3. half-life estimation for overall, fast, slow phases

#' calculate_half_life
#'
#' @param F fast phase span
#' @param S slow phase span
#' @param k_f fast phase exp decay par
#' @param k_s slow phase exp decay par
#'
#' @return a list with overall, fast phase, slow phase half-life estimates
#' @export
#'
#' @examples
#' half_life <- calculate_half_life(40,0,0.1,0.01)
calculate_half_life <- function(F, S, k_f, k_s) {
  half_life_eq <- function(t_half, F, S, k_f, k_s) {
    F * exp(-k_f * t_half) + S * exp(-k_s * t_half) - (F + S) / 2 }

  half_life_eq_fast <- function(t_half, k) {
    exp(-k * t_half) - 0.5 }

  half_life_eq_slow <- function(t_half, k) {
    exp(-k * t_half) - 0.5 }

  initial_interval <- c(0, 360)
  slow_phase_result <- 360

  # Use tryCatch to handle potential errors
  tryCatch({

    root <- uniroot(half_life_eq, interval = initial_interval, F = F, S = S, k_f = k_f, k_s = k_s)
    root_fast <- uniroot(half_life_eq_fast, interval = initial_interval, k = k_f)
    root_slow <- uniroot(half_life_eq_slow, interval = initial_interval, k = k_s)

    # Check if slow phase root was successful
    if (!is.null(root_slow$root)) {
      slow_phase_result <- root_slow$root
    }
  }, error = function(e) {
    # Do nothing or handle the error as needed
    # Omitting error output in this case
  })

  return(list(overall = ifelse(is.null(root$root), NA, root$root),
              fast_phase = ifelse(is.null(root_fast$root), NA, root_fast$root),
              slow_phase = slow_phase_result))
}
