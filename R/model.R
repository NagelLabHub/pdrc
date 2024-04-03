# 1. biphasic exponential decay model without plateau

#' model
#'
#' @param t time
#' @param F fast phase span
#' @param S slow phase span
#' @param k_f fast phase par
#' @param k_s slow phase par
#'
#' @return to be used internally
#' @export
#'
model <- function(t, F, S, k_f, k_s) {
  F * exp(-k_f * t) + S * exp(-k_s * t)
}
