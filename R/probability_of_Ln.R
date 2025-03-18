#' Calculate Probability of Final Sequence Length After BFB Cycles
#'
#' @param Ln Final sequence length after n BFB cycles
#' @param L0 Initial sequence length before BFB cycles
#' @param n Number of BFB cycles
#'
#' @details Calculates theoretical probability of observing length Ln after n BFB cycles
#'   starting from initial length L0. Returns 0 if Ln exceeds maximum possible length.
#'
#' @return Probability value between 0 and 1
probability_of_Ln = function(Ln, L0, n) {
  if (Ln > 2^n * L0) {
    p = 0  # Impossible to exceed max theoretical length (2^n * L0)
  } else {
    p = 1 / (2^n * factorial(n-1) * L0) * log(2^(n+1) * L0 / Ln)^(n-1)
  }
  return(p)
}

expected_Ln_sd = function(L0, n) {
  L0 * sqrt((4/3)^n - 1)
}
