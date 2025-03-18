#' Simulate Breakage-Fusion-Bridge (BFB) Cycles
#'
#' @param L Initial sequence length (number of genomic segments)
#' @param n Number of BFB cycles to simulate
#' @param support Distribution type for breakpoint selection ("uniform" or "beta")
#' @param alpha Shape parameter for beta distribution (only used if support="beta")
#' @param beta Shape parameter for beta distribution (only used if support="beta")
#'
#' @details For each cycle:
#'   1) Random breakpoint is selected based on specified distribution
#'   2) Segment left of breakpoint is duplicated and inverted
#'   3) Inverted segment is appended
#'
#' @return Vector containing the final sequence after n BFB cycles
sim_n_bfb = function(L, n, support = "uniform", alpha = NULL, beta = NULL) {
  # Initialize sequence
  sequence = 1:L

  # Perform n BFB cycles
  for (i in 1:n) {
    Ln = length(sequence)
    # Select random breakpoint based on specified distribution
    if (support == "uniform") {
      bp_idx = sample(1:(2*Ln), 1)
    } else if (support == "beta") {
      if (is.null(alpha) || is.null(beta)) {
        stop("For beta distribution, both alpha and beta parameters must be provided")
      }
      tau = stats::rbeta(1, alpha, beta)
      bp_idx = max(1, round(tau * 2 * Ln))  # Ensure bp_idx is at least 1
    } else {
      stop("Unsupported distribution type. Use 'uniform' or 'beta'.")
    }

    # Fusiong
    fused_seq = c(sequence, rev(sequence))
    sequence = fused_seq[1:bp_idx]
  }

  return(sequence)
}
