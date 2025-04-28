#' Simulate Breakage-Fusion-Bridge (BFB) Cycle for both daughters
#'
#' @param sequence Input sequence
#' @param support Distribution type for breakpoint selection ("uniform" or "beta")
#' @param alpha Shape parameter for beta distribution (only used if support="beta")
#' @param beta Shape parameter for beta distribution (only used if support="beta")
#'
#' @details Simulate left and right children from a BFB cycle using specified
#'   breakpoint selection distribution
#'
#' @return List containing left and right sequences
sim_bfb_left_and_right_sequences <- function(sequence, support = "uniform", alpha = NULL, beta = NULL) {
  # Calculate the total number of elements in the sequence, similar to the vectorized version
  L = get_seq_length(sequence)
  vec = seq2vec(sequence)
  bps = vec[diff(vec) == 0]

  # Select random breakpoint based on specified distribution
  # Ensure that bp_idx is different from L to obtain a proper bfb cycle
  bp_idx = L
  while (bp_idx %in% c(L, bps)) {
    if (support == "uniform") {
      bp_idx = sample(1:(2*L), 1)
    } else if (support == "beta") {
      if (is.null(alpha) || is.null(beta)) {
        stop("For beta distribution, both alpha and beta parameters must be provided")
      }
      tau = stats::rbeta(1, alpha, beta)
      bp_idx = max(1, round(tau * 2*L))  # Ensure bp_idx is at least 1
    } else {
      stop("Unsupported distribution type. Use 'uniform' or 'beta'.")
    }
  }

  # Initialize left and right sequences
  cut_seqs = cut_sequence(fuse_sequence(sequence), bp_idx)
  l_seq = cut_seqs$left_seq
  r_seq = reverse_sequence(cut_seqs$right_seq)

  # Return the left and right sequences
  return(list(l_seq = l_seq, r_seq = r_seq))
}

reverse_sequence <- function(sequence) {
  # Reverse the order of the intervals and swap start and end, flip direction
  reversed_seq = lapply(seq_along(sequence), function(i) {
    interval = sequence[[length(sequence) - i + 1]]
    list(start = interval$end,
         end = interval$start,
         direction = -interval$direction)
  })

  return(reversed_seq)
}

fuse_sequence <- function(sequence) {
  reversed_seq = reverse_sequence(sequence)
  fused_seq = c(sequence, reversed_seq)
  return(fused_seq)
}

cut_sequence <- function(sequence, cut_index) {
  # Initialize output sequences

  left_seq <- list()
  right_seq <- list()

  # Track the current position across the entire sequence
  current_length <- 0

  for (interval in sequence) {
    # Calculate the length of the current interval
    interval_length <- abs(interval$end - interval$start) + 1

    # If the cut point is before this interval, everything goes to the right
    if (current_length >= cut_index) {
      right_seq <- c(right_seq, list(interval))

      # If the cut point is after this interval, everything goes to the left
    } else if (current_length + interval_length <= cut_index) {
      left_seq <- c(left_seq, list(interval))

      # Otherwise, we split the interval
    } else {
      # How far into the current interval is the cut?
      cut_within <- cut_index - current_length

      # Handle different directions
      if (interval$direction == 1) {  # Increasing interval
        left_seq <- c(left_seq, list(list(
          start = interval$start,
          end = interval$start + cut_within - 1,
          direction = 1
        )))
        right_seq <- c(right_seq, list(list(
          start = interval$start + cut_within,
          end = interval$end,
          direction = 1
        )))
      } else if (interval$direction == -1) {  # Decreasing interval
        left_seq <- c(left_seq, list(list(
          start = interval$start,
          end = interval$start - cut_within + 1,
          direction = -1
        )))
        right_seq <- c(right_seq, list(list(
          start = interval$start - cut_within,
          end = interval$end,
          direction = -1
        )))
      } else {  # Constant interval
        left_seq <- c(left_seq, list(interval))
        right_seq <- c(right_seq, list(interval))
      }
    }

    # Update the position tracker
    current_length <- current_length + interval_length
  }

  return(list(left_seq = left_seq, right_seq = right_seq))
}
