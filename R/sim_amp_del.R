
sim_amp_del <- function(sequence, operation = "dup") {
  # Get sequence properties
  L <- get_seq_length(sequence)

  # Ensure there's a sequence to modify
  if (L == 0) {
    return(sequence)
  }

  vec <- seq2vec(sequence)

  find_consecutive_subvectors <- function(vec) {
    diffs <- diff(vec)
    runs <- rle(abs(diffs) == 1)
    lengths <- runs$lengths
    values <- runs$values

    start_indices <- cumsum(c(1, lengths[-length(lengths)]))
    good_starts <- start_indices[values]
    good_lengths <- lengths[values] + 1  # +1 because diff loses one element

    # Build a list of consecutive subvectors
    subvecs <- lapply(seq_along(good_starts), function(i) {
      idx_start <- good_starts[i]
      idx_end <- idx_start + good_lengths[i] - 1
      vec[idx_start:idx_end]
    })

    return(list(subvecs = subvecs, good_starts = good_starts))
  }



  subvectors <- find_consecutive_subvectors(vec)
  subvecs <- subvectors$subvecs
  good_starts <- subvectors$good_starts

  if (length(subvecs) == 0) {
    subvecs = lapply(vec, function(x){x})
    good_starts = 1:length(vec)
    #stop("No consecutive subsequence found!")
  }

  idx <- sample(seq_along(subvecs), 1)
  subvec <- subvecs[[idx]]
  start <- good_starts[idx]

  # Sample a section inside the selected subvector
  if (length(subvec) == 1) {
    event_lims <- c(1, 1)
  } else {
    event_lims <- sort(sample(1:length(subvec), size = 2))
  }

  s <- subvec[event_lims[1]:event_lims[2]]

  absolute_start <- start + event_lims[1] - 1
  absolute_end <- start + event_lims[2] - 1

  # Safe slicing
  before <- if (absolute_start > 1) vec[1:(absolute_start - 1)] else integer(0)
  after <- if (absolute_end < length(vec)) vec[(absolute_end + 1):length(vec)] else integer(0)

  # Apply operation
  if (operation == "dup") {
    new_vec <- c(before, s, s, after)
  } else if (operation == "del") {
    new_vec <- c(before, after)
    if (length(new_vec) == 0) {
      message("skipping deletion because of length zero")
      return(sequence)
    }
  } else {
    stop("operation not recognized!")
  }

  return(vec2seq(new_vec))
}

sim_wgd = function(sequence) {
  vec = seq2vec(sequence)
  wgd_vec = c(vec, vec)
  vec2seq(wgd_vec)
}
