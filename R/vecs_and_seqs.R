# Helper function to convert seq back to a vector
seq2vec <- function(seq) {
  # Pre-compute each interval's elements into a list, then concatenate once
  pieces <- vector("list", length(seq))
  for (j in seq_along(seq)) {
    interval <- seq[[j]]
    start <- interval$start
    end <- interval$end
    direction <- interval$direction

    if (direction == 0) {
      pieces[[j]] <- start
    } else {
      # seq() handles both increasing (direction=1) and decreasing (direction=-1)
      pieces[[j]] <- base::seq.int(start, end)
    }
  }
  unlist(pieces, use.names = FALSE)
}


# Function to convert a vector to an interval representation
# Each sequence is represented as a list with:
# - start: starting value
# - end: ending value
# - direction: 1 for increasing, -1 for decreasing, 0 for constant
# Returns a list of intervals representing the sequence
#
# O(n) implementation: track only the current direction instead of a
# growing `current_seen` set, which was O(n²) due to `%in%` + `c()`.
# A new interval is started whenever the step changes direction or is
# non-unit (gap), which correctly handles palindromic BFB sequences.
vec2seq <- function(vector) {
  n <- length(vector)
  if (n == 0L) return(list())
  if (n == 1L) return(list(list(start = vector[1L], end = vector[1L], direction = 0L)))

  intervals  <- vector("list", n)
  k          <- 0L
  seg_start  <- vector[1L]
  seg_end    <- vector[1L]
  seg_dir    <- NA_integer_

  for (i in 2L:n) {
    step <- vector[i] - vector[i - 1L]

    if (step == 1L || step == -1L) {
      if (is.na(seg_dir)) {
        # First directional step in this run: establish direction
        seg_dir <- step
        seg_end <- vector[i]
      } else if (step == seg_dir) {
        # Continuing the same monotone run
        seg_end <- vector[i]
      } else {
        # Direction reversal (e.g. BFB palindrome fold-back):
        # close the current interval, start a new one at vector[i]
        k <- k + 1L
        intervals[[k]] <- list(start = seg_start, end = seg_end, direction = seg_dir)
        seg_start <- vector[i]
        seg_end   <- vector[i]
        seg_dir   <- NA_integer_
      }
    } else {
      # Non-unit step (gap): close and restart
      k <- k + 1L
      intervals[[k]] <- list(
        start     = seg_start,
        end       = seg_end,
        direction = if (is.na(seg_dir)) 0L else seg_dir
      )
      seg_start <- vector[i]
      seg_end   <- vector[i]
      seg_dir   <- NA_integer_
    }
  }

  # Close final interval
  k <- k + 1L
  intervals[[k]] <- list(
    start     = seg_start,
    end       = seg_end,
    direction = if (is.na(seg_dir)) 0L else seg_dir
  )

  intervals[seq_len(k)]
}
