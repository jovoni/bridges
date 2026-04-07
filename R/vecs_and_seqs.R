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
vec2seq <- function(vector) {
  if (length(vector) == 0) {
    return(list())
  }

  if (length(vector) == 1) {
    return(list(list(start = vector[1], end = vector[1], direction = 0)))
  }

  # Pre-allocate to avoid O(n²) append growth; trim at the end
  intervals <- vector("list", length(vector))
  n_intervals <- 0L

  current_start <- vector[1]
  current_value <- vector[1]
  current_seen <- c(vector[1])
  current_direction <- NA

  for (i in 2:length(vector)) {
    # Check if next value is contiguous
    contiguous <- (vector[i] == vector[i-1] + 1) || (vector[i] == vector[i-1] - 1)
    # Check if value already seen
    already_seen <- vector[i] %in% current_seen

    if (!contiguous || already_seen) {
      # Close current interval
      n_intervals <- n_intervals + 1L
      intervals[[n_intervals]] <- list(
        start = current_start,
        end = current_value,
        direction = if(is.na(current_direction)) 0 else current_direction
      )
      # Start new interval
      current_start <- vector[i]
      current_value <- vector[i]
      current_seen <- c(vector[i])
      current_direction <- NA
      next
    }

    # Update seen values
    current_seen <- c(current_seen, vector[i])

    # Set direction if not set
    if (is.na(current_direction)) {
      if (vector[i] == vector[i-1] + 1) {
        current_direction <- 1
      } else if (vector[i] == vector[i-1] - 1) {
        current_direction <- -1
      }
    }

    current_value <- vector[i]
  }

  # Add the last interval
  n_intervals <- n_intervals + 1L
  intervals[[n_intervals]] <- list(
    start = current_start,
    end = current_value,
    direction = if(is.na(current_direction)) 0 else current_direction
  )

  intervals[seq_len(n_intervals)]
}
