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

  intervals <- list()
  current_start <- vector[1]
  current_value <- vector[1]
  current_direction <- NA

  # Set initial direction based on first two elements
  if (vector[2] == vector[1] + 1) {
    current_direction <- 1  # Increasing sequence
  } else if (vector[2] == vector[1] - 1) {
    current_direction <- -1  # Decreasing sequence
  } else {
    # If first two elements aren't contiguous, start a new segment
    intervals <- append(intervals, list(list(
      start = current_start,
      end = current_value,
      direction = 0
    )))
    current_start <- vector[2]
    current_value <- vector[2]
    current_direction <- NA
  }

  for (i in 2:length(vector)) {
    # Check if current elements are contiguous
    if (vector[i] == vector[i-1] + 1) {
      expected_direction <- 1
    } else if (vector[i] == vector[i-1] - 1) {
      expected_direction <- -1
    } else {
      # Non-contiguous values, end current segment
      if (!is.na(current_direction)) {
        intervals <- append(intervals, list(list(
          start = current_start,
          end = current_value,
          direction = current_direction
        )))
      }
      # Start a new segment with the current value
      current_start <- vector[i]
      current_value <- vector[i]
      current_direction <- NA
      next
    }

    # First element of a new segment
    if (is.na(current_direction)) {
      current_direction <- expected_direction
    }

    # If direction changes or numbers aren't contiguous, create a new interval
    if (expected_direction != current_direction) {
      intervals <- append(intervals, list(list(
        start = current_start,
        end = current_value,
        direction = current_direction
      )))
      current_start <- vector[i-1]
      current_direction <- expected_direction
    }

    current_value <- vector[i]
  }

  # Add the final interval
  intervals <- append(intervals, list(list(
    start = current_start,
    end = current_value,
    direction = if(is.na(current_direction)) 0 else current_direction
  )))

  return(intervals)
}
