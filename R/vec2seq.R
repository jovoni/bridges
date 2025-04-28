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
  current_seen <- c(vector[1])
  current_direction <- NA

  for (i in 2:length(vector)) {
    # Check if next value is contiguous
    contiguous <- (vector[i] == vector[i-1] + 1) || (vector[i] == vector[i-1] - 1)
    # Check if value already seen
    already_seen <- vector[i] %in% current_seen

    if (!contiguous || already_seen) {
      # Close current interval
      intervals <- append(intervals, list(list(
        start = current_start,
        end = current_value,
        direction = if(is.na(current_direction)) 0 else current_direction
      )))
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
  intervals <- append(intervals, list(list(
    start = current_start,
    end = current_value,
    direction = if(is.na(current_direction)) 0 else current_direction
  )))

  return(intervals)
}
