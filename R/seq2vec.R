# Helper function to convert seq back to a vector
seq2vec <- function(seq) {
  result <- c()

  for (interval in seq) {
    start <- interval$start
    end <- interval$end
    direction <- interval$direction

    if (direction == 0) {
      # Single element
      result <- c(result, start)
    } else if (direction == 1) {
      # Increasing sequence
      result <- c(result, seq(start, end))
    } else if (direction == -1) {
      # Decreasing sequence
      result <- c(result, seq(start, end))
    }
  }

  return(result)
}
