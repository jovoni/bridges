#' Compute the Length of a Sequence
#'
#' This function calculates the length of a sequence by summing the lengths of its intervals.
#'
#' @param sequence A list of intervals, where each interval is expected to have \code{start}, \code{end}, and \code{direction} components.
#'
#' @return An integer representing the total length of the sequence.
get_seq_length = function(sequence) {
  L = sum(sapply(sequence, function(interval) {
    if (interval$direction == 0) {
      return(1)  # A constant interval contributes 1 element
    } else {
      return(abs(interval$end - interval$start) + 1)  # Increasing or decreasing interval length
    }
  }))
  L
}
