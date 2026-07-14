#' Compute the Length of a Sequence
#'
#' This function calculates the length of a sequence by summing the lengths of its intervals.
#'
#' @param sequence A list of intervals, where each interval is expected to have \code{start}, \code{end}, and \code{direction} components.
#'
#' @return An integer representing the total length of the sequence.
#' @keywords internal
get_seq_length = function(sequence) {
  seq_length_cpp(sequence)
}
