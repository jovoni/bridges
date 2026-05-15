# Helper function to convert seq back to a vector
seq2vec <- function(seq) {
  if (length(seq) == 0L) return(integer(0L))
  seq2vec_cpp(seq)
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
  vec2seq_cpp(as.integer(vector))
}
