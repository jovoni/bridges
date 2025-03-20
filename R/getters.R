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
