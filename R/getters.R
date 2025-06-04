
get_seq_length = function(sequence) {
  if (length(sequence) == 0) return(0)
  L = sum(sapply(sequence, function(interval) {
    return(abs(interval$end - interval$start) + 1)
    # if (interval$direction == 0) {
    #   return(1)  # A constant interval contributes 1 element
    # } else {
    #   return(abs(interval$end - interval$start) + 1)  # Increasing or decreasing interval length
    # }
  }))
  L
}


get_colors = function(set) {
  if (set == "CN") {
    colors <- structure(
      c(
        "#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
        "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA"
      ),
      names = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+")
    )
  } else if (set == "gainloss") {
    colors = structure(
      c("firebrick", "steelblue"),
      names = c("Gain", "Loss")
    )
  } else if (set == "copy") {
    colors <- structure(
      c(
        "#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
        "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA"
      ),
      names = 0:11
    )
  } else {
    stop("Error: color set not known. Available sets are 'gainloss' and 'CN'")
  }

  colors
}
