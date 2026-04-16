
get_seq_length <- function(sequence) {
  n <- length(sequence)
  if (n == 0L) return(0L)
  total <- 0L
  for (j in seq_len(n)) {
    iv    <- sequence[[j]]
    total <- total + abs(iv$end - iv$start) + 1L
  }
  total
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
