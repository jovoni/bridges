#' Compute Distance Between Two Sets of Cells
#'
#' This function calculates the Euclidean distance between two sets of cells, optionally ordering them based on the proportion of the genome in a single-copy state.
#'
#' @param x1 A list representing the first set of cells. It must contain:
#'    - cells : A list of integer vectors, where each vector represents a cell sequence.
#'    - input_parameters : A list containing at least the field initial_sequence_length, specifying the expected length of each cell sequence.
#' @param x2 A list representing the second set of cells, structured identically to \code{x1}.
#' @param order Logical value indicating whether to reorder cells based on the proportion of genomic regions with exactly one copy.
#'
#' @return A numeric value representing the normalized Euclidean distance between the two sets of cells.
cells_distance <- function(x1, x2, order) {
  if (x1$input_parameters$initial_sequence_length != x2$input_parameters$initial_sequence_length) {
    stop("Input should have cells with the same initial sequence length.")
  }

  L <- x1$input_parameters$initial_sequence_length
  cells1 <- x1$cells
  cells2 <- x2$cells

  matrices <- lapply(list(cells1, cells2), function(cells) {
    cells2mat(cells, L, order)
  })

  distance <- sqrt(sum((matrices[[1]] - matrices[[2]])^2))
  normalized_distance <- distance / prod(dim(matrices[[1]]))

  return(normalized_distance)
}
