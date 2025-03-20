#' Compute Distance Between Two Sets of Cells
#'
#' This function calculates the distance between two sets of cells using a specified distance method.
#'
#' @param cells1 A list of cell sequences representing the first set of cells.
#' @param cells2 A list of cell sequences representing the second set of cells.
#' @param L An integer specifying the length of the vector representation for each cell.
#' @param order Boolean indicating whether cells should be ordered by proportion of genome in 1 copy
#'
#' @return A list containing:
#'   \item{distance}{The Euclidean distance between the two sets of cells.}
#'   \item{normalized_distance}{The normalized Euclidean distance, scaled by the dimensions of the matrices.}
#'
#' @details
#' The function converts each set of cells into a matrix representation, where each row corresponds to a cell, and columns represent bins.
#' Rows are reordered based on the proportion of features with exactly one copy in descending order. The Euclidean distance between the
#' two matrices is computed and normalized by the product of the matrix dimensions.
#'
#' @export
cells_distance <- function(cells1, cells2, L, order) {

  matrices = lapply(list(cells1, cells2), function(cells) {
    cells2mat(cells, L, order)
  })

  distance = sqrt(sum((matrices[[1]] - matrices[[2]])^2))
  normalized_distance = distance / prod(dim(matrices[[1]]))
  normalized_distance

  # return(list(
  #   distance = distance,
  #   normalized_distance = normalized_distance
  # ))
}
