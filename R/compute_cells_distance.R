#' Compute Distance Between Two Sets of Cells
#'
#' This function calculates the distance between two sets of cells using a specified distance method.
#'
#' @param cells1 A list of cell sequences representing the first set of cells.
#' @param cells2 A list of cell sequences representing the second set of cells.
#' @param L An integer specifying the length of the vector representation for each cell.
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
cells_distance <- function(cells1, cells2, L) {

  matrices = lapply(list(cells1, cells2), function(cells) {
    m = matrix(0, nrow = length(cells), ncol = L)
    idx = 1
    for (idx in 1:length(cells)) {
      tab_vec = table(seq2vec(cells[[idx]]))
      m[idx, as.numeric(names(tab_vec))] = as.numeric(tab_vec)
    }

    one_copy_percentage = rowSums(m == 1) / L

    # Get the order of the values in the ordering vector
    row_order <- order(one_copy_percentage, decreasing = TRUE)

    # Reorder the matrix rows
    m <- m[row_order, , drop = FALSE]
    m

  })

  distance = sqrt(sum((matrices[[1]] - matrices[[2]])^2))
  normalized_distance = distance / prod(dim(matrices[[1]]))

  return(list(
    distance = distance,
    normalized_distance = normalized_distance
  ))
}
