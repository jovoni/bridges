
build_cell_history <- function(tree, parent, parent_id = NA) {
  # Create an initial tibble entry for the current parent
  phylo <- dplyr::tibble(cell_id = parent, parent_id = parent_id)

  # Check if the parent exists in the tree
  if (!parent %in% names(tree)) {
    return(phylo)  # Return just the single node if no children
  }

  # Get the children of the current parent
  children <- tree[[parent]]

  # Recursively process each child
  for (child in children) {
    phylo <- dplyr::bind_rows(phylo, build_cell_history(tree, child, parent))
  }

  return(phylo)
}
