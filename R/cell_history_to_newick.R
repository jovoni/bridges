
cell_history_to_newick <- function(cell_history) {
  # Check there is a root and rename it
  root_name = cell_history$cell_id[is.na(cell_history$parent_id)]
  cell_history$cell_id[cell_history$cell_id == root_name] = "root"
  cell_history$parent_id[cell_history$parent_id == root_name] = "root"

  # Helper function to recursively build the tree
  build_tree <- function(node) {
    # Find children of the current node
    node_data <- cell_history %>%
      dplyr::filter(.data$cell_id == node) %>%
      dplyr::select(.data$bfb_event)

    # Find children of the current node
    children <- cell_history %>% dplyr::filter(.data$parent_id == node) %>% dplyr::pull(.data$cell_id)

    if (length(children) == 0) {
      # If no children, return the node itself with BFB annotation
      return(node)
    } else {
      # Recursively build subtrees for each child
      subtree <- paste(sapply(children, build_tree), collapse = ",")

      # Add BFB annotation to the current node
      return(paste0("(", subtree, ")", node))
    }
  }

  # Identify the root node (cells with no parent)
  root <- cell_history %>%
    dplyr::filter(is.na(.data$parent_id)) %>%
    dplyr::pull(.data$cell_id)

  if (length(root) != 1) {
    stop("Error: There must be exactly one root node.")
  }

  # Build the tree starting from the root
  newick_tree <- paste0(build_tree(root), ";")

  return(newick_tree)
}
