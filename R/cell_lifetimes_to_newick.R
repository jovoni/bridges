
cell_lifetimes_to_newick <- function(cell_data) {
  # Check there is a root and rename it
  root_name = cell_data$cell_id[is.na(cell_data$parent_id)]
  cell_data$cell_id[cell_data$cell_id == root_name] = "root"
  cell_data$parent_id[cell_data$parent_id == root_name] = "root"

  # # Add root if not present
  # if (!"root" %in% cell_data$cell_id) {
  #   root_row <- dplyr::tibble(
  #     cell_id = "root",
  #     parent_id = NA,
  #     bfb_event = FALSE,
  #     birth_time = 0,
  #     death_time = 0
  #   )
  #   cell_data <- dplyr::bind_rows(root_row, cell_data)
  # }


  # Helper function to recursively build the tree
  build_tree <- function(node) {
    # Find children of the current node
    node_data <- cell_data %>%
      dplyr::filter(cell_id == node) %>%
      dplyr::select(bfb_event)

    # Find children of the current node
    children <- cell_data %>% dplyr::filter(parent_id == node) %>% dplyr::pull(cell_id)

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
  root <- cell_data %>% dplyr::filter(is.na(parent_id)) %>% dplyr::pull(cell_id)

  if (length(root) != 1) {
    stop("Error: There must be exactly one root node.")
  }

  # Build the tree starting from the root
  newick_tree <- paste0(build_tree(root), ";")

  return(newick_tree)
}
