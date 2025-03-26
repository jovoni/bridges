#' Plot a Cell Lineage Tree with Detailed Annotations
#'
#' @description
#' Creates a comprehensive visualization of cell lineage with flexible annotation options
#' for tracking biological features like BFB events and hotspot amplifications.
#'
#' @param cell_data A data frame containing cell information with required columns:
#'   - `cell_id`: Unique identifier for each cell
#'   - `parent_id`: Identifier of the parent cell
#'   - `birth_time`: Timestamp of cell birth
#'   - `bfb_event`: Logical indicating breakage-fusion-bridge event
#'   - `hotspot_gained`: Logical indicating hotspot amplification
#'
#' @param legend.position Position of the legend. Defaults to "bottom".
#'   Accepts standard ggplot2 legend positions: "bottom", "top", "left", "right", "none"
#'
#' @param annotate_bfb_events Logical. If TRUE, color tree branches based on BFB events
#'
#' @param annotate_hotspot_amplifications Logical. If TRUE, add point markers for hotspot amplifications
#'
#' @param bfb_event_colors Named vector specifying colors for BFB event status
#'   Default: c("FALSE" = "black", "TRUE" = "indianred")
#'
#' @param hotspot_colors Named vector specifying colors for hotspot amplification status
#'   Default: c("FALSE" = "white", "TRUE" = "indianred")
#'
#' @param tip_label_args List of arguments to pass to geom_tiplab() for customizing tip labels
#'   Default provides alignment and offset
#'
#' @param use_computed_branch_lengths Logical. If TRUE, use custom computed branch lengths.
#'   If FALSE, use ggtree's default branch length handling.
#'
#' @return A ggtree plot object with cell lineage visualization
#' @export
plot_tree <- function(
    cell_data,
    full_tree = TRUE,
    legend.position = "bottom",
    annotate_bfb_events = TRUE,
    annotate_hotspot_amplifications = TRUE,
    bfb_event_colors = c("FALSE" = "black", "TRUE" = "indianred"),
    hotspot_colors = c("FALSE" = "white", "TRUE" = "indianred"),
    tip_label_args = list(align = TRUE, offset = 0.5, hjust = 0.5),
    use_computed_branch_lengths = TRUE
) {
  # Validate input data
  required_cols <- c("cell_id", "parent_id", "birth_time", "bfb_event", "hotspot_gained")
  if (!all(required_cols %in% names(cell_data))) {
    stop("Missing required columns in cell_data. Needed: ",
         paste(required_cols, collapse = ", "))
  }

  # Generate Newick string
  if (!full_tree) {
    cell_data = find_bfb_only_cell_data(cell_data)
  }
  newick_str <- cell_lifetimes_to_newick(cell_data)

  # Read the Newick tree
  tree <- ape::read.tree(text = newick_str)

  # Prepare annotation data
  node_data <- dplyr::as_tibble(cell_data) %>%
    dplyr::mutate(
      node = match(cell_id, c(tree$tip.label, tree$node.label)),
      is_tip = cell_id %in% tree$tip.label
    )

  # Compute branch lengths if specified
  if (use_computed_branch_lengths) {
    branch_lengths <- dplyr::mutate(
      node_data,
      branch.length = ifelse(
        is.na(parent_id),  # Root cell (cell_0) has no parent
        0,                 # Root branch length = 0
        birth_time - cell_data$birth_time[match(parent_id, cell_data$cell_id)]
      )
    )$branch.length

    tree$edge.length <- branch_lengths[match(tree$edge[,2], node_data$node)]
    tree$edge.length[is.na(tree$edge.length)] <- 1
  }

  # Create base plot
  p <- ggtree::ggtree(tree) %<+% node_data

  # Add tip labels with flexible arguments
  p <- p + do.call(ggtree::geom_tiplab, tip_label_args)

  # Conditionally add BFB event coloring
  if (annotate_bfb_events) {
    p <- p +
      ggtree::geom_tree(aes(color = bfb_event)) +
      ggplot2::scale_color_manual(values = bfb_event_colors, name = "BFB Replication")
  }

  # Conditionally add hotspot amplification markers
  if (annotate_hotspot_amplifications) {
    p <- p +
      ggtree::geom_tippoint(aes(fill = hotspot_gained), shape = 21, size = 3) +
      ggplot2::scale_fill_manual(values = hotspot_colors, name = "Hotspot Amplification")
  }

  # Add final theming
  p <- p +
    ggplot2::theme(legend.position = legend.position)

  return(p)
}
