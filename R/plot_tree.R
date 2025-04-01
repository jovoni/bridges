#' Plot a Cell Lineage Tree with Detailed Annotations
#'
#' @description
#' Creates a comprehensive visualization of cell lineage with flexible annotation options
#' for tracking biological features like BFB events and hotspot amplifications.
#'
#' @param x A list which has an element named cell_history, which is a data frame with required columns:
#'   - `cell_id`: Unique identifier for each cell
#'   - `parent_id`: Identifier of the parent cell
#'   - `birth_time`: Timestamp of cell birth
#'   - `bfb_event`: Logical indicating breakage-fusion-bridge event
#'   - `hotspot_gained`: Logical indicating hotspot amplification
#'
#' @param full_tree Logical. If TRUE, final tree will show all duplications as branches.
#'  If FALSE, only BFB replications will be shown on the tree.
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
#' @param tip_align Logical. If TRUE, align tip labels. Default is TRUE.
#'
#' @param tip_offset Numeric. Offset distance for tip labels. Default is 0.5.
#'
#' @param tip_hjust Numeric. Horizontal justification for tip labels. Default is 0.5.
#'
#' @param use_computed_branch_lengths Logical. If TRUE, use custom computed branch lengths.
#'   If FALSE, use ggtree's default branch length handling.
#'
#' @param hotspot_shape Numeric. Shape to use for hotspot amplification markers. Default is 21.
#'
#' @param hotspot_size Numeric. Size for hotspot amplification markers. Default is 2.
#'
#' @return A ggtree plot object with cell lineage visualization
#' @export
plot_tree <- function(
    x,
    full_tree = TRUE,
    legend.position = "bottom",
    annotate_bfb_events = TRUE,
    annotate_hotspot_amplifications = TRUE,
    bfb_event_colors = c("FALSE" = "black", "TRUE" = "indianred"),
    hotspot_colors = c("FALSE" = "white", "TRUE" = "indianred"),
    tip_align = TRUE,
    tip_offset = 0.5,
    tip_hjust = 0.5,
    use_computed_branch_lengths = TRUE,
    hotspot_shape = 21,
    hotspot_size = 2
) {
  cell_history = x$cell_history

  # Validate input data
  required_cols <- c("cell_id", "parent_id", "birth_time", "bfb_event", "hotspot_gained")
  if (!all(required_cols %in% names(cell_history))) {
    stop("Missing required columns in cell_history. Needed: ",
         paste(required_cols, collapse = ", "))
  }

  # Generate Newick string
  if (!full_tree) {
    cell_history = keep_only_bfb_branches(cell_history)
  }
  newick_str <- cell_history_to_newick(cell_history)

  # Read the Newick tree
  tree <- ape::read.tree(text = newick_str)

  # Prepare annotation data
  node_data <- dplyr::as_tibble(cell_history) %>%
    dplyr::mutate(
      node = match(.data$cell_id, c(tree$tip.label, tree$node.label)),
      is_tip = .data$cell_id %in% tree$tip.label
    )

  # Compute branch lengths if specified
  if (use_computed_branch_lengths) {
    branch_lengths <- dplyr::mutate(
      node_data,
      branch.length = ifelse(
        is.na(.data$parent_id),  # Root cell (root) has no parent
        0,                 # Root branch length = 0
        .data$birth_time - cell_history$birth_time[match(.data$parent_id, cell_history$cell_id)]
      )
    )$branch.length

    tree$edge.length <- branch_lengths[match(tree$edge[,2], node_data$node)]
    tree$edge.length[is.na(tree$edge.length)] <- 1
  }

  # Create base plot
  p <- ggtree::ggtree(tree) %<+% node_data

  # Add tip labels with flexible arguments
  p = p + ggtree::geom_tiplab(hjust = tip_hjust, align = tip_align, offset = tip_offset)

  # Conditionally add BFB event coloring
  if (annotate_bfb_events) {
    p <- p +
      ggtree::geom_tree(ggplot2::aes(color = .data$bfb_event)) +
      ggplot2::scale_color_manual(values = bfb_event_colors, name = "BFB Replication")
  }

  # Conditionally add hotspot amplification markers
  if (annotate_hotspot_amplifications) {
    p <- p +
      ggtree::geom_tippoint(ggplot2::aes(fill = .data$hotspot_gained), shape = hotspot_shape, size = hotspot_size) +
      ggplot2::scale_fill_manual(values = hotspot_colors, name = "Hotspot Amplification")
  }

  # Add final theming
  p <- p +
    ggplot2::theme(legend.position = legend.position)

  return(p)
}
