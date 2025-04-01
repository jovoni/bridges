#' Plot a Cell Lineage Tree with Detailed Annotations
#'
#' @description
#' Generates a cell lineage tree visualization with options to annotate BFB events and hotspot amplifications.
#'
#' @param tree A list containing:
#'   - `tree`: A phylogenetic tree object of class `phylo`.
#'   - `node_data`: A data frame with node-specific information including `cell_id`, `bfb_event`, and `hotspot_gained`.
#'
#' @param ladderize Logical. If TRUE, ladderizes the tree for better visualization. Default is TRUE.
#'
#' @param legend.position Character. Position of the legend, one of "none", "bottom", "top", "left", "right". Default is "none".
#'
#' @param annotate_tip Logical. If TRUE, displays tip labels. Default is FALSE.
#'
#' @param annotate_bfb_events Logical. If TRUE, colors tree branches based on BFB events. Default is TRUE.
#'
#' @param annotate_hotspot_amplifications Logical. If TRUE, adds point markers for hotspot amplifications. Default is TRUE.
#'
#' @param bfb_event_colors Named vector specifying colors for BFB event status. Default: c("FALSE" = "black", "TRUE" = "firebrick3").
#'
#' @param hotspot_colors Named vector specifying colors for hotspot amplification status. Default: c("FALSE" = "white", "TRUE" = "indianred").
#'
#' @param tip_align Logical. If TRUE, aligns tip labels. Default is TRUE.
#'
#' @param tip_offset Numeric. Offset distance for tip labels. Default is 0.5.
#'
#' @param size Numeric. Line width of the tree branches. Default is 0.25.
#'
#' @param tip_hjust Numeric. Horizontal justification for tip labels. Default is 0.5.
#'
#' @param expand Logical. If TRUE, expands the plot limits. Default is FALSE.
#'
#' @param hotspot_shape Numeric. Shape for hotspot amplification markers. Default is 21.
#'
#' @param hotspot_size Numeric. Size for hotspot amplification markers. Default is 0.25.
#'
#' @return A `ggtree` plot object representing the cell lineage tree.
#' @export
plot_tree = function(tree,
                     ladderize = TRUE,
                     legend.position = "none",
                     annotate_tip = FALSE,
                     annotate_bfb_events = TRUE,
                     annotate_hotspot_amplifications = TRUE,
                     bfb_event_colors = c("FALSE" = "black", "TRUE" = "firebrick3"),
                     hotspot_colors = c("FALSE" = "white", "TRUE" = "indianred"),
                     tip_align = TRUE,
                     tip_offset = 0.5,
                     size = 0.25,
                     tip_hjust = 0.5,
                     expand = FALSE,
                     hotspot_shape = 21,
                     hotspot_size = 0.25) {

  p <- ggtree::ggtree(tree$tree, size = size, ladderize = ladderize) +
    ggplot2::coord_cartesian(expand = expand) +
    ggplot2::ylim(0.5, length(tree$tree$tip.label) + 0.5)

  p$data = p$data %>%
    dplyr::left_join(tree$node_data %>% dplyr::select(!.data$node) %>% dplyr::mutate(label = .data$cell_id), by = "label")

  if (annotate_tip) {
    p = p + ggtree::geom_tiplab(hjust = tip_hjust, align = tip_align, offset = tip_offset)
  }

  if (annotate_bfb_events) {
    p <- p +
      ggtree::geom_tree(ggplot2::aes(color = .data$bfb_event), size = size) +
      ggplot2::scale_color_manual(values = bfb_event_colors, name = "BFB Replication")
  }

  if (annotate_hotspot_amplifications) {
    p <- p +
      ggtree::geom_tippoint(ggplot2::aes(fill = .data$hotspot_gained), shape = hotspot_shape, size = hotspot_size) +
      ggplot2::scale_fill_manual(values = hotspot_colors, name = "Hotspot Amplification")
  }

  p <- p +
    ggplot2::theme(legend.position = legend.position)

  return(p)
}

#' Generate Tree Data from Cell History
#'
#' @description
#' Processes a cell history dataset into a phylogenetic tree with optional root addition and branch length computation.
#'
#' @param x A list containing `cell_history`, a data frame with columns:
#'   - `cell_id`: Unique identifier for each cell.
#'   - `parent_id`: Identifier of the parent cell.
#'   - `birth_time`: Timestamp of cell birth.
#'   - `bfb_event`: Logical, indicating breakage-fusion-bridge event.
#'   - `hotspot_gained`: Logical, indicating hotspot amplification.
#'
#' @param full_tree Logical. If TRUE, keeps all branches; if FALSE, keeps only BFB-related branches. Default is TRUE.
#'
#' @param add_root Logical. If TRUE, adds a root node to cells without parents. Default is FALSE.
#'
#' @param branch_lengths Logical. If TRUE, computes custom branch lengths based on birth times. Default is TRUE.
#'
#' @return A list containing:
#'   - `tree`: A phylogenetic tree object of class `phylo`.
#'   - `node_data`: A data frame with node-specific information.
#' @export
get_tree_data = function(x, full_tree = TRUE, add_root = FALSE, branch_lengths = TRUE) {
  cell_history = x$cell_history

  required_cols <- c("cell_id", "parent_id", "birth_time", "bfb_event", "hotspot_gained")
  if (!all(required_cols %in% names(cell_history))) {
    stop("Missing required columns in cell_history. Needed: ",
         paste(required_cols, collapse = ", "))
  }

  if (!full_tree) {
    cell_history = keep_only_bfb_branches(cell_history)
  }

  if (add_root) {
    cell_history$parent_id[is.na(cell_history$parent_id)] = "root"
    cell_history = dplyr::bind_rows(
      dplyr::tibble(
        cell_id = "root",
        birth_time = -1,
        death_time = -1,
        lifetime = 0,
        is_alive = FALSE,
        parent_id = NA,
        bfb_event = FALSE,
        hotspot_gained = FALSE
      ),
      cell_history
    )
  }

  if (sum(is.na(cell_history$parent_id)) > 1) stop("Multiple cells don't have a parent. Run with 'add_root=TRUE'")
  newick_str <- cell_history_to_newick(cell_history)
  tree = ape::read.tree(text = newick_str)

  node_data <- dplyr::as_tibble(cell_history) %>%
    dplyr::mutate(
      node = match(.data$cell_id, c(tree$tip.label, tree$node.label)),
      is_tip = .data$cell_id %in% tree$tip.label
    )

  if (branch_lengths) {
    branch_lengths <- dplyr::mutate(
      node_data,
      branch.length = ifelse(
        is.na(.data$parent_id),
        0,
        .data$birth_time - cell_history$birth_time[match(.data$parent_id, cell_history$cell_id)]
      )
    )$branch.length

    tree$edge.length <- branch_lengths[match(tree$edge[,2], node_data$node)]
    tree$edge.length[is.na(tree$edge.length)] <- 1
  }

  return(list(tree = tree, node_data = node_data))
}
