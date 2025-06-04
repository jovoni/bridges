
get_ordered_cell_ids <- function(tree_plot_dat) {
  rev(dplyr::arrange(tree_plot_dat[tree_plot_dat$isTip, ], .data$y)$label)
}

#' Validate input parameters for heatmap plotting
#' @param data Input data frame
#' @param tree Phylogenetic tree object (optional)
#' @param chromosomes_to_plot Vector of chromosomes to include
#' @param to_plot Vector of features to plot
validate_heatmap_inputs <- function(data, tree, chromosomes_to_plot, to_plot) {
  # Validate data structure
  required_cols <- c("cell_id", "chr", "start")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns in data: %s", paste(missing_cols, collapse = ", ")))
  }

  # Validate to_plot columns exist
  missing_features <- setdiff(to_plot, colnames(data))
  if (length(missing_features) > 0) {
    stop(sprintf("Feature columns not found in data: %s", paste(missing_features, collapse = ", ")))
  }

  # Validate chromosomes
  available_chrs <- unique(data$chr)
  invalid_chrs <- setdiff(chromosomes_to_plot, available_chrs)
  if (length(invalid_chrs) > 0) {
    warning(sprintf("Chromosomes not found in data: %s", paste(invalid_chrs, collapse = ", ")))
  }

  # Validate tree if provided
  if (!is.null(tree)) {
    if (!inherits(tree, "phylo")) {
      stop("tree must be a phylo object")
    }

    data_cells <- unique(data$cell_id)
    tree_cells <- tree$tip.label

    if (!all(tree_cells %in% data_cells)) {
      missing_cells <- setdiff(tree_cells, data_cells)
      stop(sprintf("Tree tip labels not found in data: %s", paste(utils::head(missing_cells, 5), collapse = ", ")))
    }
  }
}

#' Prepare row annotations from a single data frame
#' @param annotations Data frame with cell_id and annotation columns
#' @param ordered_cell_ids Vector of cell IDs in desired order
#' @return ComplexHeatmap annotation object or NULL
prepare_row_annotations <- function(annotations, ordered_cell_ids) {
  if (is.null(annotations)) {
    return(NULL)
  }

  if (!is.data.frame(annotations)) {
    stop("annotations must be a data frame")
  }

  if (!"cell_id" %in% colnames(annotations)) {
    stop("annotations data frame must have a 'cell_id' column")
  }

  # Get all annotation columns (everything except cell_id)
  annotation_cols <- setdiff(colnames(annotations), "cell_id")

  if (length(annotation_cols) == 0) {
    warning("No annotation columns found (only cell_id column)")
    return(NULL)
  }

  # Validate all required cell IDs are present
  missing_cells <- setdiff(ordered_cell_ids, annotations$cell_id)
  if (length(missing_cells) > 0) {
    stop(sprintf("Cell IDs not found in annotations: %s",
                 paste(utils::head(missing_cells, 5), collapse = ", ")))
  }

  # Match and reorder to match ordered_cell_ids
  matched_indices <- match(ordered_cell_ids, annotations$cell_id)
  ordered_annotations <- annotations[matched_indices, , drop = FALSE]

  # Prepare annotation values and colors for each column
  annotation_values <- list()
  annotation_colors <- list()

  for (col_name in annotation_cols) {
    values <- ordered_annotations[[col_name]]

    if (is.numeric(values)) {
      # Continuous annotation
      annotation_values[[col_name]] <- values

      if (all(values >= 0 & values <= 1, na.rm = TRUE)) {
        # Looks like proportions/probabilities - use white to red
        annotation_colors[[col_name]] <- circlize::colorRamp2(c(0, 1), c("white", "red"))
      } else {
        # General continuous data - use blue-white-red
        value_range <- range(values, na.rm = TRUE)
        annotation_colors[[col_name]] <- circlize::colorRamp2(
          c(value_range[1], mean(value_range), value_range[2]),
          c("blue", "white", "red")
        )
      }
    } else {
      # Categorical annotation
      values <- factor(values)
      annotation_values[[col_name]] <- values

      levels_count <- length(levels(values))

      if (levels_count <= 12) {
        # Use RColorBrewer for better colors
        colors <- RColorBrewer::brewer.pal(max(3, levels_count), "Set3")[1:levels_count]
      } else {
        # Fall back to rainbow
        colors <- grDevices::rainbow(levels_count)
      }

      names(colors) <- levels(values)
      annotation_colors[[col_name]] <- colors
    }
  }

  # Create the row annotation using do.call
  annotation_args <- c(
    annotation_values,
    list(
      col = annotation_colors,
      annotation_name_side = "top",
      simple_anno_size = ggplot2::unit(0.5, "cm"),
      gap = ggplot2::unit(1, "mm")
    )
  )

  row_annotation <- do.call(ComplexHeatmap::rowAnnotation, annotation_args)

  row_annotation
}

#' Optimize tree tip ordering based on distance matrix
#' @param tree Phylogenetic tree
#' @param distance_matrix Distance matrix for optimization
#' @return Reordered phylogenetic tree
optimize_tree_ordering <- function(tree, distance_matrix) {
  if (!requireNamespace("seriation", quietly = TRUE)) {
    warning("seriation package not available, skipping tree optimization")
    return(tree)
  }

  tryCatch({
    # Use seriation to find optimal ordering
    ser_order <- seriation::seriate(distance_matrix, method = "GW")
    optimal_order <- rownames(as.matrix(distance_matrix))[seriation::get_order(ser_order)]

    # Rotate tree nodes to match optimal ordering
    tree_optimized <- ape::rotateConstr(tree, constraint = optimal_order)

    return(tree_optimized)
  }, error = function(e) {
    warning(sprintf("Tree optimization failed: %s", e$message))
    tree
  })
}

#' Process tree for heatmap visualization
#' @param tree Phylogenetic tree object
#' @param branch_length Uniform branch length (optional)
#' @param ladderize Whether to ladderize the tree
#' @param reorder_tree Whether to optimize tree ordering
#' @param distance_matrix Distance matrix for reordering (reorder_tree = TRUE)
#' @return List containing processed tree and ggplot object
process_tree <- function(tree, branch_length = NULL, ladderize = FALSE,
                         reorder_tree = FALSE, distance_matrix = NULL) {

  # Set uniform branch lengths if specified
  if (!is.null(branch_length) && !is.null(tree$edge.length)) {
    tree$edge.length <- rep(branch_length, length(tree$edge.length))
  }

  # Optimize tree ordering if requested
  if (reorder_tree && !is.null(distance_matrix)) {
    if (ladderize) {
      warning("'ladderize=TRUE' will disrupt reordering from 'reorder_tree=TRUE'")
    }
    tree <- optimize_tree_ordering(tree, distance_matrix)
  }

  # Create ggplot tree
  tree_ggplot <- ggtree::ggtree(tree, size = 0.25, ladderize = ladderize) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::ylim(0.5, length(tree$tip.label) + 0.5)

  return(list(tree = tree, ggplot = tree_ggplot))
}

#' Prepare matrix data for heatmap plotting
#' @param data Input data frame
#' @param chromosomes_to_plot Chromosomes to include
#' @param feature_name Name of feature column to extract
#' @param ordered_cell_ids Cell IDs in desired order
#' @return List containing matrix and column information
prepare_heatmap_matrix <- function(data, chromosomes_to_plot, feature_name, ordered_cell_ids) {
  # Filter and prepare data
  processed_data <- data %>%
    dplyr::filter(.data$chr %in% chromosomes_to_plot) %>%
    dplyr::group_by(.data$cell_id, .data$chr) %>%
    dplyr::arrange(.data$start, .by_group = TRUE) %>%
    dplyr::mutate(
      bin_idx = dplyr::row_number(),
      chr_bin = paste(.data$chr, .data$bin_idx, sep = "_")
    ) %>%
    dplyr::select(.data$cell_id, .data$chr_bin, .data$chr, !!rlang::sym(feature_name)) %>%
    dplyr::ungroup()

  # Create matrix
  mat_data <- processed_data %>%
    dplyr::select(.data$cell_id, .data$chr_bin, !!rlang::sym(feature_name)) %>%
    tidyr::pivot_wider(
      names_from = .data$chr_bin,
      values_from = !!rlang::sym(feature_name),
      values_fill = NA
    ) %>%
    tibble::column_to_rownames("cell_id")

  # Prepare column information
  col_info <- data.frame(
    chr_bin = colnames(mat_data),
    stringsAsFactors = FALSE
  )
  col_info$chr <- sapply(strsplit(col_info$chr_bin, "_"), `[`, 1)

  # Convert to matrix and handle extreme values
  mat <- as.matrix(mat_data)

  # Cap extreme values (customize threshold as needed)
  extreme_threshold <- 10
  extreme_mask <- !is.na(mat) & mat > extreme_threshold
  mat[extreme_mask] <- paste0(extreme_threshold + 1, "+")

  # Ensure proper ordering
  common_cells <- intersect(ordered_cell_ids, rownames(mat))
  if (length(common_cells) == 0) {
    stop("No common cell IDs found between data and ordered list")
  }

  mat <- mat[common_cells, , drop = FALSE]

  list(matrix = mat, col_info = col_info, cells_used = common_cells)
}

#' Create tree annotation for heatmap
#' @param tree_ggplot ggplot tree object
#' @param tree_width Width of tree annotation
#' @param n_cells Number of cells in heatmap
#' @return HeatmapAnnotation object
create_tree_annotation <- function(tree_ggplot, tree_width, n_cells) {
  tree_annot_func <- ComplexHeatmap::AnnotationFunction(
    fun = function(index) {
      grid::pushViewport(grid::viewport(height = 1))
      grid::grid.draw(ggplot2::ggplotGrob(tree_ggplot)$grobs[[6]])
      grid::popViewport()
    },
    var_import = list(tree_ggplot = tree_ggplot),
    width = ggplot2::unit(tree_width, "cm"),
    which = "row"
  )

  tree_annot <- ComplexHeatmap::HeatmapAnnotation(
    tree = tree_annot_func,
    which = "row",
    show_annotation_name = FALSE
  )

  tree_hm <- ComplexHeatmap::Heatmap(
    matrix(ncol = 0, nrow = n_cells),
    left_annotation = tree_annot,
    show_heatmap_legend = FALSE
  )

  tree_hm
}

#' Create comprehensive heatmap with optional phylogenetic tree
#'
#' This function creates a heatmap visualization with optional phylogenetic
#' tree display, annotations, and support for multiple features. It can display
#' genomic data across chromosomes with customizable tree coloring based on
#' reconstruction data.
#'
#' @param data Data frame containing cell_id, chr, start, and feature columns.
#'   Must include the columns specified in `to_plot` parameter.
#' @param tree Phylogenetic tree object (optional). Should be compatible with
#'   ape or phylo class objects.
#' @param chromosomes_to_plot Vector of chromosomes to include in the plot.
#'   Default: c(1:22, "X", "Y") for human chromosomes.
#' @param to_plot Vector of feature names to plot. Default: "CN" (copy number).
#'   Must correspond to column names in the data.
#' @param tree_width Width of tree display in centimeters. Default: 3.5.
#' @param show_tree Logical indicating whether to display the phylogenetic tree.
#'   Default: TRUE.
#' @param branch_length Numeric value for uniform branch length transformation
#'   of the tree. If NULL (default), original branch lengths are preserved.
#' @param annotations Data frame with cell_id column and additional annotation
#'   columns. Each non-cell_id column becomes a row annotation in the heatmap.
#'   Optional parameter.
#' @param ladderize Logical indicating whether to ladderize (sort) the tree
#'   branches. Default: FALSE.
#' @param reorder_tree Logical indicating whether to optimize tree tip ordering
#'   for better visualization. Default: TRUE.
#' @param distance_matrix Distance matrix for tree reordering optimization.
#'   Optional parameter used when reorder_tree is TRUE.
#' @param use_raster Logical indicating whether to use raster graphics for
#'   the heatmap. Default: TRUE for better performance with large datasets.
#' @param raster_quality Numeric value controlling raster image quality.
#'   Higher values produce better quality. Default: 10.
#' @param reconstruction Reconstruction data for coloring tree branches based
#'   on ancestral state reconstruction. Optional parameter.
#' @param chr_for_coloring Chromosome identifier used for tree branch coloring
#'   when reconstruction data is provided. Optional parameter.
#' @param allele_for_coloring Specific allele used for tree branch coloring
#'   when reconstruction data is provided. Optional parameter.
#' @param B_dist Distance parameter B for tree coloring algorithm. Optional.
#' @param G_dist Distance parameter G for tree coloring algorithm. Optional.
#' @param tree_colors Vector of colors for tree branch coloring. Default:
#'   c("black", "goldenrod"). First color for one state, second for another.
#' @param node_size Size of internal nodes in the tree. Default: 1.
#' @param tip_size Size of tip nodes in the tree visualization. Default: 0
#'   (tips not displayed).
#'
#' @return A ComplexHeatmap object that can be displayed or further customized.
#'   The object contains the combined heatmap(s) and optional tree annotation.
#'
#' @details
#' The function performs the following main steps:
#' \itemize{
#'   \item Validates input parameters and data structure
#'   \item Processes the phylogenetic tree (if provided) with optional coloring
#'   \item Prepares row annotations from the annotations data frame
#'   \item Creates individual heatmaps for each feature in to_plot
#'   \item Combines multiple heatmaps horizontally
#'   \item Adds tree annotation to the left side (if requested)
#' }
#'
#' The tree can be colored based on reconstruction data, which is useful for
#' visualizing ancestral state reconstructions or other phylogenetic analyses.
#' When reconstruction data is provided, the tree branches are colored according
#' to the specified parameters.
#'
#' @export
plot_heatmap <- function(data,
                         tree = NULL,
                         chromosomes_to_plot = c(1:22, "X", "Y"),
                         to_plot = "CN",
                         tree_width = 3.5,
                         show_tree = TRUE,
                         branch_length = NULL,
                         annotations = NULL,
                         ladderize = FALSE,
                         reorder_tree = TRUE,
                         distance_matrix = NULL,
                         use_raster = TRUE,
                         raster_quality = 10,
                         # New parameters for colored tree
                         reconstruction = NULL,
                         chr_for_coloring = NULL,
                         allele_for_coloring = NULL,
                         B_dist = NULL,
                         G_dist = NULL,
                         tree_colors = c("black", "goldenrod"),
                         node_size = 1,
                         tip_size = 0) {

  # Validate inputs
  validate_heatmap_inputs(data, tree, chromosomes_to_plot, to_plot)

  # Initialize variables
  to_plot_vec <- as.character(to_plot)
  all_heatmaps <- list()
  ordered_cell_ids <- unique(data$cell_id)

  # Process tree if provided
  tree_processed <- NULL
  if (!is.null(tree)) {
    # Check if we should use colored tree or regular tree
    if (!is.null(reconstruction)) {
      # Use colored tree based on reconstruction
      tree_processed <- process_tree_with_coloring(
        tree = tree,
        reconstruction = reconstruction,
        branch_length = branch_length,
        ladderize = ladderize,
        reorder_tree = reorder_tree,
        distance_matrix = distance_matrix,
        colors = tree_colors,
        node_size = node_size,
        tip_size = tip_size
      )
    } else {
      # Use regular tree processing
      tree_processed <- process_tree(
        tree = tree,
        branch_length = branch_length,
        ladderize = ladderize,
        reorder_tree = reorder_tree,
        distance_matrix = distance_matrix
      )
    }

    # Extract ordered cell IDs from tree
    ordered_cell_ids <- get_ordered_cell_ids(tree_processed$ggplot$data)
  }

  # Prepare row annotations once (shared across all heatmaps)
  row_annotations <- prepare_row_annotations(annotations, ordered_cell_ids)

  # Create heatmaps for each feature

  for (i in seq_along(to_plot_vec)) {
    feature_name <- to_plot_vec[i]
    message(sprintf("Processing feature: %s", feature_name))

    # Prepare matrix data
    matrix_data <- prepare_heatmap_matrix(
      data = data,
      chromosomes_to_plot = chromosomes_to_plot,
      feature_name = feature_name,
      ordered_cell_ids = ordered_cell_ids
    )

    # Define chromosome order
    chrs_order <- as.character(c(1:22, "X", "Y"))
    valid_chrs <- intersect(chrs_order, unique(matrix_data$col_info$chr))

    # Get color scheme
    color_scheme <- get_colors("CN")

    # Create heatmap - only add annotations to the first heatmap
    heatmap_obj <- ComplexHeatmap::Heatmap(
      matrix_data$matrix,
      name = feature_name,
      left_annotation = if (i == 1) row_annotations else NULL,
      col = color_scheme,
      cluster_rows = FALSE,
      show_row_names = FALSE,
      row_names_gp = grid::gpar(fontsize = 8),
      cluster_columns = FALSE,
      column_split = factor(matrix_data$col_info$chr, levels = valid_chrs),
      column_title_gp = grid::gpar(fontsize = 10),
      column_title_rot = 0,
      show_column_names = FALSE,
      heatmap_legend_param = list(
        title = feature_name,
        title_gp = grid::gpar(fontsize = 10),
        labels_gp = grid::gpar(fontsize = 8)
      ),
      column_gap = ggplot2::unit(0.5, "mm"),
      border = TRUE,
      raster_quality = raster_quality,
      use_raster = use_raster
    )

    all_heatmaps[[feature_name]] <- heatmap_obj
  }

  # Combine heatmaps
  combined_heatmaps <- Reduce(`+`, all_heatmaps)

  # Add tree if requested
  if (!is.null(tree_processed) && show_tree) {
    tree_annotation <- create_tree_annotation(
      tree_ggplot = tree_processed$ggplot,
      tree_width = tree_width,
      n_cells = nrow(all_heatmaps[[1]]@matrix)
    )

    final_plot <- tree_annotation + combined_heatmaps
  } else {
    final_plot <- combined_heatmaps
  }

  return(final_plot)
}

process_tree_with_coloring <- function(tree,
                                       reconstruction,
                                       branch_length = NULL,
                                       ladderize = FALSE,
                                       reorder_tree = FALSE,
                                       distance_matrix = NULL,
                                       colors = c("black", "goldenrod"),
                                       node_size = 3, tip_size = 2) {

  # First apply basic tree processing
  if (!is.null(branch_length) && !is.null(tree$edge.length)) {
    tree$edge.length <- rep(branch_length, length(tree$edge.length))
  }

  # Optimize tree ordering if requested
  if (reorder_tree && !is.null(distance_matrix)) {
    if (ladderize) {
      warning("'ladderize=TRUE' will disrupt reordering from 'reorder_tree=TRUE'")
    }
    tree <- optimize_tree_ordering(tree, distance_matrix)
  }

  # Get reconstruction results
  if (is.list(reconstruction) && "deltas" %in% names(reconstruction)) {
    recon_result <- reconstruction
  }

  # Create node-level dataframe
  internal_nodes <- as.integer(names(recon_result$deltas))
  delta_values <- unlist(recon_result$deltas)
  decision_labels <- ifelse(delta_values > 0, "Likely BFB", "Other")

  node_data <- data.frame(
    node = internal_nodes,
    delta = delta_values,
    decision = decision_labels,
    stringsAsFactors = FALSE
  )

  # Create base ggtree object
  p <- ggtree::ggtree(tree, size = 0.25, ladderize = ladderize)

  # Extract ggtree data (edge data)
  tree_data <- p$data

  # Merge decision info to internal nodes only
  tree_data$decision <- "Other"
  tree_data$decision[tree_data$node %in% node_data$node] <- node_data$decision[match(tree_data$node[tree_data$node %in% node_data$node], node_data$node)]

  # For edge coloring: assign decision to the *parent* of edge
  tree_data$edge_decision <- "Other"
  bfb_nodes <- node_data$node[node_data$decision == "Likely BFB"]
  tree_data$edge_decision[tree_data$parent %in% bfb_nodes] <- "Likely BFB"

  # Update tree plot with colored branches (edges)
  color_map <- stats::setNames(colors, c("Other", "Likely BFB"))

  p <- ggtree::ggtree(tree, size = 0.25, ladderize = ladderize) %<+% tree_data +
    ggtree::geom_tree(ggplot2::aes(color = .data$edge_decision)) +
    ggtree::geom_tippoint(size = tip_size, color = "black", alpha = 0.7) +
    ggplot2::scale_color_manual(
      values = color_map,
      name = "CN event"
    ) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::ylim(0.5, length(tree$tip.label) + 0.5)

  return(list(tree = tree, ggplot = p))
}
