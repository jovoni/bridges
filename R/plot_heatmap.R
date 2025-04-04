
get_ordered_cell_ids <- function(tree_plot_dat) {
  return(rev(dplyr::arrange(tree_plot_dat[tree_plot_dat$isTip, ], y)$label))
}

make_corrupt_tree_heatmap <- function(tree_ggplot, tree_width) {
  tree_annot_func <- ComplexHeatmap::AnnotationFunction(
    fun = function(index) {
      grid::pushViewport(grid::viewport(height = 1))
      grid::grid.draw(ggplot2::ggplotGrob(tree_ggplot)$grobs[[5]])
      grid::popViewport()
    },
    var_import = list(tree_ggplot = tree_ggplot),
    width = grid::unit(tree_width, "cm"),
    which = "row"
  )
  tree_annot <- ComplexHeatmap::HeatmapAnnotation(
    tree = tree_annot_func, which = "row", show_annotation_name = FALSE
  )

  n_cells <- sum(tree_ggplot$data$isTip)
  tree_hm <- ComplexHeatmap::Heatmap(matrix(ncol = 0, nrow = n_cells), left_annotation = tree_annot)

  return(tree_hm)
}

plot_heatmap = function(x,
                        order_heatmap = TRUE,
                        add_gain_loss_profile = FALSE,
                        add_tree = FALSE,
                        tree_width = 1,
                        full_tree = TRUE,
                        add_root = FALSE,
                        use_raster = TRUE,
                        raster_quality = 15,
                        branch_lengths = FALSE,
                        ladderize = TRUE,
                        legend.position = "none",
                        annotate_tip = FALSE,
                        annotate_bfb_events = TRUE,
                        annotate_hotspot_amplifications = TRUE,
                        bfb_event_colors = c("FALSE" = "black", "TRUE" = "firebrick3"),
                        hotspot_colors = c("FALSE" = "white", "TRUE" = "indianred"),
                        tip_align = TRUE,
                        plot_col = "state",
                        tip_offset = 1,
                        tip_hjust = 1,
                        hotspot_shape = 21,
                        hotspot_size = .25) {

  if (add_tree) {
    tree_data = get_tree_data(x, full_tree, add_root, branch_lengths)
    tree_ggplot = plot_tree(
      tree = tree_data,
      ladderize = ladderize,
      legend.position = legend.position,
      annotate_tip = annotate_tip,
      annotate_bfb_events = annotate_bfb_events,
      annotate_hotspot_amplifications = annotate_hotspot_amplifications,
      bfb_event_colors = bfb_event_colors,
      hotspot_colors = hotspot_colors,
      tip_align = tip_align,
      tip_offset = tip_offset,
      tip_hjust = tip_hjust,
      hotspot_shape = hotspot_shape,
      hotspot_size = hotspot_size
    ) + ggplot2::theme_void()
    ordered_cell_ids <- get_ordered_cell_ids(tree_ggplot$data)
    tree_hm <- make_corrupt_tree_heatmap(tree_ggplot, tree_width = tree_width)
  }

  # Create copy number matrix
  mat = cells2mat(x$cells, x$input_parameters$initial_sequence_length, order = order_heatmap)

  if (add_gain_loss_profile) {
    gain_prop <- colMeans(mat > 1) # Proportion of amplified bins
    loss_prop <- colMeans(mat < 1) # Proportion of lost bins
    # Create a stacked bar annotation
    top_anno <- ComplexHeatmap::HeatmapAnnotation(
      Gain = ComplexHeatmap::anno_barplot(
        gain_prop,
        bar_width = 1,
        beside = FALSE,
        gp = grid::gpar(fill = "indianred", col="indianred"),  # Colors for below and above x
        border = FALSE, ylim = c(0,1)
      ),
      Loss = ComplexHeatmap::anno_barplot(
        -loss_prop,
        beside = FALSE,
        bar_width = 1,
        gp = grid::gpar(fill = "steelblue", col="steelblue"),  # Colors for below and above x
        border = FALSE, ylim = c(-1,0)
      ),
      annotation_name_side = "right"
    )
  }

  if (plot_col == "state") {
    high_cn_mask = mat > 10
    mat_row_names = rownames(mat)
    mat = matrix(as.character(mat), ncol = x$input_parameters$initial_sequence_length, nrow = length(x$cells))
    mat[high_cn_mask] = "11+"
    mat = data.frame(mat)
    rownames(mat) = mat_row_names
    colvals <- get_colors("CN")
  } else if (plot_col == "copy") {
    mat = data.frame(mat)
    colvals <- circlize::colorRamp2(seq(0, 11, 1), get_colors("copy"))
  } else {
    stop("plot_col not recognized. Must be either 'copy' or 'state'")
  }

  if (add_tree) {
    mat = mat[ordered_cell_ids, , drop = FALSE]
  }

  if (add_gain_loss_profile) {
    copynumber_hm = ComplexHeatmap::Heatmap(
      name = "Copy Number",
      as.matrix(mat),
      col = colvals,
      show_column_names = FALSE,
      show_row_names = FALSE,
      use_raster = use_raster,
      raster_quality = raster_quality,
      top_annotation = top_anno,
      cluster_rows = F,
      cluster_columns = F
    )
  } else {
    copynumber_hm = ComplexHeatmap::Heatmap(
      name = "Copy Number",
      as.matrix(mat),
      col = colvals,
      show_column_names = FALSE,
      show_row_names = FALSE,
      use_raster = use_raster,
      raster_quality = raster_quality,
      cluster_rows = F,
      cluster_columns = F
    )
  }

  if (add_tree) {
    p = tree_hm + copynumber_hm
  } else {
    p = copynumber_hm
  }
  p
}

