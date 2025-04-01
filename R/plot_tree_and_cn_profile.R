
plot_tree_and_cn_profile = function(
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
    hotspot_size = 2,
    matrix_width = 3,
    matrix_offset = 2,
    plot_gain_loss = FALSE
    ) {
  # Prepare tree plot
  tree_plot = plot_tree(x,
                        full_tree,
                        legend.position,
                        annotate_bfb_events,
                        annotate_hotspot_amplifications,
                        bfb_event_colors,
                        hotspot_colors,
                        tip_align,
                        tip_offset,
                        tip_hjust,
                        use_computed_branch_lengths,
                        hotspot_shape,
                        hotspot_size
                        )
  tree_plot = tree_plot + ggnewscale::new_scale_fill()

  # Get matrix of CN profile
  mat = cells2mat(x$cells, x$input_parameters$initial_sequence_length, order = F)
  high_cn_mask = mat > 10
  mat = matrix(as.character(mat), ncol = x$input_parameters$initial_sequence_length, nrow = length(x$cells))
  mat[high_cn_mask] = "11+"
  rownames(mat) = x$cell_history %>% dplyr::filter(.data$is_alive) %>% dplyr::pull(.data$cell_id)
  mat = data.frame(mat)

  p = ggtree::gheatmap(tree_plot, mat, offset=matrix_offset, width=matrix_width, colnames = FALSE) +
    ggplot2::scale_fill_manual(values = get_colors(set = "CN")) +
    ggplot2::labs(fill = "Copy Number") +
    ggplot2::theme(legend.position = legend.position)

  if (plot_gain_loss) {
    p_upper = patchwork::plot_spacer() + plot_gain_loss_fraction_profile(x) +
      patchwork::plot_layout(ncol = 2, widths = c(1, matrix_width))
    p = p_upper / p + patchwork::plot_layout(ncol = 1, heights = c(1, 4))
  }

  p
}
