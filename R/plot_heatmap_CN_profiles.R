#
# plot_heatmap_CN_profiles = function(x, order = TRUE) {
#   cells = x$cells
#   L = x$input_parameters$initial_sequence_length
#   CNruns = multi_seqs2CNruns(cells, L)
#
#   CNruns = CNruns %>%
#     dplyr::group_by(.data$cell_idx) %>%
#     dplyr::mutate(diploid_percentage = (.data$copy_number == 1) * (.data$segment_end - .data$segment_start) / L)
#
#   if (order) {CNruns = CNruns %>% dplyr::arrange(-.data$diploid_percentage)}
#
#   CNruns$cell_idx = factor(CNruns$cell_idx, unique(CNruns$cell_idx))
#
#   colors = get_colors("CN")
#   CNruns %>%
#     dplyr::mutate(copy_number = ifelse(.data$copy_number > 10, "11+", as.character(.data$copy_number))) %>%
#     dplyr::mutate(copy_number = factor(.data$copy_number, levels = names(colors))) %>%
#     dplyr::mutate(cell_idx = as.numeric(.data$cell_idx)) %>%
#     ggplot2::ggplot(
#       mapping = ggplot2::aes(
#         xmin=.data$segment_start-0.5,
#         xmax=.data$segment_end+0.5,
#         ymin=.data$cell_idx-0.5,
#         ymax=.data$cell_idx+0.5,
#         fill=.data$copy_number,
#         col=.data$copy_number)
#       ) +
#     ggplot2::geom_rect() +
#     ggplot2::scale_fill_manual(values = colors) +
#     ggplot2::scale_color_manual(values = colors) +
#     ggplot2::theme_minimal() +
#     ggplot2::labs(x = "Bins", y="Cell")
# }
