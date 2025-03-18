
plot_heatmap_CN_profiles = function(cells, L, order = TRUE) {
  CNruns = multi_seqs2CNruns(cells, L)

  CNruns = CNruns %>%
    dplyr::group_by(cell_idx) %>%
    dplyr::mutate(diploid_percentage = (copy_number == 1) * (segment_end - segment_start) / L)

  if (order) {CNruns = CNruns %>% dplyr::arrange(-diploid_percentage)}

  CNruns$cell_idx = factor(CNruns$cell_idx, unique(CNruns$cell_idx))

  colors = get_colors("CN")
  CNruns %>%
    dplyr::mutate(copy_number = ifelse(copy_number > 10, "11+", as.character(copy_number))) %>%
    dplyr::mutate(copy_number = factor(copy_number, levels = names(colors))) %>%
    dplyr::mutate(cell_idx = as.numeric(cell_idx)) %>%
    ggplot2::ggplot(
      mapping = ggplot2::aes(
        xmin=segment_start-0.5,
        xmax=segment_end+0.5,
        ymin=cell_idx-0.5,
        ymax=cell_idx+0.5,
        fill=copy_number,
        col=copy_number)
      ) +
    ggplot2::geom_rect() +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Bins", y="Cell")
}
