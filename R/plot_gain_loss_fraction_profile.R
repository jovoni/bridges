
plot_gain_loss_fraction_profile = function(x) {
  cells = x$cells
  L = x$input_parameters$initial_sequence_length

  CN_runs = multi_seqs2CNruns(cells, L)

  colors = get_colors("gainloss")

  CN_runs %>%
    dplyr::rowwise() %>%
    dplyr::do({
      # Get the range of values for this interval
      x_values <- seq(from = .data$segment_start, to = .data$segment_end)

      # Create the output data frame with repeated values based on count
      dplyr::tibble(
        bin = x_values,
        count = .data$copy_number,
        cell_idx = .data$cell_idx
      )
    }) %>%
    dplyr::group_by(.data$bin) %>%
    dplyr::summarise(
      Loss = -sum(.data$count < 1) / length(.data$count),
      Gain = sum(.data$count > 1) / length(.data$count)
    ) %>%
    tidyr::pivot_longer(!(.data$bin)) %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x=.data$bin, y=.data$value, fill=.data$name)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Chromosome bin", y="Fraction", fill="") +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_y_continuous(limits = c(-1,1))
}


