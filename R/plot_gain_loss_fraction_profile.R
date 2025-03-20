
plot_gain_loss_fraction_profile = function(cells, L) {
  CN_runs = multi_seqs2CNruns(cells, L)

  colors = get_colors("gainloss")

  CN_runs %>%
    dplyr::rowwise() %>%
    dplyr::do({
      # Get the range of values for this interval
      x_values <- seq(from = .$segment_start, to = .$segment_end)

      # Create the output data frame with repeated values based on count
      dplyr::tibble(
        bin = x_values,
        count = .$copy_number,
        cell_idx = .$cell_idx
      )
    }) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      Loss = -sum(count < 1) / length(count),
      Gain = sum(count > 1) / length(count)
    ) %>%
    tidyr::pivot_longer(!(bin)) %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x=bin, y=value, fill=name)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Chromosome bin", y="Fraction", fill="") +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_y_continuous(limits = c(-1,1))
}


