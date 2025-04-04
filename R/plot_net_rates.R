
# Plot net birth rate
plot_net_rates = function(x, max_ploidy = 20) {
  Ns = 0:max_ploidy
  positive_selection_rates = lapply(Ns, function(n) {
    x$input_parameters$positive_selection_function(x$input_parameters$positive_selection_rate, n)
  }) %>% unlist()
  net_birth_rates = x$input_parameters$birth_rate * (1 + positive_selection_rates)

  negative_selection_rates = lapply(Ns, function(n) {
    x$input_parameters$negative_selection_function(x$input_parameters$negative_selection_rate, n)
  }) %>% unlist()
  net_death_rates = x$input_parameters$death_rate * (1 + negative_selection_rates)

  net_growth_rates = net_birth_rates - net_death_rates

  cols = c(
    "Growth rate" = "steelblue",
    "Birth rate" = "forestgreen",
    "Death rate" = "indianred"
  )

  # dplyr::tibble(Ns, "Birth rate"=net_birth_rates, "Death rate"=net_death_rates, "Growth rate"=net_growth_rates) %>%
  #   tidyr::pivot_longer(!Ns) %>%
  #   ggplot2::ggplot(mapping = ggplot2::aes(x=factor(Ns), y=value, col=name, fill=name)) +
  #   ggplot2::geom_col(width = .7, position = "dodge") +
  #   ggplot2::labs(x = "Hotspot copies", y = "Net birth rate") +
  #   ggplot2::scale_fill_manual(values = cols) +
  #   ggplot2::scale_color_manual(values = cols) +
  #   ggplot2::theme_bw()

  dplyr::tibble(Ns, "Birth rate"=net_birth_rates, "Death rate"=net_death_rates, "Growth rate"=net_growth_rates) %>%
    tidyr::pivot_longer(!Ns) %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x=Ns, y=value, col=name, fill=name)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Hotspot CN", y = "Rate", col = "", fill="") +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::theme_bw()
}
