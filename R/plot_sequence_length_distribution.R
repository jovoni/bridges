plot_sequence_length_distribution = function(L0, n, support, Nrep = 1000, bins=100, alpha=NULL, beta=NULL) {
  # Calculate the distribution of sequence lengths after n breakage-fusion-bridge cycles
  #
  # Args:
  #   L0: Initial sequence length
  #   n: Number of breakage-fusion-bridge cycles to simulate
  #   support: Type of distribution to use ("uniform" or "beta")
  #   Nrep: Number of repetitions for the simulation (default: 1000)
  #   ...: Additional parameters (e.g., alpha and beta for beta distribution)
  #
  # Returns:
  #   A vector of sequence lengths from Nrep simulations
  xs = lapply(1:Nrep, function(i) {
    sim_n_bfb(L0, n, support, alpha, beta) %>% length()
  }) %>% unlist()

  dplyr::tibble(L = xs) %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x=L)) +
    ggplot2::geom_histogram(bins = bins) +
    ggplot2::theme_bw()
}
