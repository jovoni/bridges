#' Plot the distribution of hotspot copies in a dataset
#'
#' This function generates a histogram of hotspot copy numbers from a given dataset.
#' It filters the dataset to include only living cells and plots the distribution
#' of hotspot copy numbers.
#'
#' @param x A list containing simulation data. It must have elements `input_parameters`
#' and `cell_history`, where `cell_history` includes `hotspot_copies` and `is_alive`.
#'
#' @return A ggplot object showing the histogram of hotspot copy numbers, or NULL if the dataset
#' does not include a hotspot.

#' @export
plot_hotspot_copies = function(x) {
  if (is.null(x$input_parameters$hotspot)) {
    message("Input dataset does not have a hotspot. Returning empty plot")
    return(NULL)
  }

  CN = x$cell_history %>% dplyr::filter(.data$is_alive) %>% dplyr::pull(.data$hotspot_copies)
  max_CN = max(CN)
  dplyr::tibble(CN) %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x=.data$CN)) +
    ggplot2::geom_histogram(bins = max_CN + 1, binwidth = .5) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_continuous(breaks = 0:max_CN) +
    ggplot2::labs(x = "Hotspot CN", y = "Count") +
    ggplot2::ylim(0, NA)

}
