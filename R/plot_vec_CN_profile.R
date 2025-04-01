# Plot vector CN profile
plot_vec_CN_profile = function(vec) {
  vec_table = table(vec)
  bins = as.numeric(names(vec_table))
  CN = as.numeric(vec_table)

  df = dplyr::tibble(bin=bins, CN=CN)
  df %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x=.data$bin, y=.data$CN, col=factor(.data$CN), group=.data$CN)) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = get_colors("CN")) +
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::labs(x = "bin", y="CN", col="CN")
}
