# Plot summary statistics across clonal simulation conditions

Produces a faceted box-plot grid: one panel per summary statistic, one
box per condition. If `observed_stats` is supplied (output of
[`compute_observed_stats()`](https://jovoni.github.io/bridges/reference/compute_observed_stats.md)),
a horizontal dashed line is added to each panel showing the observed
value.

## Usage

``` r
plot_clonal_comparison(
  results,
  observed_stats = NULL,
  stats = c("hotspot_mean_cn", "hotspot_max_cn", "hotspot_fraction_amp",
    "hotspot_cn_var", "max_cn_global", "mean_breakpoints", "n_alive"),
  fill_var = "condition"
)
```

## Arguments

- results:

  Output of
  [`compare_clonal_params()`](https://jovoni.github.io/bridges/reference/compare_clonal_params.md).

- observed_stats:

  Optional named numeric vector from
  [`compute_observed_stats()`](https://jovoni.github.io/bridges/reference/compute_observed_stats.md).
  A horizontal reference line is drawn in each panel for statistics
  present in both `results` and `observed_stats`.

- stats:

  Character vector of statistic names to plot.

- fill_var:

  Column name to use for box fill colour. Default `"condition"`.

## Value

A `ggplot` object.
