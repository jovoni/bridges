# Compare summary statistics across multiple parameter sets

Runs
[`run_clonal_replicates()`](https://jovoni.github.io/bridges/reference/run_clonal_replicates.md)
for each entry in `param_list` and combines the results into a single
tidy data frame. The name of each list entry becomes the `condition`
column, making it easy to facet or colour by condition in downstream
plots.

## Usage

``` r
compare_clonal_params(
  param_list,
  chromosome,
  allele,
  hotspot_pos,
  N_replicates = 100,
  max_cells = 1000,
  max_time = 300,
  sample_cells = NULL,
  bin_length = 1e+06,
  n_cores = 1
)
```

## Arguments

- param_list:

  Named list of parameter lists, each from
  [`clonal_params()`](https://jovoni.github.io/bridges/reference/clonal_params.md).
  Names become the `condition` column.

- chromosome:

  Chromosome to simulate.

- allele:

  Allele to track.

- hotspot_pos:

  1-based bin index of the position of interest.

- N_replicates:

  Number of independent replicates per condition.

- max_cells:

  Stopping criterion: maximum alive cells.

- max_time:

  Stopping criterion: simulation time limit.

- sample_cells:

  Number of cells to sample after simulation. NULL keeps all.

- bin_length:

  Bin size in bp.

- n_cores:

  Cores for parallelism (applies within each condition).

## Value

A `tibble` with all replicates across all conditions. Columns:
`condition`, `replicate`, summary statistics, `n_alive`, `failed`.

## Examples

``` r
if (FALSE) { # \dontrun{
results = compare_clonal_params(
  param_list = list(
    neutral    = clonal_params(positive_selection_rate = 0),
    selection  = clonal_params(positive_selection_rate = 2),
    bfb_only   = clonal_params(bfb_prob = 0.9, positive_selection_rate = 0),
    bfb_select = clonal_params(bfb_prob = 0.9, positive_selection_rate = 2)
  ),
  chromosome   = "7",
  allele       = "A",
  hotspot_pos  = 60,
  N_replicates = 100,
  max_time     = 20,
  sample_cells = 50
)
plot_clonal_comparison(results)
} # }
```
