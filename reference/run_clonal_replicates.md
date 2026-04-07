# Run N independent clonal simulations and summarise each

Each replicate is an independent run of
[`simulate_clone()`](https://jovoni.github.io/bridges/reference/simulate_clone.md).
The simulation stops when `max_cells` are alive or `max_time` is
reached, then `sample_cells` cells are drawn using `subsample_sim()`.
Failed simulations (e.g. all cells die) are recorded with NA summary
statistics and flagged with `failed = TRUE`.

## Usage

``` r
run_clonal_replicates(
  chromosome,
  allele,
  hotspot_pos,
  params = clonal_params(),
  N_replicates = 100,
  max_cells = 1000,
  max_time = 300,
  sample_cells = NULL,
  bin_length = 1e+06,
  n_cores = 1
)
```

## Arguments

- chromosome:

  Chromosome to simulate (character, e.g. `"7"`).

- allele:

  Allele to track (`"A"`, `"B"`, or `"CN"`).

- hotspot_pos:

  1-based bin index of the position of interest.

- params:

  Parameter list from
  [`clonal_params()`](https://jovoni.github.io/bridges/reference/clonal_params.md).

- N_replicates:

  Number of independent clonal evolutions to simulate.

- max_cells:

  Stopping criterion: maximum alive cells.

- max_time:

  Stopping criterion: simulation time limit.

- sample_cells:

  Number of cells to sample after simulation. NULL keeps all.

- bin_length:

  Bin size in bp. Default 1 Mb.

- n_cores:

  Number of parallel cores (uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) when
  \> 1).

## Value

A `tibble` with one row per replicate and columns for each summary
statistic plus `replicate`, `n_alive`, and `failed`.
