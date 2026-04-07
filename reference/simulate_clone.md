# Simulate one clonal evolution and return the CN matrix

Runs
[`bridge_sim()`](https://jovoni.github.io/bridges/reference/bridge_sim.md)
until `max_cells` are alive OR `max_time` is reached (whichever fires
first). If `sample_cells` is set, a random subsample of that many cells
is drawn from the alive population using `subsample_sim()`.

## Usage

``` r
simulate_clone(
  chromosome,
  allele,
  hotspot_pos,
  params = clonal_params(),
  max_cells = 1000,
  max_time = 300,
  sample_cells = NULL,
  bin_length = 1e+06
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

- max_cells:

  Stop when this many cells are alive. Default 1000 (effectively no
  limit when the stopping criterion is time-based).

- max_time:

  Stop at this simulation time. Default 300.

- sample_cells:

  Number of cells to randomly sample from the alive population at the
  end of the simulation. NULL keeps all alive cells.

- bin_length:

  Bin size in bp. Default 1 Mb.

## Value

A list with:

- cna_matrix:

  Integer matrix (cells x bins) for the target allele, after
  subsampling.

- n_alive:

  Number of cells alive at end of simulation, before sampling.

- n_cells:

  Number of cells in `cna_matrix` (after sampling).

- sim:

  Full
  [`bridge_sim()`](https://jovoni.github.io/bridges/reference/bridge_sim.md)
  output (for downstream use).
