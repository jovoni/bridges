# Build a parameter list for clonal simulations

Returns a complete parameter list with sensible defaults. Any argument
can be overridden. Pass the result to
[`run_clonal_replicates()`](https://jovoni.github.io/bridges/reference/run_clonal_replicates.md)
or
[`compare_clonal_params()`](https://jovoni.github.io/bridges/reference/compare_clonal_params.md).

## Usage

``` r
clonal_params(
  bfb_prob = 0.5,
  amp_rate = 0.3,
  del_rate = 0.2,
  positive_selection_rate = 0,
  negative_selection_rate = 0,
  birth_rate = 1,
  death_rate = 0.1,
  lambda = 2,
  rate = 20,
  first_round_of_bfb = TRUE
)
```

## Arguments

- bfb_prob:

  Relative probability of a BFB event per division.

- amp_rate:

  Relative probability of a focal amplification.

- del_rate:

  Relative probability of a focal deletion.

- positive_selection_rate:

  Multiplicative birth-rate boost for cells that have gained a copy at
  `hotspot_pos`. 0 = neutral.

- negative_selection_rate:

  Multiplicative death-rate boost for cells that have gained a copy at
  `hotspot_pos`. 0 = neutral.

- birth_rate:

  Base cell birth rate.

- death_rate:

  Base cell death rate.

- lambda:

  Mean number of genomic events per daughter cell per division
  (Poisson).

- rate:

  Scale of focal amp/del events: segment length is drawn from
  Exp(1/rate).

- first_round_of_bfb:

  If TRUE the founding cell already carries one BFB event, mimicking a
  clone that initiated from a single BFB.

## Value

Named list of simulation parameters.
