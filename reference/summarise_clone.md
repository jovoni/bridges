# Compute summary statistics for a single clonal CN matrix

Compute summary statistics for a single clonal CN matrix

## Usage

``` r
summarise_clone(cna_matrix, hotspot_col, base_value = 1, n_alive = NA)
```

## Arguments

- cna_matrix:

  Integer matrix (cells x bins) from
  [`simulate_clone()`](https://jovoni.github.io/bridges/reference/simulate_clone.md).

- hotspot_col:

  Column index of the hotspot bin in `cna_matrix`.

- base_value:

  Baseline (diploid) copy number: 1 for allele-specific, 2 for CN.

- n_alive:

  Number of alive cells before subsampling (from
  `simulate_clone()$n_alive`). Stored as-is for downstream analysis.

## Value

Named numeric vector of summary statistics.
