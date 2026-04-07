# Run a single ABC simulation and compute distance to observed data

Run a single ABC simulation and compute distance to observed data

## Usage

``` r
run_ABC(cna_data, allele, chromosome, params, pos = NULL, bin_length = NULL)
```

## Arguments

- cna_data:

  Observed CNA data frame (cell_id, chr, start, end, CN, A, B).

- allele:

  Allele column to use ("A", "B", or "CN").

- chromosome:

  Chromosome name as it appears in cna_data\$chr (e.g. "7").

- params:

  Named list of simulation parameters from
  [`sample_priors()`](https://jovoni.github.io/bridges/reference/sample_priors.md).

- pos:

  Optional integer bin index for the BFB hotspot. Derived automatically
  from the observed data if NULL (default).

- bin_length:

  Optional bin size in bp. Inferred from observed data if NULL.

## Value

Named list with distance components and the parameter set used.
