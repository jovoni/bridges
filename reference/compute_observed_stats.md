# Compute summary statistics from a real CNA sample

Applies the same statistics used in
[`summarise_clone()`](https://jovoni.github.io/bridges/reference/summarise_clone.md)
to an observed data frame so the result can be overlaid on simulation
comparison plots.

## Usage

``` r
compute_observed_stats(data, chromosome, allele, hotspot_pos)
```

## Arguments

- data:

  Data frame with columns `cell_id`, `chr`, `start`/`from`, `end`/`to`,
  and the allele column.

- chromosome:

  Chromosome to analyse (character, e.g. `"7"`).

- allele:

  Allele column to use (`"A"`, `"B"`, or `"CN"`).

- hotspot_pos:

  1-based bin index of the position of interest.

## Value

Named numeric vector of summary statistics, compatible with
`plot_clonal_comparison(observed_stats = ...)`.

## Details

Accepts both `start`/`end` and `from`/`to` column names.
