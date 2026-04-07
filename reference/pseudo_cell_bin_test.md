# Perform Binomial Test on Pseudo-Cell Delta Values

This function reconstructs a phylogenetic tree and performs a binomial
test to assess whether the proportion of positive delta values exceeds a
given threshold. Delta values represent the difference between G
distance and B distance costs.

## Usage

``` r
pseudo_cell_bin_test(fit, chr, allele, threshold)
```

## Arguments

- fit:

  A fitted model object containing tree structure and distance functions

- chr:

  Character string specifying the chromosome

- allele:

  Character string specifying the allele type

- threshold:

  Numeric threshold probability for the binomial test

## Value

A tibble containing:

- mean: Estimated proportion of positive delta values

- N_trials: Total number of delta values

- N_successes: Number of positive delta values

- chr: Chromosome identifier

- allele: Allele type

- p.value: P-value from binomial test

- threshold: Threshold used for testing

## Details

The binomial test uses alternative = "greater" to test if the proportion
of positive deltas significantly exceeds the threshold.
