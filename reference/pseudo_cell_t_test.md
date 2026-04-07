# Perform T-Test on Pseudo-Cell Delta Values

This function reconstructs a phylogenetic tree for a specific chromosome
and allele, then performs a one-sample t-test on the delta values
(difference between G and B distances).

## Usage

``` r
pseudo_cell_t_test(fit, chr, allele, mu = 0, sign = F)
```

## Arguments

- fit:

  A fitted model object containing tree structure and distance functions

- chr:

  Character string specifying the chromosome

- allele:

  Character string specifying the allele type

- mu:

  Numeric value for the null hypothesis mean (default: 0)

- sign:

  Logical indicating whether to use sign of delta values (default:
  FALSE)

## Value

A tibble containing:

- mean: Sample mean of delta values

- chr: Chromosome identifier

- allele: Allele type

- p.value: P-value from t-test (NA if all values are identical)

- greater: Logical indicating if mean \> mu

## Details

If all delta values are identical, no statistical test is performed and
p.value is set to NA.
