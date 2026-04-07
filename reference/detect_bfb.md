# Detect Breakpoint-Free Branches (BFB) in Phylogenetic Trees

This function performs binomial tests across all chromosomes and alleles
to detect breakpoint-free branches in a phylogenetic tree
reconstruction. It uses pseudo-cell analysis to identify branches where
copy number changes occur without breakpoints.

## Usage

``` r
detect_bfb(fit, threshold = 0.005)
```

## Arguments

- fit:

  A fitted model object containing:

  - all_input_Xs: List of chromosome data with allele information

  - b_dist_func: B distance function

  - g_dist_func: G distance function

  - tree: Phylogenetic tree structure

- threshold:

  Numeric threshold for binomial test probability (default: 0.005)

## Value

A data frame containing binomial test results for each chromosome-allele
combination:

- mean: Estimated proportion of successes

- N_trials: Total number of trials

- N_successes: Number of successes

- chr: Chromosome identifier

- allele: Allele type

- p.value: P-value from binomial test

- threshold: Threshold used for testing

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming 'fitted_model' is a properly fitted model object
bfb_results <- detect_bfb(fitted_model, threshold = 0.01)
} # }
```
