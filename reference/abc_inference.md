# ABC inference of BFB simulation parameters

## Usage

``` r
abc_inference(
  cna_data,
  allele,
  chromosome,
  n_simulations = 10000,
  tolerance_quantile = 0.01,
  n_cores = 1,
  pos = NULL,
  bin_length = NULL
)
```

## Arguments

- cna_data:

  Data frame with columns cell_id, chr, start, end, CN, A, B.

- allele:

  Allele to match ("A", "B", or "CN").

- chromosome:

  Chromosome to focus on (character, e.g. "7").

- n_simulations:

  Total number of simulations to run. Default: 10000.

- tolerance_quantile:

  Fraction of simulations to accept. Default: 0.01.

- n_cores:

  Number of parallel cores (uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) when
  \> 1). Default: 1.

- pos:

  Optional integer bin index for the BFB hotspot in the simulation.
  Automatically derived from the observed data as the bin with the
  highest CN variance if NULL (default).

- bin_length:

  Optional bin size in bp. Inferred from observed data if NULL.

## Value

A list containing:

- accepted_params:

  Data frame of accepted parameter draws.

- param_summary:

  Per-parameter mean / median / SD / 95\\ n_simulationsTotal simulations
  attempted. n_acceptedNumber of accepted simulations.
  acceptance_rateFraction accepted out of valid (non-error) simulations.
  tolerance_thresholdDistance cutoff used for acceptance.

Uses rejection ABC to infer simulation parameters that reproduce the
copy number distribution of an observed single-cell dataset. The
simulator (`bridge_sim`) is run `n_simulations` times with parameters
drawn from
[`sample_priors()`](https://jovoni.github.io/bridges/reference/sample_priors.md).
The `tolerance_quantile` fraction with the smallest distance to the
observed summary statistics is retained as the approximate posterior.
