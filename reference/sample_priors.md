# Sample parameters from prior distributions

Draws one candidate parameter set:

- amp_rate / del_rate / bfb_prob as relative proportions from a
  symmetric Dirichlet(1,1,1) — uniform on the 2-simplex.

- death_rate as a fraction of birth_rate, uniform on (0, 0.8).

- lambda (mean genomic events per daughter cell) uniform on (0.5, 5).

## Usage

``` r
sample_priors()
```

## Value

Named list of parameter values.
