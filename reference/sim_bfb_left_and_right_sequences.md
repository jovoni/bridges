# Simulate Breakage-Fusion-Bridge (BFB) Cycle for both daughters

Simulate Breakage-Fusion-Bridge (BFB) Cycle for both daughters

## Usage

``` r
sim_bfb_left_and_right_sequences(
  sequence,
  support = "uniform",
  alpha = NULL,
  beta = NULL,
  custom_breakpoints = NULL
)
```

## Arguments

- sequence:

  Input sequence

- support:

  Distribution type for breakpoint selection ("uniform" or "beta")

- alpha:

  Shape parameter for beta distribution (only used if support="beta")

- beta:

  Shape parameter for beta distribution (only used if support="beta")

- custom_breakpoints:

  .

## Value

List containing left and right sequences

## Details

Simulate left and right children from a BFB cycle using specified
breakpoint selection distribution
