# Process a birth event (diploid version)

This function handles cell division by creating two daughter cells from
a parent cell. During division, various genomic events (amplifications,
deletions, BFB) can occur based on the specified rates and lambda
parameter.

## Usage

``` r
process_birth_event(state, current_cell_id, cell_idx, lambda, rate)
```

## Arguments

- state:

  The simulation state containing cell information, sequences, and
  parameters

- current_cell_id:

  ID of the cell undergoing birth/division

- cell_idx:

  Index of the cell in the alive arrays (from get_next_event)

- lambda:

  Rate parameter for Poisson distribution used to sample the number of
  genomic events per daughter cell

- rate:

  Rate parameter used in amplification/deletion simulations. Length of
  event is sample from exponential distribution with parameter 1 / rate.

## Value

Updated simulation state with new daughter cells and updated history
