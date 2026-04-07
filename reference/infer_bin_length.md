# Infer bin length from a CNA data frame

Computes the median (end - start) for the given chromosome so the
simulator uses the same resolution as the observed data.

## Usage

``` r
infer_bin_length(data, chromosome)
```

## Arguments

- data:

  CNA data frame with columns cell_id, chr, start, end.

- chromosome:

  Chromosome to use (character, e.g. "7").

## Value

Numeric bin length in base pairs.
