# Derive hotspot bin position from a CN matrix

Returns the 1-based bin index of the genomic region most likely to be
the BFB hotspot, defined as the bin with the highest CN variance across
cells. This is data-driven and removes the need for the caller to supply
`pos`.

## Usage

``` r
derive_hotspot_pos(cna_matrix)
```

## Arguments

- cna_matrix:

  Numeric matrix (cells x bins).

## Value

Integer bin index.
