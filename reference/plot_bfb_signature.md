# Plot BFB-like Signatures Along the Genome

This function visualizes the distribution of BFB-like events (asymmetric
CN differences) across genomic segments for a given chromosome and
allele. The plot highlights which genomic regions are most frequently
affected by potential BFB events in the reconstructed phylogeny.

## Usage

``` r
plot_bfb_signature(res, chr_of_interest, allele_of_interest)
```

## Arguments

- res:

  A fitted object returned by the
  [`fit()`](https://jovoni.github.io/bridges/reference/fit.md) function.
  It must contain the `reconstructions` field populated via
  `compute_reconstructions()`.

- chr_of_interest:

  A character string indicating the chromosome to plot (e.g., `"8"` or
  `"X"`).

- allele_of_interest:

  A character string indicating the allele to analyze, typically `"A"`
  or `"B"`.

## Value

A `ggplot` object displaying a bar plot where each bar corresponds to a
genomic segment, with its height and color intensity proportional to the
number of BFB-like events affecting that segment.

## Details

The function extracts the reconstructed merged profiles and identifies
the internal nodes in the phylogeny that show asymmetric copy number
changes (i.e., delta \> 0). It computes the difference between the left
and right child profiles for these nodes, identifies affected segments,
and counts the number of times each segment is involved in a BFB-like
event. The output is a genomic bar plot showing the number of BFB-like
events per segment.
