# Plot Chromosome-wide Heatmap of Copy Number Alterations

This function plots a heatmap of copy number alterations (CNA) for a
specified chromosome and allele across multiple cells. It supports
optional ordering of cells, as well as annotations showing average copy
number or gain/loss profiles.

## Usage

``` r
plot_chr_all_heatmap(
  cna_data,
  chr,
  allele,
  order_heatmap = TRUE,
  add_avg_CN_profile = TRUE,
  add_gain_loss_profile = TRUE,
  use_raster = TRUE,
  raster_quality = 15
)
```

## Arguments

- cna_data:

  A data frame containing copy number data. Must include columns:
  `cell_id`, chromosome, bin positions, and a value column specified by
  `allele`.

- chr:

  Chromosome to plot (e.g., `"1"`, `"X"`).

- allele:

  Column name in `cna_data` containing the copy number values to be
  plotted (e.g., `"CN"`, `"cn_a"`).

- order_heatmap:

  Logical, whether to order the heatmap rows by the number of bins
  matching a reference value (2 for total CN, 1 for alleles).

- add_avg_CN_profile:

  Logical, whether to add a barplot annotation showing the average CN
  profile across cells.

- add_gain_loss_profile:

  Logical, whether to add gain/loss annotation bars (proportion of cells
  with gain/loss).

- use_raster:

  Logical, whether to rasterize the heatmap for faster rendering in
  large datasets.

- raster_quality:

  Integer, quality of the raster image (relevant only if
  `use_raster = TRUE`).

## Value

A ComplexHeatmap object representing the CNA heatmap, which can be
plotted or combined with other heatmaps.
