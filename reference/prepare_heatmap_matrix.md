# Prepare matrix data for heatmap plotting

Prepare matrix data for heatmap plotting

## Usage

``` r
prepare_heatmap_matrix(
  data,
  chromosomes_to_plot,
  feature_name,
  ordered_cell_ids
)
```

## Arguments

- data:

  Input data frame

- chromosomes_to_plot:

  Chromosomes to include

- feature_name:

  Name of feature column to extract

- ordered_cell_ids:

  Cell IDs in desired order

## Value

List containing matrix and column information
