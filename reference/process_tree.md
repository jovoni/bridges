# Process tree for heatmap visualization

Process tree for heatmap visualization

## Usage

``` r
process_tree(
  tree,
  branch_length = NULL,
  ladderize = FALSE,
  reorder_tree = FALSE,
  distance_matrix = NULL
)
```

## Arguments

- tree:

  Phylogenetic tree object

- branch_length:

  Uniform branch length (optional)

- ladderize:

  Whether to ladderize the tree

- reorder_tree:

  Whether to optimize tree ordering

- distance_matrix:

  Distance matrix for reordering (reorder_tree = TRUE)

## Value

List containing processed tree and ggplot object
