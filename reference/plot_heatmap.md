# Create comprehensive heatmap with optional phylogenetic tree

This function creates a heatmap visualization with optional phylogenetic
tree display, annotations, and support for multiple features. It can
display genomic data across chromosomes with customizable tree coloring
based on reconstruction data.

## Usage

``` r
plot_heatmap(
  data,
  tree = NULL,
  chromosomes_to_plot = c(1:22, "X", "Y"),
  to_plot = "CN",
  tree_width = 3.5,
  show_tree = TRUE,
  branch_length = NULL,
  annotations = NULL,
  ladderize = FALSE,
  reorder_tree = TRUE,
  distance_matrix = NULL,
  use_raster = TRUE,
  raster_quality = 10,
  reconstruction = NULL,
  chr_for_coloring = NULL,
  allele_for_coloring = NULL,
  tree_colors = c("black", "goldenrod"),
  node_size = 1,
  tip_size = 0
)
```

## Arguments

- data:

  Data frame containing cell_id, chr, start, and feature columns. Must
  include the columns specified in `to_plot` parameter.

- tree:

  Phylogenetic tree object (optional). Should be compatible with ape or
  phylo class objects.

- chromosomes_to_plot:

  Vector of chromosomes to include in the plot. Default: c(1:22, "X",
  "Y") for human chromosomes.

- to_plot:

  Vector of feature names to plot. Default: "CN" (copy number). Must
  correspond to column names in the data.

- tree_width:

  Width of tree display in centimeters. Default: 3.5.

- show_tree:

  Logical indicating whether to display the phylogenetic tree. Default:
  TRUE.

- branch_length:

  Numeric value for uniform branch length transformation of the tree. If
  NULL (default), original branch lengths are preserved.

- annotations:

  Data frame with cell_id column and additional annotation columns. Each
  non-cell_id column becomes a row annotation in the heatmap. Optional
  parameter.

- ladderize:

  Logical indicating whether to ladderize (sort) the tree branches.
  Default: FALSE.

- reorder_tree:

  Logical indicating whether to optimize tree tip ordering for better
  visualization. Default: TRUE.

- distance_matrix:

  Distance matrix for tree reordering optimization. Optional parameter
  used when reorder_tree is TRUE.

- use_raster:

  Logical indicating whether to use raster graphics for the heatmap.
  Default: TRUE for better performance with large datasets.

- raster_quality:

  Numeric value controlling raster image quality. Higher values produce
  better quality. Default: 10.

- reconstruction:

  Reconstruction data for coloring tree branches based on ancestral
  state reconstruction. Optional parameter.

- chr_for_coloring:

  Chromosome identifier used for tree branch coloring when
  reconstruction data is provided. Optional parameter.

- allele_for_coloring:

  Specific allele used for tree branch coloring when reconstruction data
  is provided. Optional parameter.

- tree_colors:

  Vector of colors for tree branch coloring. Default: c("black",
  "goldenrod"). First color for one state, second for another.

- node_size:

  Size of internal nodes in the tree. Default: 1.

- tip_size:

  Size of tip nodes in the tree visualization. Default: 0 (tips not
  displayed).

## Value

A ComplexHeatmap object that can be displayed or further customized. The
object contains the combined heatmap(s) and optional tree annotation.

## Details

The function performs the following main steps:

- Validates input parameters and data structure

- Processes the phylogenetic tree (if provided) with optional coloring

- Prepares row annotations from the annotations data frame

- Creates individual heatmaps for each feature in to_plot

- Combines multiple heatmaps horizontally

- Adds tree annotation to the left side (if requested)

The tree can be colored based on reconstruction data, which is useful
for visualizing ancestral state reconstructions or other phylogenetic
analyses. When reconstruction data is provided, the tree branches are
colored according to the specified parameters.
