# Main Function for Genomic Distance Calculation and Phylogenetic Tree Construction

This function computes genomic distances between samples with flexible
distance-function selection, processes the data, and constructs a
phylogenetic tree. It supports distinct metrics for greedy (G) and
breakage–fusion–bridge (BFB; B) distance calculations.

## Usage

``` r
fit(
  data,
  chromosomes = c(1:22, "X", "Y"),
  alleles = c("A", "B"),
  k_jitter_fix = 0,
  bfb_penalty = 0,
  tree_func = ape::nj,
  fillna = 0,
  g_dist_func = "greedy_fast",
  b_dist_func = "bfb_fast",
  ...
)
```

## Arguments

- data:

  Input data frame with columns:

  - `cell_id` — Cell identifier

  - `chr` — Chromosome name

  - `start` — Start of genomic bin

  - `end` — End of genomic bin

  - `CN` — Total copy number

  - `A` — Copy number of allele A

  - `B` — Copy number of allele B

- chromosomes:

  Vector of chromosomes to include (default: `c(1:22, "X", "Y")`).

- alleles:

  Alleles to consider (default: `c("A","B")`; alternative: `"CN"`).

- k_jitter_fix:

  Numeric jitter factor for numerical stability (default: `0`).

- bfb_penalty:

  Penalty to apply to BFB events (default: `0`). Ignored if
  `b_dist_func = NULL`.

- tree_func:

  Function used to build the tree from the final distance matrix
  (default: [`ape::nj`](https://rdrr.io/pkg/ape/man/nj.html)).

- fillna:

  Value to fill `NA` entries in pre-processing (default: `0`).

- g_dist_func:

  Name of the greedy (G) distance function. Must be one of
  `names(G_DISTS)` (default: `"greedy_fast"`).

- b_dist_func:

  Name of the BFB (B) distance function. Must be one of
  `names(B_DISTS)`. **Set to `NULL` to disable the BFB stage** and run a
  greedy-only analysis (default: `"bfb_fast"`).

- ...:

  Additional arguments passed to downstream helpers.

## Value

A list with:

- `tree` — The constructed phylogenetic tree (rooted, diploid removed)

- `all_input_Xs` — Processed input data

- `D` — Final distance matrix

- `greedy_Ds` — Per-chromosome/allele greedy distance matrices

- `avg_Ds` — Per-chromosome/allele balanced (G vs B) distance matrices;
  **if `b_dist_func = NULL`, then `avg_Ds` equals `greedy_Ds`**

- `g_dist_func` — Name of the G distance function used

- `b_dist_func` — Name of the B distance function used (or `NULL`)

## Details

The pipeline performs:

1.  Validation of selected distance functions

2.  Input pre-processing and diploid augmentation

3.  Computation of greedy (G) distances

4.  **Optional BFB stage**: if `b_dist_func` is provided, BFB (B)
    distances are computed and merged with G; if `b_dist_func = NULL`,
    the BFB stage is *skipped* and `avg_Ds <- greedy_Ds` (greedy-only
    analysis). In the `NULL` case, `bfb_penalty` is ignored.

5.  Merging to minimal distances and optional allele summation

6.  Tree construction via `tree_func`

Available functions can be inspected with `names(G_DISTS)` and
`names(B_DISTS)`. Using `b_dist_func = NULL` is useful for ablation
studies, speed-ups, or when BFB modelling is not desired; results reduce
to the greedy metric.
