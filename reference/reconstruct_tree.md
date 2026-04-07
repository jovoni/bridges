# Reconstruct Phylogenetic Tree with Copy Number Analysis

This function reconstructs a phylogenetic tree by traversing internal
nodes and computing distances between cell profiles. For each internal
node, it calculates both B distance (breakpoint-aware) and G distance
(general) between child nodes, then selects the reconstruction method
that minimizes cost.

## Usage

``` r
reconstruct_tree(fit, chr, allele, store_profiles = TRUE)
```

## Arguments

- fit:

  A fitted model object containing:

  - all_input_Xs: List of input data matrices by chromosome and allele

  - b_dist_func: B distance function for breakpoint-aware reconstruction

  - g_dist_func: G distance function for general reconstruction

  - tree: Phylogenetic tree structure (ape format)

- chr:

  Character string specifying the chromosome to analyze

- allele:

  Character string specifying the allele type ("CN" for copy number,
  others for allelic)

## Value

A list containing:

- internal_nodes: List of reconstructed internal node profiles

- deltas: List of delta values (G cost - B cost) for each internal node

- merged_profiles: List of merged profiles when B distance is chosen

- final_input: Final working input matrix after all reconstructions

- processed_nodes: Number of internal nodes processed

- chr: Chromosome identifier

- allele: Allele type

## Details

The function:

1.  Sets target values based on allele type (2 for "CN", 1 for others)

2.  Processes internal nodes in post-order traversal

3.  For each internal node, computes B and G distances between children

4.  Calculates delta = G_cost - B_cost

5.  Chooses B reconstruction if delta \> 0, otherwise G reconstruction

6.  Creates pseudo-cells for internal nodes and updates working dataset

## Note

The tree edge lengths are set to 1 for uniform weighting. Node labels
are automatically generated if not present.
