# TODO

## Performance: remaining slow parts of `fit()`

Measured on simulated data, 1,000 cells x 22 chromosomes, alleles A and B
(after the distance rewrite in commits 40eaea7, d4c62d6, 2be07ce). Re-measure
with `Rscript dev/benchmark_fit.R`, on the installed package:
`devtools::load_all()` compiles the C++ code without optimisation (-O0),
which makes it several times slower and distorts benchmarks.

| Stage                    | Time   | Share |
|--------------------------|--------|-------|
| Pre-processing           | 3.1 s  | 54%   |
| Ancestral reconstruction | 1.9 s  | 34%   |
| Neighbour-joining tree   | 0.33 s | 6%    |
| Distances                | 0.21 s | 4%    |
| **Total `fit()`**        | 5.7 s  |       |

The returned object takes **844 MB** (see section 5).

At 400 cells x 6 chromosomes: total 0.75 s = pre-processing 0.58 s,
reconstruction 0.15 s, distances 0.02 s, tree 0.02 s; result 138 MB.

Every change below should keep results identical: compare against the
per-pair R path as in `tests/testthat/test-fit.R`.

### 1. Pre-processing (`process_input_data()`, `R/preprocess.R`)

- [ ] `fit()` appends the synthetic diploid cell with
      `dplyr::bind_rows(data, create_diploid_data(data))`, copying the whole
      input table. Add the diploid row to each matrix instead.
- [ ] `tibble_to_matrix()` runs once per chromosome and per allele, and each
      call copies the table, builds `"start-end"` strings with `paste0()` and
      matches them. Build all chromosome/allele matrices in one
      `data.table` pass, using integer bin indices instead of strings.
- [ ] `reduce_X()` merges columns with `sapply()` + `rowMeans()` per
      segment. Without jitter correction, merged columns are identical, so
      selecting the first column of each segment is enough.
- [ ] `jitter_fix()` (only with `k_jitter_fix > 0`) loops over cells in R.
      Not measured yet.

### 2. Ancestral reconstruction (`reconstruct_tree()`, `R/bfb_mapping_and_detection.R`)

Runs once per chromosome and allele, and loops over internal nodes in R.

- [ ] Each node calls `greedy_distance_cpp()` and `greedy_bfb_distance_cpp()`
      separately, and each call allocates a result list. Move the whole
      post-order loop into C++, one call per chromosome/allele.
- [ ] `get_post_order_internal_nodes()` grows its stack and output with
      `c()` (quadratic). Use `ape::postorder()` /
      `reorder(tree, "postorder")`.
- [ ] Runs are independent across chromosomes/alleles: could run in
      parallel.

### 3. Distances (`distance_matrices_cpp()`, `src/fast_distances.cpp`)

- [ ] Parallelise the pair loop (RcppParallel, or OpenMP). Pairs are
      independent.
- [ ] `compute_fast_distances()` adds each matrix to the total in R
      (`D = D + m$min`), so peak memory is about 3 cells x cells matrices.
      Could accumulate in C++.

### 4. Tree building

- [ ] `ape::nj` grows with the cube of the number of cells: 0.3 s at 1,000
      cells, about 2.7 s at 2,000. Document `tree_func = ape::fastme.bal`
      (or similar) for large data sets, or switch automatically above some
      size. Unlike the items above, this can change the tree.

### 5. Memory of the `fit()` result (high priority)

The result takes 138 MB at 400 cells x 6 chromosomes and 844 MB at
1,000 cells x 22 chromosomes, mostly
`reconstructions` (per-node `internal_nodes` and `merged_profiles` for every
chromosome and allele) and `all_input_Xs`.

- [ ] Expose `reconstruct_tree()`'s existing `store_profiles` argument in
      `fit()`, or store only what `detect_bfb_branches()`,
      `plot_heatmap()` and `plot_bfb_signature()` need.

### Not profiled yet

- [ ] `bridge_sim()` post-processing: `sequences_to_cndata()` loops over
      cells and alleles in R, and `build_phylo_from_lineage()` builds the
      lineage tree in R.

## WGD (`fit(wgd = TRUE)`)

Known limits, from the prototype studies in `dev/wgd/README.md` (branch
`feature/wgd_develop`):

- [ ] Replicating cells or doublets have high ploidy and get a WGD count of
      1, which creates tiny spurious WGD clades (6-13 per tree on OV2295).
      Options: minimum WGD clade size, requiring the WGD explanation to win
      on that branch, or filtering such cells upstream.
- [ ] At most two WGDs per lineage; events between two WGDs are folded into
      the before/after comparison.
- [ ] Whether two doubled clades share one WGD or doubled separately is
      decided by the tree only; with near-identical sister lineages the
      data cannot tell them apart.
- [ ] Speed: WGD costs are computed only between cells whose WGD counts
      differ (`wgd_all_pairs = FALSE`). `fit(wgd = TRUE)` then costs
      1.0-2x a fit without WGD (all pairs: 1.0-14x); see
      `Rscript dev/benchmark_wgd.R`. Next step if needed: parallel pair loop.
- [ ] On noisy real data a WGD explanation sometimes wins between cells
      with the same estimated count (OV2295: 375 of 396,495 pairs at 2 Mb);
      the restricted default ignores these. Large WGD clades are unchanged,
      but at 500 kb the fine tree structure is sensitive to such tiny
      distance changes (13% of splits differ, only tiny WGD clades change).
- [ ] Noise degrades trees much more than WGD does (see the performance
      section): smoothing or segment-level input should come first on real
      data.
- [ ] The BFB distance is not used inside the WGD explanation (BFB after a
      WGD is counted as gain/loss events).
