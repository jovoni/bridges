# bridges 0.2.0

## Breaking changes

* `split_chr_arms()` and `bfb_detect_batch()` no longer accept a manual
  `centromere` argument, and no longer fall back to guessing a centromere
  position from the midpoint of observed data. Arm-awareness is now enabled
  via a required `genome` argument (`"hg19"` or `"hg38"`), which looks up
  real centromere positions from the newly bundled `centromeres` dataset.
* `bridge_sim()`'s chromosome-length table (used to size simulated
  chromosomes into bins) was previously a set of unlabeled, approximate
  values that turned out to be hg18/NCBI36 lengths. It has been replaced
  with real lengths from the same `centromeres` dataset, selectable via a
  new `genome` argument (default `"hg19"`). **This changes `bridge_sim()`'s
  simulated output bit-for-bit for the same seed** (bin counts per
  chromosome differ from the old, mislabeled values). There is no
  backwards-compatibility option to reproduce the old numbers.

## New features

* `fit()` can model whole-genome duplication (WGD) with `wgd = TRUE`. The
  distance between two cells becomes the cheaper of "no WGD" (the usual
  G/B distance) and "one cell is the other's lineage doubled once or twice,
  plus gain/loss events", with one `wgd_penalty` per WGD for the whole
  genome. WGD counts per cell (estimated from ploidy, or given with
  `wgd_counts`) are placed on the tree by parsimony; WGD branches are
  returned in `fit$wgd$edges` and as `"wgd"` events in `fit$history`, and
  their nodes are excluded from the BFB test of `detect_bfb_branches()`.
  WGD costs are only computed between cells whose WGD counts differ, which
  keeps `fit(wgd = TRUE)` within about 1-2x the time of a fit without WGD;
  `wgd_all_pairs = TRUE` computes every pair.
  The default (`wgd = FALSE`) gives exactly the same results as before.

* `fit()` is much faster with the default distances (`"greedy_fast"`,
  `"bfb_fast"`): all pairwise distances of a chromosome and allele are
  computed in one C++ call, once per pair of distinct copy-number profiles.
  On 1,000 cells x 22 chromosomes the distance step drops from about 6
  minutes to about 1 second; results are unchanged.
* `fit()` gains `keep_matrices` (default `FALSE`). The per-chromosome
  `greedy_Ds` and `avg_Ds` matrices are now only returned when it is `TRUE`.

* New bundled dataset `centromeres`: real centromere positions, arm
  boundaries, and chromosome lengths for chromosomes 1-22, X, Y under both
  hg19 (GRCh37) and hg38 (GRCh38), sourced from UCSC Genome Browser tracks.
  See `?centromeres` and `data-raw/centromeres.R`.

## Bug fixes

* `breakpoint_support = "custom"` now works throughout the simulation (it
  previously errored as soon as a BFB event occurred after the founding
  one), and the cut now falls right after one of `custom_breakpoints`
  (previously `sample()` on a length-1 vector drew a random value in
  `1:x`, and a bin value was used as a position).
* `bridge_sim(initial_sequences = ...)` with `return_phylo = TRUE` now
  returns a valid tree (founder cells were missing from the lineage).
* `bridge_sim()` no longer reports an `elapsed_time` beyond `max_time`.
* `simulate_serial_passages(n_cores > 1)` no longer errors, and
  `clonal_params()` gains `saturation_K`, which is passed through.
* `fit()` no longer warns about infinite distances on every call (the
  diagonal was Inf by construction), rejects misspelled `g_dist_func` /
  `b_dist_func` names up front, and gains `verbose`.
* Internal `tibble_to_matrix()` no longer adds columns to a caller's
  `data.table` by reference.
* `bfb_detect_batch(genome = ...)` no longer errors when a chromosome arm
  has no bins (e.g. acrocentric p-arms).
* `plot_chr_all_heatmap()`: the average-CN and gain/loss annotations are
  now shown together, restricted to the plotted chromosome, and use a
  baseline of 2 for `"CN"`.
* `plot_bfb_signature()` gives a clear error when no node is BFB-like.
* `bridge_sim()` now validates `lambda`, `rate`, `wgd_*`, `subsample` and
  `custom_breakpoints`, and warns when `hotspot$pos` lies beyond its
  chromosome. Documented defaults now match the actual defaults.
