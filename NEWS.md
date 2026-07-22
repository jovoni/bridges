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

* New bundled dataset `centromeres`: real centromere positions, arm
  boundaries, and chromosome lengths for chromosomes 1-22, X, Y under both
  hg19 (GRCh37) and hg38 (GRCh38), sourced from UCSC Genome Browser tracks.
  See `?centromeres` and `data-raw/centromeres.R`.
