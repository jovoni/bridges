
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bridges

<!-- badges: start -->

[![R-CMD-check](https://github.com/jovoni/bridges/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jovoni/bridges/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Simulate, fit, and detect breakage-fusion-bridge (BFB) signal in
single-cell copy-number data. Combines a Gillespie simulation engine
(`bridge_sim()`), a phylogenetic inference pipeline (`fit()`), two
independent BFB detection methods, and visualization tools.

**BFB in brief.** BFB cycles begin when a chromosome loses its telomere
protection. During replication the two sister chromatids fuse at their
broken ends, forming a dicentric chromosome that is pulled to both poles
during mitosis, causing a new break. Repeated cycles amplify the genomic
region near the original break, producing a characteristic “fold-back”
copy number signature at one allele.

## Installation

`bridges` depends on [`bfbtools`](https://github.com/jovoni/bfbtools)
for the underlying linear-time BFB count-vector algorithms. Because
neither package is on CRAN, installation goes through `remotes`/`pak`,
and `bfbtools` is declared in `bridges`’s `DESCRIPTION` under `Remotes:`
so it resolves automatically:

``` r
# install.packages("remotes")

# ggtree and ComplexHeatmap are Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ggtree", "ComplexHeatmap"))

remotes::install_github("jovoni/bridges")   # pulls in bfbtools automatically
```

If you’d rather install `bfbtools` explicitly first (makes failures
easier to diagnose), that works too and is equivalent:

``` r
remotes::install_github("jovoni/bfbtools")
remotes::install_github("jovoni/bridges")
```

## Usage

``` r
library(bridges)

sim <- bridge_sim(chromosomes = c(1:8), bfb_allele = "8:A", max_cells = 5000, lambda = 2)
res <- fit(data = sim$cna_data, alleles = c("A", "B"))
detect_bfb_branches(res, threshold = 0.005)
```

See `vignette("bridges")` for the full end-to-end walkthrough, and
`vignette("simulator")` / `vignette("inference")` /
`vignette("batch-detection")` / `vignette("passages-and-clones")` for
deep dives into each part of the pipeline.

## Two independent BFB detection methods

- `detect_bfb_branches()` – works on a *fitted tree* (`fit()`’s output):
  for each branch, asks whether the copy-number transition looks
  BFB-like (non-contiguous) via a binomial test.
- `bfb_detect_batch()` / `bfb_detect_vectors()` – work directly on *raw
  per-cell profiles*, no tree required, using the exact BFB
  combinatorics in `bfbtools`. Includes chromosome-arm-aware handling
  (`split_chr_arms()`) for real genomic data spanning both arms of a
  chromosome.

Neither supersedes the other – see `vignette("batch-detection")` for
when to use which.

## Part of a two-package family

- **`bridges`** (this package) – simulation, fitting, detection glue,
  and plotting. MIT licensed.
- [**`bfbtools`**](https://github.com/jovoni/bfbtools) – the core
  linear-time BFB detection algorithms. GPL-3 licensed; kept as a
  separate package for exactly that reason (see below).

## Why is `bfbtools` a separate dependency instead of being bundled in?

`bfbtools`’s core algorithm is a port of GPL-3 licensed reference code.
GPL is viral for *merged source*: folding that code into this package
would force all of `bridges` (including code that owes nothing to it) to
become GPL-3 too. A normal package *dependency* – what `Remotes:` and
`Imports:` set up here – does not have that effect; only merging source
does. This is why the two packages exist separately even though they’re
developed and released together.

## License

MIT. Depends on `bfbtools` (GPL-3) – a normal cross-package dependency
that does not change this package’s own license.
