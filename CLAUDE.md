# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`bridges` is an R package for inferring phylogenies from single-cell copy number variation (CNV) data, with robustness to BFB (Breakage-Fusion-Bridge) mutational processes common in cancer genomics.

## Development Commands

```r
devtools::load_all()    # Load package for interactive development
devtools::document()    # Regenerate Roxygen docs (NAMESPACE + man/)
devtools::check()       # Full R CMD CHECK
devtools::build()       # Build .tar.gz archive
devtools::install()     # Install locally
```

There are currently no unit tests. The package uses `testthat >= 3.0.0` in Suggests but the `tests/` directory does not exist yet.

The vignette serves as the primary end-to-end integration test:
```r
devtools::build_vignettes()
```

## Architecture

The package has three functional layers:

### 1. Simulation (`bridge_sim.R`, `utils_bridge_sim.R`)
`bridge_sim()` runs a Gillespie continuous-time simulation of cell evolution, producing synthetic copy number profiles and a ground-truth phylogeny. Supports BFB cycles, amplifications, deletions, and whole-genome duplication.

### 2. Inference (`fit.R`, `distances.R`, `bfb_mapping_and_detection.R`)
`fit()` is the main pipeline:
- Preprocesses CN data via `preprocess.R` (jitter correction, matrix conversion, diploid reference insertion)
- Computes **Greedy (G) distances** — minimum edit cost to transform one CN profile to another
- Computes **BFB/contiguity (B) distances** — penalizes non-contiguous copy number changes (BFB signature)
- Merges distances across chromosomes/alleles via `find_minimal_distances()`
- Builds a phylogenetic tree using `ape::nj()` (neighbor-joining), rooted at a synthetic diploid cell
- Reconstructs ancestral CN states per branch via `compute_reconstructions()`

`detect_bfb()` performs binomial/t-tests on the delta values (difference between G and B costs per segment) to identify genomic regions with BFB signatures.

### 3. Visualization (`plot.R`, `plot_bfb_signature.R`, `plot_chr_all.R`)
- `plot_heatmap()` — primary output: ComplexHeatmap with embedded phylogenetic tree, chromosome annotations, and CN color scales
- `plot_bfb_signature()` — bar plot of BFB event frequency per genomic segment
- `plot_chr_all_heatmap()` — detailed per-chromosome view

## Input Data Format

All main functions expect a data frame with columns:
```
cell_id | chr | start | end | CN | A | B
```
where `A` and `B` are allele-specific copy numbers (haploid), and `CN = A + B`.

## Key Design Decisions

- Distance functions are passed as arguments to `fit()` (e.g., `G_with_steps`, `A_contig`, `B_contig`), making the distance metric pluggable.
- The diploid reference cell is synthetically added as `"diploid"` before tree construction and removed/handled as root afterward.
- Reconstructions are stored as nested lists indexed by `[[chromosome]][[allele]]`.
- The `delta` values (stored in reconstructions) encode the BFB signal: high delta on a branch/segment means BFB-like events explain the CN change better than simple edits.
