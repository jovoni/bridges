# Root-to-tip length vs true number of events

In simulations every division records its events (`cell_history$cn_event`:
bfb, wgd, amp, del), so each sampled cell's true event count is the sum along
its lineage. This compares that count with the root-to-tip length of
`fit()$tree` (reconstructed events per edge, plus `wgd_penalty` per WGD).

- `events_vs_depth.R`: 4 scenarios (CNA only, BFB, WGD, BFB + WGD) ×
  lambda 0.5 / 1 / 2 × 5 replicates, 10 chromosomes, 150 cells, A + B;
  each fitted with `wgd = FALSE` and `wgd = TRUE` (~1.5 min).
- `ov2295_bfb_check.R`: OV2295 (500 kb, `wgd = TRUE`, jitter 0 and 2):
  BFB-like nodes, both BFB tests, WGD edge lengths. Used to compare
  reconstructions before and after a change.
- `analyse.R`: correlations, scale (mean depth / mean truth), per-type
  coefficients of `depth ~ amp + del + bfb + wgd`, and a scatter figure.

## Results (lambda = 1, mean over 5 replicates)

r = Pearson correlation between true events and root-to-tip length;
ratio = mean depth / mean true events; "per WGD" = depth added per true WGD.

| scenario | fit | r | ratio | per WGD |
|---|---|---|---|---|
| CNA only | either | 0.89 | 1.06 | - |
| BFB | either | 0.85 | 1.20 | - |
| WGD | `wgd = FALSE` | 0.31 | 5.15 | 99 |
| WGD | `wgd = TRUE`, plain target (before 8b06268) | 0.35 | 4.44 | 80 |
| WGD | `wgd = TRUE`, doubled target (now) | **0.92** | 1.24 | 3.1 |
| BFB + WGD | `wgd = FALSE` | 0.60 | 4.61 | 88 |
| BFB + WGD | `wgd = TRUE`, plain target (before 8b06268) | 0.62 | 4.06 | 71 |
| BFB + WGD | `wgd = TRUE`, doubled target (now) | **0.84** | 1.48 | 9.6 |

Lambda 0.5 and 2 give the same picture. Without WGD, depth tracks the truth
closely (slope about 1 per amp, del and BFB cycle; deletions slightly
undercounted, as later events can hide them).

The distance itself is fine: with `wgd = TRUE`, the correlation of
`D[cell, "diploid"]` with the truth is 0.97-0.99 in the WGD scenarios.
The problem is in the reconstruction's edge lengths.

## Cause

Since commit 96b6774, the greedy ancestor of two same-level children is
chosen closest to the diploid target (1 per allele) even inside doubled
clades. When two doubled sister lineages carry different post-WGD deletions,
the point on the greedy path closest to 1 holds *both* sets of deletions.
Those deletions are charged again on the WGD edge, and each child must
"re-gain" the deletions it never had. This repeats up the doubled clade, so
each WGD edge costs ~50 instead of 1 + a few events, and edges inside the
clade are inflated (mean 4.3 vs 1.5 outside, one WGD replicate).

## Fix (commit 8b06268)

`reconstruct_tree()` again chooses ancestors closest to `target_val * 2^k`
inside clades with k WGDs, and the BFB deltas come from that same
reconstruction. The "doubled target (now)" rows above are the package
default since then.

Effect on BFB detection (see also `dev/bfb_wgd/README.md`):
- Simulations of BFB after a WGD: unchanged (6/6, 5/6, no false calls).
- OV2295, jitter 0: chr 6 and chr 22 stay significant in both tests.
- OV2295, jitter 2: chr 6 drops from 8 to 5 BFB-like nodes and is no longer
  called by `detect_bfb_clades()` (adj. p 0.044 → 1).

OV2295 WGD branch lengths (jitter 0):
- the two real WGDs: G/H/I 68 → 71, E/F 112 → 95;
- the other 11 WGD branches, 1-6 near-diploid cells each with mean CN
  2.9-6.5, have lengths 20-363. These cells are called doubled by ploidy
  but don't look like a doubling of their neighbours, so the WGD branch has
  to explain them with many extra events. That is a WGD-count estimation
  problem, not a reconstruction one.

## Remaining limitation (not WGD-specific)

The ancestor of two children is always a state on the greedy path from one
to the other. Two lineages with different losses give a path that passes
through "both losses", never "neither", so the ancestor ends up equal to one
child. The cost of that merge is still correct, but one loss is charged on
the wrong edge and moves up the tree. Picking the ancestor bin by bin
between the two children (closest to the target) would avoid this, but
would need checking that it never increases the cost.
