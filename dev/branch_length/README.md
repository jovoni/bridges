# Root-to-tip length vs true number of events

In simulations every division records its events (`cell_history$cn_event`:
bfb, wgd, amp, del), so each sampled cell's true event count is the sum along
its lineage. This compares that count with the root-to-tip length of
`fit()$tree` (reconstructed events per edge, plus `wgd_penalty` per WGD).

- `events_vs_depth.R`: 4 scenarios (CNA only, BFB, WGD, BFB + WGD) ×
  lambda 0.5 / 1 / 2 × 5 replicates, 10 chromosomes, 150 cells, A + B;
  each fitted with `wgd = FALSE` and `wgd = TRUE` (~1.5 min).
  `REC_TARGET=doubled` reruns it with the greedy-ancestor target doubled
  inside doubled clades.
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
| WGD | `wgd = TRUE` (current) | 0.35 | 4.44 | 80 |
| WGD | `wgd = TRUE`, doubled target | **0.92** | 1.24 | 3.1 |
| BFB + WGD | `wgd = FALSE` | 0.60 | 4.61 | 88 |
| BFB + WGD | `wgd = TRUE` (current) | 0.62 | 4.06 | 71 |
| BFB + WGD | `wgd = TRUE`, doubled target | **0.84** | 1.48 | 9.6 |

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

Doubling the target (`target_val * 2^k`) removes this. It was removed in
96b6774 because, on OV2295 chr 6, the doubled target picked poor ancestors
and hid BFB signal. So the target choice trades branch-length accuracy for
BFB sensitivity.
