# Tree-aware BFB test (`detect_bfb_clades()`)

Scripts (installed package; OV2295 data in `dev/data`, not in git):

```sh
Rscript dev/bfb_clades/clade_test_sim.R out.rds   # calibration and power, simulations
Rscript dev/bfb_clades/clade_test_ov2295.R        # OV2295, 500 kb, k_jitter_fix = 2
```

## Method

**Candidates.** Every clade of the inferred tree with at least `min_nodes`
(5) tested internal nodes. No clones or clades are given.

**Clade score.** k BFB-like nodes among the clade's n tested internal
nodes:

    p_c = P(Bin(n, threshold) >= k)

The best clade has the smallest p_c. The whole tree is one candidate, so
the genome-wide test of `detect_bfb_branches()` is included.

**Correction.** The smallest p_c is compared with 999 simulations in which
every tested internal node is BFB-like independently with probability
`threshold`, on the same tree. p-values are BH-adjusted across
chromosome/alleles. WGD nodes are not tested.

## Simulations

Setup: 8 chromosomes, 150 cells, BFB only on 1:A, 10 replicates per
scenario. In the "subpopulation" scenario:
- 8 founders are grown without BFB;
- 2 founders then evolve with BFB (30 cells sampled), 6 without (120
  cells);
- the two populations are sequenced together.

| Scenario | 1:A detected, genome-wide test | 1:A detected, clade test | False calls (both) | Best clade vs BFB cells (Jaccard) |
|---|---|---|---|---|
| no BFB | – | – | 0 | – |
| no BFB, with WGD (`wgd = TRUE`) | – | – | 0 | – |
| BFB in the whole population | 10/10 | 10/10 | 0 | 0.87 |
| BFB in a subpopulation (20% of cells) | 0/10 | 2/10 | 0 | 0.42 |

**Calibration:** no false calls in any scenario.

**Power with subpopulation BFB** is limited by the signal, not the test.
- The reconstruction finds only 0–4 BFB-like nodes on 1:A per replicate,
  and all of them lie inside the BFB subpopulation.
- The uncorrected scan p is below 0.05 in 7/10 replicates, but after BH
  across 16 chromosome/alleles most are just above 0.05.
- In replicate 3 the best clade has 3 BFB-like nodes and covers exactly
  the 30 BFB cells.

**Speed:** about 0.1 s for 150 cells and 16 chromosome/alleles.

## OV2295 (500 kb, jitter 2)

| | Genome-wide test | Clade test |
|---|---|---|
| `wgd = FALSE` | chr 6 | chr 6 (best clade G/H/I, 12 of 227 nodes), chr 10 and chr 3 (both the same 8-cell near-diploid clade, 3 of 7 nodes) |
| `wgd = TRUE` | none | **chr 6** (adj. p = 0.044; best clade = the G/H/I WGD clade, 244 cells, 8 of 242 nodes BFB-like) |

With WGD the clade test recovers chr 6 inside G/H/I, consistent with the
heterogeneous chr 6 gains after the WGD (`dev/bfb_wgd/README.md`).

The chr 10 / chr 3 calls without WGD come from a small near-diploid clade
that likely sits where E/F branches off, among the near-diploid cells next
to the E/F WGD (cf. the junction artefacts in `dev/bfb_wgd`). They
disappear with WGD.

## Notes

- `min_nodes` controls the smallest clade considered. Small clades with 2–3
  BFB-like nodes can reach significance (as for chr 10/3 above), so raising
  `min_nodes` makes the test more conservative.
- With `n_sim = 999` the smallest possible `p.scan` is 0.001. More
  simulations increase power after multiple-testing correction, at
  proportional cost.
