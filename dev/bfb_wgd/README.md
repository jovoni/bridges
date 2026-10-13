# How does WGD affect BFB detection?

Scripts (installed package; OV2295 data in `dev/data`, not in git):

```sh
Rscript dev/bfb_wgd/bfb_vs_wgd_sim.R out.rds   # simulations, BFB only on 1:A
Rscript dev/bfb_wgd/bfb_vs_wgd_ov2295.R        # OV2295, 500 kb, k_jitter_fix = 2
```

Simulations: 10 chromosomes, 150 cells, 5 replicates per scenario (BFB +
WGD; WGD only, no BFB at all; BFB only). A cell counts as doubled from the
simulation's true WGD history.

## Tree-based detection (`detect_bfb_branches()`)

Results:

| Scenario | `wgd` | True BFB (1:A) called | False calls (of 19) |
|---|---|---|---|
| BFB + WGD | FALSE / TRUE | 5/5 / 5/5 | 0 / 0 |
| WGD only | FALSE / TRUE | – | 0 / 0 |
| BFB only | FALSE / TRUE | 5/5 / 5/5 | 0 / 0 |

- **Simulations:** WGD does not affect it. The true BFB is always found, no
  false calls appear, and WGD alone never produces a BFB call.
- **OV2295:**
  - without WGD: 61 BFB-like nodes over all chromosomes, 26 of which join
    cells with different WGD counts; chromosome 6 is significant
    (14 of 890 nodes, adj. p = 0.005);
  - with WGD: 33 BFB-like nodes (10 joining different WGD counts); nothing
    is significant (chromosome 6: 5 of 876 nodes).

**Correction: chromosome 6 is not a junction artefact.** An earlier
version of this note said the chromosome 6 signal came from nodes at
junctions between doubled and non-doubled cells. That was generalised from
the all-chromosome count and is wrong for chromosome 6.
- 12 of its 14 BFB-like nodes are inside the G/H/I clade, where all cells
  have the same WGD count (1). Several are large internal nodes (37, 79,
  83, 88 and 228 cells). Only 2 are in A–D.
- There, chromosome 6 shows a heterogeneous high-level gain (CN 6–9),
  which looks BFB-like on the heatmap.

The WGD option removes these nodes through two changes, separated by
swapping trees and reconstructions:

| Tree from | Reconstruction | BFB-like nodes on chr 6 | Of which ≥ 10 cells |
|---|---|---|---|
| `wgd = FALSE` | plain | 14 | 6 |
| `wgd = FALSE` | WGD-aware | 7 | 1 |
| `wgd = TRUE` | plain | 11 | 2 |
| `wgd = TRUE` | WGD-aware (`fit()` output) | 5 | 0 |

The WGD-aware reconstruction rebuilds ancestors inside doubled clades
around a baseline of 4 instead of 2, which changes the profiles compared
at deeper nodes. The different tree accounts for most of the rest. In the
simulations neither effect removed the true BFB signal, so on OV2295 it is
unknown whether the WGD option removes an artefact or hides real BFB on
chromosome 6.

Figures (`dev/data/figures/`, not in git), with BFB-like nodes on chr 6 in
red: `ov2295_500kb_jitter2_no_wgd_bfb_chr6.png`,
`ov2295_500kb_jitter2_wgd_bfb_chr6.png`.

## Per-cell detection (`bfb_detect_batch()`)

Fraction of informative count-vectors that admit BFB:

| Scenario | Chromosome | Cells | Raw | Doubled cells halved first |
|---|---|---|---|---|
| BFB + WGD | no BFB | not doubled | 0.012 | 0.012 |
| BFB + WGD | no BFB | doubled | **0.326** | 0.000 |
| BFB + WGD | BFB (1:A) | not doubled | 0.662 | 0.662 |
| BFB + WGD | BFB (1:A) | doubled | 0.545 | 0.266 |
| WGD only | no BFB | not doubled | 0.039 | 0.039 |
| WGD only | no BFB | doubled | **0.257** | 0.067 |
| BFB only | no BFB | – | 0.017 | 0.017 |
| BFB only | BFB (1:A) | – | 0.792 | 0.792 |

On OV2295 (total CN), the same pattern:

| Cells | Raw | Halved |
|---|---|---|
| Doubled | 0.410 | 0.130 |
| Near-diploid | 0.123 | 0.123 |

- **Per-cell detection is strongly biased by WGD.** Doubled cells admit BFB
  on 26–33% of chromosomes that never had BFB, against 1–4% for other
  cells. A doubled profile has a baseline of 2 and many even counts, so it
  passes the `min_max_count = 2` filter and often fits a BFB pattern.
- **Halving doubled cells first removes the false positives,** but it also
  halves sensitivity to real BFB in doubled cells (0.545 → 0.266). Rounding
  down loses the odd counts left by BFB cycles that happened after the WGD.

## Conclusions

- Tree-based detection is safe with or without WGD in simulations. On
  OV2295 the WGD option removes the chromosome 6 signal inside the doubled
  G/H/I clade, through both the tree and the doubled-baseline
  reconstruction; whether that signal is real BFB is open.
- Per-cell detection needs to know about WGD. Raw calls on doubled cells
  are unreliable; halving fixes the false positives but costs sensitivity.
  A WGD-aware version is needed (see the options in the report).
