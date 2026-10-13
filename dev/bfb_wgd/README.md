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

## Where are the WGDs, and what happened before and after? (`wgd_placement_ov2295.R`)

On `fit(wgd = TRUE)` (500 kb, jitter 2) there are 14 WGD branches.

**Two real WGD branches:**

| WGD branch | Sister clade | Note |
|---|---|---|
| G+H+I (243 cells) | a single G cell estimated as not doubled | The WGD sits at the base of G/H/I; the pre-WGD ancestor is constrained mainly by that one cell, so it is weakly determined. |
| E+F (223 cells; 1 A and 2 B cells included) | 10 near-diploid cells (A–D) | E/F branches off inside the near-diploid group. |

**Twelve 1–4-cell WGD branches**, mostly single A–D cells with high
ploidy: likely replicating cells or doublets, not real WGDs.

**Gain/loss events on each part of the lineage** (greedy events summed over
chromosomes):

| WGD | Before the WGD (diploid → pre-WGD ancestor) | On the WGD branch, after doubling | Within the clade |
|---|---|---|---|
| G/H/I | 40 | 46 | 47 |
| E/F | 95 | 124 | 74 |

Each cell–to–ancestor comparison here carries the single-cell noise (about
29 events per cell), so these totals are upper bounds.

**Chromosome 6 in G/H/I:**

| Region (Mb) | Pre-WGD ancestor (doubled) | Post-WGD ancestor | Clade consensus | Cells with CN ≥ 6 |
|---|---|---|---|---|
| 16–22.5 | 3 (6) | 6 | 7 | 98% |
| 0–16 | 3 (6) | 4 | 5 | 46–48% |
| 37.5–62.5 | 3 / 2 (6 / 4) | 4 | 6 | 52–56% |
| 110–127 | 2 (4) | 6 | 5 | 7% |
| rest of 6q | 1 (2) | 3–4 | 2–3 | 0% |

- The 16–22.5 Mb gain is clonal. It is already present (doubled) in the
  post-WGD ancestor and in 98% of cells, so it arose before G/H/I
  diversified.
- The gains at 0–16 and 37.5–62.5 Mb are heterogeneous: about half the
  cells reach CN ≥ 6, above the post-WGD ancestor (4). They arose **after
  the WGD, within G/H/I**, exactly where the BFB-like nodes of the
  `wgd = FALSE` tree are.

**Implication.** The chromosome 6 data are consistent with BFB-like
activity after the WGD. The WGD-aware reconstruction removing those nodes
therefore looks more like lost sensitivity inside doubled clades than an
artefact being removed. This needs checking with simulations of BFB after
a WGD.

Outputs (`dev/data/figures/`):
- `ov2295_wgd_ancestral_profiles.png`: genome-wide profiles around each
  WGD;
- `ov2295_wgd_ancestral_profiles.rds`: profile table.
