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
  - without WGD: 61 BFB-like nodes, 26 of which join cells with different
    WGD counts; chromosome 6 is significant (adj. p = 0.005);
  - with WGD: 33 BFB-like nodes (10 joining different WGD counts); nothing
    is significant.

  The chromosome 6 signal was largely made of nodes at junctions between
  doubled and non-doubled cells, which the WGD option excludes from the
  test. That is a likely artefact, but it cannot be confirmed without
  ground truth.

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
  OV2295 the WGD option removes BFB-like nodes at WGD junctions, which
  changes the result (chromosome 6).
- Per-cell detection needs to know about WGD. Raw calls on doubled cells
  are unreliable; halving fixes the false positives but costs sensitivity.
  A WGD-aware version is needed (see the options in the report).
