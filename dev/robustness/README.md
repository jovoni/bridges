# Robustness on OV2295 (500 kb): why are cells placed in the wrong clade?

Scripts (run on the installed package; data in `dev/data`, not in git):

```sh
Rscript dev/robustness/ov2295_robustness.R out.rds   # jitter, support filter, tree methods
Rscript dev/robustness/weighted_events.R             # length-weighted event costs
```

## Diagnosis (fit(wgd = TRUE), jitter 0)

**Metric.** A cell is "misplaced" when most cells of its local clade (the
smallest clade with at least 10 other cells) belong to another clone.

**Result.** 279 of 891 cells (31%) are misplaced. Almost all go to a
**sister clone** within the same group:

| Own clone → majority of its clade | Cells |
|---|---|
| H → G | 45 |
| B → A | 44 |
| C → D | 43 |
| D → C | 40 |
| G → H | 21 |

Cells almost never cross between the near-diploid, E/F and G/H/I groups.

**What the misplaced cells are not:** noisier than the others (median 84
breakpoints vs 83), or cells with an unusual WGD count (12 of 279).

**Why it happens: weak signal against noise.**

| Comparison | Greedy events |
|---|---|
| Consensus of sister clones C vs D | 11 |
| Consensus of sister clones G vs H | 11 |
| Consensus of sister clones A vs B | 6 |
| A single cell vs its own clone's consensus (noise + private changes) | about 29 |

**Distances or tree?** About half of the misplaced cells (140 of 279) have
mostly their own clone among their 10 nearest cells by distance, so part
of the error comes from tree building.

## What was tried

Metrics:
- misplaced %;
- sister-clone mix-ups (C/D, A/B, G/H);
- share of each cell's 10 nearest cells that belong to its clone;
- clone coherence (best F1 per clone, weighted).

The two large WGD clades (G+H+I, E+F) are the same in every configuration.

| Configuration | Misplaced | C/D | A/B | G/H | Nearest cells same clone | Clone coherence |
|---|---|---|---|---|---|---|
| jitter 0 (default) | 31.3% | 83 | 51 | 66 | 0.72 | 0.734 |
| jitter 1 | 27.3% | 74 | 44 | 55 | 0.77 | 0.741 |
| jitter 2 | 26.9% | 70 | 16 | 59 | 0.78 | 0.759 |
| **jitter 3** | **23.6%** | **63** | 22 | **56** | 0.79 | **0.778** |
| jitter 4 | 26.0% | 67 | 25 | 61 | 0.80 | 0.746 |
| support filter (≥ 3/5/10 cells) | 27.7–32.1% | 70–106 | 29–31 | 62–77 | 0.72 | 0.699–0.746 |
| support filter + jitter 1 | 26.4–28.6% | 76–82 | 35–39 | 55–61 | 0.77 | 0.728–0.755 |
| jitter 2 + BIONJ | 29.5% | 68 | 49 | 63 | 0.78 | 0.742 |
| jitter 2 + FastME (balanced) | 25.9% | 75 | 28 | 57 | 0.78 | 0.780 |
| jitter 2 + FastME (OLS) | 22.0% | 47 | 24 | 72 | 0.79 | 0.759 |
| length-weighted events, L = 20 bins | 25.4% | 52 | 5 | 47 | 0.81 | 0.691 |
| reference: plain bin-by-bin L1 distance + NJ | 19.2% | 25 | 26 | 27 | 0.90 | 0.788 |

Length-weighted events (`weighted_events.cpp`): each greedy/BFB event costs
min(1, length / L), so short noise events weigh less.

## The clone labels are not ground truth

The bin-by-bin L1 distance matches the clone labels best. The labels are
most likely produced by clustering bin-level copy-number profiles, which
favours exactly that kind of distance. On simulations with a known tree
(6 replicates, one WGD, 22 chromosomes, mean RF to the true tree; lower is
better):

| Method | No noise | Low noise |
|---|---|---|
| `fit(wgd = TRUE)`, jitter 0 | **0.302** | 0.568 |
| `fit(wgd = TRUE)`, jitter 2 | 0.400 | **0.413** |
| length-weighted events, L = 20 | 0.410 | 0.497 |
| L1 bin distance + NJ | 0.595 | 0.615 |

So the event-based distance builds much better trees than the bin-by-bin
distance; its lower agreement with OV2295's clone labels is partly a
property of how those labels were made.

## Conclusions

- Most "wrong spots" are swaps between sister clones whose consensus
  profiles differ by only 6–11 events, below the noise of a single cell
  (about 29 events). Some are expected, since the labels come from a
  different (bin-level) method.
- **`k_jitter_fix` is the effective lever on noisy data:**
  - OV2295: jitter 3 reduces misplaced cells from 31% to 24% and raises
    clone coherence from 0.734 to 0.778;
  - noisy simulations: jitter 2 reduces tree error from 0.568 to 0.413.
- **But jitter hurts clean data** (simulations without noise: 0.302 →
  0.400). It should be recommended for real data, not made the default
  for everything.
- **Did not help consistently:** the support filter, BIONJ, FastME, and
  length-weighted event costs. FastME (OLS) reduced C/D swaps (47) but
  increased G/H swaps (72) and produced more tiny WGD edges (19).
