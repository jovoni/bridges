# WGD-aware distance: prototype

Prototype only; nothing here is used by the package. Run from the package
root:

```sh
Rscript dev/wgd/check_kernel.R                        # kernel checks
Rscript dev/wgd/run_experiments.R results.rds         # first simulation study (~1 min)
Rscript dev/wgd/run_experiments_hard.R results.rds    # harder study (~40 min of compute)
Rscript dev/wgd/real_data_analysis.R out.rds figdir   # OV2295 sample in dev/data (~3 min)
```

## Idea

For cells a and b, with S0 the current `fit()` distance (sum over
chromosomes/alleles of min(G, B)):

    D(a,b) = min( S0,  w + S1(a→b),  w + S1(b→a),  2w + S2(a→b),  2w + S2(b→a) )
    Sk(a→b) = Σ_chr,allele min_x [ G(a, x) + G(2^k · x, b) ]

- `w` is the cost of one WGD.
- The WGD decision is genome-wide, while the ancestor x is chosen per
  chromosome/allele.
- Candidates for x: every profile on the greedy path from a to floor(b/2^k)
  and to ceil(b/2^k), both ends included.
- Two variants:
  - `"all"` evaluates every hypothesis;
  - `"ploidy"` uses only the one matching the cells' estimated WGD counts
    (round(log2(mean CN / 2)), 0–2).

Files:
- `wgd_proto.cpp`: kernels `wgd_cost()` and `wgd_ancestor()` (the latter
  also returns the best ancestor), and the per-chromosome matrices.
- `wgd_proto.R`:
  - genome-wide components and their combination;
  - tree building;
  - true WGD ancestry of a simulation;
  - WGD reconstruction (`assign_wgd_counts()`, `reconstruct_wgd()`);
  - scoring (`score_wgd_labels()`, `clone_coherence()`).
- `real_data_prep.R`: loads the OV2295 sample from `dev/data` (not in git).

## WGD reconstruction

1. Each tip gets its estimated WGD count (0–2, from ploidy).
2. A parsimony pass assigns counts to internal nodes. Along an edge the
   count can only increase, at cost w per WGD; the root's ancestor is
   diploid. Edges where the count increases are the predicted WGD edges.
3. Ancestral profiles are rebuilt bottom-up, per chromosome/allele:
   - two children with the parent's count: existing G/B logic;
   - one child doubled: the `wgd_ancestor()` kernel;
   - both children doubled: compare their halved profiles.

   This gives event counts per edge.

Checks:
- Toy trees: one WGD on the right clade; parallel WGDs on two tips; a WGD
  on the stem when all cells are doubled.
- Without WGD, it gives exactly the same total events as the existing
  reconstruction.

## Checks (`check_kernel.R`)

- Known cases give the expected costs:
  - W(x→2x) = 0;
  - one focal change after WGD = 1;
  - W(x→4x) = 0 with k = 2.
- Exhaustive search over x on 400 short random cases: the approximation is
  never below the exact minimum, is exact in 97.5% of cases, and is at most
  1 event above it.
- With no WGD term, the prototype distance is identical to `fit()$D`.

## Simulation study (`run_experiments.R`)

**Setup**
- 10 chromosomes, 150 sampled cells, BFB on 1:A, 6 replicates per scenario.
- Scenarios:
  - no WGD;
  - one WGD;
  - two WGDs in different clades (`wgd_available = 2`);
  - nested WGDs, cells with 0, 1 and 2 WGDs.

  The simulator never produced two WGDs in one lineage, so the nested
  scenario uses two stages: 8 founders, some doubled, then
  `initial_sequences` with one more WGD. Its true tree is the second
  stage's.
- Both allele-specific (A+B) and total CN data.

**Results, allele-specific data** (mean of 6 replicates, w = 1; total CN is
similar):

| Scenario | RF none → WGD | WGD clades monophyletic, none → WGD | D across WGD classes, none → WGD | D within class | WGD term wins (within class) |
|---|---|---|---|---|---|
| No WGD | 0.304 → 0.304 | – | – | 16.0 | 0% |
| One WGD | 0.314 → 0.300 | 100% → 100% | 37.6 → 18.7 | 15.4 | 0% |
| Parallel WGDs | 0.310 → 0.291 | 75% → 83% | 37.9 → 19.1 | 15.6 | 0% |
| Nested WGDs | 0.334 → 0.321 | 100% → 100% | 39.3 → 18.9 | 15.2 | 0% |

**Findings**
- The WGD term removes the inflation of distances between cells with
  different WGD counts (about 19 events here, with 10 chromosomes). It is
  never chosen for pairs with the same WGD count, and changes nothing when
  there is no WGD.
- Tree topology improves only slightly: RF −0.01 to −0.02 on average, and
  not in every replicate. Neighbour-joining already groups WGD clades. The
  main gain is in distances and branch lengths, which matter for
  reconstruction and for placing and rooting WGD clades.
- w = 1 works best; w = 3 is similar; w = 10 starts to lose WGD pairs on
  total CN.
- `"all"` and `"ploidy"` gave identical results, because ploidy classes
  were 100% correct in these simulations.
- Cost: about 10–16× the current distance computation (0.42 s vs 0.03 s for
  400 cells × 10 chromosomes).

**Limitations**
- The simulations are easy: no measurement noise, perfect ploidy classes,
  10 chromosomes, 6 replicates.
- Nested WGDs come from a two-stage stand-in.
- These limitations are addressed by the harder study and the WGD reconstruction below.

## Harder simulation study (`run_experiments_hard.R`)

**Setup**
- 22 chromosomes, 150 cells, 10 replicates per scenario.
- WGD clades down to 5% of cells; parallel WGD clades of at least 4 cells.
- Noise added to each cell's profile, per chromosome and allele:

  | Level | Bins changed by ±1 | Chance of a 1–5-bin spurious segment |
  |---|---|---|
  | none | 0 | 0 |
  | low | 1% | 10% |
  | high | 3% | 30% |

- Allele-specific and total CN data.
- Compared: the current distance ("none") vs the WGD-aware distance
  ("all", w = 1).

**Results, allele-specific data** (total CN is similar)

| Scenario | Noise | Tree error (RF), none → WGD | WGD label match (Jaccard), none → WGD | WGD edges predicted (true), none → WGD | Events in reconstruction, standard → WGD-aware |
|---|---|---|---|---|---|
| One WGD | none | 0.30 → 0.30 | 1.00 → 1.00 | 1.0 → 1.0 (1) | 677 → 470 |
| One WGD | low | 0.51 → 0.50 | 1.00 → 1.00 | 1.0 → 1.0 (1) | |
| One WGD | high | 0.88 → 0.87 | 1.00 → 1.00 | 1.0 → 1.0 (1) | |
| Parallel WGDs | none | 0.30 → 0.28 | 0.70 → 0.95 | 1.4 → 1.9 (2.4) | 789 → 503 |
| Parallel WGDs | low | 0.52 → 0.49 | 0.70 → 0.90 | 1.4 → 1.8 (2.4) | |
| Parallel WGDs | high | 0.86 → 0.86 | 0.70 → 0.70 | 1.4 → 1.4 (2.4) | |
| Nested WGDs | none | 0.32 → 0.30 | 1.00 → 1.00 | 2.0 → 2.0 (2) | 1095 → 490 |
| No WGD | all | identical | – | 0 → 0 | identical |

Notes on the table:
- "Parallel" replicates sometimes also contain cells with 2 WGDs, hence 2.4
  true WGD sets on average.
- "Nested" true trees have the stage-1 founders as a star at the root, so
  RF is less meaningful there.

**Findings**
- **No false positives.** With no WGD, the WGD term never wins, no WGD edge
  is predicted, and trees are identical, at every noise level.
- **The inflation is removed.** With 22 chromosomes, the distance between
  cells with different WGD counts drops from 59 to 16.6 (within-class
  distance is 13.4).
- **Parallel WGDs are where the term matters.** Without it, the two WGD
  clades are often merged under one WGD (1.4 predicted WGD edges for 2.4
  true). With it, they are mostly kept apart (1.9). The benefit disappears
  at high noise.
- **WGD-aware reconstruction explains the data with 30–55% fewer events**
  on the same tree (677 → 470 for one WGD; 1095 → 490 for nested WGDs).
- **Noise, not WGD, is the dominant problem.** Even without any WGD, tree
  error rises from 0.29 to 0.50 (low noise) and 0.88 (high noise). The
  pipeline does no smoothing by default (`k_jitter_fix = 0`).
- The `"ploidy"` and `"all"` variants gave identical results again. Ploidy
  classes were 100% correct even at high noise.

## Real data: OV2295 (`real_data_analysis.R`)

**Data**
- 891 clone-labelled cells (clones A–I), autosomes only, total CN.
- 500 kb bins, and also merged to 2 Mb (median state).
- Estimated WGD count per clone:
  - A–D near-diploid: 410 cells with k = 0, 17 with k = 1;
  - E–I: k = 1, except 2 cells of E with k = 2 and 1 cell of G with k = 0.

**WGD story**
- **Without the WGD term**, all doubled clones (E, F, G, H, I) form one
  clade, so the reconstruction places a single WGD above them.
- **With the WGD term**, E+F and G+H+I are separate clades, each with its
  own WGD: two independent WGDs.
- The deciding margin is small. With the WGD term, G/H/I is slightly closer
  to the near-diploid clones (87.2) than to E/F (90.2). Without it, G/H/I is
  closer to E/F (90.2 vs 108.0).
- Under the WGD-aware reconstruction, the two-WGD tree explains the data
  with fewer events: 32,860 vs 33,526 at 2 Mb; 45,399 vs 45,640 at 500 kb.
  These totals are dominated by noise.
- This matches the parallel-WGD simulations, where the current distance
  merges independent WGD clades.

**Clone coherence** (best F1 between each clone and any clade, weighted by
clone size)
- About 0.75 in all trees (500 kb: 0.751 without the WGD term, 0.732 with
  it).
- No single clone is monophyletic: noise scatters cells.
- Groups of clones are recovered almost perfectly (F1 ≥ 0.98): near-diploid
  A–D, E+F, and G+H+I.

**Spurious small WGD clades.** The 17 near-diploid cells with ploidy ≥ 3
(probably replicating cells or doublets) create extra WGD edges on 1–4
cells: 6–13 per tree.

Figures are in `dev/data/figures/` (git-ignored, since they show the real
sample): heatmaps with tree and clone annotation at 2 Mb, without and with
the WGD term.

## Next steps suggested by these results

1. **Handle noise first.** It degrades trees far more than WGD does. Test
   `k_jitter_fix`, bin merging, or segment-level input on the simulations
   and on OV2295.
2. **Suppress spurious small WGD clades.** Options:
   - require a minimum clade size;
   - require profile evidence (W term wins) on the WGD edge;
   - treat replicating cells and doublets upstream.
3. **Integrate into `fit()`** as an option (default off), with the WGD
   reconstruction reporting `"wgd"` events in `history`. The current
   `"ploidy"` and `"all"` variants agree; keep `"all"` (no reliance on the
   ploidy estimate) unless runtime matters.

## Alternatives tested and dropped

Two cheaper alternatives were prototyped in commit `ade6b59` and then
removed (see that commit for code and full tables):

1. **Ploidy-normalised distance**: halve the more-doubled cell, then use
   today's distance.
   - As fast as today's distance and nearly as good as the full search on
     clean data (parallel-WGD label match 0.90 vs 0.95).
   - With low noise it was no better than having no WGD term (0.70), and
     worsened trees (RF 0.59 vs 0.52). Noisy bins make many values odd,
     which count as extra changes after the doubling.
2. **Clone-level shared vs parallel test**: compare one shared WGD with two
   separate WGDs on consensus profiles of three groups.
   - Two separate WGDs were called correctly only 40–50% of the time
     (median difference 0). Pooling all non-doubled cells into one
     consensus loses the information that separates the two stories.
   - On OV2295, E+F vs G+H+I was a tie (+7 events at 2 Mb, +1 at 500 kb).

The full search is kept.

## Event model (`event_model.cpp`, `event_model.R`)

**Model.** Every event changes the copy number of one contiguous segment by
exactly one copy (+1 or −1), as in copy-number-transformation models such
as MEDICC2. There is no restriction at copy number 0. WGD multiplies the
whole profile by 2.

**Distance without WGD** (exact, linear time):

    E(a, b) = TV(b − a) / 2,   TV(v) = |v₁| + Σ |vᵢ − vᵢ₋₁| + |vₙ|

That is, half the total size of the jumps in the difference, padded with
zeros at both ends.

**Ancestor problem** (exact):

    min_x  E(f₁·x, p) + E(f₂·x, q)

Solved by dynamic programming over bins, with x's integer value per bin as
the state; it covers every case we need:
- the WGD distance (f = 1 and 2 or 4);
- same-level ancestors (f₁ = f₂ = 1);
- one or two doubled children in the reconstruction.

**Options.** `event_matrices(X, use_bfb = TRUE)` keeps, per
chromosome/allele, the smaller of E and the existing BFB distance B, as
`fit()` does with G and B. `reconstruct_wgd_event()` rebuilds ancestors
with one exact call per node.

**Checks (`check_event_model.R`)**
- E equals the true minimum number of ±1 segment events, by breadth-first
  search, on 150 random pairs.
- The ancestor program equals an exhaustive search over x on 300 random
  cases.
- Same-level ancestors lie on a shortest path.
- Known cases:
  - x → 2x = 0 events;
  - +1 gain after doubling = 1;
  - x → 4x = 0;
  - 2 → 8 on one segment = 6 events (G counts 1).

### Results (`run_event_model.R`, `real_event_model.R`)

**Simulations** (same scenarios and seeds as the harder study; mean tree
error RF over all scenarios and both data types):

| Method | No noise | Low noise |
|---|---|---|
| Greedy G/B, no WGD term | 0.308 | 0.579 |
| **Greedy G/B + full WGD search** | **0.293** | **0.566** |
| Event model, exact WGD | 0.405 | 0.598 |
| Event model + BFB distance, exact WGD | 0.309 | 0.578 |

- **Event model + BFB vs the greedy full search:** worse in 108 of 160
  paired runs, better in 38.
- **WGD labels:**
  - Event model + BFB is similar to the greedy full search: label match on
    parallel WGDs 0.90 vs 0.95 with no noise, 0.95 vs 0.90 with low noise,
    allele-specific data.
  - The pure event model predicts too many WGD edges (2.9 for 2.4 true).
- **No false WGDs:** no method predicts a WGD in simulations without one.
- **Speed:** similar to the greedy full search; somewhat faster under
  noise (about 2 s vs 3 s per data set).

**Why the pure event model loses: BFB.** One BFB cycle duplicates a whole
stretch of an already rearranged chromosome, so an amplified bin can gain
several copies in one event, while the unit model counts one event per
copy. Without BFB (`bfb_prob = 0`, 6 replicates) it is as good as greedy
(RF 0.246 vs 0.253). With chromosome 1 (the BFB chromosome) excluded, the
two are equal (0.381 vs 0.379).

**OV2295 at 2 Mb**

| | Greedy full search | Event model + BFB, with WGD |
|---|---|---|
| Clone coherence | 0.748 | 0.743 |
| WGD clades | G+H+I; E+F | G+H+I; E+F |
| Distance time | 13 s | 9 s |

- **Even without a WGD term**, the event model keeps E+F and G+H+I in
  separate WGD clades (the greedy model merges them without its WGD term).
  This supports the two-WGD reading.

**Conclusions**
- Exact counting under a unit-event model does not give better trees than
  the greedy model on these simulations, because BFB amplification is not a
  unit-event process.
- With the BFB distance added, it is a close second, with the advantage
  that the WGD distance and ancestors are exact rather than heuristic.
- Options to make the event model competitive:
  - make the event cost depend less on magnitude, e.g. cost
    1 + λ·(copies − 1) per event, which is still solvable exactly by the
    same program with modified jump costs;
  - or add BFB cycles as an explicit event type.
