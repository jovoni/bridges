# WGD-aware distance: prototype

Prototype only; nothing here is used by the package. Run from the package
root:

```sh
Rscript dev/wgd/check_kernel.R                    # kernel checks
Rscript dev/wgd/run_experiments.R results.rds     # simulation study (~1 min)
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
- `wgd_proto.cpp`: kernel `wgd_cost()` and the per-chromosome matrices.
- `wgd_proto.R`: genome-wide components, combination, tree, and the true
  WGD ancestry of a simulation.

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
- Reconstruction (labelling WGD branches) is not prototyped yet.
