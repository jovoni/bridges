# Gillespie Simulation for Break-Fusion-Bridge (BFB) Processes Modified to support diploid chromosomes (alleles A and B)

Simulates the evolution of cells undergoing Break-Fusion-Bridge cycles
using a continuous-time Gillespie algorithm. This function models cell
birth and death processes with the possibility of BFB events occurring
during replication. Cells with amplified hotspots have an increased
birth rate. Now supports modeling both alleles (A and B) of each
chromosome.

## Usage

``` r
bridge_sim(
  initial_cells = 1,
  chromosomes = c(1:22, "X", "Y"),
  bin_length = 1e+06,
  birth_rate = 0.1,
  death_rate = 0.001,
  bfb_allele = "1:A",
  normal_dup_rate = 0,
  bfb_prob = 0.5,
  amp_rate = 1,
  del_rate = 1,
  wgd_available = 0,
  wgd_probability = 0.05,
  lambda = 2,
  rate = 20,
  positive_selection_rate = 0,
  negative_selection_rate = 0,
  max_time = 300,
  max_cells = 256,
  first_round_of_bfb = TRUE,
  breakpoint_support = "uniform",
  hotspot = list(chr = "1:A", pos = 100),
  alpha = NULL,
  beta = NULL,
  custom_breakpoints = NULL
)
```

## Arguments

- initial_cells:

  Numeric. Number of cells at the start of simulation. Default: 1

- chromosomes:

  Character vector. Chromosomes to model (e.g., c("1", "2", "X")).
  Default: c(1:22, "X", "Y")

- bin_length:

  Numeric. Length of each genomic bin in base pairs. Default: 5e5

- birth_rate:

  Numeric. Base rate at which cells replicate per unit time. Default:
  0.1

- death_rate:

  Numeric. Rate at which cells die per unit time. Default: 0.001

- bfb_allele:

  Character. Allele which will be affected by BFB events. Defaul : "1:A"

- normal_dup_rate:

  Numeric. Rate of normal duplication events. Default: 0.5

- bfb_prob:

  Numeric. Probability of BFB event occurring during replication.
  Default: 0.01

- amp_rate:

  Numeric. Rate of amplification events. Default: 0.1

- del_rate:

  Numeric. Rate of deletion events. Default: 0.1

- wgd_available:

  Numeric. Number of whole-genome duplication events available. Default:
  0

- wgd_probability:

  Numeric. Probability of whole-genome duplication event occurring.
  Default: 0.05

- lambda:

  Rate parameter for Poisson distribution used to sample the number of
  genomic events per daughter cell

- rate:

  Rate parameter used in amplification/deletion simulations. Length of
  event is sample from exponential distribution with parameter 1 / rate.

- positive_selection_rate:

  Numeric. Selection advantage for cells with amplified hotspot.
  Default: 0

- negative_selection_rate:

  Numeric. Selection disadvantage for cells without amplified hotspot.
  Default: 0

- max_time:

  Numeric. Maximum simulation time. Default: 50

- max_cells:

  Numeric. Maximum number of cells allowed before simulation stops.
  Default: 100

- first_round_of_bfb:

  Logical. Whether to apply BFB to initial cells. Default: TRUE

- breakpoint_support:

  Character. Distribution used for breakpoint selection ("uniform",
  "beta", etc.). Default: "uniform"

- hotspot:

  Named list. Hotspot positions for each chromosome allele (e.g.,
  list(chr = "1:A", pos = 100), which is default)

- alpha:

  Numeric. First parameter for beta distribution if used for breakpoint
  selection. Default: NULL

- beta:

  Numeric. Second parameter for beta distribution if used for breakpoint
  selection. Default: NULL

## Value

A list containing:

- cells:

  List of cell chromosome sequences at the end of simulation

- cell_history:

  Data frame with cell birth/death times and lineage information

- input_parameters:

  List of input parameters used in the simulation
