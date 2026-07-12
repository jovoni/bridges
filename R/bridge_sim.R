#' Gillespie Simulation for Break-Fusion-Bridge (BFB) Processes
#' Modified to support diploid chromosomes (alleles A and B)
#'
#' @description
#' Simulates the evolution of cells undergoing Break-Fusion-Bridge cycles using
#' a continuous-time Gillespie algorithm. This function models cell birth and death
#' processes with the possibility of BFB events occurring during replication.
#' Cells with amplified hotspots have an increased birth rate.
#' Now supports modeling both alleles (A and B) of each chromosome.
#'
#' @param initial_cells Numeric. Number of cells at the start of simulation. Default: 1
#' @param chromosomes Character vector. Chromosomes to model (e.g., c("1", "2", "X")). Default: c(1:22, "X", "Y")
#' @param bin_length Numeric. Length of each genomic bin in base pairs. Default: 5e5
#' @param birth_rate Numeric. Base rate at which cells replicate per unit time. Default: 0.1
#' @param death_rate Numeric. Rate at which cells die per unit time. Default: 0.001
#' @param bfb_allele Character. Allele which will be affected by BFB events. Defaul : "1:A"
#' @param normal_dup_rate Numeric. Rate of normal duplication events. Default: 0.5
#' @param bfb_prob Numeric. Probability of BFB event occurring during replication. Default: 0.01
#' @param amp_rate Numeric. Rate of amplification events. Default: 0.1
#' @param del_rate Numeric. Rate of deletion events. Default: 0.1
#' @param wgd_available Numeric. Number of whole-genome duplication events available. Default: 0
#' @param wgd_probability Numeric. Probability of whole-genome duplication event occurring. Default: 0.05
#' @param lambda Rate parameter for Poisson distribution used to sample
#' the number of genomic events per daughter cell
#' @param rate Rate parameter used in amplification/deletion simulations.
#'  Length of event is sample from exponential distribution with parameter 1 / rate.
#' @param positive_selection_rate Numeric. Selection advantage for cells with amplified hotspot. Default: 0
#' @param negative_selection_rate Numeric. Selection disadvantage for cells without amplified hotspot. Default: 0
#' @param max_time Numeric. Maximum simulation time. Default: 50
#' @param max_cells Numeric. Maximum number of cells allowed before simulation stops. Default: 100
#' @param first_round_of_bfb Logical. Whether to apply BFB to initial cells. Default: TRUE
#' @param return_phylo Logical. Whether to build and return the phylogenetic tree.
#'   Set to FALSE when only CNA data is needed (e.g. ABC or clonal comparison
#'   workflows) to skip the expensive tree-building step. Default: TRUE
#' @param breakpoint_support Character. Distribution used for breakpoint
#'  selection ("uniform", "beta", etc.). Default: "beta"
#' @param hotspot Named list. Hotspot positions for each chromosome allele
#'  (e.g., list(chr = "1:A", pos = 100), which is default)
#' @param alpha Numeric. First parameter for beta distribution if used for
#'  breakpoint selection. Default: 50
#' @param beta Numeric. Second parameter for beta distribution if used for
#'  breakpoint selection. Default: 50
#' @param custom_breakpoints Integer vector of bin indices to use as the breakpoint
#'  position distribution when \code{breakpoint_support = "custom"}.  Each draw
#'  during the simulation is sampled uniformly from this vector.  Ignored for
#'  other breakpoint support modes.  Default: NULL
#' @param subsample Integer or NULL.  If set, this many alive cells are randomly
#'  sampled \emph{before} \code{prepare_results()} and \code{sequences_to_cndata()}
#'  are called, so the phylogeny and CNA tibble are built only once on the
#'  subsampled cells.  The \code{n_alive} field in the return value still records
#'  the pre-subsample count.  Default: NULL (keep all alive cells).
#'
#' @return A named list containing:
#' \describe{
#'   \item{cells}{Named list of alive cell chromosome sequences (internal interval
#'     representation).  Length equals \code{subsample} when set, otherwise the
#'     number of cells alive at simulation end.}
#'   \item{cell_history}{Tibble recording every birth/death event: \code{cell_id},
#'     \code{parent_id}, \code{bfb_event}, \code{wgd_event}, \code{cn_event},
#'     \code{chr_allele}, \code{is_alive}.}
#'   \item{tree}{A \code{phylo} object (from \pkg{ape}) pruned to alive cells.
#'     \code{NULL} when \code{return_phylo = FALSE}.}
#'   \item{cna_data}{Long-format tibble of copy number profiles with columns
#'     \code{cell_id}, \code{chr}, \code{start}, \code{end}, \code{CN}, \code{A},
#'     \code{B}.}
#'   \item{elapsed_time}{Numeric.  Simulation clock value at stop — proportional
#'     to evolutionary depth.}
#'   \item{n_alive}{Integer.  Number of alive cells before any subsampling.
#'     Equals \code{length(cells)} unless \code{subsample} was set.}
#'   \item{input_parameters}{List of all input parameters used.}
#' }
#'
#' @examples
#' \dontrun{
#' # Minimal simulation: chromosome 8, BFB on allele A, 128 cells
#' sim <- bridge_sim(
#'   chromosomes = "8",
#'   bfb_allele  = "8:A",
#'   max_cells   = 128,
#'   lambda      = 2
#' )
#' head(sim$cna_data)
#' cat("Elapsed time:", sim$elapsed_time, "\n")
#'
#' # Simulate 1000 cells but only post-process 100 (faster)
#' sim2 <- bridge_sim(
#'   chromosomes = "1",
#'   max_cells   = 1000,
#'   subsample   = 100,
#'   return_phylo = FALSE
#' )
#' cat("Alive:", sim2$n_alive, "  Sequenced:", length(sim2$cells), "\n")
#' }
#'
#' @export
bridge_sim <- function(
  initial_cells = 1,
  initial_sequences = NULL,
  chromosomes = c(1:22, "X", "Y"),
  bin_length = 1e6,
  birth_rate = 0.1,
  death_rate = 0.001,
  bfb_allele = "1:A",
  normal_dup_rate = 0,
  bfb_prob = 0.5,
  amp_rate = 1,
  del_rate = 1,
  wgd_available = 0,
  wgd_probability = .05,
  lambda = 2,
  rate = 20,
  positive_selection_rate = 0,
  negative_selection_rate = 0,
  max_time = 300,
  max_cells = 256,
  subsample = NULL,
  first_round_of_bfb = TRUE,
  return_phylo = TRUE,
  return_cna_data = TRUE,
  breakpoint_support = "beta",
  hotspot = list(chr = "1:A", pos = 100),
  alpha = 50,
  beta = 50,
  custom_breakpoints = NULL,
  selection_type = "constant",
  saturation_K = 10
) {
  # When seeding from pre-evolved sequences, derive initial_cells from them.
  if (!is.null(initial_sequences)) {
    initial_cells <- length(initial_sequences)
  }

  validate_bridge_sim_params(
    initial_cells,
    chromosomes,
    bin_length,
    birth_rate,
    death_rate,
    bfb_allele,
    normal_dup_rate,
    bfb_prob,
    amp_rate,
    del_rate,
    positive_selection_rate,
    negative_selection_rate,
    max_time,
    max_cells,
    first_round_of_bfb,
    breakpoint_support,
    hotspot,
    alpha,
    beta,
    selection_type,
    saturation_K
  )

  # Default human chromosome lengths (approximate, in base pairs)
  default_chr_lengths <- c(
    "1" = 247249719,
    "2" = 242193529,
    "3" = 198295559,
    "4" = 190214555,
    "5" = 181538259,
    "6" = 170805979,
    "7" = 159345973,
    "8" = 145138636,
    "9" = 138394717,
    "10" = 133797422,
    "11" = 135086622,
    "12" = 133275309,
    "13" = 114364328,
    "14" = 107043718,
    "15" = 101991189,
    "16" = 90338345,
    "17" = 83257441,
    "18" = 80373285,
    "19" = 58617616,
    "20" = 64444167,
    "21" = 46709983,
    "22" = 50818468,
    "X" = 156040895,
    "Y" = 57227415
  )
  chr_lengths <- default_chr_lengths[as.character(chromosomes)]

  # Calculate sequence lengths for each chromosome
  chr_seq_lengths <- round(chr_lengths / bin_length)
  names(chr_seq_lengths) <- names(chr_lengths)

  # Create chromosome allele names (A and B for each chromosome)
  chr_alleles <- paste0(rep(names(chr_seq_lengths), each = 2), ":", c("A", "B"))

  # Normalize rates
  sum_rates <- sum(normal_dup_rate, bfb_prob, amp_rate, del_rate)
  normal_dup_rate <- normal_dup_rate / sum_rates
  bfb_prob <- bfb_prob / sum_rates
  amp_rate <- amp_rate / sum_rates
  del_rate <- del_rate / sum_rates

  rates <- list(
    normal = normal_dup_rate,
    bfb = bfb_prob,
    amp = amp_rate,
    del = del_rate
  )

  # Init state with parameters
  input_parameters <- list(
    initial_cells = initial_cells,
    chromosomes = chromosomes,
    chr_alleles = chr_alleles,
    chr_seq_lengths = chr_seq_lengths,
    bin_length = bin_length,
    chr_lengths = chr_lengths,
    birth_rate = birth_rate,
    death_rate = death_rate,
    bfb_allele = bfb_allele,
    rates = rates,
    wgd_available = wgd_available,
    wgd_probability = wgd_probability,
    positive_selection_rate = positive_selection_rate,
    negative_selection_rate = negative_selection_rate,
    max_time = max_time,
    max_cells = max_cells,
    first_round_of_bfb = first_round_of_bfb,
    breakpoint_support = breakpoint_support,
    hotspot = hotspot,
    alpha = alpha,
    beta = beta,
    custom_breakpoints = custom_breakpoints,
    selection_type = selection_type,
    saturation_K = saturation_K
  )

  # Initialize simulation state
  if (is.null(initial_sequences)) {
    sim_state <- initialize_simulation(input_parameters)
  } else {
    sim_state <- initialize_from_sequences(initial_sequences, input_parameters)
  }

  # Main simulation loop (C++ Gillespie engine)
  sim_state <- bridge_sim_loop_cpp(sim_state, lambda = lambda, rate = rate)

  # Optional pre-processing subsample: discard surplus alive cells before
  # prepare_results() so phylogeny and CNA tibble are built only once on k cells
  n_alive_total <- length(sim_state$cell_ids)
  if (!is.null(subsample) && subsample < n_alive_total) {
    keep_ids    <- sample(sim_state$cell_ids, subsample)
    discard_ids <- setdiff(sim_state$cell_ids, keep_ids)

    sim_state$cell_ids       <- keep_ids
    sim_state$cell_sequences <- sim_state$cell_sequences[keep_ids]

    keep_mask <- !(sim_state$h_cell_id %in% discard_ids)
    sim_state$h_cell_id    <- sim_state$h_cell_id[keep_mask]
    sim_state$h_parent_id  <- sim_state$h_parent_id[keep_mask]
    sim_state$h_bfb_event  <- sim_state$h_bfb_event[keep_mask]
    sim_state$h_wgd_event  <- sim_state$h_wgd_event[keep_mask]
    sim_state$h_cn_event   <- sim_state$h_cn_event[keep_mask]
    sim_state$h_chr_allele <- sim_state$h_chr_allele[keep_mask]
    sim_state$h_n          <- sum(keep_mask)
  }
  sim_state$n_alive_total <- n_alive_total

  # Finalize and prepare results
  sim_state <- prepare_results(sim_state, return_phylo = return_phylo)

  # Prune tree to exactly the alive cells: build_phylo_from_lineage() includes
  # dead-cell leaf nodes; ape::keep.tip() removes them so tree tips == cna_data cells.
  if (!is.null(sim_state$tree)) {
    kept_ids <- names(sim_state$cells)
    sim_state$tree <- ape::keep.tip(sim_state$tree, kept_ids)
  }

  sim_state$cna_data <- if (return_cna_data) {
    sequences_to_cndata(
      sequences       = sim_state$cells,
      chr_seq_lengths = chr_seq_lengths,
      bin_length      = bin_length
    )
  } else {
    NULL
  }

  return(sim_state)
}
