#' Simulate single-cell copy-number evolution with breakage-fusion-bridge cycles
#'
#' @description
#' Simulates a growing population of diploid cells (alleles A and B for every
#' chromosome) using a continuous-time Gillespie birth-death process. At each
#' division, each daughter receives a Poisson(\code{lambda}) number of genomic
#' events drawn from: breakage-fusion-bridge (BFB) cycles on \code{bfb_allele},
#' focal amplifications, focal deletions, neutral "no-op" events and,
#' optionally, whole-genome duplication (WGD). Cells carrying extra copies of
#' a \code{hotspot} bin can be given a birth (and/or death) rate advantage.
#'
#' @param initial_cells Positive integer. Number of founder cells. Ignored
#'   when \code{initial_sequences} is supplied. Default: 1
#' @param initial_sequences Optional list of pre-evolved cell sequences to
#'   seed the simulation from (the \code{cells} element of a previous
#'   \code{bridge_sim()} result), e.g. for serial-passage workflows. When
#'   supplied, the number of founders is \code{length(initial_sequences)}.
#' @param chromosomes Character vector of chromosomes to model
#'   (any of \code{1:22}, \code{"X"}, \code{"Y"}). Default: \code{c(1:22, "X", "Y")}
#' @param genome Genome build used to look up chromosome lengths: \code{"hg19"}
#'   or \code{"hg38"}. See \code{\link{centromeres}}. Default: \code{"hg19"}
#' @param bin_length Bin size in base pairs. Default: 1e6
#' @param birth_rate Base per-cell division rate. Default: 0.1
#' @param death_rate Base per-cell death rate. Default: 0.001
#' @param bfb_allele The allele undergoing BFB, as \code{"chr:allele"}
#'   (e.g. \code{"8:A"}). Default: \code{"1:A"}
#' @param normal_dup_rate,bfb_prob,amp_rate,del_rate Non-negative
#'   \emph{relative weights} of the four event types (neutral no-op, BFB cycle,
#'   focal amplification, focal deletion). They are normalised to sum to one,
#'   so only their ratios matter. Defaults: 0, 0.5, 1, 1
#' @param wgd_available Non-negative integer. Maximum number of whole-genome
#'   duplications allowed in the whole simulation. Default: 0
#' @param wgd_probability Probability that a division attempts a WGD (while
#'   \code{wgd_available > 0}). Default: 0.05
#' @param lambda Mean of the Poisson number of genomic events per daughter
#'   cell. Default: 2
#' @param rate Mean length, in bins, of focal amplifications and deletions
#'   (lengths are exponentially distributed). Default: 20
#' @param positive_selection_rate Relative birth-rate increase for cells with
#'   an amplified hotspot (see \code{selection_type}). Default: 0
#' @param negative_selection_rate Relative death-rate increase for cells with
#'   an amplified hotspot (see \code{selection_type}). Default: 0
#' @param max_time Simulation stops when the clock reaches this time.
#'   Default: 300
#' @param max_cells Simulation stops when this many cells are alive.
#'   Default: 256
#' @param subsample \code{NULL} or a positive integer. If set, this many alive
#'   cells are sampled at random before building the tree and copy-number
#'   table (faster for large populations). \code{n_alive} still reports the
#'   pre-subsampling count. Default: \code{NULL}
#' @param first_round_of_bfb Logical. If \code{TRUE}, each founder cell
#'   undergoes one BFB cycle before the simulation starts. Default: \code{TRUE}
#' @param return_phylo Logical. Whether to build and return the true lineage
#'   tree (\code{tree}). Set to \code{FALSE} to skip this step. Default: \code{TRUE}
#' @param return_cna_data Logical. Whether to build and return the long-format
#'   copy-number table (\code{cna_data}). Set to \code{FALSE} when only the
#'   \code{cells} sequences are needed (e.g. intermediate serial passages).
#'   Default: \code{TRUE}
#' @param breakpoint_support Distribution of BFB breakpoints along the fused
#'   (dicentric) chromosome: \code{"beta"}, \code{"uniform"} or
#'   \code{"custom"}. Default: \code{"beta"}
#' @param hotspot List with elements \code{chr} (\code{"chr:allele"}) and
#'   \code{pos} (bin index) defining the bin whose copy number drives
#'   selection. Its chromosome must be in \code{chromosomes}.
#'   Default: \code{list(chr = "1:A", pos = 100)}
#' @param alpha,beta Shape parameters of the beta breakpoint distribution
#'   (\code{breakpoint_support = "beta"}). Default: 10, 10
#' @param custom_breakpoints Integer vector of bin indices eligible as BFB
#'   breakpoints when \code{breakpoint_support = "custom"}; at each BFB cycle
#'   the cut is drawn uniformly among the positions of the fused chromosome
#'   carrying one of these bins. A cell whose BFB allele no longer carries
#'   any of these bins cannot undergo further BFB cycles. Default: \code{NULL}
#' @param selection_type How hotspot copy number \eqn{n} maps to the selection
#'   multiplier \eqn{s} applied to \code{positive_selection_rate} and
#'   \code{negative_selection_rate}: \code{"constant"} (\eqn{s = 1} if
#'   \eqn{n > 1}, else 0), \code{"linear"} (\eqn{s = n - 1}) or
#'   \code{"saturation"} (\eqn{s = (n-1)/((n-1) + K)}). Default: \code{"constant"}
#' @param saturation_K Saturation constant \eqn{K} for
#'   \code{selection_type = "saturation"}. Default: 10
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
#'   \item{elapsed_time}{Numeric. Simulation clock value at stop (at most
#'     \code{max_time}).}
#'   \item{n_alive}{Integer.  Number of alive cells before any subsampling.
#'     Equals \code{length(cells)} unless \code{subsample} was set.}
#'   \item{input_parameters}{List of all input parameters used.}
#' }
#' If the population goes extinct (no cells alive), a warning is given and
#' \code{cells} is empty, \code{tree} is \code{NULL} and \code{cna_data} has
#' no rows.
#'
#' @examples
#' set.seed(1)
#' # BFB on chromosome 8, allele A; the selection hotspot must lie on a
#' # simulated chromosome too
#' sim <- bridge_sim(
#'   chromosomes = "8",
#'   bfb_allele  = "8:A",
#'   hotspot     = list(chr = "8:A", pos = 100),
#'   max_cells   = 32
#' )
#' head(sim$cna_data)
#' sim$tree
#'
#' # Grow 200 cells but only post-process 50 of them (faster)
#' sim2 <- bridge_sim(
#'   chromosomes  = "1",
#'   max_cells    = 200,
#'   subsample    = 50,
#'   return_phylo = FALSE
#' )
#' c(alive = sim2$n_alive, sampled = length(sim2$cells))
#'
#' @export
bridge_sim <- function(
  initial_cells = 1,
  initial_sequences = NULL,
  chromosomes = c(1:22, "X", "Y"),
  genome = "hg19",
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
  alpha = 10,
  beta = 10,
  custom_breakpoints = NULL,
  selection_type = "constant",
  saturation_K = 10
) {
  # When seeding from pre-evolved sequences, derive initial_cells from them.
  if (!is.null(initial_sequences)) {
    initial_cells <- length(initial_sequences)
  }

  # Real chromosome lengths for the requested genome build, from the
  # bundled `centromeres` dataset (see ?centromeres).
  chr_lengths <- .chr_lengths_for(genome, chromosomes)
  chr_seq_lengths <- round(chr_lengths / bin_length)
  names(chr_seq_lengths) <- names(chr_lengths)

  validate_bridge_sim_params(
    initial_cells,
    chr_seq_lengths,
    lambda,
    rate,
    wgd_available,
    wgd_probability,
    subsample,
    custom_breakpoints,
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
    genome = genome,
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

  # All cells died: nothing to sequence. Return an empty result instead of
  # failing while building the tree and copy-number table.
  extinct <- length(sim_state$cell_ids) == 0
  if (extinct) {
    warning("The population went extinct (no cells alive at time ",
            signif(sim_state$time, 3), "); returning no cells.", call. = FALSE)
  }

  # Finalize and prepare results
  sim_state <- prepare_results(sim_state, return_phylo = return_phylo && !extinct)

  # Prune tree to exactly the alive cells: build_phylo_from_lineage() includes
  # dead-cell leaf nodes; ape::keep.tip() removes them so tree tips == cna_data cells.
  if (!is.null(sim_state$tree)) {
    kept_ids <- names(sim_state$cells)
    sim_state$tree <- ape::keep.tip(sim_state$tree, kept_ids)
  }

  sim_state$cna_data <- if (return_cna_data && extinct) {
    dplyr::tibble(cell_id = character(0), chr = character(0), bin_idx = integer(0),
                  A = numeric(0), B = numeric(0), CN = numeric(0),
                  start = numeric(0), end = numeric(0))
  } else if (return_cna_data) {
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
