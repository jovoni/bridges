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
#' @param wgd_available Numeric. Number of whole-genome duplication events available. Default: 1
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
#' @param breakpoint_support Character. Distribution used for breakpoint
#'  selection ("uniform", "beta", etc.). Default: "uniform"
#' @param hotspot Named list. Hotspot positions for each chromosome allele
#'  (e.g., list(chr = "1:A", pos = 100), which is default)
#' @param alpha Numeric. First parameter for beta distribution if used for
#'  breakpoint selection. Default: NULL
#' @param beta Numeric. Second parameter for beta distribution if used for
#'  breakpoint selection. Default: NULL
#'
#' @return A list containing:
#'   \item{cells}{List of cell chromosome sequences at the end of simulation}
#'   \item{cell_history}{Data frame with cell birth/death times and lineage information}
#'   \item{input_parameters}{List of input parameters used in the simulation}
#'
#' @export
bridge_sim <- function(
  initial_cells = 1,
  chromosomes = c(1:22, "X", "Y"),
  bin_length = 1e6,
  birth_rate = 0.1,
  death_rate = 0.001,
  bfb_allele = "1:A",
  normal_dup_rate = 0,
  bfb_prob = 0.5,
  amp_rate = 1,
  del_rate = 1,
  wgd_available = 1,
  wgd_probability = .05,
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
) {
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
    beta
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
    custom_breakpoints = custom_breakpoints
  )

  # Initialize simulation state
  sim_state <- initialize_simulation(input_parameters)

  # Main simulation loop
  while (continue_simulation(sim_state)) {
    # Find the next cell to have an event
    next_event_info <- get_next_event(sim_state)
    sim_state$time <- next_event_info$time
    current_cell_id <- next_event_info$cell_id
    current_event_type <- next_event_info$event_type

    if (current_event_type == "birth") {
      # Handle birth event
      sim_state <- process_birth_event(
        state = sim_state,
        current_cell_id = current_cell_id,
        lambda = lambda,
        rate = rate
      )
    } else {
      # Handle death event
      sim_state <- process_death_event(sim_state, current_cell_id)
    }
  }

  # Finalize and prepare results (without subsampling)
  sim_state <- prepare_results(sim_state)
  sim_state$cna_data <- sequences_to_cndata(
    sequences = sim_state$cells,
    chr_seq_lengths = chr_seq_lengths,
    bin_length = bin_length
  )

  return(sim_state)
}
