#' Gillespie Simulation for Break-Fusion-Bridge (BFB) Processes with Hotspot Selection
#'
#' @description
#' Simulates the evolution of cells undergoing Break-Fusion-Bridge cycles using
#' a continuous-time Gillespie algorithm. This function models cell birth and death
#' processes with the possibility of BFB events occurring during replication.
#' Cells with amplified hotspots have an increased birth rate.
#'
#' @param initial_cells Numeric. Number of cells at the start of simulation. Default: 1
#' @param initial_sequence_length Numeric. Length of the initial genomic sequence. Default: 100
#' @param birth_rate Numeric. Base rate at which cells replicate per unit time. Default: 0.1
#' @param death_rate Numeric. Rate at which cells die per unit time. Default: 0.001
#' @param bfb_prob Numeric. Probability of BFB event occurring during replication. Default: 0.01
#' @param positive_selection_rate Numeric. Selection advantage for cells with amplified hotspot. Default: 0
#' @param negative_selection_rate Numeric. Selection disadvantage for cells without amplified hotspot. Default: 0
#' @param positive_selection_function Function. A function that takes as input a positive selection rate and a number
#'  of copies and returns a value between -1 and Infinity. This determines how selection advantage scales with hotspot amplification.
#' @param negative_selection_function Function. A function that takes as input a negative selection rate and a number
#'  of copies and returns a value between -1 and Infinity. This determines how selection disadvantage scales with hotspot amplification.
#' @param max_time Numeric. Maximum simulation time. Default: 50
#' @param max_cells Numeric. Maximum number of cells allowed before simulation stops. Default: 100
#' @param first_round_of_bfb Logical. Whether to apply BFB to initial cells. Default: TRUE
#' @param breakpoint_support Character. Distribution used for breakpoint selection ("uniform", "beta", etc.). Default: "uniform"
#' @param hotspot Numeric vector. Genomic positions considered as hotspots. Default: NULL
#' @param alpha Numeric. First parameter for beta distribution if used for breakpoint selection. Default: NULL
#' @param beta Numeric. Second parameter for beta distribution if used for breakpoint selection. Default: NULL
#'
#' @return A list containing:
#'   \item{final_cells}{List of cell sequences at the end of simulation}
#'   \item{sequence_distribution}{Frequency table of unique sequences}
#'   \item{cell_history}{Data frame with cell birth/death times and lineage information}
#'   \item{final_time}{Time at which simulation ended}
#'   \item{birth_count}{Total number of birth events}
#'   \item{death_count}{Total number of death events}
#'   \item{newick_tree}{Newick format representation of the phylogenetic tree}
#'   \item{phylogenetic_info}{List containing BFB event history for each cell, tracking gain/loss events}
#'   \item{hotspot_info}{Data frame tracking hotspot amplification status of cells}
#'   \item{parameters}{List of input parameters used in the simulation}
#'
#' @export
gillespie_sim <- function(
    initial_cells = 1,
    initial_sequence_length = 100,
    birth_rate = 0.1,
    death_rate = 0.001,
    normal_dup_rate = .5,
    bfb_prob = 0.01,
    amp_rate = .1,
    del_rate = .1,
    allow_wgd = TRUE,
    positive_selection_rate = 0,
    negative_selection_rate = 0,
    positive_selection_function = positive_selection_function,
    negative_selection_function = negative_selection_function,
    max_time = 50,
    max_cells = 100,
    #first_n_bfb_cycles = 0,
    first_round_of_bfb = TRUE,
    breakpoint_support = "uniform",
    hotspot = NULL,
    alpha = NULL,
    beta = NULL
) {
  # Normalize rates
  sum_rates = sum(normal_dup_rate, bfb_prob, amp_rate, del_rate)
  normal_dup_rate = normal_dup_rate / sum_rates
  bfb_prob = bfb_prob / sum_rates
  amp_rate = amp_rate / sum_rates
  del_rate = del_rate / sum_rates
  rates = list(
    normal=normal_dup_rate,
    bfb=bfb_prob,
    amp=amp_rate,
    del=del_rate
  )

  # Init state with parameters
  input_parameters = list(
    initial_cells = initial_cells,
    initial_sequence_length = initial_sequence_length,
    birth_rate = birth_rate,
    death_rate = death_rate,
    # normal_dup_rate = normal_dup_rate,
    # bfb_prob = bfb_prob,
    # amp_rate = amp_rate,
    # del_rate = del_rate,
    rates = rates,
    allow_wgd = allow_wgd,
    positive_selection_rate = positive_selection_rate,
    negative_selection_rate = negative_selection_rate,
    positive_selection_function = positive_selection_function,
    negative_selection_function = negative_selection_function,
    max_time = max_time,
    max_cells = max_cells,
    #first_n_bfb_cycles = first_n_bfb_cycles,
    first_round_of_bfb = first_round_of_bfb,
    breakpoint_support = breakpoint_support,
    hotspot = hotspot,
    alpha = alpha,
    beta = beta
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
      sim_state <- process_birth_event(sim_state, current_cell_id)
    } else {
      # Handle death event
      sim_state <- process_death_event(sim_state, current_cell_id)
    }

  }

  sim_state$cell_replication_history %>% table()


  # Finalize and prepare results
  sim_state = prepare_results(sim_state)
  return(sim_state)
}
