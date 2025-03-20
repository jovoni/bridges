
#' Gillespie Simulation for Break-Fusion-Bridge (BFB) Processes
#'
#' @description
#' Simulates the evolution of cells undergoing Break-Fusion-Bridge cycles using
#' a continuous-time Gillespie algorithm. This function models cell birth and death
#' processes with the possibility of BFB events occurring during replication.
#'
#' @param initial_cells Numeric. Number of cells at the start of simulation. Default: 1
#' @param initial_sequence_length Numeric. Length of the initial genomic sequence. Default: 100
#' @param birth_rate Numeric. Rate at which cells replicate per unit time. Default: 0.1
#' @param death_rate Numeric. Rate at which cells die per unit time. Default: 0.001
#' @param bfb_prob Numeric. Probability of BFB event occurring during replication. Default: 0.01
#' @param max_time Numeric. Maximum simulation time. Default: 50
#' @param max_cells Numeric. Maximum number of cells allowed before simulation stops. Default: 100
#' @param first_n_bfb_cycles Numeric. Number of initial birth events that will force BFB events. Default: 0
#' @param first_round_of_bfb Logical. Whether to apply BFB to initial cells. Default: TRUE
#' @param breakpoint_support Character. Distribution used for breakpoint selection ("uniform", "beta", etc.). Default: "uniform"
#' @param alpha Numeric. First parameter for beta distribution if used for breakpoint selection. Default: NULL
#' @param beta Numeric. Second parameter for beta distribution if used for breakpoint selection. Default: NULL
#'
#' @return A list containing:
#'   \item{final_cells}{List of cell sequences at the end of simulation}
#'   \item{sequence_distribution}{Frequency table of unique sequences}
#'   \item{cell_lifetimes}{Data frame with cell birth/death times and lineage information}
#'   \item{final_time}{Time at which simulation ended}
#'   \item{birth_count}{Total number of birth events}
#'   \item{death_count}{Total number of death events}
#'   \item{newick_tree}{Newick format representation of the phylogenetic tree}
#'   \item{phylogenetic_info}{List containing BFB event history for each cell, tracking gain/loss events}
#'   \item{parameters}{List of input parameters used in the simulation}
#'
#' @details
#' The function implements a Gillespie algorithm to simulate stochastic birth and death
#' processes with possible BFB events. Events occur at rates proportional to the number of cells
#' and depend on time since the last event. The simulation tracks cell lineages, creating
#' a phylogenetic history of BFB events.
#'
#' BFB (Break-Fusion-Bridge) cycles are a mechanism of genomic instability where
#' chromosome breakage followed by fusion of broken ends creates genomic rearrangements.
#' Each BFB event produces two daughter cells with different rearranged sequences.
#'
#' The phylogenetic information tracks whether each BFB event resulted in a sequence
#' length gain or loss for each affected cell lineage.
#'
#' @export
gillespie_sim <- function(
    initial_cells = 1,
    initial_sequence_length = 100,
    birth_rate = 0.1,
    death_rate = 0.001,
    bfb_prob = 0.01,
    max_time = 50,
    max_cells = 100,
    first_n_bfb_cycles = 0,
    first_round_of_bfb = TRUE,
    breakpoint_support = "uniform",
    alpha = NULL,
    beta = NULL
) {
  # Initialize
  time <- 0
  cells <- list()
  birth_count <- 0  # Track number of birth events
  death_count <- 0

  # Track cell creation times, IDs, parents, whether they resulted from BFB, and last event time
  cell_ids <- integer(0)
  cell_birth_times <- numeric(0)
  cell_death_times <- numeric(0)
  cell_parents <- integer(0)  # Track parent of each cell
  cell_bfb <- logical(0)  # Track whether the cell resulted from a BFB event
  cell_last_event <- numeric(0)  # Track the last time each cell underwent an event
  next_cell_id <- 1

  # Track BFB event history for phylogenetic information
  bfb_events <- list()
  next_bfb_id <- 1
  cell_bfb_history <- list()  # Store BFB event history for each cell

  # Create the initial cell with its sequence
  initial_sequence <- vec2seq(1:initial_sequence_length)
  initial_length <- get_seq_length(initial_sequence)

  # Apply BFB to all initial cells if first_round_of_bfb is TRUE
  if (first_round_of_bfb) {
    # Create temporary storage for new cells after BFB
    new_cells <- list()
    new_cell_ids <- integer(0)
    new_birth_times <- numeric(0)
    new_death_times <- numeric(0)
    new_parents <- integer(0)
    new_bfb <- logical(0)
    new_last_event <- numeric(0)
    new_bfb_history <- list()

    for (i in 1:initial_cells) {
      # Apply BFB to create daughter sequences
      daughter_seqs <- sim_bfb_left_and_right_sequences(initial_sequence, breakpoint_support, alpha, beta)

      # Error checking
      if (any(is.na(daughter_seqs$l_seq))) {
        stop("NA values in left daughter sequence")
      }

      if (any(is.na(daughter_seqs$r_seq))) {
        stop("NA values in right daughter sequence")
      }

      # Record BFB event
      bfb_events[[next_bfb_id]] <- list(
        id = next_bfb_id,
        time = time,
        parent_id = 0,
        parent_cell_id = 0
      )

      # Determine if each daughter experienced gain or loss
      l_length <- get_seq_length(daughter_seqs$l_seq)
      r_length <- get_seq_length(daughter_seqs$r_seq)

      l_effect <- ifelse(l_length > initial_length, "gain", "loss")
      r_effect <- ifelse(r_length > initial_length, "gain", "loss")

      # Create BFB history for left and right daughters
      left_history <- list(paste(next_bfb_id, l_effect, sep = "-"))
      right_history <- list(paste(next_bfb_id, r_effect, sep = "-"))

      # Add both daughters to new cells list
      new_cells[[length(new_cells) + 1]] <- daughter_seqs$l_seq
      new_cells[[length(new_cells) + 1]] <- daughter_seqs$r_seq

      # Add BFB history for both daughters
      new_bfb_history[[length(new_bfb_history) + 1]] <- left_history
      new_bfb_history[[length(new_bfb_history) + 1]] <- right_history

      # Record cell creation data for both new cells
      new_cell_ids <- c(new_cell_ids, next_cell_id, next_cell_id + 1)
      new_birth_times <- c(new_birth_times, time, time)
      new_death_times <- c(new_death_times, NA, NA)  # Not dead yet
      new_parents <- c(new_parents, 0, 0)  # Initial cells have no parent
      new_bfb <- c(new_bfb, TRUE, TRUE)  # Both daughters resulted from BFB
      new_last_event <- c(new_last_event, time, time)  # Set last event time to current time
      next_cell_id <- next_cell_id + 2
      next_bfb_id <- next_bfb_id + 1
    }

    # Replace the original cells with new BFB-processed cells
    cells <- new_cells
    cell_ids <- new_cell_ids
    cell_birth_times <- new_birth_times
    cell_death_times <- new_death_times
    cell_parents <- new_parents
    cell_bfb <- new_bfb
    cell_last_event <- new_last_event
    cell_bfb_history <- new_bfb_history
  } else {
    # Original initialization without BFB
    for (i in 1:initial_cells) {
      cells[[i]] <- initial_sequence

      # Record cell creation data
      cell_ids <- c(cell_ids, next_cell_id)
      cell_birth_times <- c(cell_birth_times, time)
      cell_death_times <- c(cell_death_times, NA)  # Not dead yet
      cell_parents <- c(cell_parents, 0)  # Initial cells have no parent
      cell_bfb <- c(cell_bfb, FALSE)  # Initial cells did not result from BFB
      cell_last_event <- c(cell_last_event, time)  # Set last event time to current time
      cell_bfb_history[[i]] <- list()  # Empty BFB history for initial cells
      next_cell_id <- next_cell_id + 1
    }
  }

  # Main simulation loop
  while (time < max_time && length(cells) > 0 && length(cells) < max_cells) {
    # Calculate propensities
    total_birth_prop <- birth_rate * length(cells)
    total_death_prop <- death_rate * length(cells)
    total_prop <- total_birth_prop + total_death_prop

    # If all cells died or total propensity is 0, exit
    if (length(cells) == 0 || total_prop == 0) {
      break
    }

    # Time until next event
    dt <- stats::rexp(1, total_prop)
    time <- time + dt

    # Determine event type
    if (stats::runif(1) < total_birth_prop / total_prop) {
      # Birth event
      birth_count <- birth_count + 1

      # Select parent cell based on time since last event
      time_since_last_event <- time - cell_last_event
      selection_probs <- time_since_last_event / sum(time_since_last_event)
      parent_idx <- sample(1:length(cells), 1, prob = selection_probs)
      parent_seq <- cells[[parent_idx]]
      parent_id <- cell_ids[parent_idx]
      parent_bfb_history <- cell_bfb_history[[parent_idx]]
      parent_length <- get_seq_length(parent_seq)

      # Update last event time for the parent cell
      cell_last_event[parent_idx] <- time

      # Determine if BFB should occur
      force_bfb <- birth_count <= first_n_bfb_cycles
      random_bfb <- stats::runif(1) < bfb_prob

      if (force_bfb || random_bfb) {
        # Apply the mutation function to get two daughter sequences
        daughter_seqs <- sim_bfb_left_and_right_sequences(parent_seq, breakpoint_support, alpha, beta)

        if (any(is.na(daughter_seqs$l_seq))) {
          stop("NA values in left daughter sequence")
        }

        if (any(is.na(daughter_seqs$r_seq))) {
          stop("NA values in right daughter sequence")
        }

        # Record BFB event
        bfb_events[[next_bfb_id]] <- list(
          id = next_bfb_id,
          time = time,
          parent_id = parent_id,
          parent_cell_id = parent_id
        )

        # Determine if each daughter experienced gain or loss compared to parent
        l_length <- get_seq_length(daughter_seqs$l_seq)
        r_length <- get_seq_length(daughter_seqs$r_seq)

        l_effect <- ifelse(l_length > parent_length, "gain", "loss")
        r_effect <- ifelse(r_length > parent_length, "gain", "loss")

        # Create BFB histories for both daughters
        left_history <- c(parent_bfb_history, list(paste(next_bfb_id, l_effect, sep = "-")))
        right_history <- c(parent_bfb_history, list(paste(next_bfb_id, r_effect, sep = "-")))

        # Replace parent with left daughter
        cells[[parent_idx]] <- daughter_seqs$l_seq
        cell_bfb_history[[parent_idx]] <- left_history

        # Add right daughter as new cell
        cells[[length(cells) + 1]] <- daughter_seqs$r_seq
        cell_bfb_history[[length(cell_bfb_history) + 1]] <- right_history

        # Record cell creation data for the new cell
        cell_ids <- c(cell_ids, next_cell_id)
        cell_birth_times <- c(cell_birth_times, time)
        cell_death_times <- c(cell_death_times, NA)  # Not dead yet
        cell_parents <- c(cell_parents, parent_id)  # Record parent
        cell_bfb <- c(cell_bfb, TRUE)  # New cell resulted from BFB
        cell_last_event <- c(cell_last_event, time)  # Set last event time for new cell
        next_cell_id <- next_cell_id + 1
        next_bfb_id <- next_bfb_id + 1
      } else {
        # Normal replication without mutation
        # Parent cell remains unchanged
        # Add identical daughter cell
        cells[[length(cells) + 1]] <- parent_seq
        cell_bfb_history[[length(cell_bfb_history) + 1]] <- parent_bfb_history

        # Record cell creation data for the new cell
        cell_ids <- c(cell_ids, next_cell_id)
        cell_birth_times <- c(cell_birth_times, time)
        cell_death_times <- c(cell_death_times, NA)  # Not dead yet
        cell_parents <- c(cell_parents, parent_id)  # Record parent
        cell_bfb <- c(cell_bfb, FALSE)  # New cell did not result from BFB
        cell_last_event <- c(cell_last_event, time)  # Set last event time for new cell
        next_cell_id <- next_cell_id + 1
      }

    } else {
      # Death event - remove a random cell
      death_count <- death_count + 1

      if (length(cells) > 0) {
        # Select cell to remove based on time since last event
        time_since_last_event <- time - cell_last_event
        selection_probs <- time_since_last_event / sum(time_since_last_event)
        cell_to_remove <- sample(1:length(cells), 1, prob = selection_probs)

        # Record death time
        cell_death_times[cell_to_remove] <- time

        # Remove the cell from the list of alive cells
        cells <- cells[-cell_to_remove]
        cell_bfb_history <- cell_bfb_history[-cell_to_remove]
        cell_last_event <- cell_last_event[-cell_to_remove]
      }
    }
  }

  # Create cell lifetime data
  cell_death_times[is.na(cell_death_times)] <- time

  cell_lifetimes <- data.frame(
    cell_id = cell_ids,
    birth_time = cell_birth_times,
    death_time = cell_death_times,
    lifetime = cell_death_times - cell_birth_times,
    is_alive = cell_death_times == max(cell_death_times),
    parent_id = cell_parents,
    bfb_event = cell_bfb
  )

  # Function to recursively build Newick string
  build_newick <- function(node_id) {
    children <- cell_lifetimes$cell_id[cell_lifetimes$parent_id == node_id & cell_lifetimes$bfb_event]
    if (length(children) == 0) {
      return(as.character(node_id))
    } else {
      child_newick <- sapply(children, build_newick)
      return(paste0("(", paste(child_newick, collapse = ","), ")", node_id))
    }
  }

  # Build the Newick tree starting from the root (parent_id = 0)
  newick_tree <- build_newick(0)
  if (!grepl(";$", newick_tree)) {
    newick_tree <- paste0(newick_tree, ";")
  }

  # Prepare result list
  result <- list(
    final_cells = cells,
    sequence_distribution = table(sapply(cells, function(seq) paste(unlist(seq), collapse = ","))),
    cell_lifetimes = cell_lifetimes,
    final_time = time,
    birth_count = birth_count,
    death_count = death_count,
    newick_tree = newick_tree,
    phylogenetic_info = list(
      bfb_events = bfb_events,
      cell_bfb_history = cell_bfb_history
    )
  )

  result$parameters <- list(
    initial_cells = initial_cells,
    initial_sequence_length = initial_sequence_length,
    birth_rate = birth_rate,
    death_rate = death_rate,
    bfb_prob = bfb_prob,
    max_time = max_time,
    max_cells = max_cells,
    first_n_bfb_cycles = first_n_bfb_cycles,
    first_round_of_bfb = first_round_of_bfb,
    breakpoint_support = breakpoint_support,
    alpha = alpha,
    beta = beta
  )

  return(result)
}
