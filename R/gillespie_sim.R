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
#' @param selection_rate Numeric. Selection advantage for cells with amplified hotspot. Default: 0
#' @param max_time Numeric. Maximum simulation time. Default: 50
#' @param max_cells Numeric. Maximum number of cells allowed before simulation stops. Default: 100
#' @param first_n_bfb_cycles Numeric. Number of initial birth events that will force BFB events. Default: 0
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
#' @details
#' The function implements a Gillespie algorithm to simulate stochastic birth and death
#' processes with possible BFB events. Events occur at rates proportional to the number of cells
#' and depend on time since the last event. The simulation tracks cell lineages, creating
#' a phylogenetic history of BFB events.
#'
#' Cells with amplified hotspots replicate at an increased rate of birth_rate * (1 + selection_rate).
#'
#' BFB (Break-Fusion-Bridge) cycles are a mechanism of genomic instability where
#' chromosome breakage followed by fusion of broken ends creates genomic rearrangements.
#' Each BFB event produces two daughter cells with different rearranged sequences.
#'
#' @export
gillespie_sim <- function(
  initial_cells = 1,
  initial_sequence_length = 100,
  birth_rate = 0.1,
  death_rate = 0.001,
  bfb_prob = 0.01,
  selection_rate = 0,
  max_time = 50,
  max_cells = 100,
  first_n_bfb_cycles = 0,
  first_round_of_bfb = TRUE,
  breakpoint_support = "uniform",
  hotspot = NULL,
  alpha = NULL,
  beta = NULL
  ) {
  # Initialize
  time <- 0
  cells <- list()
  birth_count <- 0  # Track number of birth events
  death_count <- 0

  # Track cell information with character IDs
  cell_ids <- character(0)
  cell_sequences <- list()
  cell_birth_times <- numeric(0)
  cell_death_times <- numeric(0)
  cell_parents <- character(0)
  cell_bfb <- logical(0)
  cell_hotspot_gained <- logical(0)
  cell_is_alive <- logical(0)
  cell_next_event_times <- numeric(0)
  cell_bfb_history <- list()

  # Track BFB event history
  bfb_events <- list()
  next_bfb_id <- 1

  # Create initial sequences and cells
  initial_sequence <- vec2seq(1:initial_sequence_length)
  initial_length <- get_seq_length(initial_sequence)
  initial_hotspot_gained <- is_hotspot_gained(initial_sequence, hotspot)

  if (first_round_of_bfb) {
    for (i in 1:initial_cells) {
      # Create character ID
      cell_id <- paste0("cell_", length(cell_ids)+1)

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
        parent_id = NA,
        parent_cell_id = NA
      )

      # Determine effects and hotspot status
      l_length <- get_seq_length(daughter_seqs$l_seq)
      r_length <- get_seq_length(daughter_seqs$r_seq)

      l_effect <- ifelse(l_length > initial_length, "gain", "loss")
      r_effect <- ifelse(r_length > initial_length, "gain", "loss")

      l_hotspot_gained <- is_hotspot_gained(daughter_seqs$l_seq, hotspot)
      r_hotspot_gained <- is_hotspot_gained(daughter_seqs$r_seq, hotspot)

      # Create cell IDs for daughters
      l_cell_id <- paste0("cell_", length(cell_ids)+1)
      r_cell_id <- paste0("cell_", length(cell_ids)+2)

      # Store cell information
      cell_ids <- c(cell_ids, l_cell_id, r_cell_id)
      cell_sequences[[l_cell_id]] <- daughter_seqs$l_seq
      cell_sequences[[r_cell_id]] <- daughter_seqs$r_seq

      cell_birth_times <- c(cell_birth_times, time, time)
      cell_death_times <- c(cell_death_times, NA, NA)
      cell_parents <- c(cell_parents, NA, NA)
      cell_bfb <- c(cell_bfb, TRUE, TRUE)
      cell_hotspot_gained <- c(cell_hotspot_gained, l_hotspot_gained, r_hotspot_gained)
      cell_is_alive <- c(cell_is_alive, TRUE, TRUE)

      # Initialize next event times with future birth events
      # Modify birth rate based on hotspot gain
      l_birth_rate <- birth_rate * (1 + selection_rate * l_hotspot_gained)
      r_birth_rate <- birth_rate * (1 + selection_rate * r_hotspot_gained)
      cell_next_event_times <- c(cell_next_event_times, time + stats::rexp(1, l_birth_rate), time + stats::rexp(1, r_birth_rate))

      # Track BFB history
      cell_bfb_history[[l_cell_id]] <- list(paste(next_bfb_id, l_effect, sep = "-"))
      cell_bfb_history[[r_cell_id]] <- list(paste(next_bfb_id, r_effect, sep = "-"))

      next_bfb_id <- next_bfb_id + 1
    }
  } else {
    # Initial setup without BFB (similar to previous version but with character IDs)
    for (i in 1:initial_cells) {
      cell_id <- paste0("cell_", length(cell_ids)+1)

      cell_ids <- c(cell_ids, cell_id)
      cell_sequences[[cell_id]] <- initial_sequence
      cell_birth_times <- c(cell_birth_times, time)
      cell_death_times <- c(cell_death_times, NA)
      cell_parents <- c(cell_parents, NA)
      cell_bfb <- c(cell_bfb, FALSE)
      cell_hotspot_gained <- c(cell_hotspot_gained, initial_hotspot_gained)
      cell_is_alive <- c(cell_is_alive, TRUE)

      # Initialize next event time with modified birth rate
      initial_cell_birth_rate <- birth_rate * (1 + selection_rate * initial_hotspot_gained)
      cell_next_event_times <- c(cell_next_event_times,
                                 time + stats::rexp(1, initial_cell_birth_rate))

      cell_bfb_history[[cell_id]] <- list()
    }
  }

  # Hotspot gain events tracking
  hotspot_gain_events <- data.frame(
    time = numeric(0),
    cell_id = character(0)
  )

  # Main simulation loop
  while (time < max_time && sum(cell_is_alive) > 0 && sum(cell_is_alive) < max_cells) {
    # Find the next cell to have an event
    alive_indices <- which(cell_is_alive)
    next_event_cells <- cell_next_event_times[alive_indices]

    # If no more events possible, break
    if (length(next_event_cells) == 0) break

    # Determine next event cell and time
    next_event_idx <- alive_indices[which.min(next_event_cells)]
    time <- cell_next_event_times[next_event_idx]
    current_cell_id <- cell_ids[next_event_idx]

    # Calculate total event rates
    alive_birth_rates <- birth_rate * (1 + selection_rate * cell_hotspot_gained[alive_indices])
    total_birth_prop <- sum(alive_birth_rates)
    total_death_prop <- death_rate * length(alive_indices)
    total_prop <- total_birth_prop + total_death_prop

    # Determine event type (birth or death)
    current_event_type <- if (stats::runif(1) < total_birth_prop / total_prop) "birth" else "death"

    if (current_event_type == "birth") {
      birth_count <- birth_count + 1

      # Determine if BFB occurs
      force_bfb <- birth_count <= first_n_bfb_cycles
      random_bfb <- stats::runif(1) < bfb_prob

      parent_seq <- cell_sequences[[current_cell_id]]
      parent_length <- get_seq_length(parent_seq)
      parent_bfb_history <- cell_bfb_history[[current_cell_id]]
      parent_hotspot_gained <- cell_hotspot_gained[next_event_idx]

      if (force_bfb || random_bfb) {
        # BFB event: simulate left and right sequences
        daughter_seqs <- sim_bfb_left_and_right_sequences(parent_seq, breakpoint_support, alpha, beta)
      } else {
        # Normal duplication: create identical sequences
        daughter_seqs <- list(l_seq = parent_seq, r_seq = parent_seq)
      }

      # Create new cell IDs
      l_cell_id <- paste0("cell_", length(cell_ids) + 1)
      r_cell_id <- paste0("cell_", length(cell_ids) + 2)

      # Mark parent cell as dead
      cell_is_alive[next_event_idx] <- FALSE
      cell_death_times[next_event_idx] <- time

      # Record BFB event if applicable
      if (force_bfb || random_bfb) {
        bfb_events[[next_bfb_id]] <- list(
          id = next_bfb_id,
          time = time,
          parent_id = current_cell_id,
          parent_cell_id = current_cell_id
        )
        next_bfb_id <- next_bfb_id + 1
      }

      # Determine effects and hotspot status
      l_length <- get_seq_length(daughter_seqs$l_seq)
      r_length <- get_seq_length(daughter_seqs$r_seq)

      l_effect <- ifelse(l_length > parent_length, "gain", "loss")
      r_effect <- ifelse(r_length > parent_length, "gain", "loss")

      l_hotspot_gained <- is_hotspot_gained(daughter_seqs$l_seq, hotspot)
      r_hotspot_gained <- is_hotspot_gained(daughter_seqs$r_seq, hotspot)

      # Record hotspot gain events
      if (l_hotspot_gained && !parent_hotspot_gained) {
        hotspot_gain_events <- rbind(hotspot_gain_events, data.frame(time = time, cell_id = l_cell_id))
      }
      if (r_hotspot_gained && !parent_hotspot_gained) {
        hotspot_gain_events <- rbind(hotspot_gain_events, data.frame(time = time, cell_id = r_cell_id))
      }

      # Store new cell information
      cell_ids <- c(cell_ids, l_cell_id, r_cell_id)
      cell_sequences[[l_cell_id]] <- daughter_seqs$l_seq
      cell_sequences[[r_cell_id]] <- daughter_seqs$r_seq

      cell_birth_times <- c(cell_birth_times, time, time)
      cell_death_times <- c(cell_death_times, NA, NA)
      cell_parents <- c(cell_parents, current_cell_id, current_cell_id)
      cell_bfb <- c(cell_bfb, force_bfb || random_bfb, force_bfb || random_bfb)
      cell_hotspot_gained <- c(cell_hotspot_gained, l_hotspot_gained, r_hotspot_gained)
      cell_is_alive <- c(cell_is_alive, TRUE, TRUE)

      # Initialize next event times for new cells
      l_birth_rate <- birth_rate * (1 + selection_rate * l_hotspot_gained)
      r_birth_rate <- birth_rate * (1 + selection_rate * r_hotspot_gained)
      cell_next_event_times <- c(cell_next_event_times, time + stats::rexp(1, l_birth_rate), time + stats::rexp(1, r_birth_rate))

      # Update BFB history
      l_bfb_history <- c(parent_bfb_history, list(paste(next_bfb_id, l_effect, sep = "-")))
      r_bfb_history <- c(parent_bfb_history, list(paste(next_bfb_id, r_effect, sep = "-")))

      cell_bfb_history[[l_cell_id]] <- l_bfb_history
      cell_bfb_history[[r_cell_id]] <- r_bfb_history

    } else {
      # Death event: mark the cell as dead
      death_count <- death_count + 1
      cell_is_alive[next_event_idx] <- FALSE
      cell_death_times[next_event_idx] <- time
    }

    # Regenerate next event times for living cells
    alive_indices <- which(cell_is_alive)
    if (length(alive_indices) > 0) {
      alive_birth_rates <- birth_rate * (1 + selection_rate * cell_hotspot_gained[alive_indices])
      combined_rates <- alive_birth_rates + death_rate
      cell_next_event_times[alive_indices] <- time + stats::rexp(length(alive_indices), combined_rates)
    }
  }

  # Set death time for any remaining alive cells to the final simulation time
  cell_death_times[is.na(cell_death_times)] <- time

  # Create cell lifetime data
  cell_history <- data.frame(
    cell_id = cell_ids,
    birth_time = cell_birth_times,
    death_time = cell_death_times,
    lifetime = cell_death_times - cell_birth_times,
    is_alive = cell_is_alive,
    parent_id = cell_parents,
    bfb_event = cell_bfb,
    hotspot_gained = cell_hotspot_gained
  )

  # Create summary of hotspot status
  hotspot_summary <- data.frame(
    time_points = c(0, sort(unique(c(cell_birth_times, cell_death_times)))),
    cells_with_hotspot = NA,
    total_cells = NA
  )

  # Calculate number of cells with hotspot at each time point
  for (i in 1:nrow(hotspot_summary)) {
    t <- hotspot_summary$time_points[i]
    cells_alive <- cell_birth_times <= t & (cell_death_times > t | is.na(cell_death_times))
    hotspot_summary$total_cells[i] <- sum(cells_alive)
    hotspot_summary$cells_with_hotspot[i] <- sum(cells_alive & cell_hotspot_gained)
  }

  # Get only alive cells for final_cells output
  alive_indices <- which(cell_is_alive)
  final_cells <- lapply(alive_indices, function(i) {
    cell_sequences[[paste0("cell_",i)]]
  })

  cell_history = cell_history %>%
    dplyr::mutate(parent_id = ifelse(is.na(.data$parent_id), "root", .data$parent_id))

  cell_history = dplyr::bind_rows(
    dplyr::tibble(
      cell_id="root",
      birth_time = -1,
      death_time=-1,
      lifetime = 0,
      is_alive=FALSE,
      parent_id=NA,
      bfb_event=FALSE,
      hotspot_gained=FALSE),
    cell_history
  )

  # Prepare result list
  result <- list(
    cells = final_cells,
    #sequence_distribution = table(sapply(final_cells, function(seq) paste(unlist(seq), collapse = ","))),
    cell_history = cell_history#,
    #final_time = time,
    #birth_count = birth_count,
    #death_count = death_count,
    # phylogenetic_info = list(
    #   bfb_events = bfb_events,
    #   cell_bfb_history = cell_bfb_history
    # ),
    # hotspot_info = list(
    #   hotspot_summary = hotspot_summary,
    #   hotspot_gain_events = hotspot_gain_events,
    #   hotspot_positions = hotspot
    # )
  )

  result$input_parameters <- list(
    initial_cells = initial_cells,
    initial_sequence_length = initial_sequence_length,
    birth_rate = birth_rate,
    death_rate = death_rate,
    bfb_prob = bfb_prob,
    selection_rate = selection_rate,
    max_time = max_time,
    max_cells = max_cells,
    first_n_bfb_cycles = first_n_bfb_cycles,
    first_round_of_bfb = first_round_of_bfb,
    breakpoint_support = breakpoint_support,
    hotspot = hotspot,
    alpha = alpha,
    beta = beta
  )

  return(result)
}
