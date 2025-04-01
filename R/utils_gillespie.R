#' Determine the next event to occur
#'
#' @param state The simulation state
#' @param birth_rate Base birth rate
#' @param death_rate Base death rate
#' @param positive_selection_rate Selection advantage for hotspot cells
#' @param negative_selection_rate Selection disadvantage for non hotspot cells
#'
#' @return List with time, cell ID, and event type of next event
get_next_event <- function(state, birth_rate, death_rate, positive_selection_rate, negative_selection_rate) {
  # Find the next cell to have an event
  alive_indices <- which(state$cell_is_alive)

  # If no more events possible, return NULL
  if (length(alive_indices) == 0) {
    return(NULL)
  }

  # Determine next event cell and time
  next_event_times <- state$cell_next_event_times[alive_indices]
  min_event_idx <- which.min(next_event_times)
  next_event_idx <- alive_indices[min_event_idx]
  next_event_time <- next_event_times[min_event_idx]
  current_cell_id <- state$cell_ids[next_event_idx]

  # Calculate modified birth and death rates for this specific cell
  cell_hotspot_gained <- state$cell_hotspot_gained[next_event_idx]
  cell_birth_rate <- birth_rate * (1 + positive_selection_rate * cell_hotspot_gained)
  cell_death_rate <- death_rate * (1 - negative_selection_rate * (!cell_hotspot_gained))

  # Determine event type (birth or death) based on relative rates
  event_probability <- cell_birth_rate / (cell_birth_rate + cell_death_rate)
  event_type <- if (stats::runif(1) < event_probability) "birth" else "death"

  return(list(
    time = next_event_time,
    cell_id = current_cell_id,
    event_type = event_type,
    cell_idx = next_event_idx
  ))
}

#' Process a birth event
#'
#' @param state The simulation state
#' @param current_cell_id ID of the cell undergoing birth
#' @param birth_rate Base birth rate
#' @param death_rate Base death rate
#' @param bfb_prob Probability of BFB event
#' @param positive_selection_rate Selection advantage for hotspot cells
#' @param negative_selection_rate Selection disadvantage for non hotspot cells
#' @param breakpoint_support Distribution for breakpoint selection
#' @param first_n_bfb_cycles Number of initial birth events with forced BFB
#' @param hotspot Positions considered hotspots
#' @param alpha Beta distribution parameter
#' @param beta Beta distribution parameter
#'
#' @return Updated simulation state
process_birth_event <- function(
    state,
    current_cell_id,
    birth_rate,
    death_rate,
    bfb_prob,
    positive_selection_rate,
    negative_selection_rate,
    breakpoint_support,
    first_n_bfb_cycles,
    hotspot,
    alpha,
    beta
) {
  # Increment birth count
  state$birth_count <- state$birth_count + 1

  # Find index of current cell
  cell_idx <- which(state$cell_ids == current_cell_id)

  # Determine if BFB occurs
  force_bfb <- state$birth_count <= first_n_bfb_cycles
  random_bfb <- stats::runif(1) < bfb_prob
  bfb_occurs <- force_bfb || random_bfb

  # Get parent cell information
  parent_seq <- state$cell_sequences[[current_cell_id]]
  parent_length <- get_seq_length(parent_seq)
  parent_bfb_history <- state$cell_bfb_history[[current_cell_id]]
  parent_hotspot_gained <- state$cell_hotspot_gained[cell_idx]

  # Create daughter sequences
  if (bfb_occurs) {
    # BFB event: simulate left and right sequences
    daughter_seqs <- sim_bfb_left_and_right_sequences(parent_seq, breakpoint_support, alpha, beta)
  } else {
    # Normal duplication: create identical sequences
    daughter_seqs <- list(l_seq = parent_seq, r_seq = parent_seq)
  }

  # Create new cell IDs
  l_cell_id <- paste0("cell_", length(state$cell_ids) + 1)
  r_cell_id <- paste0("cell_", length(state$cell_ids) + 2)

  # Mark parent cell as dead
  state$cell_is_alive[cell_idx] <- FALSE
  state$cell_death_times[cell_idx] <- state$time

  # Record BFB event if applicable
  if (bfb_occurs) {
    state$bfb_events[[state$next_bfb_id]] <- list(
      id = state$next_bfb_id,
      time = state$time,
      parent_id = current_cell_id,
      parent_cell_id = current_cell_id
    )
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
    state$hotspot_gain_events <- rbind(
      state$hotspot_gain_events,
      data.frame(time = state$time, cell_id = l_cell_id)
    )
  }
  if (r_hotspot_gained && !parent_hotspot_gained) {
    state$hotspot_gain_events <- rbind(
      state$hotspot_gain_events,
      data.frame(time = state$time, cell_id = r_cell_id)
    )
  }

  # Store new cell information
  state$cell_ids <- c(state$cell_ids, l_cell_id, r_cell_id)
  state$cell_sequences[[l_cell_id]] <- daughter_seqs$l_seq
  state$cell_sequences[[r_cell_id]] <- daughter_seqs$r_seq

  state$cell_birth_times <- c(state$cell_birth_times, state$time, state$time)
  state$cell_death_times <- c(state$cell_death_times, NA, NA)
  state$cell_parents <- c(state$cell_parents, current_cell_id, current_cell_id)
  state$cell_bfb <- c(state$cell_bfb, bfb_occurs, bfb_occurs)
  state$cell_hotspot_gained <- c(state$cell_hotspot_gained, l_hotspot_gained, r_hotspot_gained)
  state$cell_is_alive <- c(state$cell_is_alive, TRUE, TRUE)

  # Calculate modified birth and death rates for new cells
  l_birth_rate <- birth_rate * (1 + positive_selection_rate * l_hotspot_gained)
  r_birth_rate <- birth_rate * (1 + positive_selection_rate * r_hotspot_gained)

  l_death_rate <- death_rate * (1 - negative_selection_rate * (!l_hotspot_gained))
  r_death_rate <- death_rate * (1 - negative_selection_rate * (!r_hotspot_gained))

  # Calculate next event times for new cells based on combined rates
  l_combined_rate <- l_birth_rate + l_death_rate
  r_combined_rate <- r_birth_rate + r_death_rate

  # Initialize next event times ONLY for the new cells
  state$cell_next_event_times <- c(
    state$cell_next_event_times,
    state$time + stats::rexp(1, l_combined_rate),
    state$time + stats::rexp(1, r_combined_rate)
  )

  # Update BFB history
  if (bfb_occurs) {
    l_bfb_history <- c(parent_bfb_history, list(paste(state$next_bfb_id, l_effect, sep = "-")))
    r_bfb_history <- c(parent_bfb_history, list(paste(state$next_bfb_id, r_effect, sep = "-")))

    state$cell_bfb_history[[l_cell_id]] <- l_bfb_history
    state$cell_bfb_history[[r_cell_id]] <- r_bfb_history

    state$next_bfb_id <- state$next_bfb_id + 1
  } else {
    # For normal division, copy parent history
    state$cell_bfb_history[[l_cell_id]] <- parent_bfb_history
    state$cell_bfb_history[[r_cell_id]] <- parent_bfb_history
  }

  return(state)
}

#' Process a death event
#'
#' @param state The simulation state
#' @param current_cell_id ID of the cell undergoing death
#'
#' @return Updated simulation state
process_death_event <- function(state, current_cell_id) {
  # Find index of current cell
  cell_idx <- which(state$cell_ids == current_cell_id)

  # Increment death count
  state$death_count <- state$death_count + 1

  # Mark the cell as dead
  state$cell_is_alive[cell_idx] <- FALSE
  state$cell_death_times[cell_idx] <- state$time

  return(state)
}

#' Initialize simulation state
#'
#' @param initial_cells Number of cells to start with
#' @param initial_sequence_length Length of initial genomic sequence
#' @param birth_rate Base birth rate
#' @param death_rate Base death rate
#' @param positive_selection_rate Selection advantage for hotspot cells
#' @param negative_selection_rate Selection disadvantage for non hotspot cells
#' @param first_round_of_bfb Whether to apply BFB to initial cells
#' @param breakpoint_support Distribution for breakpoint selection
#' @param hotspot Positions considered hotspots
#' @param alpha Beta distribution parameter
#' @param beta Beta distribution parameter
#'
#' @return A list containing the initialized simulation state
initialize_simulation <- function(
    initial_cells,
    initial_sequence_length,
    birth_rate,
    death_rate,
    positive_selection_rate,
    negative_selection_rate,
    first_round_of_bfb,
    breakpoint_support,
    hotspot,
    alpha,
    beta
) {
  # Initialize state variables
  state <- list(
    time = 0,
    birth_count = 0,
    death_count = 0,
    cell_ids = character(0),
    cell_sequences = list(),
    cell_birth_times = numeric(0),
    cell_death_times = numeric(0),
    cell_parents = character(0),
    cell_bfb = logical(0),
    cell_hotspot_gained = logical(0),
    cell_is_alive = logical(0),
    cell_next_event_times = numeric(0),
    cell_bfb_history = list(),
    bfb_events = list(),
    next_bfb_id = 1,
    hotspot_gain_events = data.frame(time = numeric(0), cell_id = character(0))
  )

  # Create initial sequence
  initial_sequence <- vec2seq(1:initial_sequence_length)
  initial_length <- get_seq_length(initial_sequence)
  initial_hotspot_gained <- is_hotspot_gained(initial_sequence, hotspot)

  if (first_round_of_bfb) {
    state <- initialize_with_bfb(
      state,
      initial_cells,
      initial_sequence,
      initial_length,
      birth_rate,
      death_rate,
      positive_selection_rate,
      negative_selection_rate,
      breakpoint_support,
      hotspot,
      alpha,
      beta
    )
  } else {
    state <- initialize_without_bfb(
      state,
      initial_cells,
      initial_sequence,
      initial_hotspot_gained,
      birth_rate,
      death_rate,
      positive_selection_rate,
      negative_selection_rate
    )
  }

  return(state)
}

#' Initialize simulation with BFB events for initial cells
#'
#' @param state The simulation state
#' @param initial_cells Number of cells to start with
#' @param initial_sequence Initial genomic sequence
#' @param initial_length Length of initial sequence
#' @param birth_rate Base birth rate
#' @param death_rate Base death rate
#' @param positive_selection_rate Selection advantage for hotspot cells
#' @param negative_selection_rate Selection disadvantage for non hotspot cells
#' @param breakpoint_support Distribution for breakpoint selection
#' @param hotspot Positions considered hotspots
#' @param alpha Beta distribution parameter
#' @param beta Beta distribution parameter
#'
#' @return Updated simulation state
initialize_with_bfb <- function(
    state,
    initial_cells,
    initial_sequence,
    initial_length,
    birth_rate,
    death_rate,
    positive_selection_rate,
    negative_selection_rate,
    breakpoint_support,
    hotspot,
    alpha,
    beta
) {
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
    state$bfb_events[[state$next_bfb_id]] <- list(
      id = state$next_bfb_id,
      time = state$time,
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
    l_cell_id <- paste0("cell_", length(state$cell_ids) + 1)
    r_cell_id <- paste0("cell_", length(state$cell_ids) + 2)

    # Store cell information
    state$cell_ids <- c(state$cell_ids, l_cell_id, r_cell_id)
    state$cell_sequences[[l_cell_id]] <- daughter_seqs$l_seq
    state$cell_sequences[[r_cell_id]] <- daughter_seqs$r_seq

    state$cell_birth_times <- c(state$cell_birth_times, state$time, state$time)
    state$cell_death_times <- c(state$cell_death_times, NA, NA)
    state$cell_parents <- c(state$cell_parents, NA, NA)
    state$cell_bfb <- c(state$cell_bfb, TRUE, TRUE)
    state$cell_hotspot_gained <- c(state$cell_hotspot_gained, l_hotspot_gained, r_hotspot_gained)
    state$cell_is_alive <- c(state$cell_is_alive, TRUE, TRUE)

    # Calculate modified birth and death rates
    l_birth_rate <- birth_rate * (1 + positive_selection_rate * l_hotspot_gained)
    r_birth_rate <- birth_rate * (1 + positive_selection_rate * r_hotspot_gained)

    l_death_rate <- death_rate * (1 - negative_selection_rate * (!l_hotspot_gained))
    r_death_rate <- death_rate * (1 - negative_selection_rate * (!r_hotspot_gained))

    # Calculate combined rates
    l_combined_rate <- l_birth_rate + l_death_rate
    r_combined_rate <- r_birth_rate + r_death_rate

    # Initialize next event times with future events based on combined rates
    state$cell_next_event_times <- c(
      state$cell_next_event_times,
      state$time + stats::rexp(1, l_combined_rate),
      state$time + stats::rexp(1, r_combined_rate)
    )

    # Track BFB history
    state$cell_bfb_history[[l_cell_id]] <- list(paste(state$next_bfb_id, l_effect, sep = "-"))
    state$cell_bfb_history[[r_cell_id]] <- list(paste(state$next_bfb_id, r_effect, sep = "-"))

    state$next_bfb_id <- state$next_bfb_id + 1
  }

  return(state)
}

#' Initialize simulation without BFB for initial cells
#'
#' @param state The simulation state
#' @param initial_cells Number of cells to start with
#' @param initial_sequence Initial genomic sequence
#' @param initial_hotspot_gained Whether initial cells have hotspot gain
#' @param birth_rate Base birth rate
#' @param death_rate Base death rate
#' @param positive_selection_rate Selection advantage for hotspot cells
#' @param negative_selection_rate Selection disadvantage for non hotspot cells
#'
#' @return Updated simulation state
initialize_without_bfb <- function(
    state,
    initial_cells,
    initial_sequence,
    initial_hotspot_gained,
    birth_rate,
    death_rate,
    positive_selection_rate,
    negative_selection_rate
) {
  for (i in 1:initial_cells) {
    cell_id <- paste0("cell_", length(state$cell_ids) + 1)

    state$cell_ids <- c(state$cell_ids, cell_id)
    state$cell_sequences[[cell_id]] <- initial_sequence
    state$cell_birth_times <- c(state$cell_birth_times, state$time)
    state$cell_death_times <- c(state$cell_death_times, NA)
    state$cell_parents <- c(state$cell_parents, NA)
    state$cell_bfb <- c(state$cell_bfb, FALSE)
    state$cell_hotspot_gained <- c(state$cell_hotspot_gained, initial_hotspot_gained)
    state$cell_is_alive <- c(state$cell_is_alive, TRUE)

    # Calculate modified birth and death rates
    cell_birth_rate <- birth_rate * (1 + positive_selection_rate * initial_hotspot_gained)
    cell_death_rate <- death_rate * (1 - negative_selection_rate * (!initial_hotspot_gained))

    # Calculate combined rate
    combined_rate <- cell_birth_rate + cell_death_rate

    # Initialize next event time based on combined rate
    state$cell_next_event_times <- c(
      state$cell_next_event_times,
      state$time + stats::rexp(1, combined_rate)
    )

    state$cell_bfb_history[[cell_id]] <- list()
  }

  return(state)
}

#' Check if simulation should continue
#'
#' @param state The simulation state
#' @param max_time Maximum simulation time
#' @param max_cells Maximum number of cells allowed
#'
#' @return Logical indicating whether simulation should continue
continue_simulation <- function(state, max_time, max_cells) {
  alive_count <- sum(state$cell_is_alive)
  return(state$time < max_time && alive_count > 0 && alive_count < max_cells)
}

#' Prepare final results from simulation
#'
#' @param state The simulation state
#'
#' @return A list containing simulation results
prepare_results <- function(state) {
  # Set death time for any remaining alive cells to the final simulation time
  state$cell_death_times[is.na(state$cell_death_times)] <- state$time

  # Create cell lifetime data
  cell_history <- data.frame(
    cell_id = state$cell_ids,
    birth_time = state$cell_birth_times,
    death_time = state$cell_death_times,
    lifetime = state$cell_death_times - state$cell_birth_times,
    is_alive = state$cell_is_alive,
    parent_id = state$cell_parents,
    bfb_event = state$cell_bfb,
    hotspot_gained = state$cell_hotspot_gained
  )

  # Get only alive cells for final_cells output
  alive_indices <- which(state$cell_is_alive)
  final_cells <- lapply(state$cell_ids[alive_indices], function(id) {
    state$cell_sequences[[id]]
  })

  # Fix parent IDs for display
  cell_history <- cell_history %>%
    dplyr::mutate(parent_id = ifelse(is.na(.data$parent_id), "root", .data$parent_id))

  # Add root node
  cell_history <- dplyr::bind_rows(
    dplyr::tibble(
      cell_id = "root",
      birth_time = -1,
      death_time = -1,
      lifetime = 0,
      is_alive = FALSE,
      parent_id = NA,
      bfb_event = FALSE,
      hotspot_gained = FALSE
    ),
    cell_history
  )

  # Prepare and return results
  result <- list(
    cells = final_cells,
    cell_history = cell_history
  )

  return(result)
}
