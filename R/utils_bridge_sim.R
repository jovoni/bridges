#' Helper functions modified for diploid chromosome modeling

#' Initialize simulation state for diploid chromosomes
#'
#' @param input_parameters List of parameters given as input to bridge_sim function.
#'
#' @return A list containing the initialized simulation state
initialize_simulation <- function(input_parameters) {
  # Pre-allocate history arrays. Upper bound: initial cells create 2 cells each,
  # then each birth creates 2 daughters. Total entries ≤ 2*initial + 2*max_cells births.
  max_history <- 6L * as.integer(input_parameters$max_cells) +
                 4L * as.integer(input_parameters$initial_cells) + 20L

  state <- list(
    time = 0,
    # Active arrays — contain ONLY currently alive cells (no dead cells kept).
    # This keeps which.min and length fast throughout the simulation.
    cell_ids               = character(0),
    cell_sequences         = list(),
    cell_next_event_times  = numeric(0),
    hotspot_status         = logical(0),
    # Pre-allocated history vectors with integer counter h_n.
    # Avoids O(n) dplyr::bind_rows growth in the hot loop.
    h_cell_id    = character(max_history),
    h_parent_id  = character(max_history),
    h_bfb_event  = logical(max_history),
    h_wgd_event  = logical(max_history),
    h_cn_event   = character(max_history),
    h_chr_allele = character(max_history),
    h_n          = 0L
  )
  state$input_parameters <- input_parameters
  state$next_cell_id <- 1L

  initial_sequences <- create_initial_chromosome_sequences(input_parameters)

  if (state$input_parameters$first_round_of_bfb) {
    state <- initialize_with_bfb(state, initial_sequences)
  } else {
    state <- initialize_without_bfb(state, initial_sequences)
  }

  return(state)
}

#' Create initial chromosome sequences for all alleles
#'
#' @param input_parameters List of input parameters
#'
#' @return Named list of initial sequences for each chromosome allele
create_initial_chromosome_sequences <- function(input_parameters) {
  initial_sequences <- list()

  for (chr_allele in input_parameters$chr_alleles) {
    # Extract chromosome name (remove _A or _B suffix)
    chr_name <- sub(":[AB]$", "", chr_allele)
    seq_length <- input_parameters$chr_seq_lengths[chr_name]

    # Create sequence as a vector from 1 to sequence length
    initial_sequences[[chr_allele]] <- vec2seq(1:seq_length)
  }

  return(initial_sequences)
}

#' Initialize simulation with BFB events for initial cells (diploid version)
#'
#' @param state The simulation state
#' @param initial_sequences List of initial chromosome sequences
#'
#' @return Updated simulation state
initialize_with_bfb <- function(state, initial_sequences) {
  for (i in 1:state$input_parameters$initial_cells) {
    # Get the BFB allele and generate daughter sequences
    bfb_allele <- state$input_parameters$bfb_allele
    daughter_seqs <- sim_bfb_left_and_right_sequences(
      sequence = initial_sequences[[bfb_allele]],
      support = state$input_parameters$breakpoint_support,
      alpha = state$input_parameters$alpha,
      beta = state$input_parameters$beta,
      custom_breakpoints = state$input_parameters$custom_breakpoints
    )

    # Left child
    left_cell <- initial_sequences
    left_cell[[bfb_allele]] <- daughter_seqs$l_seq
    left_cell_id <- paste0("cell_", state$next_cell_id)
    state$next_cell_id <- state$next_cell_id + 1L

    # Right child
    right_cell <- initial_sequences
    right_cell[[bfb_allele]] <- daughter_seqs$r_seq
    right_cell_id <- paste0("cell_", state$next_cell_id)
    state$next_cell_id <- state$next_cell_id + 1L

    if (!is.null(state$input_parameters$hotspot)) {
      hotspot <- state$input_parameters$hotspot
      left_hotspot <- is_hotspot_gained(left_cell[[hotspot$chr]], hotspot = hotspot$pos)
      right_hotspot <- is_hotspot_gained(right_cell[[hotspot$chr]], hotspot = hotspot$pos)
    } else {
      left_hotspot <- right_hotspot <- FALSE
    }

    l_birth_rate <- state$input_parameters$birth_rate *
      (1 + state$input_parameters$positive_selection_rate * left_hotspot)
    l_death_rate <- state$input_parameters$death_rate *
      (1 + state$input_parameters$negative_selection_rate * left_hotspot)

    r_birth_rate <- state$input_parameters$birth_rate *
      (1 + state$input_parameters$positive_selection_rate * right_hotspot)
    r_death_rate <- state$input_parameters$death_rate *
      (1 + state$input_parameters$negative_selection_rate * right_hotspot)

    l_combined_rate <- l_birth_rate + l_death_rate
    r_combined_rate <- r_birth_rate + r_death_rate

    # Update active arrays (only alive cells)
    state$cell_ids               <- c(state$cell_ids, left_cell_id, right_cell_id)
    state$cell_sequences[[left_cell_id]]  <- left_cell
    state$cell_sequences[[right_cell_id]] <- right_cell
    state$hotspot_status         <- c(state$hotspot_status, left_hotspot, right_hotspot)
    state$cell_next_event_times  <- c(
      state$cell_next_event_times,
      state$time + stats::rexp(1, l_combined_rate),
      state$time + stats::rexp(1, r_combined_rate)
    )

    # Append to pre-allocated history
    i_l <- state$h_n + 1L
    i_r <- state$h_n + 2L
    state$h_n <- i_r

    state$h_cell_id[i_l]    <- left_cell_id
    state$h_parent_id[i_l]  <- "root"
    state$h_bfb_event[i_l]  <- TRUE
    state$h_wgd_event[i_l]  <- FALSE
    state$h_cn_event[i_l]   <- "bfb"
    state$h_chr_allele[i_l] <- bfb_allele

    state$h_cell_id[i_r]    <- right_cell_id
    state$h_parent_id[i_r]  <- "root"
    state$h_bfb_event[i_r]  <- TRUE
    state$h_wgd_event[i_r]  <- FALSE
    state$h_cn_event[i_r]   <- "bfb"
    state$h_chr_allele[i_r] <- bfb_allele
  }
  return(state)
}


#' Initialize simulation without BFB for initial cells (diploid version)
#'
#' @param state The simulation state
#' @param initial_sequences List of initial chromosome sequences
#'
#' @return Updated simulation state
initialize_without_bfb <- function(state, initial_sequences) {
  for (i in 1:state$input_parameters$initial_cells) {
    cell_id <- paste0("cell_", state$next_cell_id)
    state$next_cell_id <- state$next_cell_id + 1L

    # Hotspot logic
    if (!is.null(state$input_parameters$hotspot)) {
      hotspot <- state$input_parameters$hotspot
      hotspot_gained <- is_hotspot_gained(
        initial_sequences[[hotspot$chr]],
        hotspot = hotspot$pos
      )
    } else {
      hotspot_gained <- FALSE
    }

    # Birth and death rate with selection
    birth_rate <- state$input_parameters$birth_rate *
      (1 + state$input_parameters$positive_selection_rate * hotspot_gained)
    death_rate <- state$input_parameters$death_rate *
      (1 + state$input_parameters$negative_selection_rate * hotspot_gained)

    combined_rate  <- birth_rate + death_rate
    next_event_time <- state$time + stats::rexp(1, combined_rate)

    # Update active arrays
    state$cell_ids                    <- c(state$cell_ids, cell_id)
    state$cell_sequences[[cell_id]]   <- initial_sequences
    state$cell_next_event_times       <- c(state$cell_next_event_times, next_event_time)
    state$hotspot_status              <- c(state$hotspot_status, hotspot_gained)

    # Append to pre-allocated history
    idx <- state$h_n + 1L
    state$h_n <- idx

    state$h_cell_id[idx]    <- cell_id
    state$h_parent_id[idx]  <- "root"
    state$h_bfb_event[idx]  <- FALSE
    state$h_wgd_event[idx]  <- FALSE
    state$h_cn_event[idx]   <- "none"
    state$h_chr_allele[idx] <- NA_character_
  }

  return(state)
}

#' Process a birth event (diploid version)
#'
#' This function handles cell division by creating two daughter cells from a parent cell.
#' During division, various genomic events (amplifications, deletions, BFB) can occur
#' based on the specified rates and lambda parameter.
#'
#' @param state The simulation state containing cell information,
#'  sequences, and parameters
#' @param current_cell_id ID of the cell undergoing birth/division
#' @param cell_idx Index of the cell in the alive arrays (from get_next_event)
#' @param lambda Rate parameter for Poisson distribution used to sample
#' the number of genomic events per daughter cell
#' @param rate Rate parameter used in amplification/deletion simulations.
#'  Length of event is sample from exponential distribution with parameter 1 / rate.
#'
#' @return Updated simulation state with new daughter cells and updated history
process_birth_event <- function(state, current_cell_id, cell_idx, lambda, rate) {
  # Sample number of events for each daughter cell
  n_left  <- stats::rpois(1, lambda)
  n_right <- stats::rpois(1, lambda)

  # Get parent cell information
  parent_sequences <- state$cell_sequences[[current_cell_id]]

  # Initialize daughter sequences with parent sequences
  l_sequences <- parent_sequences
  r_sequences <- parent_sequences

  # Track if special events have occurred and which events happened
  bfb_occurred <- FALSE
  wgd_occurred <- FALSE
  l_events      <- character(0)
  r_events      <- character(0)
  l_chr_alleles <- character(0)
  r_chr_alleles <- character(0)

  # Helper function to apply a single event to sequences
  apply_single_event <- function(
    sequences,
    event_name,
    selected_chr_allele = NULL
  ) {
    if (event_name == "normal") {
      return(list(sequences = sequences, chr_allele = selected_chr_allele))
    } else if (event_name == "amp") {
      if (is.null(selected_chr_allele)) {
        selected_chr_allele <- sample(names(sequences), 1)
      }
      sequences[[selected_chr_allele]] <- sim_amp_del(
        sequences[[selected_chr_allele]],
        operation = "dup",
        rate = rate
      )
    } else if (event_name == "del") {
      if (is.null(selected_chr_allele)) {
        selected_chr_allele <- sample(names(sequences), 1)
      }
      sequences[[selected_chr_allele]] <- sim_amp_del(
        sequences[[selected_chr_allele]],
        operation = "del",
        rate = rate
      )
    } else if (event_name == "wgd") {
      # WGD affects all chromosomes
      for (chr_name in names(sequences)) {
        sequences[[chr_name]] <- sim_wgd(sequences[[chr_name]])
      }
      return(list(sequences = sequences, chr_allele = "all"))
    } else if (event_name == "bfb") {
      selected_chr_allele <- state$input_parameters$bfb_allele
      bfb_result <- sim_bfb_left_and_right_sequences(
        sequences[[selected_chr_allele]],
        state$input_parameters$breakpoint_support,
        state$input_parameters$alpha,
        state$input_parameters$beta,
        state$input_parameters$custom_breakpoints
      )
      # For BFB, we return both left and right sequences
      return(list(
        l_seq = bfb_result$l_seq,
        r_seq = bfb_result$r_seq,
        chr_allele = selected_chr_allele
      ))
    }

    return(list(sequences = sequences, chr_allele = selected_chr_allele))
  }

  # Check if WGD should occur first (randomly and if available)
  wgd_will_occur <- FALSE
  if (state$input_parameters$wgd_available > 0) {
    wgd_probability <- state$input_parameters$wgd_probability
    wgd_will_occur  <- stats::runif(1) < wgd_probability
  }

  # Check if we need to do special events first (BFB or WGD)
  total_events <- n_left + n_right
  if (total_events > 0) {
    # Sample all events that will occur (excluding WGD from rates)
    event_rates <- state$input_parameters$rates
    if ("wgd" %in% names(event_rates)) {
      event_rates <- event_rates[names(event_rates) != "wgd"]
    }

    all_event_names <- sample(
      names(event_rates),
      size = total_events,
      prob = unlist(event_rates),
      replace = TRUE
    )

    # Check for BFB events only (WGD is handled separately)
    bfb_indices <- which(all_event_names == "bfb")

    # Handle special events (only one can occur per division)
    special_event_occurred <- FALSE

    if (wgd_will_occur && length(bfb_indices) > 0) {
      # If both BFB and WGD are scheduled, randomly choose one
      chosen_event <- sample(c("bfb", "wgd"), 1)
      if (chosen_event == "wgd") {
        bfb_indices <- integer(0)  # Cancel BFB
      } else {
        wgd_will_occur <- FALSE  # Cancel WGD
      }
    }

    if (length(bfb_indices) > 0) {
      # Apply BFB (affects both daughter cells differently)
      bfb_occurred <- TRUE
      special_event_occurred <- TRUE

      bfb_result <- apply_single_event(parent_sequences, "bfb")
      l_sequences[[bfb_result$chr_allele]] <- bfb_result$l_seq
      r_sequences[[bfb_result$chr_allele]] <- bfb_result$r_seq

      # Record BFB event for both cells
      l_events      <- c(l_events, "bfb")
      r_events      <- c(r_events, "bfb")
      l_chr_alleles <- c(l_chr_alleles, bfb_result$chr_allele)
      r_chr_alleles <- c(r_chr_alleles, bfb_result$chr_allele)

      # Remove BFB events from the list
      remaining_events <- all_event_names[-bfb_indices]

    } else if (wgd_will_occur) {
      # Apply WGD (affects both daughter cells identically)
      wgd_occurred <- TRUE
      special_event_occurred <- TRUE

      wgd_result <- apply_single_event(parent_sequences, "wgd")
      l_sequences <- wgd_result$sequences
      r_sequences <- wgd_result$sequences

      # Record WGD event for both cells
      l_events      <- c(l_events, "wgd")
      r_events      <- c(r_events, "wgd")
      l_chr_alleles <- c(l_chr_alleles, wgd_result$chr_allele)
      r_chr_alleles <- c(r_chr_alleles, wgd_result$chr_allele)

      # Decrement available WGD count
      state$input_parameters$wgd_available <- state$input_parameters$wgd_available - 1

      # All originally sampled events remain
      remaining_events <- all_event_names

    } else {
      remaining_events <- all_event_names
    }

    # Handle remaining events after special events
    if (special_event_occurred) {
      total_remaining <- length(remaining_events)

      # Adjust n_left and n_right to account for special event
      n_left  <- max(0, n_left  - 1)
      n_right <- max(0, n_right - 1)

      # Redistribute remaining events
      if (total_remaining > 0 && (n_left + n_right) > 0) {
        # Randomly assign remaining events to left and right
        left_additional  <- min(n_left,  total_remaining)
        right_additional <- min(n_right, total_remaining - left_additional)

        if (left_additional > 0) {
          left_events_idx        <- sample(total_remaining, left_additional)
          left_additional_events <- remaining_events[left_events_idx]
          remaining_events       <- remaining_events[-left_events_idx]
          total_remaining        <- total_remaining - left_additional
        } else {
          left_additional_events <- character(0)
        }

        if (right_additional > 0 && total_remaining > 0) {
          right_additional_events <- remaining_events[
            1:min(right_additional, total_remaining)
          ]
        } else {
          right_additional_events <- character(0)
        }

        # Apply additional events to left cell
        for (event in left_additional_events) {
          result      <- apply_single_event(l_sequences, event)
          l_sequences <- result$sequences
          l_events      <- c(l_events, event)
          l_chr_alleles <- c(l_chr_alleles, result$chr_allele)
        }

        # Apply additional events to right cell
        for (event in right_additional_events) {
          result      <- apply_single_event(r_sequences, event)
          r_sequences <- result$sequences
          r_events      <- c(r_events, event)
          r_chr_alleles <- c(r_chr_alleles, result$chr_allele)
        }
      }
    } else {
      # No special events - distribute events normally
      left_events_count  <- min(n_left,  total_events)
      right_events_count <- min(n_right, total_events - left_events_count)

      # Apply events to left cell
      if (left_events_count > 0) {
        left_event_names <- all_event_names[1:left_events_count]
        for (event in left_event_names) {
          result      <- apply_single_event(l_sequences, event)
          l_sequences <- result$sequences
          l_events      <- c(l_events, event)
          l_chr_alleles <- c(l_chr_alleles, result$chr_allele)
        }
      }

      # Apply events to right cell
      if (right_events_count > 0) {
        right_event_names <- all_event_names[
          (left_events_count + 1):(left_events_count + right_events_count)
        ]
        for (event in right_event_names) {
          result      <- apply_single_event(r_sequences, event)
          r_sequences <- result$sequences
          r_events      <- c(r_events, event)
          r_chr_alleles <- c(r_chr_alleles, result$chr_allele)
        }
      }
    }
    # Handle WGD separately if no other events are sampled but WGD should occur
  } else if (wgd_will_occur) {
    # Apply WGD even when no other events are sampled
    wgd_occurred <- TRUE

    wgd_result <- apply_single_event(parent_sequences, "wgd")
    l_sequences <- wgd_result$sequences
    r_sequences <- wgd_result$sequences

    # Record WGD event for both cells
    l_events      <- c(l_events, "wgd")
    r_events      <- c(r_events, "wgd")
    l_chr_alleles <- c(l_chr_alleles, wgd_result$chr_allele)
    r_chr_alleles <- c(r_chr_alleles, wgd_result$chr_allele)

    # Decrement available WGD count
    state$input_parameters$wgd_available <- state$input_parameters$wgd_available - 1
  }

  # Create new cell IDs
  l_cell_id <- paste0("cell_", state$next_cell_id)
  state$next_cell_id <- state$next_cell_id + 1L
  r_cell_id <- paste0("cell_", state$next_cell_id)
  state$next_cell_id <- state$next_cell_id + 1L

  # Hotspot status for daughters
  hotspot   <- state$input_parameters$hotspot
  l_hotspot <- if (!is.null(hotspot))
    is_hotspot_gained(l_sequences[[hotspot$chr]], hotspot = hotspot$pos) else FALSE
  r_hotspot <- if (!is.null(hotspot))
    is_hotspot_gained(r_sequences[[hotspot$chr]], hotspot = hotspot$pos) else FALSE

  # Selection-adjusted rates
  l_birth_rate <- state$input_parameters$birth_rate *
    (1 + state$input_parameters$positive_selection_rate * l_hotspot)
  l_death_rate <- state$input_parameters$death_rate *
    (1 + state$input_parameters$negative_selection_rate * l_hotspot)

  r_birth_rate <- state$input_parameters$birth_rate *
    (1 + state$input_parameters$positive_selection_rate * r_hotspot)
  r_death_rate <- state$input_parameters$death_rate *
    (1 + state$input_parameters$negative_selection_rate * r_hotspot)

  # ---- Update active arrays (only alive cells) ----
  # Remove parent by swap-with-last then shrink — avoids shifting the entire vector.
  n_alive <- length(state$cell_ids)
  if (cell_idx < n_alive) {
    state$cell_ids[cell_idx]              <- state$cell_ids[n_alive]
    state$cell_next_event_times[cell_idx] <- state$cell_next_event_times[n_alive]
    state$hotspot_status[cell_idx]        <- state$hotspot_status[n_alive]
  }
  state$cell_ids              <- state$cell_ids[-n_alive]
  state$cell_next_event_times <- state$cell_next_event_times[-n_alive]
  state$hotspot_status        <- state$hotspot_status[-n_alive]

  # Free parent sequence memory
  state$cell_sequences[[current_cell_id]] <- NULL

  # Add daughters
  state$cell_ids              <- c(state$cell_ids, l_cell_id, r_cell_id)
  state$cell_next_event_times <- c(
    state$cell_next_event_times,
    state$time + stats::rexp(1, l_birth_rate + l_death_rate),
    state$time + stats::rexp(1, r_birth_rate + r_death_rate)
  )
  state$hotspot_status           <- c(state$hotspot_status, l_hotspot, r_hotspot)
  state$cell_sequences[[l_cell_id]] <- l_sequences
  state$cell_sequences[[r_cell_id]] <- r_sequences

  # ---- Append to pre-allocated history (O(1), no bind_rows) ----
  i_l <- state$h_n + 1L
  i_r <- state$h_n + 2L
  state$h_n <- i_r

  state$h_cell_id[i_l]    <- l_cell_id
  state$h_parent_id[i_l]  <- current_cell_id
  state$h_bfb_event[i_l]  <- bfb_occurred
  state$h_wgd_event[i_l]  <- wgd_occurred
  state$h_cn_event[i_l]   <- if (length(l_events) == 0) "none" else paste(l_events, collapse = ",")
  state$h_chr_allele[i_l] <- if (length(l_chr_alleles) == 0) NA_character_ else paste(l_chr_alleles, collapse = ",")

  state$h_cell_id[i_r]    <- r_cell_id
  state$h_parent_id[i_r]  <- current_cell_id
  state$h_bfb_event[i_r]  <- bfb_occurred
  state$h_wgd_event[i_r]  <- wgd_occurred
  state$h_cn_event[i_r]   <- if (length(r_events) == 0) "none" else paste(r_events, collapse = ",")
  state$h_chr_allele[i_r] <- if (length(r_chr_alleles) == 0) NA_character_ else paste(r_chr_alleles, collapse = ",")

  state
}

#' Prepare final results from simulation
#'
#' @param state The simulation state
#' @param return_phylo Logical. Whether to build the phylogenetic tree. Default TRUE.
#'
#' @return A list containing simulation results
prepare_results <- function(state, return_phylo = TRUE) {
  # Alive cells = whatever remains in the active arrays at simulation end
  alive_cell_ids <- state$cell_ids

  # Build history from pre-allocated vectors (trim to used length)
  n <- state$h_n
  if (n > 0L) {
    cell_history <- dplyr::tibble(
      cell_id    = state$h_cell_id[1:n],
      parent_id  = state$h_parent_id[1:n],
      bfb_event  = state$h_bfb_event[1:n],
      wgd_event  = state$h_wgd_event[1:n],
      cn_event   = state$h_cn_event[1:n],
      chr_allele = state$h_chr_allele[1:n],
      # is_alive computed in O(n) once at the end — no per-step mutate needed
      is_alive   = state$h_cell_id[1:n] %in% alive_cell_ids
    )
  } else {
    cell_history <- dplyr::tibble(
      cell_id    = character(0),
      parent_id  = character(0),
      bfb_event  = logical(0),
      wgd_event  = logical(0),
      cn_event   = character(0),
      chr_allele = character(0),
      is_alive   = logical(0)
    )
  }

  # Add root node
  root_row <- tibble::tibble(
    cell_id    = "root",
    parent_id  = NA_character_,
    bfb_event  = FALSE,
    wgd_event  = FALSE,
    cn_event   = "none",
    chr_allele = NA_character_,
    is_alive   = FALSE
  )
  cell_history <- dplyr::bind_rows(root_row, cell_history)

  # Final cells: only alive sequences (already cleaned during simulation)
  final_cells <- state$cell_sequences[alive_cell_ids]

  result <- list(
    cells            = final_cells,
    cell_history     = cell_history,
    tree             = if (return_phylo) build_phylo_from_lineage(cell_history) else NULL,
    input_parameters = state$input_parameters
  )

  return(result)
}

#' Process a death event
#'
#' @param state The simulation state
#' @param current_cell_id ID of the cell undergoing death
#' @param cell_idx Index of the cell in the alive arrays (from get_next_event)
#'
#' @return Updated simulation state
process_death_event <- function(state, current_cell_id, cell_idx) {
  # Remove dead cell from active arrays by swap-with-last then shrink.
  # This avoids shifting the entire vector and keeps arrays compact.
  n_alive <- length(state$cell_ids)
  if (cell_idx < n_alive) {
    state$cell_ids[cell_idx]              <- state$cell_ids[n_alive]
    state$cell_next_event_times[cell_idx] <- state$cell_next_event_times[n_alive]
    state$hotspot_status[cell_idx]        <- state$hotspot_status[n_alive]
  }
  state$cell_ids              <- state$cell_ids[-n_alive]
  state$cell_next_event_times <- state$cell_next_event_times[-n_alive]
  state$hotspot_status        <- state$hotspot_status[-n_alive]

  # Free sequence memory
  state$cell_sequences[[current_cell_id]] <- NULL

  # No history entry needed for explicit deaths:
  # is_alive is computed at the end as cell_id %in% alive_cell_ids.
  # Cells removed here will simply not appear in alive_cell_ids.

  state
}


#' Check if simulation should continue
#'
#' @param state The simulation state
#'
#' @return Logical indicating whether simulation should continue
continue_simulation <- function(state) {
  # Active arrays contain only alive cells, so length() is O(1).
  alive_count <- length(state$cell_ids)
  state$time < state$input_parameters$max_time &&
    alive_count > 0 &&
    alive_count < state$input_parameters$max_cells
}

#' Determine the next event to occur
#'
#' @param state The simulation state
#'
#' @return List with time, cell ID, index in alive arrays, and event type
get_next_event <- function(state) {
  # All entries in cell_ids are alive — no alive filtering needed.
  n_alive <- length(state$cell_ids)
  if (n_alive == 0L) return(NULL)

  # which.min over the compact alive vector only
  min_idx        <- which.min(state$cell_next_event_times)
  current_cell_id <- state$cell_ids[min_idx]
  hotspot_status  <- state$hotspot_status[min_idx]

  # Calculate modified birth and death rates for this specific cell
  cell_birth_rate <- state$input_parameters$birth_rate *
    (1 + state$input_parameters$positive_selection_rate * hotspot_status)
  cell_death_rate <- state$input_parameters$death_rate *
    (1 + state$input_parameters$negative_selection_rate * hotspot_status)

  event_probability <- cell_birth_rate / (cell_birth_rate + cell_death_rate)
  event_type <- if (stats::runif(1) < event_probability) "birth" else "death"

  list(
    time      = state$cell_next_event_times[min_idx],
    cell_id   = current_cell_id,
    cell_idx  = min_idx,   # index in alive arrays — passed to event handlers
    event_type = event_type
  )
}


# Cell to cn copy data
# sequences_to_cndata <- function(sequences, chr_seq_lengths, bin_length) {
#   cndata <- lapply(names(sequences), function(cell_id) {
#     seqs <- sequences[[cell_id]]
#     lapply(names(seqs), function(chr_allele) {
#       s <- seqs[[chr_allele]]
#
#       chr <- strsplit(chr_allele, ":")[[1]][1]
#       allele <- strsplit(chr_allele, ":")[[1]][2]
#
#       bin_idxs <- 1:chr_seq_lengths[[chr]]
#       t <- table(seq2vec(s))[bin_idxs]
#       d <- dplyr::tibble(
#         cell_id = cell_id,
#         bin_idx = bin_idxs,
#         allele = allele,
#         chr = chr,
#         state = as.integer(t[as.character(bin_idxs)])
#       )
#       d$state[is.na(d$state)] <- 0
#       d
#     }) %>%
#       do.call(dplyr::bind_rows, .)
#   }) %>%
#     do.call(dplyr::bind_rows, .)
#
#   cndata <- cndata %>%
#     dplyr::group_by(.data$cell_id, .data$chr) %>%
#     tidyr::pivot_wider(values_from = .data$state, names_from = .data$allele) %>%
#     dplyr::mutate(CN = .data$A + .data$B) %>%
#     dplyr::mutate(
#       start = (.data$bin_idx - 1) * bin_length + 1,
#       end = .data$bin_idx * bin_length
#     )
#
#   cndata
# }


sim_amp_del <- function(sequence, operation = "dup", rate = 1e7) {
  # Work directly on the interval representation — avoids the expensive
  # seq2vec (expand) + find_consecutive_subvectors + vec2seq (recompress) cycle.
  # Each stored interval is already a consecutive monotone run, which is
  # exactly the kind of sub-sequence the original code searched for.

  n_iv <- length(sequence)
  if (n_iv == 0L) return(sequence)

  # Compute interval lengths in O(n_intervals)
  iv_lengths <- integer(n_iv)
  for (j in seq_len(n_iv)) {
    iv <- sequence[[j]]
    iv_lengths[j] <- abs(iv$end - iv$start) + 1L
  }
  total_len <- sum(iv_lengths)
  if (total_len == 0L) return(sequence)

  # Only operate on directional intervals (length >= 2); fall back to all
  # intervals (including single-element ones) if none qualify.
  eligible <- which(iv_lengths >= 2L)
  if (length(eligible) == 0L) eligible <- seq_len(n_iv)

  # Sample one interval uniformly (mirrors the original's uniform-over-runs)
  idx    <- eligible[sample.int(length(eligible), 1L)]
  iv     <- sequence[[idx]]
  iv_len <- iv_lengths[idx]

  # Sample event length
  event_length  <- max(1L, round(stats::rexp(1, rate = 1 / rate)))
  actual_length <- min(event_length, iv_len)

  # Sample starting offset within the chosen interval
  max_offset <- iv_len - actual_length
  offset     <- if (max_offset <= 0L) 0L else sample.int(max_offset + 1L, 1L) - 1L

  # Derive the sub-segment and the before/after fragments of the interval
  dir <- iv$direction
  if (dir == 1L) {
    seg_s    <- iv$start + offset
    seg_e    <- seg_s + actual_length - 1L
    before_iv <- if (offset > 0L)
      list(start = iv$start, end = seg_s - 1L, direction = 1L) else NULL
    after_iv  <- if (seg_e < iv$end)
      list(start = seg_e + 1L, end = iv$end,   direction = 1L) else NULL
  } else if (dir == -1L) {
    seg_s    <- iv$start - offset
    seg_e    <- seg_s - actual_length + 1L
    before_iv <- if (offset > 0L)
      list(start = iv$start, end = seg_s + 1L, direction = -1L) else NULL
    after_iv  <- if (seg_e > iv$end)
      list(start = seg_e - 1L, end = iv$end,   direction = -1L) else NULL
  } else {
    # direction == 0: single-element interval
    seg_s <- iv$start; seg_e <- iv$start
    before_iv <- NULL;  after_iv  <- NULL
  }
  seg <- list(start = seg_s, end = seg_e, direction = dir)

  # Prefix and suffix of the interval list surrounding the chosen interval
  head_ivs <- if (idx > 1L)    sequence[seq_len(idx - 1L)]      else list()
  tail_ivs <- if (idx < n_iv)  sequence[(idx + 1L):n_iv]        else list()

  if (operation == "dup") {
    mid <- c(
      if (!is.null(before_iv)) list(before_iv),
      list(seg, seg),
      if (!is.null(after_iv)) list(after_iv)
    )
  } else if (operation == "del") {
    mid <- c(
      if (!is.null(before_iv)) list(before_iv),
      if (!is.null(after_iv)) list(after_iv)
    )
    if (length(head_ivs) == 0L && length(mid) == 0L && length(tail_ivs) == 0L) {
      message("skipping deletion because of length zero")
      return(sequence)
    }
  } else {
    stop("operation not recognized!")
  }

  c(head_ivs, mid, tail_ivs)
}

sim_wgd = function(sequence) {
  vec = seq2vec(sequence)
  wgd_vec = c(vec, vec)
  vec2seq(wgd_vec)
}

#' Simulate Breakage-Fusion-Bridge (BFB) Cycle for both daughters
#'
#' @param sequence Input sequence
#' @param support Distribution type for breakpoint selection ("uniform" or "beta")
#' @param alpha Shape parameter for beta distribution (only used if support="beta")
#' @param beta Shape parameter for beta distribution (only used if support="beta")
#' @param custom_breakpoints .
#'
#' @details Simulate left and right children from a BFB cycle using specified
#'   breakpoint selection distribution
#'
#' @return List containing left and right sequences
sim_bfb_left_and_right_sequences <- function(
    sequence,
    support = "uniform",
    alpha = NULL,
    beta = NULL,
    custom_breakpoints = NULL
) {
  # Calculate the total number of elements in the sequence
  L <- get_seq_length(sequence)

  # Find "fusion values": bin values that appear at the junction between two
  # consecutive intervals (i.e. interval[k]$end == interval[k+1]$start).
  # These correspond to existing BFB fold-back points (where diff(vec)==0 in
  # the expanded representation).  Computing directly from the interval list
  # is O(n_intervals) and avoids the expensive seq2vec expansion.
  n_iv <- length(sequence)
  bps  <- integer(0L)
  if (n_iv >= 2L) {
    for (k in seq_len(n_iv - 1L)) {
      if (sequence[[k]]$end == sequence[[k + 1L]]$start) {
        bps <- c(bps, sequence[[k]]$end)
      }
    }
  }

  # Select random breakpoint based on specified distribution
  # Ensure that bp_idx is different from L to obtain a proper bfb cycle
  bp_idx = L
  attempts = 0
  max_attempts = 10

  while (bp_idx %in% c(L, bps) && attempts < max_attempts) {
    attempts = attempts + 1

    if (support == "uniform") {
      bp_idx = sample(1:(2 * L), 1)
    } else if (support == "beta") {
      if (is.null(alpha) || is.null(beta)) {
        stop(
          "For beta distribution, both alpha and beta parameters must be provided"
        )
      }
      tau = stats::rbeta(1, alpha, beta)
      bp_idx = max(1, round(tau * 2 * L)) # Ensure bp_idx is at least 1
    } else if (support == "custom") {
      if (is.null(custom_breakpoints) || length(custom_breakpoints) == 0) {
        stop(
          "For custom distribution, custom_breakpoints vector must be provided and non-empty"
        )
      }
      # Find all indices in vec that match the custom breakpoints
      # This creates a vector of potential bp_idx values
      valid_indices <- which(vec %in% custom_breakpoints)
      if (length(valid_indices) == 0) {
        stop(
          "None of the custom breakpoints are present in the sequence"
        )
      }
      valid_idx = sample(valid_indices, size = 1)
      breakpoint_positions <- vec[valid_idx]
      bp_idx = sample(breakpoint_positions, size = 1)
    } else {
      stop("Unsupported distribution type. Use 'uniform', 'beta', or 'custom'.")
    }
  }

  # Check if we exceeded maximum attempts
  if (attempts >= max_attempts && bp_idx %in% c(L, bps)) {
    warning("BFB breakpoint selection failed after ", max_attempts, " attempts (sequence may be too fragmented). Returning unchanged sequence as fallback.")
    return(list(l_seq = sequence, r_seq = sequence))
  }

  # Initialize left and right sequences
  cut_seqs = cut_sequence(fuse_sequence(sequence), bp_idx)
  l_seq = cut_seqs$left_seq
  r_seq = reverse_sequence(cut_seqs$right_seq)

  # Return the left and right sequences
  if (stats::runif(1, 0, 1) > .5) {
    return(list(l_seq = l_seq, r_seq = r_seq))
  } else {
    return(list(l_seq = r_seq, r_seq = l_seq))
  }
}

reverse_sequence <- function(sequence) {
  # Reverse the order of the intervals and swap start and end, flip direction
  reversed_seq = lapply(seq_along(sequence), function(i) {
    interval = sequence[[length(sequence) - i + 1]]
    list(
      start = interval$end,
      end = interval$start,
      direction = -interval$direction
    )
  })

  return(reversed_seq)
}

fuse_sequence <- function(sequence) {
  reversed_seq = reverse_sequence(sequence)
  fused_seq = c(sequence, reversed_seq)
  return(fused_seq)
}

cut_sequence <- function(sequence, cut_index) {
  # Initialize output sequences

  left_seq <- list()
  right_seq <- list()

  # Track the current position across the entire sequence
  current_length <- 0

  for (interval in sequence) {
    # Calculate the length of the current interval
    interval_length <- abs(interval$end - interval$start) + 1

    # If the cut point is before this interval, everything goes to the right
    if (current_length >= cut_index) {
      right_seq <- c(right_seq, list(interval))

      # If the cut point is after this interval, everything goes to the left
    } else if (current_length + interval_length <= cut_index) {
      left_seq <- c(left_seq, list(interval))

      # Otherwise, we split the interval
    } else {
      # How far into the current interval is the cut?
      cut_within <- cut_index - current_length

      # Handle different directions
      if (interval$direction == 1) {
        # Increasing interval
        left_seq <- c(
          left_seq,
          list(list(
            start = interval$start,
            end = interval$start + cut_within - 1,
            direction = 1
          ))
        )
        right_seq <- c(
          right_seq,
          list(list(
            start = interval$start + cut_within,
            end = interval$end,
            direction = 1
          ))
        )
      } else if (interval$direction == -1) {
        # Decreasing interval
        left_seq <- c(
          left_seq,
          list(list(
            start = interval$start,
            end = interval$start - cut_within + 1,
            direction = -1
          ))
        )
        right_seq <- c(
          right_seq,
          list(list(
            start = interval$start - cut_within,
            end = interval$end,
            direction = -1
          ))
        )
      } else {
        # Constant interval
        left_seq <- c(left_seq, list(interval))
        right_seq <- c(right_seq, list(interval))
      }
    }

    # Update the position tracker
    current_length <- current_length + interval_length
  }

  return(list(left_seq = left_seq, right_seq = right_seq))
}

is_hotspot_gained <- function(cell, hotspot) {
  hc <- get_hotspot_copies(cell, hotspot)
  if (is.nan(hc)) return(NA)
  hc > 1
}


# Count how many intervals in `cell` contain `hotspot` (a bin index).
# O(n_intervals) — avoids the seq2vec + table expansion used previously.
get_hotspot_copies <- function(cell, hotspot) {
  if (is.null(hotspot)) return(NaN)
  count <- 0L
  for (iv in cell) {
    lo <- if (iv$start <= iv$end) iv$start else iv$end
    hi <- if (iv$start <= iv$end) iv$end   else iv$start
    if (hotspot >= lo && hotspot <= hi) count <- count + 1L
  }
  count
}

# cell_history_to_newick <- function(cell_history) {
#   # Check there is a root and rename it
#   root_name = cell_history$cell_id[is.na(cell_history$parent_id)]
#   cell_history$cell_id[cell_history$cell_id == root_name] = "root"
#   cell_history$parent_id[cell_history$parent_id == root_name] = "root"
#
#   # Helper function to recursively build the tree
#   build_tree <- function(node) {
#     # Find children of the current node
#     node_data <- cell_history %>%
#       dplyr::filter(.data$cell_id == node) %>%
#       dplyr::select(.data$bfb_event)
#
#     # Find children of the current node
#     children <- cell_history %>% dplyr::filter(.data$parent_id == node) %>% dplyr::pull(.data$cell_id)
#
#     if (length(children) == 0) {
#       # If no children, return the node itself with BFB annotation
#       return(node)
#     } else {
#       # Recursively build subtrees for each child
#       subtree <- paste(sapply(children, build_tree), collapse = ",")
#
#       # Add BFB annotation to the current node
#       return(paste0("(", subtree, ")", node))
#     }
#   }
#
#   # Identify the root node (cells with no parent)
#   root <- cell_history %>%
#     dplyr::filter(is.na(.data$parent_id)) %>%
#     dplyr::pull(.data$cell_id)
#
#   if (length(root) != 1) {
#     stop("Error: There must be exactly one root node.")
#   }
#
#   # Build the tree starting from the root
#   newick_tree <- paste0(build_tree(root), ";")
#
#   return(newick_tree)
# }

cell_history_to_newick <- function(cell_history) {
  # Check there is a root and rename it
  root_name = cell_history$cell_id[is.na(cell_history$parent_id)]
  cell_history$cell_id[cell_history$cell_id == root_name] = "root"
  cell_history$parent_id[cell_history$parent_id == root_name] = "root"

  # Helper function to check if a node has any living descendants
  has_living_descendants <- function(node) {
    # Get node data
    node_data <- cell_history %>%
      dplyr::filter(.data$cell_id == node)

    # If node doesn't exist, return FALSE
    if (nrow(node_data) == 0) {
      return(FALSE)
    }

    # If this node is alive, return TRUE
    if (node_data$is_alive) {
      return(TRUE)
    }

    # Find children of the current node
    children <- cell_history %>%
      dplyr::filter(.data$parent_id == node) %>%
      dplyr::pull(.data$cell_id)

    # If no children, this is a leaf - return whether it's alive
    if (length(children) == 0) {
      return(node_data$is_alive)
    }

    # If has children, check if any child has living descendants
    child_results <- vapply(children, has_living_descendants, logical(1))
    return(any(child_results))
  }

  # Helper function to recursively build the tree (only for nodes with living descendants)
  build_tree <- function(node) {
    # Skip this node if it has no living descendants
    if (!has_living_descendants(node)) {
      return(NULL)
    }

    # Find children of the current node
    children <- cell_history %>%
      dplyr::filter(.data$parent_id == node) %>%
      dplyr::pull(.data$cell_id)

    # Filter children to only those with living descendants
    living_children <- children[vapply(
      children,
      has_living_descendants,
      logical(1)
    )]

    if (length(living_children) == 0) {
      # If no living children, return the node itself (this should be a living leaf)
      return(node)
    } else {
      # Recursively build subtrees for each living child
      subtrees <- sapply(living_children, build_tree)
      # Remove any NULL subtrees (shouldn't happen with our filtering, but safety check)
      subtrees <- subtrees[!sapply(subtrees, is.null)]

      if (length(subtrees) == 0) {
        return(node)
      } else {
        subtree_str <- paste(subtrees, collapse = ",")
        return(paste0("(", subtree_str, ")", node))
      }
    }
  }

  # Identify the root node (cells with no parent)
  root <- cell_history %>%
    dplyr::filter(is.na(.data$parent_id)) %>%
    dplyr::pull(.data$cell_id)

  if (length(root) != 1) {
    stop("Error: There must be exactly one root node.")
  }

  # Check if root has any living descendants
  if (!has_living_descendants(root)) {
    stop("Error: No living cells found in the tree.")
  }

  # Build the tree starting from the root
  newick_tree <- paste0(build_tree(root), ";")
  return(newick_tree)
}

validate_bridge_sim_params <- function(
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
) {
  # Valid chromosome names
  valid_chromosomes <- c(as.character(1:22), "X", "Y")

  # Valid breakpoint support distributions
  valid_breakpoint_supports <- c("uniform", "beta", "custom")

  # --- Numeric Parameter Validation ---

  # initial_cells: positive integer
  if (
    !is.numeric(initial_cells) ||
      initial_cells <= 0 ||
      initial_cells != round(initial_cells)
  ) {
    stop("initial_cells must be a positive integer")
  }

  # bin_length: positive number
  if (!is.numeric(bin_length) || bin_length <= 0) {
    stop("bin_length must be a positive number")
  }

  # birth_rate: non-negative number
  if (!is.numeric(birth_rate) || birth_rate < 0) {
    stop("birth_rate must be a non-negative number")
  }

  # death_rate: non-negative number
  if (!is.numeric(death_rate) || death_rate < 0) {
    stop("death_rate must be a non-negative number")
  }

  # Rate parameters: non-negative numbers
  rate_params <- list(
    normal_dup_rate = normal_dup_rate,
    bfb_prob = bfb_prob,
    amp_rate = amp_rate,
    del_rate = del_rate
  )

  for (param_name in names(rate_params)) {
    param_value <- rate_params[[param_name]]
    if (!is.numeric(param_value) || param_value < 0) {
      stop(paste(param_name, "must be a non-negative number"))
    }
  }

  # Check that at least one rate is positive
  if (sum(unlist(rate_params)) == 0) {
    stop(
      "At least one of normal_dup_rate, bfb_prob, amp_rate, or del_rate must be positive"
    )
  }

  # Selection rates: numeric
  if (!is.numeric(positive_selection_rate)) {
    stop("positive_selection_rate must be numeric")
  }

  if (!is.numeric(negative_selection_rate)) {
    stop("negative_selection_rate must be numeric")
  }

  # max_time: positive number
  if (!is.numeric(max_time) || max_time <= 0) {
    stop("max_time must be a positive number")
  }

  # max_cells: positive integer
  if (
    !is.numeric(max_cells) || max_cells <= 0 || max_cells != round(max_cells)
  ) {
    stop("max_cells must be a positive integer")
  }

  # --- Character/String Parameter Validation ---

  # chromosomes: valid chromosome names
  chromosomes <- as.character(chromosomes)
  invalid_chrs <- setdiff(chromosomes, valid_chromosomes)
  if (length(invalid_chrs) > 0) {
    stop(paste(
      "Invalid chromosome(s):",
      paste(invalid_chrs, collapse = ", "),
      "\nValid chromosomes are:",
      paste(valid_chromosomes, collapse = ", ")
    ))
  }

  # bfb_allele: must be in format "chr:allele"
  if (!is.character(bfb_allele) || length(bfb_allele) != 1) {
    stop("bfb_allele must be a single character string")
  }

  # Parse bfb_allele format
  bfb_parts <- strsplit(bfb_allele, ":")[[1]]
  if (length(bfb_parts) != 2) {
    stop("bfb_allele must be in format 'chromosome:allele' (e.g., '1:A')")
  }

  bfb_chr <- bfb_parts[1]
  bfb_allele_letter <- bfb_parts[2]

  if (!bfb_chr %in% valid_chromosomes) {
    stop(paste("Invalid chromosome in bfb_allele:", bfb_chr))
  }

  if (!bfb_allele_letter %in% c("A", "B")) {
    stop("Allele in bfb_allele must be 'A' or 'B'")
  }

  # breakpoint_support: valid distribution
  if (!breakpoint_support %in% valid_breakpoint_supports) {
    stop(paste(
      "breakpoint_support must be one of:",
      paste(valid_breakpoint_supports, collapse = ", ")
    ))
  }

  # --- Logical Parameter Validation ---

  logical_params <- list(
    first_round_of_bfb = first_round_of_bfb
  )

  for (param_name in names(logical_params)) {
    param_value <- logical_params[[param_name]]
    if (!is.logical(param_value) || length(param_value) != 1) {
      stop(paste(param_name, "must be a single logical value (TRUE or FALSE)"))
    }
  }

  # --- List Parameter Validation ---

  # hotspot: must be a named list with chr and pos
  if (!is.list(hotspot)) {
    stop("hotspot must be a list")
  }

  required_hotspot_names <- c("chr", "pos")
  if (!all(required_hotspot_names %in% names(hotspot))) {
    stop(paste(
      "hotspot must contain named elements:",
      paste(required_hotspot_names, collapse = ", ")
    ))
  }

  # Validate hotspot chromosome format
  if (!is.character(hotspot$chr) || length(hotspot$chr) != 1) {
    stop("hotspot$chr must be a single character string")
  }

  hotspot_parts <- strsplit(hotspot$chr, ":")[[1]]
  if (length(hotspot_parts) != 2) {
    stop("hotspot$chr must be in format 'chromosome:allele' (e.g., '1:A')")
  }

  hotspot_chr <- hotspot_parts[1]
  hotspot_allele <- hotspot_parts[2]

  if (!hotspot_chr %in% valid_chromosomes) {
    stop(paste("Invalid chromosome in hotspot$chr:", hotspot_chr))
  }

  if (!hotspot_allele %in% c("A", "B")) {
    stop("Allele in hotspot$chr must be 'A' or 'B'")
  }

  # Validate hotspot position
  if (!is.numeric(hotspot$pos) || hotspot$pos <= 0) {
    stop("hotspot$pos must be a positive number")
  }

  # --- Beta Distribution Parameter Validation ---

  if (breakpoint_support == "beta") {
    if (is.null(alpha) || is.null(beta)) {
      stop(
        "alpha and beta parameters must be provided when breakpoint_support is 'beta'"
      )
    }

    if (!is.numeric(alpha) || alpha <= 0) {
      stop("alpha must be a positive number for beta distribution")
    }

    if (!is.numeric(beta) || beta <= 0) {
      stop("beta must be a positive number for beta distribution")
    }
  }

  # --- Cross-Parameter Validation ---

  # Check that bfb_allele chromosome is included in chromosomes
  if (!bfb_chr %in% chromosomes) {
    stop(paste(
      "bfb_allele chromosome",
      bfb_chr,
      "must be included in the chromosomes parameter"
    ))
  }

  # Check that hotspot chromosome is included in chromosomes
  if (!hotspot_chr %in% chromosomes) {
    stop(paste(
      "hotspot chromosome",
      hotspot_chr,
      "must be included in the chromosomes parameter"
    ))
  }

  # Warn if birth_rate is much smaller than death_rate
  if (birth_rate > 0 && death_rate > birth_rate * 10) {
    warning(
      "death_rate is much larger than birth_rate - population may quickly go extinct"
    )
  }

  # Warn if max_cells is very large
  if (max_cells > 10000) {
    warning(
      "max_cells is very large - simulation may be slow or use excessive memory"
    )
  }

  # All validations passed
  return(invisible(NULL))
}


subsample_sim <- function(sim_result, f_subsample = 1) {
  if (f_subsample < 0 || f_subsample > 1) {
    stop("f_subsample must be between 0 and 1")
  }

  if (f_subsample == 1) {
    return(sim_result)
  }

  # Get original alive cell IDs
  original_cell_ids <- names(sim_result$cells)

  # Subsample cells
  n_subsample <- as.integer(f_subsample * length(original_cell_ids))
  if (n_subsample == 0) {
    stop("Subsampling fraction too small - results in 0 cells")
  }

  sampled_cell_ids <- sample(original_cell_ids, n_subsample)

  # Create subsampled cells
  subsampled_cells <- sim_result$cells[sampled_cell_ids]

  # Filter cell history to remove non-sampled alive cells
  alive_cells_not_sampled <- sim_result$cell_history %>%
    dplyr::filter(.data$is_alive, !.data$cell_id %in% sampled_cell_ids) %>%
    dplyr::pull(.data$cell_id)

  subsampled_cell_history <- sim_result$cell_history %>%
    dplyr::filter(!.data$cell_id %in% alive_cells_not_sampled)

  # Rebuild tree from filtered history only when the original had one.
  # When return_phylo = FALSE was passed to bridge_sim, sim_result$tree is NULL
  # and we skip this expensive step entirely.
  if (!is.null(sim_result$tree)) {
    subsampled_tree <- build_phylo_from_lineage(subsampled_cell_history)
    subsampled_tree <- ape::keep.tip(subsampled_tree, sampled_cell_ids)
  } else {
    subsampled_tree <- NULL
  }

  # Create new CNA data for subsampled cells
  chr_seq_lengths <- sim_result$input_parameters$chr_seq_lengths
  bin_length <- sim_result$input_parameters$bin_length
  subsampled_cna_data <- sequences_to_cndata(
    subsampled_cells,
    chr_seq_lengths,
    bin_length
  )

  # Return subsampled result
  subsampled_result <- list(
    cells = subsampled_cells,
    cell_history = subsampled_cell_history,
    tree = subsampled_tree,
    cna_data = subsampled_cna_data,
    input_parameters = sim_result$input_parameters
  )

  return(subsampled_result)
}

sequences_to_cndata <- function(sequences, chr_seq_lengths, bin_length) {
  # Pre-compute total rows needed
  are_null = lapply(sequences, function(s) {
    is.null(s)
  }) %>% unlist()



  total_rows <- sum(unlist(lapply(sequences, function(seqs) {
    sum(chr_seq_lengths[sapply(names(seqs), function(x) strsplit(x, ":")[[1]][1])])
  })))

  # Pre-allocate vectors for maximum efficiency
  cell_ids <- character(total_rows)
  bin_idxs <- integer(total_rows)
  alleles <- character(total_rows)
  chrs <- character(total_rows)
  states <- integer(total_rows)

  row_idx <- 1

  for (cell_id in names(sequences)) {
    seqs <- sequences[[cell_id]]

    for (chr_allele in names(seqs)) {
      s <- seqs[[chr_allele]]
      chr_parts <- strsplit(chr_allele, ":", fixed = TRUE)[[1]]
      chr <- chr_parts[1]
      allele <- chr_parts[2]

      n_bins <- chr_seq_lengths[[chr]]
      bin_range <- seq_len(n_bins)

      # Fill pre-allocated vectors
      end_idx <- row_idx + n_bins - 1
      cell_ids[row_idx:end_idx] <- cell_id
      bin_idxs[row_idx:end_idx] <- bin_range
      alleles[row_idx:end_idx] <- allele
      chrs[row_idx:end_idx] <- chr

      # Process states
      seq_vec <- seq2vec(s)
      t <- table(seq_vec)
      state_vec <- integer(n_bins)
      matching_bins <- intersect(as.integer(names(t)), bin_range)
      if (length(matching_bins) > 0) {
        state_vec[matching_bins] <- t[as.character(matching_bins)]
      }
      states[row_idx:end_idx] <- state_vec

      row_idx <- end_idx + 1
    }
  }

  # Create data.table from pre-allocated vectors
  cndata <- data.table::data.table(
    cell_id = cell_ids,
    bin_idx = bin_idxs,
    allele = alleles,
    chr = chrs,
    state = states
  )

  # Pivot and compute final columns
  cndata_wide <- data.table::dcast(cndata, cell_id + chr + bin_idx ~ allele,
                       value.var = "state", fill = 0)

  cndata_wide[, CN := A + B]
  cndata_wide[, `:=`(
    start = (bin_idx - 1) * bin_length + 1,
    end = bin_idx * bin_length
  )]

  dplyr::as_tibble(cndata_wide)
}


build_phylo_from_lineage <- function(cell_history) {
  # Handle root naming consistently
  if (any(is.na(cell_history$parent_id))) {
    root_rows <- which(is.na(cell_history$parent_id))
    if (length(root_rows) != 1) {
      stop("Error: There must be exactly one root node.")
    }

    # If root is not already named "root", rename it
    root_name <- cell_history$cell_id[root_rows]
    if (root_name != "root") {
      cell_history$cell_id[cell_history$cell_id == root_name] <- "root"
      cell_history$parent_id[cell_history$parent_id == root_name] <- "root"
    }
  }

  # Create a mapping from cell names to node numbers
  # Tips (terminal nodes) get numbers 1:n_tips
  # Internal nodes get numbers (n_tips+1):(n_tips+n_internal)

  all_cells <- unique(cell_history$cell_id)
  n_cells <- length(all_cells)

  # Identify terminal nodes (nodes that are not parents)
  parent_cells <- unique(cell_history$parent_id[!is.na(cell_history$parent_id)])
  terminal_cells <- setdiff(all_cells, parent_cells)
  internal_cells <- intersect(all_cells, parent_cells)

  n_tips <- length(terminal_cells)
  n_internal <- length(internal_cells)

  # Create node number mapping
  # Tips: 1 to n_tips
  # Internal nodes: (n_tips + 1) to (n_tips + n_internal)
  cell_to_node <- integer(n_cells)
  names(cell_to_node) <- c(terminal_cells, internal_cells)
  cell_to_node[terminal_cells] <- 1:n_tips
  cell_to_node[internal_cells] <- (n_tips + 1):(n_tips + n_internal)

  # Build the edge matrix
  # Each row represents an edge: [parent_node, child_node]
  edges <- data.frame(
    parent = character(0),
    child = character(0),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(cell_history)) {
    child <- cell_history$cell_id[i]
    parent <- cell_history$parent_id[i]

    if (!is.na(parent)) {
      edges <- rbind(edges, data.frame(parent = parent, child = child, stringsAsFactors = FALSE))
    }
  }

  # Convert to node numbers
  edge_matrix <- matrix(0, nrow = nrow(edges), ncol = 2)
  edge_matrix[, 1] <- cell_to_node[edges$parent]  # parent nodes
  edge_matrix[, 2] <- cell_to_node[edges$child]   # child nodes

  # Create edge lengths (all set to 1 for simplicity, can be modified)
  edge_lengths <- rep(1, nrow(edge_matrix))

  # Create tip labels
  tip_labels <- terminal_cells

  # Find root node number
  root_node <- cell_to_node["root"]

  # Create the phylo object
  phylo_tree <- list(
    edge = edge_matrix,
    edge.length = edge_lengths,
    tip.label = tip_labels,
    Nnode = n_internal,
    root.edge = NULL
  )

  class(phylo_tree) <- "phylo"

  # Validate the tree structure
  if (!is.null(phylo_tree)) {
    tryCatch({
      # Basic validation - check if tree is valid
      checkValidPhylo(phylo_tree)
    }, error = function(e) {
      warning("Created tree may have structural issues: ", e$message)
    })
  }

  return(phylo_tree)
}

# Helper function to check phylo validity (basic checks)
checkValidPhylo <- function(phylo_obj) {
  if (!inherits(phylo_obj, "phylo")) {
    stop("Object is not of class 'phylo'")
  }

  if (nrow(phylo_obj$edge) == 0) {
    stop("Edge matrix is empty")
  }

  if (length(phylo_obj$tip.label) == 0) {
    stop("No tip labels found")
  }

  # Check that edge matrix has valid node numbers
  max_node <- max(phylo_obj$edge)
  expected_max <- length(phylo_obj$tip.label) + phylo_obj$Nnode

  if (max_node > expected_max) {
    stop("Edge matrix contains invalid node numbers")
  }

  return(TRUE)
}
