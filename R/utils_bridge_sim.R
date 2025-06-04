#' Helper functions modified for diploid chromosome modeling

#' Initialize simulation state for diploid chromosomes
#'
#' @param input_parameters List of parameters given as input to bridge_sim function.
#'
#' @return A list containing the initialized simulation state
initialize_simulation <- function(input_parameters) {
  # Initialize state variables
  state <- list(
    time = 0,
    #birth_count = 0,
    #death_count = 0,
    cell_ids = character(0),
    cell_sequences = list(),  # Will contain chromosome sequences for each cell
    #cell_birth_times = numeric(0),
    #cell_death_times = numeric(0),
    #cell_parents = character(0),
    #cell_bfb = logical(0),
    #cell_replication_history = character(0),
    #cell_hotspot_gained = logical(0),
    #cell_hotspot_copies = numeric(0),
    cell_is_alive = logical(0),
    cell_next_event_times = numeric(0)#
    #cell_bfb_history = list()
  )
  state$input_parameters = input_parameters
  state$next_cell_id = 1

  # Create initial sequences for all chromosome alleles
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
      initial_sequences[[bfb_allele]],
      state$input_parameters$breakpoint_support,
      state$input_parameters$alpha,
      state$input_parameters$beta
    )

    # Left child
    left_cell = initial_sequences
    left_cell[[bfb_allele]] = daughter_seqs$l_seq
    left_cell_id = paste0("cell_", state$next_cell_id)
    state$next_cell_id = state$next_cell_id + 1

    # Right child
    right_cell = initial_sequences
    right_cell[[bfb_allele]] = daughter_seqs$r_seq
    right_cell_id = paste0("cell_", state$next_cell_id)
    state$next_cell_id = state$next_cell_id + 1

    if (!is.null(state$input_parameters$hotspot)) {
      hotspot = state$input_parameters$hotspot
      left_hotspot = is_hotspot_gained(left_cell[[hotspot$chr]], hotspot = hotspot$pos)
      right_hotspot = is_hotspot_gained(right_cell[[hotspot$chr]], hotspot = hotspot$pos)
    } else {
      left_hotspot = right_hotspot = FALSE
    }

    l_birth_rate <- state$input_parameters$birth_rate * (1 + state$input_parameters$positive_selection_rate * left_hotspot)
    l_death_rate <- state$input_parameters$death_rate * (1 + state$input_parameters$negative_selection_rate * left_hotspot)

    r_birth_rate <- state$input_parameters$birth_rate * (1 + state$input_parameters$positive_selection_rate * right_hotspot)
    r_death_rate <- state$input_parameters$death_rate * (1 + state$input_parameters$negative_selection_rate * right_hotspot)

    # Calculate combined rates
    l_combined_rate <- l_birth_rate + l_death_rate
    r_combined_rate <- r_birth_rate + r_death_rate

    # Update state
    state$cell_ids <- c(left_cell_id, right_cell_id)
    state$cell_sequences[[left_cell_id]] <- left_cell
    state$cell_sequences[[right_cell_id]] <- right_cell
    state$cell_is_alive <- c(state$cell_is_alive, TRUE, TRUE)
    state$hotspot_status = c(state$hotspot_status, left_hotspot, right_hotspot)

    # Initialize next event times with future events based on combined rates
    state$cell_next_event_times <- c(
      state$cell_next_event_times,
      state$time + stats::rexp(1, l_combined_rate),
      state$time + stats::rexp(1, r_combined_rate)
    )

    state$history = dplyr::bind_rows(
      state$history,
      dplyr::tibble(cell_id = state$cell_ids,
                    parent_id = "root",
                    bfb_event = TRUE,
                    cn_event = "bfb",
                    chr_allele = bfb_allele,
                    is_alive = state$cell_is_alive)
    )
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
    state$next_cell_id <- state$next_cell_id + 1

    # Hotspot logic
    if (!is.null(state$input_parameters$hotspot)) {
      hotspot <- state$input_parameters$hotspot
      hotspot_gained <- is_hotspot_gained(initial_sequences[[hotspot$chr]], hotspot = hotspot$pos)
    } else {
      hotspot_gained <- FALSE
    }

    # Birth and death rate with selection
    birth_rate <- state$input_parameters$birth_rate * (1 + state$input_parameters$positive_selection_rate * hotspot_gained)
    death_rate <- state$input_parameters$death_rate * (1 + state$input_parameters$negative_selection_rate * hotspot_gained)

    # Combined rate and next event
    combined_rate <- birth_rate + death_rate
    next_event_time <- state$time + stats::rexp(1, combined_rate)

    # Update state
    state$cell_ids <- c(state$cell_ids, cell_id)
    state$cell_sequences[[cell_id]] <- initial_sequences
    state$cell_is_alive <- c(state$cell_is_alive, TRUE)
    state$cell_next_event_times <- c(state$cell_next_event_times, next_event_time)
    state$hotspot_status = c(state$hotspot_status, hotspot_gained)

    # Update history
    chr_allele <- NA
    state$history <- dplyr::bind_rows(
      state$history,
      dplyr::tibble(
        cell_id = cell_id,
        parent_id = "root",
        bfb_event = FALSE,
        cn_event = "none",
        chr_allele = chr_allele,
        is_alive = TRUE
      )
    )
  }

  return(state)
}


#' Process a birth event (diploid version)
#'
#' @param state The simulation state
#' @param current_cell_id ID of the cell undergoing birth
#'
#' @return Updated simulation state
#' Process a birth event (diploid version)
#'
#' @param state The simulation state
#' @param current_cell_id ID of the cell undergoing birth
#'
#' @return Updated simulation state
process_birth_event <- function(state, current_cell_id) {
  # Find index of current cell
  cell_idx <- which(state$cell_ids == current_cell_id)

  # Decide type of event
  event_name <- sample(names(state$input_parameters$rates), size = 1, prob = unlist(state$input_parameters$rates))

  # Get parent cell information
  parent_sequences <- state$cell_sequences[[current_cell_id]]
  # parent_hotspot <- state$hotspot_status[cell_idx]

  # Create daughter cell sequences
  if (event_name == "normal") {
    l_sequences <- parent_sequences
    r_sequences <- parent_sequences
  } else {
    selected_chr_allele <- sample(names(parent_sequences), 1)

    if (event_name == "amp") {
      modified_seq <- sim_amp_del(parent_sequences[[selected_chr_allele]], operation = "dup")
    } else if (event_name == "del") {
      modified_seq <- sim_amp_del(parent_sequences[[selected_chr_allele]], operation = "del")
    } else if (event_name == "bfb") {
      selected_chr_allele = state$input_parameters$bfb_allele
      bfb_result <- sim_bfb_left_and_right_sequences(
        parent_sequences[[selected_chr_allele]],
        state$input_parameters$breakpoint_support,
        state$input_parameters$alpha,
        state$input_parameters$beta
      )
    }

    # Base daughter sequences on parent
    l_sequences <- parent_sequences
    r_sequences <- parent_sequences

    if (event_name == "bfb") {
      l_sequences[[selected_chr_allele]] <- bfb_result$l_seq
      r_sequences[[selected_chr_allele]] <- bfb_result$r_seq
    } else {
      if (stats::runif(1) < 0.5) {
        l_sequences[[selected_chr_allele]] <- modified_seq
      } else {
        r_sequences[[selected_chr_allele]] <- modified_seq
      }
    }
  }

  # Create new cell IDs
  l_cell_id <- paste0("cell_", state$next_cell_id)
  state$next_cell_id <- state$next_cell_id + 1
  r_cell_id <- paste0("cell_", state$next_cell_id)
  state$next_cell_id <- state$next_cell_id + 1

  # Mark parent cell as dead
  state$cell_is_alive[cell_idx] <- FALSE

  # Hotspot status for daughters
  hotspot <- state$input_parameters$hotspot
  l_hotspot <- if (!is.null(hotspot)) is_hotspot_gained(l_sequences[[hotspot$chr]], hotspot = hotspot$pos) else FALSE
  r_hotspot <- if (!is.null(hotspot)) is_hotspot_gained(r_sequences[[hotspot$chr]], hotspot = hotspot$pos) else FALSE

  # Store new cell data
  state$cell_ids <- c(state$cell_ids, l_cell_id, r_cell_id)
  state$cell_sequences[[l_cell_id]] <- l_sequences
  state$cell_sequences[[r_cell_id]] <- r_sequences
  state$cell_is_alive <- c(state$cell_is_alive, TRUE, TRUE)
  state$hotspot_status <- c(state$hotspot_status, l_hotspot, r_hotspot)

  # Selection-adjusted rates
  l_birth_rate <- state$input_parameters$birth_rate * (1 + state$input_parameters$positive_selection_rate * l_hotspot)
  l_death_rate <- state$input_parameters$death_rate * (1 + state$input_parameters$negative_selection_rate * l_hotspot)

  r_birth_rate <- state$input_parameters$birth_rate * (1 + state$input_parameters$positive_selection_rate * r_hotspot)
  r_death_rate <- state$input_parameters$death_rate * (1 + state$input_parameters$negative_selection_rate * r_hotspot)

  # Schedule next events
  state$cell_next_event_times <- c(
    state$cell_next_event_times,
    state$time + stats::rexp(1, l_birth_rate + l_death_rate),
    state$time + stats::rexp(1, r_birth_rate + r_death_rate)
  )

  # Update history
  state$history <- state$history %>%
    dplyr::mutate(is_alive = ifelse(.data$cell_id == current_cell_id, F, .data$is_alive))
  state$history <- dplyr::bind_rows(
    state$history,
    dplyr::tibble(
      cell_id = l_cell_id,
      parent_id = current_cell_id,
      bfb_event = (event_name == "bfb"),
      cn_event = ifelse(event_name == "normal", "none", event_name),
      chr_allele = ifelse(event_name == "normal", NA, selected_chr_allele),
      is_alive = TRUE
    ),
    dplyr::tibble(
      cell_id = r_cell_id,
      parent_id = current_cell_id,
      bfb_event = (event_name == "bfb"),
      cn_event = ifelse(event_name == "normal", "none", event_name),
      chr_allele = ifelse(event_name == "normal", NA, selected_chr_allele),
      is_alive = TRUE
    )
  )

  state
}

#' Prepare final results from simulation
#'
#' @param state The simulation state
#'
#' @return A list containing simulation results
prepare_results <- function(state) {
  # Get only alive cells for final_cells output (no subsampling)
  alive_indices <- which(state$cell_is_alive)

  final_cells <- lapply(state$cell_ids[alive_indices], function(id) {
    state$cell_sequences[[id]]
  })
  names(final_cells) <- state$cell_ids[alive_indices]

  # Extract history table
  cell_history <- state$history

  # Add root node
  root_row <- tibble::tibble(
    cell_id = "root",
    parent_id = NA,
    bfb_event = FALSE,
    cn_event = "none",
    chr_allele = NA,
    is_alive = FALSE
  )
  cell_history <- dplyr::bind_rows(root_row, cell_history)

  # Return results
  result <- list(
    cells = final_cells,
    cell_history = cell_history,
    tree = ape::read.tree(text = cell_history_to_newick(cell_history)),
    input_parameters = state$input_parameters
  )

  return(result)
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

  # Mark the cell as dead
  state$cell_is_alive[cell_idx] <- FALSE

  # Update history (optional, to track death in history table)
  state$history <- state$history %>%
    dplyr::filter(.data$cell_id != current_cell_id)

  state
}


#' Check if simulation should continue
#'
#' @param state The simulation state
#'
#' @return Logical indicating whether simulation should continue
continue_simulation <- function(state) {
  alive_count <- sum(state$cell_is_alive)
  c1 <- state$time < state$input_parameters$max_time
  c2 <- alive_count > 0
  c3 <- alive_count < state$input_parameters$max_cells
  c1 && c2 && c3
}

#' Determine the next event to occur
#'
#' @param state The simulation state
#'
#' @return List with time, cell ID, and event type of next event
get_next_event <- function(state) {

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
  hotspot_status = state$hotspot_status[min_event_idx]

  # Calculate modified birth and death rates for this specific cell
  cell_birth_rate <- state$input_parameters$birth_rate * (1 + state$input_parameters$positive_selection_rate * hotspot_status)
  cell_death_rate <- state$input_parameters$death_rate * (1 + state$input_parameters$negative_selection_rate * hotspot_status)

  # Determine event type (birth or death) based on relative rates
  event_probability <- cell_birth_rate / (cell_birth_rate + cell_death_rate)
  event_type <- if (stats::runif(1) < event_probability) "birth" else "death"

  list(
    time = next_event_time,
    cell_id = current_cell_id,
    event_type = event_type,
    cell_idx = next_event_idx
  )
}




# Cell to cn copy data
sequences_to_cndata <- function(sequences, chr_seq_lengths, bin_length) {

  cndata <- lapply(names(sequences), function(cell_id) {
    seqs <- sequences[[cell_id]]
    lapply(names(seqs), function(chr_allele) {
      s <- seqs[[chr_allele]]

      chr <- strsplit(chr_allele, ":")[[1]][1]
      allele <- strsplit(chr_allele, ":")[[1]][2]

      bin_idxs <- 1:chr_seq_lengths[[chr]]
      t <- table(seq2vec(s))[bin_idxs]
      d <- dplyr::tibble(
        cell_id = cell_id,
        bin_idx = bin_idxs,
        allele = allele,
        chr = chr,
        state = as.integer(t[as.character(bin_idxs)])
      )
      d$state[is.na(d$state)] <- 0
      d
    }) %>%
      do.call(dplyr::bind_rows, .)
  }) %>%
    do.call(dplyr::bind_rows, .)

  cndata <- cndata %>%
    dplyr::group_by(.data$cell_id, .data$chr) %>%
    tidyr::pivot_wider(values_from = .data$state, names_from = .data$allele) %>%
    dplyr::mutate(CN = .data$A + .data$B) %>%
    dplyr::mutate(start = (.data$bin_idx - 1) * bin_length + 1, end = .data$bin_idx * bin_length)

  cndata
}



sim_amp_del <- function(sequence, operation = "dup") {
  # Get sequence properties
  L <- get_seq_length(sequence)

  # Ensure there's a sequence to modify
  if (L == 0) {
    return(sequence)
  }

  vec <- seq2vec(sequence)

  find_consecutive_subvectors <- function(vec) {
    diffs <- diff(vec)
    runs <- rle(abs(diffs) == 1)
    lengths <- runs$lengths
    values <- runs$values

    start_indices <- cumsum(c(1, lengths[-length(lengths)]))
    good_starts <- start_indices[values]
    good_lengths <- lengths[values] + 1  # +1 because diff loses one element

    # Build a list of consecutive subvectors
    subvecs <- lapply(seq_along(good_starts), function(i) {
      idx_start <- good_starts[i]
      idx_end <- idx_start + good_lengths[i] - 1
      vec[idx_start:idx_end]
    })

    return(list(subvecs = subvecs, good_starts = good_starts))
  }



  subvectors <- find_consecutive_subvectors(vec)
  subvecs <- subvectors$subvecs
  good_starts <- subvectors$good_starts

  if (length(subvecs) == 0) {
    subvecs = lapply(vec, function(x){x})
    good_starts = 1:length(vec)
    #stop("No consecutive subsequence found!")
  }

  idx <- sample(seq_along(subvecs), 1)
  subvec <- subvecs[[idx]]
  start <- good_starts[idx]

  # Sample a section inside the selected subvector
  if (length(subvec) == 1) {
    event_lims <- c(1, 1)
  } else {
    event_lims <- sort(sample(1:length(subvec), size = 2))
  }

  s <- subvec[event_lims[1]:event_lims[2]]

  absolute_start <- start + event_lims[1] - 1
  absolute_end <- start + event_lims[2] - 1

  # Safe slicing
  before <- if (absolute_start > 1) vec[1:(absolute_start - 1)] else integer(0)
  after <- if (absolute_end < length(vec)) vec[(absolute_end + 1):length(vec)] else integer(0)

  # Apply operation
  if (operation == "dup") {
    new_vec <- c(before, s, s, after)
  } else if (operation == "del") {
    new_vec <- c(before, after)
    if (length(new_vec) == 0) {
      message("skipping deletion because of length zero")
      return(sequence)
    }
  } else {
    stop("operation not recognized!")
  }

  return(vec2seq(new_vec))
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
#'
#' @details Simulate left and right children from a BFB cycle using specified
#'   breakpoint selection distribution
#'
#' @return List containing left and right sequences
sim_bfb_left_and_right_sequences <- function(sequence, support = "uniform", alpha = NULL, beta = NULL) {
  # Calculate the total number of elements in the sequence, similar to the vectorized version
  L = get_seq_length(sequence)
  vec = seq2vec(sequence)
  bps = vec[diff(vec) == 0]

  # Select random breakpoint based on specified distribution
  # Ensure that bp_idx is different from L to obtain a proper bfb cycle
  bp_idx = L
  while (bp_idx %in% c(L, bps)) {
    if (support == "uniform") {
      bp_idx = sample(1:(2*L), 1)
    } else if (support == "beta") {
      if (is.null(alpha) || is.null(beta)) {
        stop("For beta distribution, both alpha and beta parameters must be provided")
      }
      tau = stats::rbeta(1, alpha, beta)
      bp_idx = max(1, round(tau * 2*L))  # Ensure bp_idx is at least 1
    } else {
      stop("Unsupported distribution type. Use 'uniform' or 'beta'.")
    }
  }

  # Initialize left and right sequences
  cut_seqs = cut_sequence(fuse_sequence(sequence), bp_idx)
  l_seq = cut_seqs$left_seq
  r_seq = reverse_sequence(cut_seqs$right_seq)

  # Return the left and right sequences
  return(list(l_seq = l_seq, r_seq = r_seq))
}

reverse_sequence <- function(sequence) {
  # Reverse the order of the intervals and swap start and end, flip direction
  reversed_seq = lapply(seq_along(sequence), function(i) {
    interval = sequence[[length(sequence) - i + 1]]
    list(start = interval$end,
         end = interval$start,
         direction = -interval$direction)
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
      if (interval$direction == 1) {  # Increasing interval
        left_seq <- c(left_seq, list(list(
          start = interval$start,
          end = interval$start + cut_within - 1,
          direction = 1
        )))
        right_seq <- c(right_seq, list(list(
          start = interval$start + cut_within,
          end = interval$end,
          direction = 1
        )))
      } else if (interval$direction == -1) {  # Decreasing interval
        left_seq <- c(left_seq, list(list(
          start = interval$start,
          end = interval$start - cut_within + 1,
          direction = -1
        )))
        right_seq <- c(right_seq, list(list(
          start = interval$start - cut_within,
          end = interval$end,
          direction = -1
        )))
      } else {  # Constant interval
        left_seq <- c(left_seq, list(interval))
        right_seq <- c(right_seq, list(interval))
      }
    }

    # Update the position tracker
    current_length <- current_length + interval_length
  }

  return(list(left_seq = left_seq, right_seq = right_seq))
}

is_hotspot_gained = function(cell, hotspot) {
  hc = get_hotspot_copies(cell, hotspot)
  if (is.nan(hc)) return(NA)
  hc > 1
}


get_hotspot_copies = function(cell, hotspot) {
  if (is.null(hotspot)) return(NaN)

  table_vec = table(seq2vec(cell))
  hotspot = c(hotspot)
  flag = names(table_vec) %in% hotspot
  if (any(flag)) {
    table_vec[names(table_vec) %in% hotspot] %>% mean()
  } else {
    0
  }
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
    living_children <- children[vapply(children, has_living_descendants, logical(1))]

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
    allow_wgd,
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
  valid_breakpoint_supports <- c("uniform", "beta")

  # --- Numeric Parameter Validation ---

  # initial_cells: positive integer
  if (!is.numeric(initial_cells) || initial_cells <= 0 || initial_cells != round(initial_cells)) {
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
    stop("At least one of normal_dup_rate, bfb_prob, amp_rate, or del_rate must be positive")
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
  if (!is.numeric(max_cells) || max_cells <= 0 || max_cells != round(max_cells)) {
    stop("max_cells must be a positive integer")
  }

  # --- Character/String Parameter Validation ---

  # chromosomes: valid chromosome names
  chromosomes <- as.character(chromosomes)
  invalid_chrs <- setdiff(chromosomes, valid_chromosomes)
  if (length(invalid_chrs) > 0) {
    stop(paste("Invalid chromosome(s):", paste(invalid_chrs, collapse = ", "),
               "\nValid chromosomes are:", paste(valid_chromosomes, collapse = ", ")))
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
    stop(paste("breakpoint_support must be one of:", paste(valid_breakpoint_supports, collapse = ", ")))
  }

  # --- Logical Parameter Validation ---

  logical_params <- list(
    allow_wgd = allow_wgd,
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
    stop(paste("hotspot must contain named elements:", paste(required_hotspot_names, collapse = ", ")))
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
      stop("alpha and beta parameters must be provided when breakpoint_support is 'beta'")
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
    stop(paste("bfb_allele chromosome", bfb_chr, "must be included in the chromosomes parameter"))
  }

  # Check that hotspot chromosome is included in chromosomes
  if (!hotspot_chr %in% chromosomes) {
    stop(paste("hotspot chromosome", hotspot_chr, "must be included in the chromosomes parameter"))
  }

  # Warn if birth_rate is much smaller than death_rate
  if (birth_rate > 0 && death_rate > birth_rate * 10) {
    warning("death_rate is much larger than birth_rate - population may quickly go extinct")
  }

  # Warn if max_cells is very large
  if (max_cells > 10000) {
    warning("max_cells is very large - simulation may be slow or use excessive memory")
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

  # Rebuild tree from filtered history
  subsampled_tree <- ape::read.tree(text = cell_history_to_newick(subsampled_cell_history))

  # Create new CNA data for subsampled cells
  chr_seq_lengths <- sim_result$input_parameters$chr_seq_lengths
  bin_length <- sim_result$input_parameters$bin_length
  subsampled_cna_data <- sequences_to_cndata(subsampled_cells, chr_seq_lengths, bin_length)

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
