
# Gillespie algorithm for birth-death process with sequence tracking
gillespie_sim <- function(
    initial_cells = 5,
    initial_sequence_length = 100,
    birth_rate = 0.1,
    death_rate = 0.001,
    bfb_prob = 0.01,
    max_time = 50,
    max_cells = 1000,
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

  # Track cell creation times and IDs
  cell_ids <- integer(0)
  cell_birth_times <- numeric(0)
  cell_death_times <- numeric(0)
  next_cell_id <- 1

  # Create the initial cell with its sequence
  initial_sequence = vec2seq(1:initial_sequence_length)

  # Apply BFB to all initial cells if first_round_of_bfb is TRUE
  if (first_round_of_bfb) {
    # Create temporary storage for new cells after BFB
    new_cells <- list()
    new_cell_ids <- integer(0)
    new_birth_times <- numeric(0)
    new_death_times <- numeric(0)

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

      # Add both daughters to new cells list
      new_cells[[length(new_cells) + 1]] <- daughter_seqs$l_seq
      new_cells[[length(new_cells) + 1]] <- daughter_seqs$r_seq

      # Record cell creation data for both new cells
      new_cell_ids <- c(new_cell_ids, next_cell_id, next_cell_id + 1)
      new_birth_times <- c(new_birth_times, time, time)
      new_death_times <- c(new_death_times, NA, NA)  # Not dead yet
      next_cell_id <- next_cell_id + 2
    }

    # Replace the original cells with new BFB-processed cells
    cells <- new_cells
    cell_ids <- new_cell_ids
    cell_birth_times <- new_birth_times
    cell_death_times <- new_death_times
  } else {
    # Original initialization without BFB
    for (i in 1:initial_cells) {
      cells[[i]] <- initial_sequence

      # Record cell creation data
      cell_ids <- c(cell_ids, next_cell_id)
      cell_birth_times <- c(cell_birth_times, time)
      cell_death_times <- c(cell_death_times, NA)  # Not dead yet
      next_cell_id <- next_cell_id + 1
    }
  }

  # Record history
  # Calculate initial sequence statistics
  seq_lengths <- sapply(cells, get_seq_length)

  seq_hashes <- sapply(cells, function(seq) paste(unlist(seq), collapse=","))
  unique_seqs <- length(unique(seq_hashes))

  history <- data.frame(
    time = time,
    cell_count = length(cells),
    unique_sequences = unique_seqs,
    avg_seq_length = mean(seq_lengths),
    birth_count = 0,
    death_count = 0
  )

  death_count <- 0  # Track number of death events

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

      # Select parent cell
      parent_idx <- sample(1:length(cells), 1)
      parent_seq <- cells[[parent_idx]]
      parent_id <- cell_ids[parent_idx]

      # Determine if BFB should occur
      # Force BFB for first n births, otherwise use probability
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

        # Replace parent with left daughter
        cells[[parent_idx]] <- daughter_seqs$l_seq

        # Add right daughter as new cell
        cells[[length(cells) + 1]] <- daughter_seqs$r_seq

        # Record cell creation data for the new cell
        cell_ids <- c(cell_ids, next_cell_id)
        cell_birth_times <- c(cell_birth_times, time)
        cell_death_times <- c(cell_death_times, NA)  # Not dead yet
        next_cell_id <- next_cell_id + 1
      } else {
        # Normal replication without mutation
        # Parent cell remains unchanged
        # Add identical daughter cell
        cells[[length(cells) + 1]] <- parent_seq

        # Record cell creation data for the new cell
        cell_ids <- c(cell_ids, next_cell_id)
        cell_birth_times <- c(cell_birth_times, time)
        cell_death_times <- c(cell_death_times, NA)  # Not dead yet
        next_cell_id <- next_cell_id + 1
      }

    } else {
      # Death event - remove a random cell
      death_count <- death_count + 1

      if (length(cells) > 0) {
        cell_to_remove <- sample(1:length(cells), 1)

        # Record death time before removing
        cell_death_times[which(cell_ids == cell_ids[cell_to_remove])] <- time

        # Remove the cell
        cells <- cells[-cell_to_remove]
        cell_ids <- cell_ids[-cell_to_remove]
      }
    }

    # Calculate sequence statistics
    # seq_lengths <- sapply(cells, length)
    # seq_hashes <- sapply(cells, function(seq) paste(seq, collapse=","))
    # unique_seqs <- length(unique(seq_hashes))

    # Record current state
    # history <- rbind(history, data.frame(
    #   time = time,
    #   cell_count = length(cells),
    #   unique_sequences = unique_seqs,
    #   avg_seq_length = ifelse(length(seq_lengths) > 0, mean(seq_lengths), 0),
    #   birth_count = birth_count,
    #   death_count = death_count
    # ))
  }

  # Prepare sequence distribution data
  seq_hashes <- sapply(cells, function(seq) paste(unlist(seq), collapse=","))
  if (length(seq_hashes) > 0) {
    seq_counts <- table(seq_hashes)
  } else {
    seq_counts = NULL
  }

  # Create cell lifetime data
  # For cells still alive, use the final simulation time
  cell_death_times[is.na(cell_death_times)] <- time

  cell_lifetimes <- data.frame(
    cell_id = seq_len(next_cell_id - 1),
    birth_time = c(cell_birth_times, cell_birth_times[match(seq_len(next_cell_id - 1), cell_ids)]),
    death_time = c(cell_death_times, cell_death_times[match(seq_len(next_cell_id - 1), cell_ids)]),
    lifetime = c(cell_death_times - cell_birth_times,
                 (cell_death_times - cell_birth_times)[match(seq_len(next_cell_id - 1), cell_ids)]),
    is_alive = c(is.na(cell_death_times), is.na(cell_death_times)[match(seq_len(next_cell_id - 1), cell_ids)])
  )

  # Fix NA values and create final status
  cell_lifetimes$is_alive <- cell_lifetimes$is_alive
  cell_lifetimes$status <- ifelse(cell_lifetimes$is_alive, "alive", "dead")
  cell_lifetimes$lifetime <- ifelse(cell_lifetimes$is_alive, time - cell_lifetimes$birth_time,
                                    cell_lifetimes$death_time - cell_lifetimes$birth_time)

  # Return results
  return(list(
    #history = history,
    final_cells = cells,
    sequence_distribution = seq_counts,
    cell_lifetimes = cell_lifetimes,
    final_time = time,
    birth_count = birth_count,
    death_count = death_count
  ))
}
