# Fixed version of Gillespie algorithm with phylogenetic tree tracking
gillespie_sim_with_phylogeny <- function(
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

  # Track all cells with their complete lifecycle information
  all_cells <- data.frame(
    cell_id = integer(0),
    birth_time = numeric(0),
    death_time = numeric(0),
    status = character(0),
    stringsAsFactors = FALSE
  )

  # Track active cells (current population)
  active_cell_ids <- integer(0)
  next_cell_id <- 1

  # Track parent-child relationships and mutation events for phylogeny
  phylogeny_edges <- data.frame(
    parent = integer(0),
    child = integer(0),
    time = numeric(0),
    event_type = character(0),
    stringsAsFactors = FALSE
  )

  # Create the initial cell with its sequence
  initial_sequence = vec2seq(1:initial_sequence_length)

  # Apply BFB to all initial cells if first_round_of_bfb is TRUE
  if (first_round_of_bfb) {
    # Create temporary storage for new cells after BFB
    new_cells <- list()
    new_cell_ids <- integer(0)

    for (i in 1:initial_cells) {
      # Apply BFB to create daughter sequences
      daughter_seqs <- sim_bfb_left_and_right_sequences(
        initial_sequence,
        breakpoint_support,
        alpha,
        beta
      )

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

      # Record cell IDs for active population
      left_id <- next_cell_id
      right_id <- next_cell_id + 1
      new_cell_ids <- c(new_cell_ids, left_id, right_id)

      # Add to all cells tracking dataframe
      all_cells <- rbind(all_cells, data.frame(
        cell_id = left_id,
        birth_time = time,
        death_time = NA,
        status = "alive",
        stringsAsFactors = FALSE
      ))

      all_cells <- rbind(all_cells, data.frame(
        cell_id = right_id,
        birth_time = time,
        death_time = NA,
        status = "alive",
        stringsAsFactors = FALSE
      ))

      # Record phylogeny - both cells originate from a virtual "root" with ID 0
      phylogeny_edges <- rbind(phylogeny_edges,
                               data.frame(
                                 parent = 0,
                                 child = left_id,
                                 time = time,
                                 event_type = "initial_bfb_left",
                                 stringsAsFactors = FALSE
                               )
      )

      phylogeny_edges <- rbind(phylogeny_edges,
                               data.frame(
                                 parent = 0,
                                 child = right_id,
                                 time = time,
                                 event_type = "initial_bfb_right",
                                 stringsAsFactors = FALSE
                               )
      )

      next_cell_id <- next_cell_id + 2
    }

    # Replace the original cells with new BFB-processed cells
    cells <- new_cells
    active_cell_ids <- new_cell_ids

  } else {
    # Original initialization without BFB
    for (i in 1:initial_cells) {
      cells[[i]] <- initial_sequence

      # Record cell ID for active population
      active_cell_ids <- c(active_cell_ids, next_cell_id)

      # Add to all cells tracking dataframe
      all_cells <- rbind(all_cells, data.frame(
        cell_id = next_cell_id,
        birth_time = time,
        death_time = NA,
        status = "alive",
        stringsAsFactors = FALSE
      ))

      # Record phylogeny - all initial cells originate from root
      phylogeny_edges <- rbind(phylogeny_edges,
                               data.frame(
                                 parent = 0,
                                 child = next_cell_id,
                                 time = time,
                                 event_type = "initial",
                                 stringsAsFactors = FALSE
                               )
      )

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
      parent_id <- active_cell_ids[parent_idx]

      # Determine if BFB should occur
      # Force BFB for first n births, otherwise use probability
      force_bfb <- birth_count <= first_n_bfb_cycles
      random_bfb <- stats::runif(1) < bfb_prob

      if (force_bfb || random_bfb) {
        # Apply the mutation function to get two daughter sequences
        daughter_seqs <- sim_bfb_left_and_right_sequences(
          parent_seq,
          breakpoint_support,
          alpha,
          beta
        )

        if (any(is.na(daughter_seqs$l_seq))) {
          stop("NA values in left daughter sequence")
        }
        if (any(is.na(daughter_seqs$r_seq))) {
          stop("NA values in right daughter sequence")
        }

        # Replace parent with left daughter (sequence changes but ID stays the same)
        cells[[parent_idx]] <- daughter_seqs$l_seq

        # Add right daughter as new cell
        cells[[length(cells) + 1]] <- daughter_seqs$r_seq

        # Add new right daughter to active cells
        active_cell_ids <- c(active_cell_ids, next_cell_id)

        # Add right daughter to all cells tracking dataframe
        all_cells <- rbind(all_cells, data.frame(
          cell_id = next_cell_id,
          birth_time = time,
          death_time = NA,
          status = "alive",
          stringsAsFactors = FALSE
        ))

        # Record phylogeny for both daughters
        # Left daughter maintains parent ID but we record the event
        phylogeny_edges <- rbind(phylogeny_edges,
                                 data.frame(
                                   parent = parent_id,
                                   child = parent_id,
                                   time = time,
                                   event_type = "bfb_left_daughter",
                                   stringsAsFactors = FALSE
                                 )
        )

        # Right daughter is completely new
        phylogeny_edges <- rbind(phylogeny_edges,
                                 data.frame(
                                   parent = parent_id,
                                   child = next_cell_id,
                                   time = time,
                                   event_type = "bfb_right_daughter",
                                   stringsAsFactors = FALSE
                                 )
        )

        next_cell_id <- next_cell_id + 1

      } else {
        # Normal replication without mutation
        # Parent cell remains unchanged
        # Add identical daughter cell
        cells[[length(cells) + 1]] <- parent_seq

        # Add new daughter to active cells
        active_cell_ids <- c(active_cell_ids, next_cell_id)

        # Add new daughter to all cells tracking dataframe
        all_cells <- rbind(all_cells, data.frame(
          cell_id = next_cell_id,
          birth_time = time,
          death_time = NA,
          status = "alive",
          stringsAsFactors = FALSE
        ))

        # Record phylogeny for simple replication
        phylogeny_edges <- rbind(phylogeny_edges,
                                 data.frame(
                                   parent = parent_id,
                                   child = next_cell_id,
                                   time = time,
                                   event_type = "mitosis",
                                   stringsAsFactors = FALSE
                                 )
        )

        next_cell_id <- next_cell_id + 1
      }

    } else {
      # Death event - remove a random cell
      death_count <- death_count + 1

      if (length(cells) > 0) {
        cell_to_remove <- sample(1:length(cells), 1)
        dying_cell_id <- active_cell_ids[cell_to_remove]

        # Record death in all cells dataframe
        all_cells$death_time[all_cells$cell_id == dying_cell_id] <- time
        all_cells$status[all_cells$cell_id == dying_cell_id] <- "dead"

        # Record death event in phylogeny
        phylogeny_edges <- rbind(phylogeny_edges,
                                 data.frame(
                                   parent = dying_cell_id,
                                   child = -1 * dying_cell_id,
                                   time = time,
                                   event_type = "death",
                                   stringsAsFactors = FALSE
                                 )
        )

        # Remove the cell from active population
        cells <- cells[-cell_to_remove]
        active_cell_ids <- active_cell_ids[-cell_to_remove]
      }
    }

    # Calculate sequence statistics for history if needed
    if (time %% 5 < dt || length(cells) %% 100 == 0) {  # Record every ~5 time units or every 100 cells
      seq_lengths <- sapply(cells, get_seq_length)
      seq_hashes <- sapply(cells, function(seq) paste(unlist(seq), collapse=","))
      unique_seqs <- length(unique(seq_hashes))

      # Record current state
      history <- rbind(history, data.frame(
        time = time,
        cell_count = length(cells),
        unique_sequences = unique_seqs,
        avg_seq_length = ifelse(length(seq_lengths) > 0, mean(seq_lengths), 0),
        birth_count = birth_count,
        death_count = death_count
      ))
    }
  }

  # Prepare final sequence distribution data
  seq_lengths <- sapply(cells, get_seq_length)
  seq_hashes <- sapply(cells, function(seq) paste(unlist(seq), collapse=","))
  if (length(seq_hashes) > 0) {
    seq_counts <- table(seq_hashes)
  } else {
    seq_counts = NULL
  }

  # Calculate lifetimes for all cells
  # For cells still alive, use the final simulation time
  all_cells$death_time[all_cells$status == "alive"] <- time
  all_cells$lifetime <- all_cells$death_time - all_cells$birth_time

  # Create cell lifetimes data frame (now using our complete tracking dataframe)
  cell_lifetimes <- all_cells

  # Return results
  return(list(
    history = history,
    final_cells = cells,
    active_cell_ids = active_cell_ids,  # IDs of cells alive at end of simulation
    all_cells = all_cells,  # Complete cell tracking dataframe
    sequence_distribution = seq_counts,
    cell_lifetimes = cell_lifetimes,
    phylogeny = phylogeny_edges,
    final_time = time,
    birth_count = birth_count,
    death_count = death_count
  ))
}
