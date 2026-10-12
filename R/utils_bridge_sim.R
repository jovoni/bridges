
#' Helper functions modified for diploid chromosome modeling

#' Initialize simulation state for diploid chromosomes
#'
#' @param input_parameters List of parameters given as input to bridge_sim function.
#'
#' @return A list containing the initialized simulation state
#' @keywords internal
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
    hotspot_counts         = integer(0),
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

#' Initialize simulation state from pre-evolved cell sequences
#'
#' Used by \code{bridge_sim()} when \code{initial_sequences} is provided (e.g.
#' for serial passage simulations).  Builds the same \code{sim_state} structure
#' as \code{initialize_simulation()} but seeds the population from an existing
#' set of heterogeneous cell sequences rather than a fresh diploid cell.
#'
#' @param sequences Named list of cell sequences in the format returned by
#'   \code{bridge_sim()$cells}: \code{list(cell_id = list("chr:allele" = seq))}.
#' @param input_parameters Full parameter list as built inside \code{bridge_sim()}.
#'
#' @return A simulation state list ready to be passed to \code{bridge_sim_loop_cpp()}.
#' @keywords internal
initialize_from_sequences <- function(sequences, input_parameters) {
  n <- length(sequences)
  if (n == 0L) stop("initialize_from_sequences: sequences list is empty")

  # Rename cells to cell_1 … cell_n for clean per-passage IDs.
  new_ids <- paste0("cell_", seq_len(n))
  sequences <- stats::setNames(sequences, new_ids)

  max_history <- 6L * as.integer(input_parameters$max_cells) + 4L * n + 20L

  state <- list(
    time                   = 0,
    cell_ids               = character(0),
    cell_sequences         = list(),
    cell_next_event_times  = numeric(0),
    hotspot_counts         = integer(0),
    h_cell_id    = character(max_history),
    h_parent_id  = character(max_history),
    h_bfb_event  = logical(max_history),
    h_wgd_event  = logical(max_history),
    h_cn_event   = character(max_history),
    h_chr_allele = character(max_history),
    h_n          = 0L
  )
  state$input_parameters <- input_parameters
  state$next_cell_id <- n + 1L

  p       <- input_parameters
  hotspot <- p$hotspot

  for (i in seq_len(n)) {
    cid  <- new_ids[i]
    seqs <- sequences[[cid]]

    # Determine hotspot copy count from actual copy number in the BFB allele.
    if (!is.null(hotspot)) {
      hc <- get_hotspot_copies(seqs[[hotspot$chr]], hotspot = hotspot$pos)
      if (is.nan(hc)) hc <- 0L
    } else {
      hc <- 0L
    }

    sel_mult      <- compute_selection_multiplier(hc, p$selection_type, p$saturation_K)
    birth_rate    <- p$birth_rate * (1 + p$positive_selection_rate * sel_mult)
    death_rate    <- p$death_rate * (1 + p$negative_selection_rate * sel_mult)
    combined_rate <- birth_rate + death_rate

    state$cell_ids                   <- c(state$cell_ids, cid)
    state$cell_sequences[[cid]]      <- seqs
    state$cell_next_event_times      <- c(state$cell_next_event_times,
                                          stats::rexp(1, combined_rate))
    state$hotspot_counts             <- c(state$hotspot_counts, as.integer(hc))

    # Record each founder as a child of "root" so the lineage (and the tree
    # built from it) is connected.
    idx <- state$h_n + 1L
    state$h_n <- idx
    state$h_cell_id[idx]    <- cid
    state$h_parent_id[idx]  <- "root"
    state$h_bfb_event[idx]  <- FALSE
    state$h_wgd_event[idx]  <- FALSE
    state$h_cn_event[idx]   <- "none"
    state$h_chr_allele[idx] <- NA_character_
  }

  state
}

#' Create initial chromosome sequences for all alleles
#'
#' @param input_parameters List of input parameters
#'
#' @return Named list of initial sequences for each chromosome allele
#' @keywords internal
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
#' @keywords internal
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
      left_hc  <- get_hotspot_copies(left_cell[[hotspot$chr]],  hotspot = hotspot$pos)
      right_hc <- get_hotspot_copies(right_cell[[hotspot$chr]], hotspot = hotspot$pos)
      if (is.nan(left_hc))  left_hc  <- 0L
      if (is.nan(right_hc)) right_hc <- 0L
    } else {
      left_hc <- right_hc <- 0L
    }

    sel_type <- state$input_parameters$selection_type
    l_mult   <- compute_selection_multiplier(left_hc,  sel_type, state$input_parameters$saturation_K)
    r_mult   <- compute_selection_multiplier(right_hc, sel_type, state$input_parameters$saturation_K)

    l_birth_rate <- state$input_parameters$birth_rate *
      (1 + state$input_parameters$positive_selection_rate * l_mult)
    l_death_rate <- state$input_parameters$death_rate *
      (1 + state$input_parameters$negative_selection_rate * l_mult)

    r_birth_rate <- state$input_parameters$birth_rate *
      (1 + state$input_parameters$positive_selection_rate * r_mult)
    r_death_rate <- state$input_parameters$death_rate *
      (1 + state$input_parameters$negative_selection_rate * r_mult)

    l_combined_rate <- l_birth_rate + l_death_rate
    r_combined_rate <- r_birth_rate + r_death_rate

    # Update active arrays (only alive cells)
    state$cell_ids               <- c(state$cell_ids, left_cell_id, right_cell_id)
    state$cell_sequences[[left_cell_id]]  <- left_cell
    state$cell_sequences[[right_cell_id]] <- right_cell
    state$hotspot_counts         <- c(state$hotspot_counts, as.integer(left_hc), as.integer(right_hc))
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
#' @keywords internal
initialize_without_bfb <- function(state, initial_sequences) {
  for (i in 1:state$input_parameters$initial_cells) {
    cell_id <- paste0("cell_", state$next_cell_id)
    state$next_cell_id <- state$next_cell_id + 1L

    # Hotspot logic
    if (!is.null(state$input_parameters$hotspot)) {
      hotspot <- state$input_parameters$hotspot
      hc <- get_hotspot_copies(
        initial_sequences[[hotspot$chr]],
        hotspot = hotspot$pos
      )
      if (is.nan(hc)) hc <- 0L
    } else {
      hc <- 0L
    }

    # Birth and death rate with selection
    sel_mult   <- compute_selection_multiplier(hc, state$input_parameters$selection_type, state$input_parameters$saturation_K)
    birth_rate <- state$input_parameters$birth_rate *
      (1 + state$input_parameters$positive_selection_rate * sel_mult)
    death_rate <- state$input_parameters$death_rate *
      (1 + state$input_parameters$negative_selection_rate * sel_mult)

    combined_rate  <- birth_rate + death_rate
    next_event_time <- state$time + stats::rexp(1, combined_rate)

    # Update active arrays
    state$cell_ids                    <- c(state$cell_ids, cell_id)
    state$cell_sequences[[cell_id]]   <- initial_sequences
    state$cell_next_event_times       <- c(state$cell_next_event_times, next_event_time)
    state$hotspot_counts              <- c(state$hotspot_counts, as.integer(hc))

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

#' Prepare final results from simulation
#'
#' Assembles the final return list from a completed (or in-progress)
#' simulation state: cell sequences, cell history, optionally the
#' reconstructed phylogeny, and summary fields.
#'
#' @param state The simulation state to finalize
#' @param return_phylo Logical; whether to reconstruct and include the
#'   phylogenetic tree (via \code{build_phylo_from_lineage}). Default
#'   TRUE.
#'
#' @return A named list with the final simulation state and derived
#'   outputs (see \code{bridge_sim()}'s return value for the full
#'   structure).
#' @keywords internal
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
    elapsed_time     = state$time,
    n_alive          = if (!is.null(state$n_alive_total)) state$n_alive_total else length(final_cells),
    input_parameters = state$input_parameters
  )

  return(result)
}

#' Simulate Breakage-Fusion-Bridge (BFB) cycle for both daughters
#'
#' Generates the left and right daughter sequences resulting from one
#' BFB cycle applied to \code{sequence}, via the C++ engine
#' (\code{sim_bfb_cpp}).
#'
#' @param sequence Interval-encoded chromosome sequence to break.
#' @param support Breakpoint distribution: \code{"uniform"},
#'   \code{"beta"}, or \code{"custom"}.
#' @param alpha,beta Beta-distribution shape parameters, used when
#'   \code{support = "beta"}.
#' @param custom_breakpoints Integer vector of bin indices eligible as
#'   breakpoints, required when \code{support = "custom"}.
#'
#' @return A list with \code{l_seq} and \code{r_seq}, the two daughter
#'   sequences.
#'
#' @keywords internal
sim_bfb_left_and_right_sequences <- function(
    sequence,
    support = "uniform",
    alpha = NULL,
    beta = NULL,
    custom_breakpoints = NULL
) {
  if (support == "custom" && length(custom_breakpoints) == 0) {
    stop("For custom distribution, custom_breakpoints must be provided and non-empty")
  }
  sim_bfb_cpp(
    seq_list           = sequence,
    support            = support,
    alpha              = if (is.null(alpha)) NA_real_ else alpha,
    beta_param         = if (is.null(beta))  NA_real_ else beta,
    custom_breakpoints = as.integer(if (is.null(custom_breakpoints)) integer(0) else custom_breakpoints)
  )
}

get_hotspot_copies <- function(cell, hotspot) {
  if (is.null(hotspot)) return(NaN)
  hotspot_copies_cpp(cell, hotspot)
}

# Translate raw hotspot copy count to a selection multiplier.
# "constant": 1 if copies > 1, else 0 (backward-compatible).
# "linear":   the raw copy count.
compute_selection_multiplier <- function(hotspot_count, selection_type, saturation_K = 10) {
  if (selection_type == "linear") {
    return(as.numeric(max(0L, as.integer(hotspot_count) - 1L)))
  }
  if (selection_type == "saturation") {
    n <- as.numeric(hotspot_count)
    if (n <= 1) return(0.0)
    return((n - 1) / ((n - 1) + saturation_K))
  }
  as.numeric(hotspot_count > 1)
}

validate_bridge_sim_params <- function(
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

  # lambda / rate: event count and event-length parameters
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) {
    stop("lambda must be a single non-negative number")
  }
  if (!is.numeric(rate) || length(rate) != 1 || rate <= 0) {
    stop("rate must be a single positive number")
  }

  # WGD parameters
  if (!is.numeric(wgd_available) || length(wgd_available) != 1 ||
      wgd_available < 0 || wgd_available != round(wgd_available)) {
    stop("wgd_available must be a single non-negative integer")
  }
  if (!is.numeric(wgd_probability) || length(wgd_probability) != 1 ||
      wgd_probability < 0 || wgd_probability > 1) {
    stop("wgd_probability must be a single number in [0, 1]")
  }

  # subsample: NULL or positive integer
  if (!is.null(subsample) &&
      (!is.numeric(subsample) || length(subsample) != 1 ||
       subsample < 1 || subsample != round(subsample))) {
    stop("subsample must be NULL or a single positive integer")
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

  # selection_type: must be "constant", "linear", or "saturation"
  if (!is.character(selection_type) || length(selection_type) != 1 ||
      !selection_type %in% c("constant", "linear", "saturation")) {
    stop('selection_type must be "constant", "linear", or "saturation"')
  }

  if (selection_type == "saturation") {
    if (!is.numeric(saturation_K) || length(saturation_K) != 1 || saturation_K <= 0) {
      stop("saturation_K must be a single positive number when selection_type = 'saturation'")
    }
  }

  if (breakpoint_support == "custom" && length(custom_breakpoints) == 0) {
    stop("custom_breakpoints must be a non-empty integer vector when breakpoint_support is 'custom'")
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

  # Warn if the hotspot lies beyond the end of its chromosome: selection
  # would then silently never act.
  if (hotspot$pos > chr_seq_lengths[[hotspot_chr]]) {
    warning(sprintf(
      "hotspot$pos (%s) is beyond the last bin of chromosome %s (%d bins at this bin_length); selection will have no effect",
      hotspot$pos, hotspot_chr, as.integer(chr_seq_lengths[[hotspot_chr]])
    ))
  }

  # Warn if birth_rate is much smaller than death_rate
  if (birth_rate > 0 && death_rate > birth_rate * 10) {
    warning(
      "death_rate is much larger than birth_rate and population may quickly go extinct"
    )
  }

  # Warn if max_cells is very large
  if (max_cells > 10000) {
    warning(
      "max_cells is very large and simulation may be slow or use excessive memory"
    )
  }

  # All validations passed
  return(invisible(NULL))
}


sequences_to_cndata <- function(sequences, chr_seq_lengths, bin_length) {
  # Pre-compute total rows needed
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

  # Build the edge matrix: one row per (parent, child) pair
  has_parent <- !is.na(cell_history$parent_id)
  edge_matrix <- cbind(
    unname(cell_to_node[cell_history$parent_id[has_parent]]),
    unname(cell_to_node[cell_history$cell_id[has_parent]])
  )
  storage.mode(edge_matrix) <- "integer"

  # Create edge lengths (all set to 1 for simplicity, can be modified)
  edge_lengths <- rep(1, nrow(edge_matrix))

  # Create tip labels
  tip_labels <- terminal_cells

  # Create the phylo object
  phylo_tree <- list(
    edge = edge_matrix,
    edge.length = edge_lengths,
    tip.label = tip_labels,
    Nnode = n_internal
  )

  class(phylo_tree) <- "phylo"

  # Validate the tree structure
  if (!is.null(phylo_tree)) {
    tryCatch({
      # Basic validation to check if tree is valid
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
