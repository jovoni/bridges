
# ── Passage-series simulation ─────────────────────────────────────────────────
#
# Self-consistency test: do simulations with more cells capture more complex
# CNA landscapes, mimicking cells that have been passaged more times?
#
# Key insight:
#   - max_cells = the depth control — to produce N cells, the simulation must
#     run until N cells are alive, accumulating CNA events throughout.
#   - elapsed_time (recorded from sim$elapsed_time) = the actual evolutionary
#     clock when the simulation stopped.
#   - sample_cells = observation window — subsample to this many cells at every
#     depth so comparisons are fair.
#   - expected_doublings = elapsed_time × birth_rate / log(2) — a biologically
#     calibrated x-axis (number of cell doublings since founding).
#
# For a self-consistent simulator, all complexity metrics should increase
# monotonically with elapsed_time (= passage depth).
#
# Workflow:
#   series <- simulate_passage_series(c("p1"=50, "p5"=500, "p10"=5000), ...)
#   stats  <- summarise_passage_series(series)
#   plot_passage_series(stats)


#' Run simulations at a grid of evolutionary depths
#'
#' For each cell-count target in \code{passage_cell_counts}, runs
#' \code{n_replicates} independent clonal evolutions (via
#' \code{simulate_clone()}) with \code{max_time = Inf} so that only
#' \code{max_cells} determines when the simulation stops.  At the end of each
#' run, \code{sample_cells} cells are randomly drawn so that all depths are
#' compared at the same observation scale.
#'
#' The actual elapsed simulation time (from \code{sim\$elapsed_time}) is
#' stored alongside each result and becomes the x-axis in
#' \code{plot_passage_series()}.
#'
#' @param passage_cell_counts Named numeric vector of \code{max_cells} targets,
#'   one per "passage depth".  E.g.
#'   \code{c("early"=50, "mid"=500, "late"=5000)}.  Names become the passage
#'   labels in plots.
#' @param chromosome    Chromosome to simulate (character, e.g. \code{"1"}).
#' @param allele        Allele to track (\code{"A"}, \code{"B"}, or
#'   \code{"CN"}).
#' @param hotspot_pos   1-based bin index of the position of interest.
#' @param params        Parameter list from \code{clonal_params()}.
#' @param sample_cells  Number of cells to draw from each simulation for fair
#'   comparison.  Replicates with fewer alive cells are kept in full and
#'   flagged with a warning.  Default 50.
#' @param n_replicates  Independent runs per depth.  Default 10.
#' @param n_cores       Parallel cores via \code{parallel::mclapply}.
#'   Default 1.
#' @param bin_length    Bin size in bp.  Default 1 Mb.
#'
#' @return A list of named lists, one per \code{(passage_label, replicate)}
#'   combination.  Each element contains \code{passage_label},
#'   \code{passage_cells}, \code{replicate}, \code{cna_matrix},
#'   \code{hotspot_pos}, \code{n_alive}, \code{elapsed_time},
#'   \code{birth_rate}, and \code{failed}.
#'
#' @examples
#' \dontrun{
#' series <- simulate_passage_series(
#'   passage_cell_counts = c(early = 50, mid = 200, late = 600),
#'   chromosome   = "7",
#'   allele       = "A",
#'   hotspot_pos  = 60,
#'   n_replicates = 5,
#'   sample_cells = 30
#' )
#' stats  <- summarise_passage_series(series)
#' plot_passage_series(stats)
#' }
#'
#' @export
simulate_passage_series <- function(
  passage_cell_counts,
  chromosome,
  allele,
  hotspot_pos,
  params       = clonal_params(),
  sample_cells = 50L,
  n_replicates = 10L,
  n_cores      = 1L,
  bin_length   = 1e6
) {
  stopifnot(is.numeric(passage_cell_counts), length(passage_cell_counts) >= 1)
  if (is.null(names(passage_cell_counts)))
    names(passage_cell_counts) <- paste0("d", seq_along(passage_cell_counts))

  jobs <- do.call(rbind, lapply(seq_along(passage_cell_counts), function(di) {
    data.frame(
      passage_label  = names(passage_cell_counts)[di],
      passage_cells  = passage_cell_counts[[di]],
      replicate      = seq_len(n_replicates),
      stringsAsFactors = FALSE
    )
  }))

  run_one <- function(j) {
    tryCatch({
      res <- simulate_clone(
        chromosome   = chromosome,
        allele       = allele,
        hotspot_pos  = hotspot_pos,
        params       = params,
        max_cells    = jobs$passage_cells[j],
        max_time     = 1e15,   # only max_cells stops the simulation
        sample_cells = sample_cells,
        bin_length   = bin_length
      )
      list(
        passage_label = jobs$passage_label[j],
        passage_cells = jobs$passage_cells[j],
        replicate     = jobs$replicate[j],
        cna_matrix    = res$cna_matrix,
        hotspot_pos   = hotspot_pos,
        n_alive       = res$n_alive,
        elapsed_time  = res$sim$elapsed_time,
        birth_rate    = params$birth_rate,
        failed        = FALSE
      )
    }, error = function(e) {
      list(
        passage_label = jobs$passage_label[j],
        passage_cells = jobs$passage_cells[j],
        replicate     = jobs$replicate[j],
        cna_matrix    = matrix(integer(0), nrow = 0, ncol = 0),
        hotspot_pos   = hotspot_pos,
        n_alive       = 0L,
        elapsed_time  = NA_real_,
        birth_rate    = params$birth_rate,
        failed        = TRUE
      )
    })
  }

  if (n_cores > 1L) {
    parallel::mclapply(seq_len(nrow(jobs)), run_one, mc.cores = n_cores)
  } else {
    lapply(seq_len(nrow(jobs)), run_one)
  }
}


#' Summarise a passage series into a tidy tibble
#'
#' Applies \code{summarise_clone()} to every element of
#' \code{simulate_passage_series()} and returns a long-format tibble suitable
#' for \code{plot_passage_series()}.
#'
#' A \code{doublings} column is added as a biologically calibrated x-axis:
#' \deqn{\text{doublings} = \text{elapsed\_time} \times \text{birth\_rate} / \ln 2}
#'
#' @param series     Output of \code{simulate_passage_series()}.
#' @param base_value Baseline copy number: 1 for allele-specific, 2 for CN.
#'
#' @return A tibble with columns \code{passage_label}, \code{passage_cells},
#'   \code{elapsed_time}, \code{doublings}, \code{replicate}, \code{stat},
#'   \code{value}, and \code{failed}.
#'
#' @examples
#' \dontrun{
#' series <- simulate_passage_series(
#'   c(early = 50, mid = 200, late = 600), chromosome = "7",
#'   allele = "A", hotspot_pos = 60, n_replicates = 5
#' )
#' stats <- summarise_passage_series(series)
#' head(stats)
#' }
#'
#' @export
summarise_passage_series <- function(series, base_value = 1) {
  rows <- lapply(series, function(el) {
    stats_vec <- if (el$failed || nrow(el$cna_matrix) == 0) {
      c(
        hotspot_mean_cn      = NA_real_,
        hotspot_max_cn       = NA_real_,
        hotspot_fraction_amp = NA_real_,
        hotspot_cn_var       = NA_real_,
        max_cn_global        = NA_real_,
        mean_breakpoints     = NA_real_,
        n_alive              = NA_real_,
        n_cells              = NA_real_
      )
    } else {
      summarise_clone(
        cna_matrix  = el$cna_matrix,
        hotspot_col = el$hotspot_pos,
        base_value  = base_value,
        n_alive     = el$n_alive
      )
    }

    doublings <- if (!is.na(el$elapsed_time))
      el$elapsed_time * el$birth_rate / log(2) else NA_real_

    dplyr::tibble(
      passage_label = el$passage_label,
      passage_cells = el$passage_cells,
      elapsed_time  = el$elapsed_time,
      doublings     = doublings,
      replicate     = el$replicate,
      stat          = names(stats_vec),
      value         = unname(as.numeric(stats_vec)),
      failed        = el$failed
    )
  })

  dplyr::bind_rows(rows)
}


# ── Serial passage simulation ──────────────────────────────────────────────────
#
# Models biological serial passaging: grow to N cells, subsample N_sub, repeat.
# Each bottleneck introduces drift and allows clonal sweeps to accumulate
# across passages — a richer model than simply running a longer simulation.
#
# Workflow:
#   df <- simulate_serial_passages(params, chromosome = "8", allele = "A",
#                                  hotspot_pos = 139, N = 500, N_sub = 50, K = 10)
#   plot_serial_passages(df)   # wide tibble — no summarise step needed


#' Simulate serial passages and return summary statistics
#'
#' Implements biological serial passaging: a single cell grows to \code{N}
#' cells, \code{N_sub} cells are randomly sampled (the bottleneck), and those
#' cells seed the next passage.  This is repeated \code{K} times.  At each
#' passage listed in \code{passages_to_keep}, summary statistics are computed
#' immediately and the CNA matrix is discarded.
#'
#' The return format mirrors \code{compare_clonal_params()}: one row per
#' \code{(passage, replicate)} with one column per summary statistic.
#'
#' Each replicate runs all \code{K} passages sequentially.  Replicates are
#' run in parallel when \code{n_cores > 1}.
#'
#' @param params          Parameter list from \code{clonal_params()}.
#' @param chromosome      Chromosome to simulate (character, e.g. \code{"8"}).
#' @param allele          Allele to track (\code{"A"} or \code{"B"}).
#' @param hotspot_pos     1-based bin index of the position of interest.
#' @param N               Target population size that ends each passage
#'   (\code{max_cells}).
#' @param N_sub           Number of cells sampled between passages (the
#'   bottleneck).  Must be \eqn{\leq N}.
#' @param K               Total number of passages to simulate.
#' @param passages_to_keep Integer vector of passage numbers to record.
#'   Defaults to all passages \code{1:K}.
#' @param n_replicates    Number of independent replicate series.  Default 1.
#' @param n_cores         Parallel cores via \code{parallel::mclapply}.
#'   Default 1.
#' @param bin_length      Bin size in bp.  Default 1 Mb.
#'
#' @return A \code{tibble} with one row per \code{(passage, replicate)}
#'   combination in \code{passages_to_keep}.  Columns: \code{passage},
#'   \code{replicate}, \code{hotspot_mean_cn}, \code{hotspot_max_cn},
#'   \code{hotspot_fraction_amp}, \code{hotspot_cn_var}, \code{max_cn_global},
#'   \code{mean_breakpoints}, \code{n_alive}, \code{n_cells}, \code{failed}.
#'
#' @examples
#' \dontrun{
#' p <- clonal_params(positive_selection_rate = 0.25, bfb_prob = 0.95, lambda = 0.1)
#' df <- simulate_serial_passages(
#'   params = p, chromosome = "8", allele = "A", hotspot_pos = 139,
#'   N = 500, N_sub = 50, K = 10, passages_to_keep = c(1, 5, 10),
#'   n_replicates = 20
#' )
#' plot_serial_passages(df)
#' }
#'
#' @export
simulate_serial_passages <- function(
    params,
    chromosome,
    allele,
    hotspot_pos,
    N,
    N_sub,
    K,
    passages_to_keep = seq_len(K),
    n_replicates     = 1L,
    n_cores          = 1L,
    bin_length       = 1e6
) {
  stopifnot(N_sub <= N, K >= 1L)
  chr_allele <- paste0(chromosome, ":", allele)
  base_value <- ifelse(allele == "CN", 2, 1)
  na_stats   <- c(
    hotspot_mean_cn = NA_real_, hotspot_max_cn = NA_real_,
    hotspot_fraction_amp = NA_real_, hotspot_cn_var = NA_real_,
    max_cn_global = NA_real_, mean_breakpoints = NA_real_,
    n_alive = NA_real_, n_cells = 0L
  )

  # ── progress setup ──────────────────────────────────────────────────────────
  total_steps <- n_replicates * K
  pb <- cli::cli_progress_bar(
    name   = "Serial passages",
    total  = total_steps,
    format = paste0(
      "{cli::pb_spin} Rep {rep_counter}/{n_replicates} | ",
      "Passage {pass_counter}/{K} | ",
      "{cli::pb_bar} {cli::pb_percent} | ",
      "ETA {cli::pb_eta}"
    ),
    .envir = environment()
  )
  rep_counter  <- 0L
  pass_counter <- 0L
  pb_lock <- if (n_cores > 1L) parallel::makeCluster(0L) else NULL  # unused; note below
  # ────────────────────────────────────────────────────────────────────────────

  run_one_replicate <- function(rep_idx) {
    current_seqs <- NULL
    rows         <- list()

    for (p in seq_len(K)) {
      # update counters visible to the format string
      pass_counter <<- p
      rep_counter  <<- rep_idx

      need_cna <- p %in% passages_to_keep
      sim <- tryCatch(
        bridge_sim(
          initial_sequences       = current_seqs,
          chromosomes             = chromosome,
          bin_length              = bin_length,
          birth_rate              = params$birth_rate,
          death_rate              = params$death_rate,
          bfb_allele              = chr_allele,
          normal_dup_rate         = 0,
          bfb_prob                = params$bfb_prob,
          amp_rate                = params$amp_rate,
          del_rate                = params$del_rate,
          lambda                  = params$lambda,
          rate                    = params$rate,
          positive_selection_rate = params$positive_selection_rate,
          negative_selection_rate = params$negative_selection_rate,
          max_cells               = N,
          max_time                = 1e15,
          subsample               = N_sub,
          first_round_of_bfb      = is.null(current_seqs) && params$first_round_of_bfb,
          return_phylo            = FALSE,
          return_cna_data         = need_cna,
          hotspot                 = list(chr = chr_allele, pos = hotspot_pos),
          selection_type          = params$selection_type
        ),
        error = function(e) NULL
      )

      # tick after each passage attempt, whether it succeeded or failed
      cli::cli_progress_update(id = pb, .envir = parent.env(environment()))

      if (is.null(sim)) {
        if (p %in% passages_to_keep) {
          rows[[length(rows) + 1L]] <- c(passage = p, replicate = rep_idx,
                                         na_stats, failed = TRUE)
        }
        # fill remaining ticks for this replicate so the bar stays accurate
        remaining <- K - p
        if (remaining > 0L)
          cli::cli_progress_update(id = pb, inc = remaining,
                                   .envir = parent.env(environment()))
        break
      }

      # current_seqs <- compress_cell_sequences_to_cn(
      #   sim$cells,
      #   sim$input_parameters$chr_seq_lengths
      # )
      current_seqs = sim$cells

      if (need_cna) {
        cna_matrix <- tibble_to_matrix(sim$cna_data, value_column = allele)
        stats      <- summarise_clone(cna_matrix, hotspot_col = hotspot_pos,
                                      base_value = base_value, n_alive = sim$n_alive)
        rows[[length(rows) + 1L]] <- c(passage = p, replicate = rep_idx,
                                       stats, failed = FALSE)
      }
    }
    rows
  }

  all_rows <- do.call(c,
                      if (n_cores > 1L) {
                        parallel::mclapply(seq_len(n_replicates), run_one_replicate, mc.cores = n_cores)
                      } else {
                        lapply(seq_len(n_replicates), run_one_replicate)
                      }
  )

  cli::cli_progress_done(id = pb)

  result           <- do.call(rbind, lapply(all_rows, function(r) as.data.frame(t(r))))
  result           <- tibble::as_tibble(result)
  result$passage   <- as.integer(result$passage)
  result$replicate <- as.integer(result$replicate)
  result$failed    <- as.logical(result$failed)
  dplyr::select(result, passage, replicate, dplyr::everything())
}


#' Plot CNA complexity across serial passages
#'
#' Takes the wide tibble from \code{simulate_serial_passages()} and displays
#' how each summary statistic evolves with passage number (median ± IQR ribbon
#' across replicates, faceted by statistic).
#'
#' @param results       Output of \code{simulate_serial_passages()}.
#' @param stats_to_show Character vector of column names to plot.
#'
#' @return A \code{ggplot} object (faceted by statistic).
#'
#' @examples
#' \dontrun{
#' p <- clonal_params(positive_selection_rate = 0.25, bfb_prob = 0.95, lambda = 0.1)
#' df <- simulate_serial_passages(
#'   params = p, chromosome = "8", allele = "A", hotspot_pos = 139,
#'   N = 500, N_sub = 50, K = 8, n_replicates = 10
#' )
#' plot_serial_passages(df)
#' }
#'
#' @export
plot_serial_passages <- function(
  results,
  stats_to_show = c("mean_breakpoints", "hotspot_cn_var",
                    "max_cn_global", "hotspot_max_cn")
) {
  stats_to_show <- intersect(stats_to_show, names(results))
  if (length(stats_to_show) == 0)
    stop("None of the requested stats are present in results.")

  df <- results |>
    dplyr::filter(!failed) |>
    tidyr::pivot_longer(cols      = dplyr::all_of(stats_to_show),
                        names_to  = "stat",
                        values_to = "value") |>
    dplyr::group_by(passage, stat) |>
    dplyr::summarise(
      med = stats::median(value, na.rm = TRUE),
      lo  = stats::quantile(value, 0.25, na.rm = TRUE),
      hi  = stats::quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot2::ggplot(df, ggplot2::aes(x = passage, y = med, ymin = lo, ymax = hi)) +
    ggplot2::geom_ribbon(alpha = 0.2) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty()) +
    ggplot2::facet_wrap(~stat, scales = "free_y") +
    ggplot2::labs(
      x     = "Passage number",
      y     = "Statistic value  (median ± IQR)",
      title = "CNA complexity across serial passages"
    ) +
    ggplot2::theme_bw()
}


#' Plot CNA complexity across evolutionary depths
#'
#' Displays how each summary statistic evolves with simulation depth.  The
#' x-axis is either \code{"doublings"} (expected cell doublings, biologically
#' calibrated) or \code{"passage_cells"} (the target max_cells value).  For a
#' self-consistent simulator, all complexity metrics should increase
#' monotonically.
#'
#' @param summary_tibble Output of \code{summarise_passage_series()}.
#' @param x             X-axis: \code{"doublings"} (default) or
#'   \code{"passage_cells"}.
#' @param stats_to_show Character vector of stat names to show.  Defaults to
#'   the four main complexity metrics.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' series <- simulate_passage_series(
#'   c(early = 50, mid = 200, late = 600), chromosome = "7",
#'   allele = "A", hotspot_pos = 60, n_replicates = 5
#' )
#' stats <- summarise_passage_series(series)
#' plot_passage_series(stats)
#' plot_passage_series(stats, x = "passage_cells")  # raw cell count on x-axis
#' }
#'
#' @export
plot_passage_series <- function(
  summary_tibble,
  x             = c("doublings", "passage_cells"),
  stats_to_show = c("mean_breakpoints", "hotspot_cn_var",
                    "max_cn_global", "hotspot_max_cn")
) {
  x <- match.arg(x)

  df <- summary_tibble |>
    dplyr::filter(!failed, stat %in% stats_to_show) |>
    dplyr::group_by(passage_label, passage_cells, stat) |>
    dplyr::summarise(
      med_doublings = stats::median(doublings,     na.rm = TRUE),
      med           = stats::median(value,         na.rm = TRUE),
      lo            = stats::quantile(value, 0.25, na.rm = TRUE),
      hi            = stats::quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  x_col <- if (x == "doublings") "med_doublings" else "passage_cells"
  x_lab <- if (x == "doublings") "Expected cell doublings (median)" else "Cells at simulation stop (max_cells)"

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[[x_col]], y = med, ymin = lo, ymax = hi)
  ) +
    ggplot2::geom_ribbon(alpha = 0.2) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ stat, scales = "free_y") +
    ggplot2::labs(
      x     = x_lab,
      y     = "Statistic value  (median ± IQR)",
      title = "CNA complexity vs evolutionary depth"
    ) +
    ggplot2::theme_bw()
}


get_serial_passages_example_data <- function(
    params,
    chromosome,
    allele,
    hotspot_pos,
    N,
    N_sub,
    K,
    passages_to_keep = seq_len(K),
    bin_length        = 1e6,
    return_phylo      = FALSE
) {
  stopifnot(N_sub <= N, K >= 1L)
  chr_allele <- paste0(chromosome, ":", allele)

  current_seqs <- NULL
  sims <- list()

  pb <- cli::cli_progress_bar(
    name   = "Serial passages (example data)",
    total  = K,
    format = paste0(
      "{cli::pb_spin} Passage {pass_counter}/{K} | ",
      "{cli::pb_bar} {cli::pb_percent} | ETA {cli::pb_eta}"
    ),
    .envir = environment()
  )
  pass_counter <- 0L

  for (p in seq_len(K)) {
    pass_counter <- p
    need_cna <- p %in% passages_to_keep

    sim <- tryCatch(
      bridge_sim(
        initial_sequences       = current_seqs,
        chromosomes             = chromosome,
        bin_length              = bin_length,
        birth_rate              = params$birth_rate,
        death_rate              = params$death_rate,
        bfb_allele              = chr_allele,
        normal_dup_rate         = 0,
        bfb_prob                = params$bfb_prob,
        amp_rate                = params$amp_rate,
        del_rate                = params$del_rate,
        lambda                  = params$lambda,
        rate                    = params$rate,
        positive_selection_rate = params$positive_selection_rate,
        negative_selection_rate = params$negative_selection_rate,
        max_cells               = N,
        max_time                = 1e15,
        subsample               = N_sub,
        first_round_of_bfb      = is.null(current_seqs) && params$first_round_of_bfb,
        return_phylo            = return_phylo,
        return_cna_data         = need_cna,
        hotspot                 = list(chr = chr_allele, pos = hotspot_pos),
        selection_type          = params$selection_type,
        breakpoint_support = "beta",
        alpha = 50,
        beta = 50
      ),
      error = function(e) NULL
    )

    cli::cli_progress_update(id = pb, .envir = environment())

    if (is.null(sim)) {
      warning(sprintf("Simulation failed at passage %d; stopping early.", p))
      break
    }

    current_seqs <- sim$cells

    if (need_cna) {
      sims[[paste0("passage_", p)]] <- sim
    }
  }

  cli::cli_progress_done(id = pb)

  sims
}
