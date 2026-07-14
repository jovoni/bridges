# Models biological serial passaging: grow to N cells, subsample N_sub, repeat.
# Each bottleneck introduces drift and allows clonal sweeps to accumulate
# across passages. It's a richer model than simply running a longer simulation.
#
# Workflow:
#   df <- simulate_serial_passages(params, chromosome = "8", allele = "A",
#                                  hotspot_pos = 139, N = 500, N_sub = 50, K = 10)
#   plot_serial_passages(df)   # wide tibble


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

  # progress setup
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
  #

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
#' how each summary statistic evolves with passage number (median and IQR ribbon
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
      y     = "Statistic value  (median +- IQR)",
      title = "CNA complexity across serial passages"
    ) +
    ggplot2::theme_bw()
}
