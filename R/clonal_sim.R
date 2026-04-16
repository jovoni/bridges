
# ── Clonal simulation study ───────────────────────────────────────────────────
#
# Tools for studying how BFB and selection shape the copy number distribution
# of an evolving clone. The core workflow is:
#
#   1.  Define one or more parameter sets (BFB rate, selection coefficients, …).
#   2.  Call run_clonal_replicates() to grow N independent clones, stopping at
#       max_cells or max_time (whichever fires first), then sample sample_cells
#       cells using subsample_sim().
#   3.  Call compare_clonal_params() to sweep over multiple parameter sets and get
#       a tidy data frame ready for plotting.
#   4.  Call plot_clonal_comparison() to visualise how summary statistics differ
#       across conditions.
#
# hotspot_pos is the 1-based bin index of the genomic position of interest
# (e.g. the bin most likely to be amplified by BFB).  It corresponds directly
# to the `pos` field inside bridge_sim's `hotspot` argument.

# ── Default parameter constructor ─────────────────────────────────────────────

#' Build a parameter list for clonal simulations
#'
#' Returns a complete parameter list with sensible defaults. Any argument can
#' be overridden.  Pass the result to \code{run_clonal_replicates()} or
#' \code{compare_clonal_params()}.
#'
#' @param bfb_prob            Relative probability of a BFB event per division.
#' @param amp_rate            Relative probability of a focal amplification.
#' @param del_rate            Relative probability of a focal deletion.
#' @param positive_selection_rate  Multiplicative birth-rate boost for cells that
#'   have gained a copy at \code{hotspot_pos}.  0 = neutral.
#' @param negative_selection_rate  Multiplicative death-rate boost for cells that
#'   have gained a copy at \code{hotspot_pos}.  0 = neutral.
#' @param birth_rate          Base cell birth rate.
#' @param death_rate          Base cell death rate.
#' @param lambda              Mean number of genomic events per daughter cell per
#'   division (Poisson).
#' @param rate                Scale of focal amp/del events: segment length is
#'   drawn from Exp(1/rate).
#' @param first_round_of_bfb  If TRUE the founding cell already carries one BFB
#'   event, mimicking a clone that initiated from a single BFB.
#'
#' @return Named list of simulation parameters.
#' @export
clonal_params = function(
  bfb_prob                 = 0.5,
  amp_rate                 = 0.3,
  del_rate                 = 0.2,
  positive_selection_rate  = 0,
  negative_selection_rate  = 0,
  birth_rate               = 1.0,
  death_rate               = 0.1,
  lambda                   = 2,
  rate                     = 20,
  first_round_of_bfb       = TRUE
) {
  list(
    bfb_prob                = bfb_prob,
    amp_rate                = amp_rate,
    del_rate                = del_rate,
    positive_selection_rate = positive_selection_rate,
    negative_selection_rate = negative_selection_rate,
    birth_rate              = birth_rate,
    death_rate              = death_rate,
    lambda                  = lambda,
    rate                    = rate,
    first_round_of_bfb      = first_round_of_bfb
  )
}


# ── Single replicate ──────────────────────────────────────────────────────────

#' Simulate one clonal evolution and return the CN matrix
#'
#' Runs \code{bridge_sim()} until \code{max_cells} are alive OR \code{max_time}
#' is reached (whichever fires first).  If \code{sample_cells} is set, a random
#' subsample of that many cells is drawn from the alive population using
#' \code{subsample_sim()}.
#'
#' @param chromosome   Chromosome to simulate (character, e.g. \code{"7"}).
#' @param allele       Allele to track (\code{"A"}, \code{"B"}, or \code{"CN"}).
#' @param hotspot_pos  1-based bin index of the position of interest.
#' @param params       Parameter list from \code{clonal_params()}.
#' @param max_cells    Stop when this many cells are alive.  Default 1000 (effectively
#'   no limit when the stopping criterion is time-based).
#' @param max_time     Stop at this simulation time.  Default 300.
#' @param sample_cells Number of cells to randomly sample from the alive population
#'   at the end of the simulation.  NULL keeps all alive cells.
#' @param bin_length   Bin size in bp.  Default 1 Mb.
#'
#' @return A list with:
#' \describe{
#'   \item{cna_matrix}{Integer matrix (cells x bins) for the target allele,
#'     after subsampling.}
#'   \item{n_alive}{Number of cells alive at end of simulation, before sampling.}
#'   \item{n_cells}{Number of cells in \code{cna_matrix} (after sampling).}
#'   \item{sim}{Full \code{bridge_sim()} output (for downstream use).}
#' }
simulate_clone = function(
  chromosome,
  allele,
  hotspot_pos,
  params       = clonal_params(),
  max_cells    = 1000,
  max_time     = 300,
  sample_cells = NULL,
  bin_length   = 1e6
) {
  chr_allele = paste0(chromosome, ":", allele)

  sim = bridge_sim(
    initial_cells           = 1,
    chromosomes             = c(chromosome),
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
    max_cells               = max_cells,
    max_time                = max_time,
    first_round_of_bfb      = params$first_round_of_bfb,
    return_phylo            = FALSE,
    hotspot                 = list(chr = chr_allele, pos = hotspot_pos),
    breakpoint_support      = "uniform"
  )

  n_alive = length(sim$cells)

  # Subsample using the existing subsample_sim() if requested
  if (!is.null(sample_cells) && n_alive > sample_cells) {
    f = sample_cells / n_alive
    sim = subsample_sim(sim, f_subsample = f)
  } else if (!is.null(sample_cells) && n_alive < sample_cells) {
    warning(sprintf(
      "Requested sample_cells=%d but only %d cells are alive. Keeping all.",
      sample_cells, n_alive
    ))
  }

  cna_matrix = tibble_to_matrix(sim$cna_data, value_column = allele)
  list(
    cna_matrix = cna_matrix,
    n_alive    = n_alive,
    n_cells    = nrow(cna_matrix),
    sim        = sim
  )
}


# ── Summary statistics for one replicate ──────────────────────────────────────

#' Compute summary statistics for a single clonal CN matrix
#'
#' @param cna_matrix  Integer matrix (cells x bins) from \code{simulate_clone()}.
#' @param hotspot_col Column index of the hotspot bin in \code{cna_matrix}.
#' @param base_value  Baseline (diploid) copy number: 1 for allele-specific, 2 for CN.
#' @param n_alive     Number of alive cells before subsampling (from
#'   \code{simulate_clone()$n_alive}).  Stored as-is for downstream analysis.
#'
#' @return Named numeric vector of summary statistics.
summarise_clone = function(cna_matrix, hotspot_col, base_value = 1, n_alive = NA) {
  if (nrow(cna_matrix) == 0 || ncol(cna_matrix) == 0) {
    return(c(
      hotspot_mean_cn      = NA_real_,
      hotspot_max_cn       = NA_real_,
      hotspot_fraction_amp = NA_real_,
      hotspot_cn_var       = NA_real_,
      max_cn_global        = NA_real_,
      mean_breakpoints     = NA_real_,
      n_alive              = n_alive,
      n_cells              = 0L
    ))
  }

  hotspot_col = min(hotspot_col, ncol(cna_matrix))  # guard against out-of-range
  hotspot_cn  = cna_matrix[, hotspot_col]

  bp_per_cell = apply(cna_matrix, 1, function(r) sum(diff(r) != 0))

  c(
    hotspot_mean_cn      = mean(hotspot_cn),
    hotspot_max_cn       = max(hotspot_cn),
    hotspot_fraction_amp = mean(hotspot_cn > base_value),
    hotspot_cn_var       = stats::var(hotspot_cn),
    max_cn_global        = max(cna_matrix),
    mean_breakpoints     = mean(bp_per_cell),
    n_alive              = n_alive,
    n_cells              = nrow(cna_matrix)
  )
}


# ── Replicated simulation ─────────────────────────────────────────────────────

#' Run N independent clonal simulations and summarise each
#'
#' Each replicate is an independent run of \code{simulate_clone()}.  The
#' simulation stops when \code{max_cells} are alive or \code{max_time} is
#' reached, then \code{sample_cells} cells are drawn using
#' \code{subsample_sim()}.  Failed simulations (e.g. all cells die) are
#' recorded with NA summary statistics and flagged with \code{failed = TRUE}.
#'
#' @param chromosome   Chromosome to simulate (character, e.g. \code{"7"}).
#' @param allele       Allele to track (\code{"A"}, \code{"B"}, or \code{"CN"}).
#' @param hotspot_pos  1-based bin index of the position of interest.
#' @param params       Parameter list from \code{clonal_params()}.
#' @param N_replicates Number of independent clonal evolutions to simulate.
#' @param max_cells    Stopping criterion: maximum alive cells.
#' @param max_time     Stopping criterion: simulation time limit.
#' @param sample_cells Number of cells to sample after simulation.  NULL keeps all.
#' @param bin_length   Bin size in bp.  Default 1 Mb.
#' @param n_cores      Number of parallel cores (uses \code{parallel::mclapply}
#'   when > 1).
#'
#' @return A \code{tibble} with one row per replicate and columns for each
#'   summary statistic plus \code{replicate}, \code{n_alive}, and \code{failed}.
#' @export
run_clonal_replicates = function(
  chromosome,
  allele,
  hotspot_pos,
  params        = clonal_params(),
  N_replicates  = 100,
  max_cells     = 1000,
  max_time      = 300,
  sample_cells  = NULL,
  bin_length    = 1e6,
  n_cores       = 1
) {
  base_value = ifelse(allele == "CN", 2, 1)

  run_one = function(i) {
    tryCatch({
      res   = simulate_clone(chromosome, allele, hotspot_pos, params,
                             max_cells, max_time, sample_cells, bin_length)
      stats = summarise_clone(res$cna_matrix,
                              hotspot_col = hotspot_pos,
                              base_value  = base_value,
                              n_alive     = res$n_alive)
      c(stats, failed = FALSE)
    }, error = function(e) {
      c(hotspot_mean_cn = NA_real_, hotspot_max_cn = NA_real_,
        hotspot_fraction_amp = NA_real_, hotspot_cn_var = NA_real_,
        max_cn_global = NA_real_, mean_breakpoints = NA_real_,
        n_alive = NA_real_, n_cells = 0L, failed = TRUE)
    })
  }

  if (n_cores > 1) {
    rows = parallel::mclapply(seq_len(N_replicates), run_one, mc.cores = n_cores)
  } else {
    rows = lapply(seq_len(N_replicates), run_one)
  }

  result           = do.call(rbind, lapply(rows, function(r) as.data.frame(t(r))))
  result$replicate = seq_len(N_replicates)
  result$failed    = as.logical(result$failed)
  tibble::as_tibble(result)
}


# ── Parameter sweep ───────────────────────────────────────────────────────────

#' Compare summary statistics across multiple parameter sets
#'
#' Runs \code{run_clonal_replicates()} for each entry in \code{param_list} and
#' combines the results into a single tidy data frame.  The name of each list
#' entry becomes the \code{condition} column, making it easy to facet or
#' colour by condition in downstream plots.
#'
#' @param param_list   Named list of parameter lists, each from
#'   \code{clonal_params()}.  Names become the \code{condition} column.
#' @param chromosome   Chromosome to simulate.
#' @param allele       Allele to track.
#' @param hotspot_pos  1-based bin index of the position of interest.
#' @param N_replicates Number of independent replicates per condition.
#' @param max_cells    Stopping criterion: maximum alive cells.
#' @param max_time     Stopping criterion: simulation time limit.
#' @param sample_cells Number of cells to sample after simulation.  NULL keeps all.
#' @param bin_length   Bin size in bp.
#' @param n_cores      Cores for parallelism (applies within each condition).
#'
#' @return A \code{tibble} with all replicates across all conditions.  Columns:
#'   \code{condition}, \code{replicate}, summary statistics, \code{n_alive},
#'   \code{failed}.
#'
#' @examples
#' \dontrun{
#' results = compare_clonal_params(
#'   param_list = list(
#'     neutral    = clonal_params(positive_selection_rate = 0),
#'     selection  = clonal_params(positive_selection_rate = 2),
#'     bfb_only   = clonal_params(bfb_prob = 0.9, positive_selection_rate = 0),
#'     bfb_select = clonal_params(bfb_prob = 0.9, positive_selection_rate = 2)
#'   ),
#'   chromosome   = "7",
#'   allele       = "A",
#'   hotspot_pos  = 60,
#'   N_replicates = 100,
#'   max_time     = 20,
#'   sample_cells = 50
#' )
#' plot_clonal_comparison(results)
#' }
#' @export
compare_clonal_params = function(
  param_list,
  chromosome,
  allele,
  hotspot_pos,
  N_replicates = 100,
  max_cells    = 1000,
  max_time     = 300,
  sample_cells = NULL,
  bin_length   = 1e6,
  n_cores      = 1
) {
  if (is.null(names(param_list)) || any(names(param_list) == "")) {
    stop("All entries in param_list must be named.")
  }

  all_results = lapply(names(param_list), function(cond) {
    message("Running condition: ", cond)
    df = run_clonal_replicates(
      chromosome   = chromosome,
      allele       = allele,
      hotspot_pos  = hotspot_pos,
      params       = param_list[[cond]],
      N_replicates = N_replicates,
      max_cells    = max_cells,
      max_time     = max_time,
      sample_cells = sample_cells,
      bin_length   = bin_length,
      n_cores      = n_cores
    )
    df$condition = cond
    df
  })

  dplyr::bind_rows(all_results) %>%
    dplyr::select(.data$condition, .data$replicate, dplyr::everything())
}


# ── Observed data summary ─────────────────────────────────────────────────────

#' Compute summary statistics from a real CNA sample
#'
#' Applies the same statistics used in \code{summarise_clone()} to an observed
#' data frame so the result can be overlaid on simulation comparison plots.
#'
#' Accepts both \code{start}/\code{end} and \code{from}/\code{to} column names.
#'
#' @param data        Data frame with columns \code{cell_id}, \code{chr},
#'   \code{start}/\code{from}, \code{end}/\code{to}, and the allele column.
#' @param chromosome  Chromosome to analyse (character, e.g. \code{"7"}).
#' @param allele      Allele column to use (\code{"A"}, \code{"B"}, or \code{"CN"}).
#' @param hotspot_pos 1-based bin index of the position of interest.
#'
#' @return Named numeric vector of summary statistics, compatible with
#'   \code{plot_clonal_comparison(observed_stats = ...)}.
#' @export
compute_observed_stats = function(data, chromosome, allele, hotspot_pos) {
  # Normalise column names: accept from/to as aliases for start/end
  if (!("start" %in% names(data)) && "from" %in% names(data)) {
    data = dplyr::rename(data, start = "from", end = "to")
  }
  if (!("chr" %in% names(data)) && "chromosome" %in% names(data)) {
    data = dplyr::rename(data, chr = "chromosome")
  }

  missing = setdiff(c("cell_id", "chr", "start", "end", allele), names(data))
  if (length(missing) > 0) {
    stop("Missing columns in data: ", paste(missing, collapse = ", "))
  }

  base_value = ifelse(allele == "CN", 2, 1)
  cna_matrix = tibble_to_matrix(data, value_column = allele)

  summarise_clone(cna_matrix, hotspot_col = hotspot_pos,
                  base_value = base_value, n_alive = nrow(cna_matrix))
}


# ── Visualisation ─────────────────────────────────────────────────────────────

#' Plot summary statistics across clonal simulation conditions
#'
#' Produces a faceted box-plot grid: one panel per summary statistic, one box
#' per condition.  If \code{observed_stats} is supplied (output of
#' \code{compute_observed_stats()}), a horizontal dashed line is added to each
#' panel showing the observed value.
#'
#' @param results        Output of \code{compare_clonal_params()}.
#' @param observed_stats Optional named numeric vector from
#'   \code{compute_observed_stats()}.  A horizontal reference line is drawn in
#'   each panel for statistics present in both \code{results} and
#'   \code{observed_stats}.
#' @param stats          Character vector of statistic names to plot.
#' @param fill_var       Column name to use for box fill colour.  Default
#'   \code{"condition"}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_clonal_comparison = function(
  results,
  observed_stats = NULL,
  stats    = c("hotspot_mean_cn", "hotspot_max_cn", "hotspot_fraction_amp",
               "hotspot_cn_var", "max_cn_global", "mean_breakpoints", "n_alive"),
  fill_var = "condition"
) {
  stats = intersect(stats, names(results))
  if (length(stats) == 0) stop("None of the requested stats are present in results.")

  plot_data = results[!results$failed, ]

  long = tidyr::pivot_longer(
    plot_data,
    cols      = dplyr::all_of(stats),
    names_to  = "statistic",
    values_to = "value"
  )
  long$condition = factor(long$condition, levels = unique(results$condition))

  stat_labels = c(
    hotspot_mean_cn      = "Hotspot mean CN",
    hotspot_max_cn       = "Hotspot max CN",
    hotspot_fraction_amp = "Fraction cells amplified at hotspot",
    hotspot_cn_var       = "Hotspot CN variance",
    max_cn_global        = "Max CN (global)",
    mean_breakpoints     = "Mean breakpoints per cell",
    n_alive              = "Alive cells at end of simulation"
  )
  long$statistic = factor(long$statistic, levels = stats,
                          labels = stat_labels[stats])

  p = ggplot2::ggplot(long, ggplot2::aes(x = .data$condition,
                                          y = .data$value,
                                          fill = .data[[fill_var]])) +
    ggplot2::geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
    ggplot2::facet_wrap(~ statistic, scales = "free_y") +
    ggplot2::labs(x = "Condition", y = NULL, fill = fill_var) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position = "none"
    )

  if (!is.null(observed_stats)) {
    # Keep only stats present in the plot and with non-NA observed values
    obs_stats_in_plot = intersect(stats, names(observed_stats))
    obs_df = data.frame(
      statistic = factor(obs_stats_in_plot, levels = stats,
                         labels = stat_labels[obs_stats_in_plot]),
      value     = as.numeric(observed_stats[obs_stats_in_plot])
    )
    obs_df = obs_df[!is.na(obs_df$value), ]

    if (nrow(obs_df) > 0) {
      p = p + ggplot2::geom_hline(
        data        = obs_df,
        ggplot2::aes(yintercept = .data$value),
        colour      = "firebrick",
        linetype    = "dashed",
        linewidth   = 0.7,
        inherit.aes = FALSE
      )
    }
  }

  p
}
