
# ── Summary statistics ────────────────────────────────────────────────────────

# Per-bin gain rate (fraction of cells above baseline)
get_gain_profile = function(cna_matrix, base_value) {
  colMeans(cna_matrix > base_value)
}

# Per-bin loss rate (fraction of cells below baseline)
get_loss_profile = function(cna_matrix, base_value) {
  colMeans(cna_matrix < base_value)
}

# Per-bin mean CN
get_avg_cn_profile = function(cna_matrix) {
  colMeans(cna_matrix)
}

# Per-bin CN variance across cells — captures clonal heterogeneity
get_cn_variance_profile = function(cna_matrix) {
  apply(cna_matrix, 2, stats::var)
}

# Mean number of CN breakpoints per cell — BFB creates many segment boundaries
get_mean_breakpoints_per_cell = function(cna_matrix) {
  bp_counts = apply(cna_matrix, 1, function(row) sum(diff(row) != 0))
  mean(bp_counts)
}

# Maximum CN value across all cells — BFB drives extreme focal amplification
get_max_cn = function(cna_matrix) {
  max(cna_matrix)
}

# Mean CN-segment length (used to calibrate the `rate` param of bridge_sim)
get_cn_length_rate = function(cna_matrix) {
  cna_lengths = lapply(1:nrow(cna_matrix), function(j) {
    idxs = which(diff(cna_matrix[j, ]) != 0)
    if (length(idxs) > 1) diff(idxs) else NULL
  })
  vals = unlist(cna_lengths)
  if (length(vals) == 0) return(ncol(cna_matrix))  # fallback: whole chromosome
  mean(vals)
}

# Breakpoint position distribution (used for custom_breakpoints in bridge_sim)
get_breakpoints_dist = function(cna_matrix) {
  lapply(1:nrow(cna_matrix), function(j) {
    which(diff(cna_matrix[j, ]) != 0)
  }) %>% unlist()
}


# ── Data-driven calibration helpers ──────────────────────────────────────────

#' Infer bin length from a CNA data frame
#'
#' Computes the median (end - start) for the given chromosome so the simulator
#' uses the same resolution as the observed data.
#'
#' @param data CNA data frame with columns cell_id, chr, start, end.
#' @param chromosome Chromosome to use (character, e.g. "7").
#' @return Numeric bin length in base pairs.
infer_bin_length = function(data, chromosome) {
  chr_data = data[data$chr == chromosome, ]
  as.integer(stats::median(chr_data$end - chr_data$start))
}

#' Derive hotspot bin position from a CN matrix
#'
#' Returns the 1-based bin index of the genomic region most likely to be the
#' BFB hotspot, defined as the bin with the highest CN variance across cells.
#' This is data-driven and removes the need for the caller to supply `pos`.
#'
#' @param cna_matrix Numeric matrix (cells x bins).
#' @return Integer bin index.
derive_hotspot_pos = function(cna_matrix) {
  which.max(apply(cna_matrix, 2, stats::var))
}


# ── Core ABC simulation ───────────────────────────────────────────────────────

#' Run a single ABC simulation and compute distance to observed data
#'
#' @param cna_data Observed CNA data frame (cell_id, chr, start, end, CN, A, B).
#' @param allele Allele column to use ("A", "B", or "CN").
#' @param chromosome Chromosome name as it appears in cna_data$chr (e.g. "7").
#' @param params Named list of simulation parameters from \code{sample_priors()}.
#' @param pos Optional integer bin index for the BFB hotspot. Derived automatically
#'   from the observed data if NULL (default).
#' @param bin_length Optional bin size in bp. Inferred from observed data if NULL.
#'
#' @return Named list with distance components and the parameter set used.
run_ABC = function(cna_data, allele, chromosome, params, pos = NULL, bin_length = NULL) {
  base_value = ifelse(allele == "CN", 2, 1)

  target_cna_matrix = tibble_to_matrix(cna_data, chromosome = chromosome, value_column = allele)
  n_target_cells    = nrow(target_cna_matrix)
  n_target_bins     = ncol(target_cna_matrix)

  # ── Infer bin_length from observed data if not supplied ──────────────────
  if (is.null(bin_length)) {
    bin_length = infer_bin_length(cna_data, chromosome)
  }

  # ── Auto-derive hotspot position if not supplied ─────────────────────────
  if (is.null(pos)) {
    pos = derive_hotspot_pos(target_cna_matrix)
  }

  # ── Compute observed summary statistics ──────────────────────────────────
  target_gain_profile     = get_gain_profile(target_cna_matrix, base_value)
  target_loss_profile     = get_loss_profile(target_cna_matrix, base_value)
  target_avg_cn_profile   = get_avg_cn_profile(target_cna_matrix)
  target_var_profile      = get_cn_variance_profile(target_cna_matrix)
  target_max_cn           = get_max_cn(target_cna_matrix)
  target_mean_bp          = get_mean_breakpoints_per_cell(target_cna_matrix)
  target_cna_rate         = get_cn_length_rate(target_cna_matrix)
  custom_breakpoint_dist  = get_breakpoints_dist(target_cna_matrix)

  # ── Run simulation ────────────────────────────────────────────────────────
  chr_allele = paste0(chromosome, ":", allele)
  sim = bridge_sim(
    initial_cells           = 1,
    bin_length              = bin_length,
    max_time                = 300,
    max_cells               = n_target_cells,
    normal_dup_rate         = 0,
    amp_rate                = params$amp_rate,
    del_rate                = params$del_rate,
    bfb_prob                = params$bfb_prob,
    birth_rate              = params$birth_rate,
    death_rate              = params$death_rate,
    rate                    = target_cna_rate,
    lambda                  = params$lambda,
    bfb_allele              = chr_allele,
    chromosomes             = c(chromosome),
    first_round_of_bfb      = TRUE,
    positive_selection_rate = 0,
    negative_selection_rate = 0,
    hotspot                 = list(chr = chr_allele, pos = pos),
    breakpoint_support      = "custom",
    custom_breakpoints      = custom_breakpoint_dist
  )

  sim_cna_matrix = tibble_to_matrix(sim$cna_data, chromosome = chromosome, value_column = allele)
  n_sim_bins     = ncol(sim_cna_matrix)

  # Guard: bin counts must match for a meaningful profile comparison
  if (n_sim_bins != n_target_bins) {
    stop(sprintf(
      "Simulated data has %d bins but observed data has %d bins (chromosome %s). ",
      n_sim_bins, n_target_bins, chromosome,
      "Check that bin_length matches the resolution of your input data."
    ))
  }

  # ── Compute predicted summary statistics ─────────────────────────────────
  pred_gain_profile   = get_gain_profile(sim_cna_matrix, base_value)
  pred_loss_profile   = get_loss_profile(sim_cna_matrix, base_value)
  pred_avg_cn_profile = get_avg_cn_profile(sim_cna_matrix)
  pred_var_profile    = get_cn_variance_profile(sim_cna_matrix)
  pred_max_cn         = get_max_cn(sim_cna_matrix)
  pred_mean_bp        = get_mean_breakpoints_per_cell(sim_cna_matrix)

  # ── Distance components ───────────────────────────────────────────────────
  rmse = function(a, b) sqrt(mean((a - b)^2))

  d_gain    = rmse(pred_gain_profile,   target_gain_profile)
  d_loss    = rmse(pred_loss_profile,   target_loss_profile)
  d_avg_cn  = rmse(pred_avg_cn_profile, target_avg_cn_profile)
  d_var     = rmse(pred_var_profile,    target_var_profile)
  # Scalars: normalise by observed value so they are scale-invariant
  d_max_cn  = abs(pred_max_cn  - target_max_cn)  / max(target_max_cn,  1)
  d_mean_bp = abs(pred_mean_bp - target_mean_bp) / max(target_mean_bp, 1)

  total_distance = d_gain + d_loss + d_avg_cn + d_var + d_max_cn + d_mean_bp

  list(
    distance  = total_distance,
    d_gain    = d_gain,
    d_loss    = d_loss,
    d_avg_cn  = d_avg_cn,
    d_var     = d_var,
    d_max_cn  = d_max_cn,
    d_mean_bp = d_mean_bp,
    params    = params
  )
}


# ── Prior sampling ────────────────────────────────────────────────────────────

#' Sample parameters from prior distributions
#'
#' Draws one candidate parameter set:
#' - amp_rate / del_rate / bfb_prob as relative proportions from a symmetric
#'   Dirichlet(1,1,1) — uniform on the 2-simplex.
#' - death_rate as a fraction of birth_rate, uniform on (0, 0.8).
#' - lambda (mean genomic events per daughter cell) uniform on (0.5, 5).
#'
#' @return Named list of parameter values.
sample_priors = function() {
  relative_rates = MCMCpack::rdirichlet(1, c(1, 1, 1))[1, ]

  birth_rate     = 1.0
  death_fraction = stats::runif(1, 0, 0.8)
  death_rate     = birth_rate * death_fraction

  list(
    amp_rate       = relative_rates[1],
    del_rate       = relative_rates[2],
    bfb_prob       = relative_rates[3],
    birth_rate     = birth_rate,
    death_rate     = death_rate,
    lambda         = stats::runif(1, 0.5, 5),
    death_fraction = death_fraction
  )
}


# ── Main ABC rejection algorithm ──────────────────────────────────────────────

#' ABC inference of BFB simulation parameters
#'
#' Uses rejection ABC to infer simulation parameters that reproduce the
#' copy number distribution of an observed single-cell dataset. The simulator
#' (\code{bridge_sim}) is run \code{n_simulations} times with parameters drawn
#' from \code{sample_priors()}. The \code{tolerance_quantile} fraction with
#' the smallest distance to the observed summary statistics is retained as the
#' approximate posterior.
#'
#' @param cna_data Data frame with columns cell_id, chr, start, end, CN, A, B.
#' @param allele Allele to match ("A", "B", or "CN").
#' @param chromosome Chromosome to focus on (character, e.g. "7").
#' @param n_simulations Total number of simulations to run. Default: 10000.
#' @param tolerance_quantile Fraction of simulations to accept. Default: 0.01.
#' @param n_cores Number of parallel cores (uses \code{parallel::mclapply} when > 1).
#'   Default: 1.
#' @param pos Optional integer bin index for the BFB hotspot in the simulation.
#'   Automatically derived from the observed data as the bin with the highest CN
#'   variance if NULL (default).
#' @param bin_length Optional bin size in bp. Inferred from observed data if NULL.
#'
#' @return A list containing:
#' \describe{
#'   \item{accepted_params}{Data frame of accepted parameter draws.}
#'   \item{param_summary}{Per-parameter mean / median / SD / 95\% CI.}
#'   \item{n_simulations}{Total simulations attempted.}
#'   \item{n_accepted}{Number of accepted simulations.}
#'   \item{acceptance_rate}{Fraction accepted out of valid (non-error) simulations.}
#'   \item{tolerance_threshold}{Distance cutoff used for acceptance.}
#' }
#'
#' @export
abc_inference = function(cna_data,
                         allele,
                         chromosome,
                         n_simulations      = 10000,
                         tolerance_quantile = 0.01,
                         n_cores            = 1,
                         pos                = NULL,
                         bin_length         = NULL) {

  # Pre-compute bin_length and pos once so every simulation reuses them
  if (is.null(bin_length)) {
    bin_length = infer_bin_length(cna_data, chromosome)
    message("Inferred bin_length from data: ", bin_length, " bp")
  }

  if (is.null(pos)) {
    obs_matrix = tibble_to_matrix(cna_data, chromosome = chromosome, value_column = allele)
    pos = derive_hotspot_pos(obs_matrix)
    message("Derived hotspot pos from data: bin ", pos)
  }

  run_one = function(i) {
    params = sample_priors()
    tryCatch(
      run_ABC(cna_data, allele, chromosome, params, pos = pos, bin_length = bin_length),
      error = function(e) list(distance = Inf, params = params, error = e$message)
    )
  }

  if (n_cores > 1) {
    results = parallel::mclapply(1:n_simulations, run_one, mc.cores = n_cores)
  } else {
    pb = utils::txtProgressBar(min = 0, max = n_simulations, style = 3)
    results = vector("list", n_simulations)
    for (i in 1:n_simulations) {
      results[[i]] = run_one(i)
      utils::setTxtProgressBar(pb, i)
    }
    close(pb)
  }

  valid_results = results[sapply(results, function(x) is.finite(x$distance))]

  if (length(valid_results) == 0) {
    stop("No valid simulations completed. Check your model parameters and data.")
  }

  distances      = sapply(valid_results, function(x) x$distance)
  n_accept       = max(1L, floor(length(valid_results) * tolerance_quantile))
  accepted       = valid_results[order(distances)[seq_len(n_accept)]]

  accepted_params = data.frame(
    amp_rate       = sapply(accepted, function(x) x$params$amp_rate),
    del_rate       = sapply(accepted, function(x) x$params$del_rate),
    bfb_prob       = sapply(accepted, function(x) x$params$bfb_prob),
    birth_rate     = sapply(accepted, function(x) x$params$birth_rate),
    death_rate     = sapply(accepted, function(x) x$params$death_rate),
    lambda         = sapply(accepted, function(x) x$params$lambda),
    death_fraction = sapply(accepted, function(x) x$params$death_fraction),
    distance       = sapply(accepted, function(x) x$distance)
  )

  param_cols   = setdiff(names(accepted_params), "distance")
  param_summary = data.frame(
    parameter = param_cols,
    mean      = sapply(accepted_params[param_cols], mean),
    median    = sapply(accepted_params[param_cols], stats::median),
    sd        = sapply(accepted_params[param_cols], stats::sd),
    q025      = sapply(accepted_params[param_cols], stats::quantile, 0.025),
    q975      = sapply(accepted_params[param_cols], stats::quantile, 0.975),
    row.names = NULL
  )

  list(
    accepted_params    = accepted_params,
    param_summary      = param_summary,
    n_simulations      = n_simulations,
    n_accepted         = n_accept,
    acceptance_rate    = n_accept / length(valid_results),
    tolerance_threshold = max(accepted_params$distance)
  )
}


# ── Diagnostics ───────────────────────────────────────────────────────────────

#' Plot ABC posterior distributions
#'
#' @param abc_results Output of \code{abc_inference()}.
#' @export
plot_abc_results = function(abc_results) {
  params = abc_results$accepted_params

  make_hist = function(df, x, fill, title) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]])) +
      ggplot2::geom_histogram(bins = 30, alpha = 0.7, fill = fill) +
      ggplot2::geom_vline(xintercept = stats::median(df[[x]]),
                          color = "red", linetype = "dashed") +
      ggplot2::labs(title = title, x = x, y = "Count") +
      ggplot2::theme_minimal()
  }

  plots = list(
    make_hist(params, "amp_rate",       "steelblue", "Relative: amp_rate"),
    make_hist(params, "del_rate",       "steelblue", "Relative: del_rate"),
    make_hist(params, "bfb_prob",       "steelblue", "Relative: bfb_prob"),
    make_hist(params, "death_fraction", "darkgreen", "Posterior: death_fraction"),
    make_hist(params, "lambda",         "darkgreen", "Posterior: lambda"),
    make_hist(params, "distance",       "grey40",    "Accepted distances")
  )

  # Amp vs del coloured by BFB proportion
  plots[[7]] = ggplot2::ggplot(params, ggplot2::aes(x = .data$amp_rate,
                                                     y = .data$del_rate,
                                                     color = .data$bfb_prob)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::scale_color_gradient(low = "blue", high = "red", name = "BFB prob") +
    ggplot2::labs(title = "Amp vs Del (coloured by BFB)", x = "amp_rate", y = "del_rate") +
    ggplot2::theme_minimal()

  gridExtra::grid.arrange(grobs = plots, ncol = 3)
}
