# Parameter constructor and summary statistics used by
# simulate_serial_passages() (see passage_sim.R).
#
# hotspot_pos is the 1-based bin index of the genomic position of interest
# (e.g. the bin most likely to be amplified by BFB). It corresponds directly
# to the `pos` field inside bridge_sim's `hotspot` argument.

#' Build a parameter list for clonal simulations
#'
#' Returns a complete parameter list with sensible defaults. Any argument can
#' be overridden. Pass the result to \code{\link{simulate_serial_passages}}.
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
#' @param rate                Mean length, in bins, of focal amplifications and
#'   deletions (exponentially distributed).
#' @param first_round_of_bfb  If TRUE the founding cell already carries one BFB
#'   event, mimicking a clone that initiated from a single BFB.
#' @param selection_type How hotspot copy number translates into a
#'   selection advantage: \code{"constant"} (any gain gives the same
#'   fixed advantage), \code{"linear"}, or \code{"saturation"} -- passed
#'   straight through to \code{bridge_sim()}. Default:
#'   \code{"constant"}.
#' @param saturation_K Saturation constant used when
#'   \code{selection_type = "saturation"}. Default: 10.
#'
#' @return Named list of simulation parameters.
#'
#' @examples
#' # Default parameters (BFB-driven, neutral selection)
#' p <- clonal_params()
#'
#' # High BFB probability with positive selection
#' p_sel <- clonal_params(bfb_prob = 0.8, positive_selection_rate = 2)
#'
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
  first_round_of_bfb       = TRUE,
  selection_type           = "constant",
  saturation_K             = 10
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
    first_round_of_bfb      = first_round_of_bfb,
    selection_type          = selection_type,
    saturation_K            = saturation_K
  )
}


#' Summarize a copy-number matrix into scalar statistics
#'
#' Computes a fixed set of summary statistics from a cells-by-bins
#' copy-number matrix, centered on one "hotspot" column of interest:
#' mean/max/variance of copy number at the hotspot, the fraction of
#' cells with a hotspot gain above \code{base_value}, the genome-wide
#' maximum copy number, and the mean number of breakpoints per cell.
#' Used internally by \code{simulate_serial_passages()} to summarize
#' each passage.
#'
#' @param cna_matrix Integer matrix (cells x bins), e.g. one allele's
#'   worth of copy-number calls for a set of cells.
#' @param hotspot_col 1-based column index of the position of interest
#'   (clamped to \code{ncol(cna_matrix)} if out of range).
#' @param base_value Baseline copy number used to define a "gain" at
#'   the hotspot (\code{hotspot_cn > base_value}). Default 1.
#' @param n_alive Optional: number of cells alive before any
#'   subsampling, stored as-is in the output for downstream reference.
#'
#' @return A named numeric vector: \code{hotspot_mean_cn},
#'   \code{hotspot_max_cn}, \code{hotspot_fraction_amp},
#'   \code{hotspot_cn_var}, \code{max_cn_global}, \code{mean_breakpoints},
#'   \code{n_alive}, \code{n_cells}. All \code{NA} (with \code{n_cells = 0})
#'   if \code{cna_matrix} is empty.
#' @examples
#' m <- rbind(cell_1 = c(1, 1, 3, 3, 1), cell_2 = c(1, 2, 2, 1, 1))
#' summarise_clone(m, hotspot_col = 3)
#' @export
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

