
#' Main Function for Genomic Distance Calculation and Phylogenetic Tree Construction
#'
#' This function computes genomic distances between samples with flexible
#' distance-function selection, processes the data, and constructs a phylogenetic
#' tree. It supports distinct metrics for greedy (G) and breakage–fusion–bridge
#' (BFB; B) distance calculations.
#'
#' @param data Input data frame with columns:
#' \itemize{
#'   \item \code{cell_id} — Cell identifier
#'   \item \code{chr} — Chromosome name
#'   \item \code{start} — Start of genomic bin
#'   \item \code{end} — End of genomic bin
#'   \item \code{CN} — Total copy number
#'   \item \code{A} — Copy number of allele A
#'   \item \code{B} — Copy number of allele B
#' }
#' @param chromosomes Vector of chromosomes to include (default: \code{c(1:22, "X", "Y")}).
#' @param alleles Alleles to consider (default: \code{c("A","B")}; alternative: \code{"CN"}).
#' @param k_jitter_fix Numeric jitter factor for numerical stability (default: \code{0}).
#' @param bfb_penalty Penalty to apply to BFB events (default: \code{0}). Ignored if
#'   \code{b_dist_func = NULL}.
#' @param tree_func Function used to build the tree from the final distance
#'   matrix (default: \code{ape::nj}).
#' @param fillna Value to fill \code{NA} entries in pre-processing (default: \code{0}).
#' @param g_dist_func Name of the greedy (G) distance function: one of
#'   \code{"greedy_fast"} (C++, default), \code{"greedy"} or \code{"G"}
#'   (pure-R reference implementations).
#' @param b_dist_func Name of the BFB (B) distance function: one of
#'   \code{"bfb_fast"} (C++, default), \code{"bfb"}, \code{"avg"}, \code{"A"}
#'   or \code{"B_contig"}. \strong{Set to \code{NULL} to disable the BFB stage}
#'   and run a greedy-only analysis.
#' @param keep_matrices Logical; also return the per-chromosome, per-allele
#'   greedy and BFB distance matrices (\code{greedy_Ds}, \code{avg_Ds}).
#'   These take \eqn{2 \times} (chromosomes \eqn{\times} alleles) matrices of
#'   size cells \eqn{\times} cells, so they are not kept by default.
#'   Default: \code{FALSE}.
#' @param wgd Logical; also consider whole-genome duplication (WGD) when
#'   comparing cells, and label WGD branches on the tree. See Details.
#'   Default: \code{FALSE}.
#' @param wgd_penalty Cost of one WGD event, in the same units as the
#'   gain/loss events. Default: 1.
#' @param wgd_counts Optional named integer vector (names = cell ids) giving
#'   the number of WGDs (0, 1 or 2) of each cell. If \code{NULL}, it is
#'   estimated from ploidy as \code{round(log2(mean CN / 2))}, clamped to 0-2.
#' @param verbose Logical; print progress messages. Default: \code{TRUE}.
#' @param ... Currently unused.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{tree} — The constructed phylogenetic tree (rooted, diploid removed)
#'   \item \code{all_input_Xs} — Processed input data
#'   \item \code{D} — Final distance matrix
#'   \item \code{greedy_Ds} — Per-chromosome/allele greedy distance matrices
#'     (\code{NULL} unless \code{keep_matrices = TRUE})
#'   \item \code{avg_Ds} — Per-chromosome/allele BFB distance matrices
#'     (\code{NULL} unless \code{keep_matrices = TRUE});
#'     \strong{if \code{b_dist_func = NULL}, then \code{avg_Ds} equals \code{greedy_Ds}}
#'   \item \code{g_dist_func} — Name of the G distance function used
#'   \item \code{b_dist_func} — Name of the B distance function used (or \code{NULL})
#'   \item \code{reconstructions} — Nested list \code{[[chr]][[allele]]} of ancestral CN
#'     state reconstructions per tree branch, as returned by
#'     \code{compute_reconstructions()}.  Each element contains \code{profiles},
#'     \code{merged_profiles}, and \code{deltas} (the per-branch BFB signal).
#'   \item \code{history} — Reconstructed events per branch, chromosome and allele
#'     (\code{"gain"}, \code{"loss"}, \code{"bfb"}, \code{"wgd"}, ...)
#'   \item \code{wgd} — Only with \code{wgd = TRUE}: a list with \code{counts}
#'     (WGDs per cell), \code{penalty}, \code{k_node} (WGD count of every tree
#'     node, ape numbering), \code{root_count} (WGDs above the root) and
#'     \code{edges} (one row per WGD branch, with its descendant \code{cells})
#' }
#'
#' @details
#' The pipeline performs:
#' \enumerate{
#'   \item Validation of selected distance functions
#'   \item Input pre-processing and diploid augmentation
#'   \item Computation of greedy (G) distances
#'   \item \strong{Optional BFB stage}: if \code{b_dist_func} is provided, BFB (B) distances
#'         are computed and merged with G; if \code{b_dist_func = NULL}, the BFB stage
#'         is \emph{skipped} and \code{avg_Ds <- greedy_Ds} (greedy-only analysis).
#'         In the \code{NULL} case, \code{bfb_penalty} is ignored.
#'   \item Merging to minimal distances and optional allele summation
#'   \item Tree construction via \code{tree_func}
#' }
#'
#' With the default distance functions, all pairwise distances for a
#' chromosome and allele are computed in a single C++ call, once per pair of
#' distinct copy-number profiles, and only their sum over chromosomes and
#' alleles is kept. The pure-R alternatives are computed pair by pair and
#' are much slower; they are kept as references.
#'
#' \strong{Whole-genome duplication.} A WGD doubles every chromosome at once,
#' but the gain/loss distance would count it as many separate events. With
#' \code{wgd = TRUE}, the distance between two cells is the cheaper of two
#' explanations: (i) no WGD, the usual distance; (ii) one cell descends from
#' the other's lineage doubled once (or twice) plus gain/loss events, costing
#' \code{wgd_penalty} per WGD plus the events before and after the doubling.
#' The WGD is counted once for the whole genome. After the tree is built,
#' each cell's WGD count (\code{wgd_counts}) is placed on the tree by
#' parsimony (doublings only, never halving); branches where the count
#' increases are reported in \code{wgd$edges} and as \code{"wgd"} events in
#' \code{history}, and their nodes are left out of the BFB test of
#' \code{\link{detect_bfb_branches}}. Only the default distance functions
#' support \code{wgd = TRUE}.
#'
#' Using \code{b_dist_func = NULL} is useful for ablation studies, speed-ups, or when
#' BFB modelling is not desired; results reduce to the greedy metric.
#'
#' @examples
#' \dontrun{
#' sim <- bridge_sim(chromosomes = "8", bfb_allele = "8:A",
#'                   max_cells = 128, lambda = 2)
#'
#' res <- fit(
#'   data    = sim$cna_data,
#'   alleles = c("A", "B")
#' )
#'
#' # Compare inferred vs true tree (Robinson-Foulds distance)
#' phangorn::RF.dist(sim$tree, res$tree, normalize = TRUE)
#'
#' # Greedy-only run (no BFB modelling, faster)
#' res_g <- fit(data = sim$cna_data, b_dist_func = NULL)
#'
#' # With whole-genome duplication
#' sim_w <- bridge_sim(chromosomes = "8", bfb_allele = "8:A", max_cells = 128,
#'                     wgd_available = 1, wgd_probability = 0.1)
#' res_w <- fit(sim_w$cna_data, chromosomes = "8", wgd = TRUE)
#' res_w$wgd$edges
#' }
#'
#' @export
fit = function(data,
               chromosomes = c(1:22, "X", "Y"),
               alleles = c("A", "B"),
               k_jitter_fix = 0,
               bfb_penalty = 0,
               tree_func = ape::nj,
               fillna = 0,
               g_dist_func = "greedy_fast",
               b_dist_func = "bfb_fast",
               keep_matrices = FALSE,
               wgd = FALSE,
               wgd_penalty = 1,
               wgd_counts = NULL,
               verbose = TRUE,
               ...) {
  say = function(...) if (isTRUE(verbose)) message(...)

  # Fail early on misspelled distance names
  get_g_dist(g_dist_func)
  if (!is.null(b_dist_func)) get_b_dist(b_dist_func)

  fast_path = identical(g_dist_func, "greedy_fast") &&
    (is.null(b_dist_func) || identical(b_dist_func, "bfb_fast"))
  if (isTRUE(wgd) && !fast_path) {
    stop("wgd = TRUE requires the default distances (g_dist_func = \"greedy_fast\", ",
         "b_dist_func = \"bfb_fast\" or NULL).", call. = FALSE)
  }
  if (!is.numeric(wgd_penalty) || length(wgd_penalty) != 1 || wgd_penalty < 0) {
    stop("wgd_penalty must be a single non-negative number.", call. = FALSE)
  }

  data = dplyr::bind_rows(data, create_diploid_data(data))
  say("Pre-processing input...")
  all_input_Xs = process_input_data(data, chromosomes, alleles, k_jitter_fix, fillna)

  if (fast_path) {
    say("Computing distances...")
    dists = compute_fast_distances(all_input_Xs, chromosomes, alleles,
                                   bfb_penalty = bfb_penalty,
                                   use_bfb = !is.null(b_dist_func),
                                   keep_matrices = keep_matrices,
                                   wgd = isTRUE(wgd), wgd_penalty = wgd_penalty)
    D = dists$D
    greedy_Ds = dists$greedy_Ds
    avg_Ds = dists$avg_Ds
  } else {
    say("Computing G matrices using ", g_dist_func, "...")
    greedy_Ds = compute_greedy_distances(all_input_Xs, chromosomes, alleles, g_dist_func)

    if (is.null(b_dist_func)) {
      avg_Ds = greedy_Ds
    } else {
      say("Computing B matrices using ", b_dist_func, "...")
      avg_Ds = compute_avg_distances(all_input_Xs, chromosomes, alleles, bfb_penalty, b_dist_func)
    }

    say("Computing final D matrix...")
    min_Ds = find_minimal_distances(greedy_Ds, avg_Ds)
    min_Ds = sum_across_alleles(min_Ds, alleles)
    D = Reduce("+", min_Ds)
    if (is.list(D)) D = D[[1]]
    if (!keep_matrices) greedy_Ds = avg_Ds = NULL
  }

  # The diagonal is Inf by construction (see compute_distance_matrix());
  # only off-diagonal infinities indicate a problem.
  diag(D) = 0
  if (any(is.infinite(D))) {
    warning("Infinite distances detected; replacing with max finite distance. Check input data for pathological cell pairs.")
    D[is.infinite(D)] = max(D[is.finite(D)], na.rm = TRUE)
  }

  say("Building tree...")
  tree = tree_func(D)
  tree = ape::root(tree, outgroup = "diploid", resolve.root = TRUE)
  tree = ape::drop.tip(tree, "diploid")

  wgd_info = NULL
  if (isTRUE(wgd)) {
    say("Placing WGD events...")
    cells = tree$tip.label
    counts = if (is.null(wgd_counts)) estimate_wgd_counts(data)[cells] else
      validate_wgd_counts(wgd_counts, cells)
    k_node = assign_wgd_counts(tree, counts, penalty = wgd_penalty)
    wgd_info = list(counts = counts, penalty = wgd_penalty, k_node = k_node,
                    root_count = k_node[length(cells) + 1],
                    edges = wgd_edges(tree, k_node))
  }

  res = list(
    tree = tree,
    all_input_Xs = all_input_Xs,
    D = D,
    greedy_Ds = greedy_Ds,
    avg_Ds = avg_Ds,
    g_dist_func = g_dist_func,
    b_dist_func = b_dist_func,
    reconstructions = list()
  )
  if (!is.null(wgd_info)) res$wgd = wgd_info

  say("Reconstructing events...")
  recon = compute_reconstructions(fit = res, chromosomes = chromosomes, alleles = alleles)
  res$reconstructions = recon$reconstructions
  res$tree$edge.length = recon$edge.length
  res$history = recon$history

  res
}
