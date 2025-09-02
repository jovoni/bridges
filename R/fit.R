
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
#' @param g_dist_func Name of the greedy (G) distance function. Must be one of
#'   \code{names(G_DISTS)} (default: \code{"greedy_fast"}).
#' @param b_dist_func Name of the BFB (B) distance function. Must be one of
#'   \code{names(B_DISTS)}. \strong{Set to \code{NULL} to disable the BFB stage}
#'   and run a greedy-only analysis (default: \code{"bfb_fast"}).
#' @param ... Additional arguments passed to downstream helpers.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{tree} — The constructed phylogenetic tree (rooted, diploid removed)
#'   \item \code{all_input_Xs} — Processed input data
#'   \item \code{D} — Final distance matrix
#'   \item \code{greedy_Ds} — Per-chromosome/allele greedy distance matrices
#'   \item \code{avg_Ds} — Per-chromosome/allele balanced (G vs B) distance matrices;
#'     \strong{if \code{b_dist_func = NULL}, then \code{avg_Ds} equals \code{greedy_Ds}}
#'   \item \code{g_dist_func} — Name of the G distance function used
#'   \item \code{b_dist_func} — Name of the B distance function used (or \code{NULL})
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
#' Available functions can be inspected with \code{names(G_DISTS)} and \code{names(B_DISTS)}.
#' Using \code{b_dist_func = NULL} is useful for ablation studies, speed-ups, or when
#' BFB modelling is not desired; results reduce to the greedy metric.
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
               ...) {

  data = dplyr::bind_rows(data, create_diploid_data(data))
  message("Pre-processing input...")
  all_input_Xs = process_input_data(data, chromosomes, alleles, k_jitter_fix, fillna)

  message(paste("Computing G matrices using", g_dist_func, "..."))
  greedy_Ds = compute_greedy_distances(all_input_Xs, chromosomes, alleles, g_dist_func)

  if (is.null(b_dist_func)) {
    avg_Ds = greedy_Ds
  } else {
    message(paste("Computing B matrices using", b_dist_func, "..."))
    avg_Ds = compute_avg_distances(all_input_Xs, chromosomes, alleles, bfb_penalty, b_dist_func)
  }

  message("Computing final D matrix...")
  min_Ds = find_minimal_distances(greedy_Ds, avg_Ds)
  min_Ds = sum_across_alleles(min_Ds, alleles)
  D = Reduce("+", min_Ds)
  if (is.list(D)) D = D[[1]]
  D[D == Inf] = 0

  message("Building tree...")
  tree = tree_func(D)
  tree = ape::root(tree, outgroup = "diploid", resolve.root = TRUE)
  tree = ape::drop.tip(tree, "diploid")

  res = list(
    tree = tree,
    all_input_Xs = all_input_Xs,
    D = D,
    greedy_Ds = greedy_Ds,
    avg_Ds = avg_Ds,
    g_dist_func = g_dist_func,
    b_dist_func = b_dist_func,
    reconstructions = list()  # new field
  )

  message("Reconstruct events...")
  recon = compute_reconstructions(fit = res, chromosomes = chromosomes, alleles = alleles)
  res$reconstructions = recon$reconstructions
  res$tree$edge.length = recon$edge.length
  res$history = recon$history

  res
}
