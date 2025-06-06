#' Main Function for Genomic Distance Calculation and Phylogenetic Tree Construction
#'
#' This is the enhanced main function that computes genomic distances between samples
#' with flexible distance function selection, processes the data, and constructs
#' a phylogenetic tree. It supports different distance metrics for greedy (G) and
#' BFB (B) distance calculations.
#'
#' @param data Input data containing genomic information. Should be in a format
#'   that can be processed by `process_input_data()`.
#' @param chromosomes Vector of chromosomes to include in analysis (default: 1:22, "X", "Y").
#' @param alleles Alleles to consider (default: c(A", "B"), alternative: "CN").
#' @param k_jitter_fix Jitter factor for numerical stability (default: 0).
#' @param target_vals Target values for distance calculation (default: c(2, 1, 1)).
#' @param bfb_penalty Penalty for breakage-fusion-bridge events (default: 0).
#' @param sum_across_minor_alleles Logical indicating whether to sum distances across
#'   minor alleles (default: TRUE).
#' @param tree_func Function to use for tree construction (default: ape::nj).
#' @param fillna Value to fill NA entries (default: 0).
#' @param g_dist_func Distance function to use for greedy calculation. Must be one
#'   of the available G functions (default: "G_with_steps").
#' @param b_dist_func Distance function to use for BFB calculation. Must be one
#'   of the available B functions (default: "A_contig").
#' @param ... Additional arguments passed to downstream functions.
#'
#' @return A list containing:
#' \itemize{
#'   \item tree - The constructed phylogenetic tree
#'   \item all_input_Xs - Processed input data
#'   \item D - Final distance matrix
#'   \item greedy_Ds - Greedy distance matrices
#'   \item avg_Ds - Balanced distance matrices
#'   \item g_dist_func - Name of the used G distance function
#'   \item b_dist_func - Name of the used B distance function
#' }
#'
#' @details
#' The function performs the following steps:
#' 1. Validates the selected distance functions
#' 2. Processes input data
#' 3. Computes distances using selected functions
#' 4. Merges distance matrices
#' 5. Optionally sums across minor alleles
#' 6. Constructs phylogenetic tree
#'
#' Available distance functions can be checked via `names(G_DISTS)` and `names(B_DISTS)`.
#'
#' @export
fit = function(data,
                chromosomes = c(1:22, "X", "Y"),
                alleles = c("CN", "A", "B"),
                k_jitter_fix = 0,
                target_vals = c(2, 1, 1),
                bfb_penalty = 0,
                sum_across_minor_alleles = TRUE,
                tree_func = ape::nj,
                fillna = 0,
                # New parameters for distance function selection
                g_dist_func = "G_with_steps",
                b_dist_func = "A_contig",
                ...) {

  # Process data
  message("Pre-processing input...")
  all_input_Xs = process_input_data(data, chromosomes, alleles, k_jitter_fix, fillna)

  # Compute distances with selected functions
  message(paste("Computing G matrices using", g_dist_func, "..."))
  greedy_Ds = compute_greedy_distances(all_input_Xs, chromosomes, alleles, target_vals, g_dist_func)

  message(paste("Computing B matrices using", b_dist_func, "..."))
  avg_Ds = compute_avg_distances(all_input_Xs, chromosomes, alleles, bfb_penalty, b_dist_func)

  # Merge distances
  message("Computing final D matrix...")
  min_Ds = find_minimal_distances(greedy_Ds, avg_Ds)

  # Sum across minor alleles
  if (sum_across_minor_alleles) {
    min_Ds = sum_across_alleles(min_Ds, alleles)
  }

  # Sum across chromosomes
  D = Reduce("+", min_Ds)
  if (is.list(D)) {D = D[[1]]}
  D[D == Inf] = 0

  # Build tree
  message("Building tree...")
  tree = tree_func(D)

  # Return results with metadata about chosen functions
  list(
    tree = tree,
    all_input_Xs = all_input_Xs,
    D = D,
    greedy_Ds = greedy_Ds,
    avg_Ds = avg_Ds,
    g_dist_func = g_dist_func,
    b_dist_func = b_dist_func
  )
}
