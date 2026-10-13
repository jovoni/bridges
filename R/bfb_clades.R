#' Find clades of the tree enriched in BFB-like nodes
#'
#' A tree-aware version of \code{\link{detect_bfb_branches}}. Instead of one
#' test over all internal nodes, every clade of the inferred tree is a
#' candidate, so a group of BFB-like nodes concentrated in one part of the
#' tree (for example BFB in a subpopulation only) is not diluted by the rest
#' of the tree. No clones or clades need to be defined in advance.
#'
#' @param fit A result of \code{\link{fit}}.
#' @param threshold Background probability that an internal node looks
#'   BFB-like by chance, as in \code{\link{detect_bfb_branches}}.
#'   Default: 0.005
#' @param min_nodes Minimum number of tested internal nodes a clade must
#'   contain to be a candidate. Default: 5
#' @param n_sim Number of simulations of the null model used to correct for
#'   scanning many clades. Default: 999
#' @param seed Optional integer seed for the simulations. The random-number
#'   state of the session is restored afterwards. Default: \code{NULL}
#'
#' @return A tibble with one row per chromosome and allele:
#' \describe{
#'   \item{chr, allele}{Chromosome and allele.}
#'   \item{p.scan}{P-value of the best clade, corrected for scanning all
#'     clades (Monte Carlo).}
#'   \item{adj.pval}{Benjamini-Hochberg adjusted \code{p.scan} across all
#'     chromosome-allele combinations.}
#'   \item{clade_node}{Node of the best clade (ape numbering of
#'     \code{fit$tree}).}
#'   \item{clade_nodes, clade_bfb_nodes}{Tested internal nodes and BFB-like
#'     nodes in the best clade.}
#'   \item{clade_p}{Uncorrected binomial p-value of the best clade.}
#'   \item{clade_n_cells, clade_cells}{Number and ids of the cells in the
#'     best clade.}
#'   \item{total_nodes, total_bfb_nodes}{The same counts over the whole tree,
#'     as used by \code{\link{detect_bfb_branches}}.}
#'   \item{threshold}{The background rate used.}
#' }
#'
#' @details
#' For each chromosome and allele, an internal node is BFB-like when its
#' children are better explained by a BFB event than by greedy gains and
#' losses (\eqn{\delta > 0}, see \code{\link{detect_bfb_branches}}). Nodes
#' where a WGD happens (\code{fit(wgd = TRUE)}) are not tested.
#'
#' For every clade with at least \code{min_nodes} tested internal nodes, the
#' number \eqn{k} of BFB-like nodes among its \eqn{n} tested internal nodes is
#' compared with the background rate:
#' \eqn{p_c = P(\mathrm{Bin}(n, \mathrm{threshold}) \ge k)}. The best clade is
#' the one with the smallest \eqn{p_c}; the whole tree is one of the
#' candidates, so the usual genome-wide test is included. To account for
#' having looked at many nested clades, the smallest \eqn{p_c} is compared
#' with its distribution when every internal node is BFB-like independently
#' with probability \code{threshold} on the same tree (\code{n_sim}
#' simulations).
#'
#' @examples
#' set.seed(1)
#' sim <- bridge_sim(chromosomes = "1", bfb_allele = "1:A", max_cells = 60)
#' res <- fit(sim$cna_data, chromosomes = "1", verbose = FALSE)
#' detect_bfb_clades(res, n_sim = 199, seed = 1)
#' @export
detect_bfb_clades = function(fit, threshold = 0.005, min_nodes = 5, n_sim = 999,
                             seed = NULL) {
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0 || threshold >= 1) {
    stop("threshold must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(min_nodes) || length(min_nodes) != 1 || min_nodes < 1) {
    stop("min_nodes must be a single number >= 1.", call. = FALSE)
  }
  if (!is.numeric(n_sim) || length(n_sim) != 1 || n_sim < 1) {
    stop("n_sim must be a single number >= 1.", call. = FALSE)
  }
  if (!is.null(seed)) {
    old_seed = if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit({
      if (is.null(old_seed)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old_seed, envir = globalenv())
    }, add = TRUE)
    set.seed(seed)
  }

  tree = fit$tree
  if (!ape::is.rooted(tree)) stop("fit$tree must be rooted.", call. = FALSE)
  n_tip = length(tree$tip.label)
  n_node = n_tip + tree$Nnode
  post = ape::reorder.phylo(tree, "postorder")$edge
  clade_tips = clade_tip_sets(tree)

  alleles = names(fit$all_input_Xs[[1]])
  chromosomes = names(fit$all_input_Xs)
  rows = list()
  for (allele in alleles) {
    for (chr in chromosomes) {
      r = fit$reconstructions[[chr]][[allele]]
      if (is.null(r)) r = reconstruct_tree(fit, chr, allele)
      rows[[length(rows) + 1]] = bfb_clade_scan(r$deltas, tree, post, n_tip, n_node, clade_tips,
                                                threshold, min_nodes, n_sim, chr, allele)
    }
  }
  out = dplyr::bind_rows(rows)
  out$adj.pval = stats::p.adjust(out$p.scan, method = "BH")
  dplyr::relocate(out, "adj.pval", .after = "p.scan")
}

# Tip labels of the clade below every internal node (index = node number)
clade_tip_sets = function(tree) {
  n_tip = length(tree$tip.label)
  out = vector("list", n_tip + tree$Nnode)
  for (i in seq_len(n_tip)) out[[i]] = i
  post = ape::reorder.phylo(tree, "postorder")$edge
  for (e in seq_len(nrow(post))) out[[post[e, 1]]] = c(out[[post[e, 1]]], out[[post[e, 2]]])
  lapply(out, function(ix) tree$tip.label[ix])
}

# Sum a per-node quantity (rows of x: nodes, columns: replicates) over each
# clade, by one post-order pass over the edges
clade_sums = function(x, post) {
  for (e in seq_len(nrow(post))) {
    x[post[e, 1], ] = x[post[e, 1], ] + x[post[e, 2], ]
  }
  x
}

bfb_clade_scan = function(deltas, tree, post, n_tip, n_node, clade_tips, threshold,
                          min_nodes, n_sim, chr, allele) {
  d = unlist(deltas)
  internal = n_tip + seq_len(tree$Nnode)
  delta = rep(NA_real_, n_node)
  delta[as.integer(names(d))] = d
  tested = !is.na(delta) & seq_len(n_node) %in% internal   # WGD nodes (NA) are not tested
  bfb = tested & delta > 0

  n_c = clade_sums(matrix(as.numeric(tested), ncol = 1), post)[, 1]
  k_c = clade_sums(matrix(as.numeric(bfb), ncol = 1), post)[, 1]
  candidates = internal[n_c[internal] >= min_nodes]
  root = n_tip + 1
  if (!length(candidates)) candidates = root

  tail_p = function(k, n) stats::pbinom(k - 1, n, threshold, lower.tail = FALSE)
  p_c = tail_p(k_c[candidates], n_c[candidates])
  # best clade: smallest p, ties resolved towards the smaller clade
  best = candidates[order(p_c, n_c[candidates])[1]]
  s_obs = min(p_c)

  # Null: each tested node BFB-like independently with probability threshold
  sims = matrix(0, n_node, n_sim)
  idx = which(tested)
  sims[idx, ] = stats::rbinom(length(idx) * n_sim, 1, threshold)
  k_sim = clade_sums(sims, post)[candidates, , drop = FALSE]
  s_sim = apply(matrix(tail_p(k_sim, n_c[candidates]), nrow = length(candidates)), 2, min)
  p_scan = (1 + sum(s_sim <= s_obs)) / (n_sim + 1)

  dplyr::tibble(
    chr = chr,
    allele = allele,
    p.scan = p_scan,
    clade_node = best,
    clade_nodes = n_c[best],
    clade_bfb_nodes = k_c[best],
    clade_p = tail_p(k_c[best], n_c[best]),
    clade_n_cells = length(clade_tips[[best]]),
    clade_cells = list(clade_tips[[best]]),
    total_nodes = n_c[root],
    total_bfb_nodes = k_c[root],
    threshold = threshold
  )
}
