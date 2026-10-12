# Prototype WGD-aware genome-wide distance (not part of the package).
# Requires the package loaded with devtools::load_all() and
# Rcpp::sourceCpp("dev/wgd/wgd_proto.cpp").

# Number of WGDs per cell estimated from ploidy: round(log2(mean CN / 2)),
# clamped to 0..2. `data` is a cna_data table (CN column, or A + B).
estimate_wgd_class <- function(data) {
  cn <- if ("CN" %in% names(data)) data$CN else data$A + data$B
  ploidy <- tapply(cn, data$cell_id, mean)
  k <- round(log2(pmax(ploidy, 1e-9) / 2))
  pmin(pmax(k, 0), 2)
}

# Per-chromosome distances summed over the genome:
#   S0 = sum of pmin(G, B) (current fit() distance)
#   W1[i, j] = sum of W(cell i -> cell j, f = 2), W2 the same with f = 4
wgd_components <- function(all_input_Xs, chromosomes, alleles) {
  alleles <- alleles[alleles %in% names(all_input_Xs[[1]])]
  chromosomes <- chromosomes[chromosomes %in% names(all_input_Xs)]
  cells <- rownames(all_input_Xs[[chromosomes[1]]][[alleles[1]]])
  n <- length(cells)
  S0 <- W1 <- W2 <- matrix(0, n, n, dimnames = list(cells, cells))
  for (chromosome in chromosomes) for (allele in alleles) {
    X <- all_input_Xs[[chromosome]][[allele]][cells, , drop = FALSE]
    storage.mode(X) <- "double"
    m <- wgd_matrices(X)
    S0 <- S0 + m$S0
    W1 <- W1 + m$W1
    W2 <- W2 + m$W2
  }
  list(S0 = S0, W1 = W1, W2 = W2)
}

# Genome-wide distance from the components, with optional WGD hypotheses.
#   mode = "none":   D = S0
#   mode = "all":    D = min(S0, w + W1(a->b), w + W1(b->a), 2w + W2(a->b), 2w + W2(b->a))
#   mode = "ploidy": only the hypothesis matching the difference in estimated
#                    WGD class, in the direction low -> high class
# Returns D and the per-pair winning hypothesis (0 = no WGD, 1 or 2 WGDs).
wgd_combine <- function(comp, w = 1, mode = c("none", "all", "ploidy"), wgd_class = NULL) {
  mode <- match.arg(mode)
  S0 <- comp$S0; W1 <- comp$W1; W2 <- comp$W2
  cells <- rownames(S0); n <- length(cells)
  if (mode == "all") {
    hyp <- list(S0, w + W1, w + t(W1), 2 * w + W2, 2 * w + t(W2))
    label <- c(0, 1, 1, 2, 2)
  } else if (mode == "ploidy") {
    k <- wgd_class[cells]
    dk <- outer(k, k, function(ki, kj) kj - ki)   # > 0: column cell has more WGDs
    up1 <- ifelse(dk == 1, w + W1, Inf)
    up2 <- ifelse(dk == 2, 2 * w + W2, Inf)
    hyp <- list(S0, pmin(up1, t(up1)), pmin(up2, t(up2)))
    label <- c(0, 1, 2)
  } else {
    hyp <- list(S0); label <- 0
  }
  D <- Reduce(pmin, hyp)
  winner <- matrix(0, n, n, dimnames = list(cells, cells))
  for (h in rev(seq_along(hyp))) winner[hyp[[h]] == D] <- label[h]
  diag(D) <- 0
  list(D = D, winner = winner)
}

# Tree from a distance matrix, rooted on the synthetic diploid cell as in fit()
tree_from_distance <- function(D) {
  tree <- ape::nj(D)
  tree <- ape::root(tree, outgroup = "diploid", resolve.root = TRUE)
  ape::drop.tip(tree, "diploid")
}

# True WGD ancestry of sampled cells from bridge_sim()$cell_history.
# Returns, per sampled cell, the number of WGDs on its lineage, and for each
# WGD event (identified by the parent cell where it happened) the set of
# sampled descendants.
true_wgd <- function(sim) {
  h <- sim$cell_history
  parent <- stats::setNames(h$parent_id, h$cell_id)
  wgd <- stats::setNames(h$wgd_event, h$cell_id)
  cells <- names(sim$cells)
  n_wgd <- integer(length(cells)); names(n_wgd) <- cells
  events <- list()
  for (cell in cells) {
    node <- cell
    while (!is.na(node) && node != "root") {
      if (isTRUE(wgd[[node]])) {
        n_wgd[cell] <- n_wgd[cell] + 1L
        ev <- parent[[node]]   # both daughters of this division are doubled
        events[[ev]] <- c(events[[ev]], cell)
      }
      node <- parent[[node]]
    }
  }
  list(n_wgd = n_wgd, events = lapply(events, unique))
}

# ---- WGD reconstruction -----------------------------------------------------

# Parsimony assignment of WGD counts (0..kmax) to every node of a rooted
# tree. Tips are fixed to `k_tip` (named by tip label). Along an edge the
# count may only increase, at cost w per WGD; the root's parent has k = 0.
# Returns integer k for nodes 1..(Ntip + Nnode), in ape numbering.
assign_wgd_counts <- function(tree, k_tip, w = 1, kmax = 2) {
  n_tip <- length(tree$tip.label)
  n_node <- n_tip + tree$Nnode
  ks <- 0:kmax
  cost <- matrix(Inf, n_node, kmax + 1)
  cost[cbind(seq_len(n_tip), k_tip[tree$tip.label] + 1)] <- 0
  # cost of a child subtree given the parent's count kp
  child_cost <- function(c, kp) min(cost[c, ks >= kp] + w * (ks[ks >= kp] - kp))
  post <- ape::reorder.phylo(tree, "postorder")$edge
  for (node in unique(post[, 1])) {
    kids <- post[post[, 1] == node, 2]
    for (k in ks) cost[node, k + 1] <- sum(vapply(kids, child_cost, numeric(1), kp = k))
  }
  root <- n_tip + 1
  k_node <- integer(n_node)
  k_node[root] <- ks[which.min(cost[root, ] + w * ks)]
  pre <- ape::reorder.phylo(tree, "cladewise")$edge
  for (i in seq_len(nrow(pre))) {
    p <- pre[i, 1]; c <- pre[i, 2]
    ok <- ks >= k_node[p]
    total <- cost[c, ok] + w * (ks[ok] - k_node[p])
    k_node[c] <- ks[ok][which.min(total)]   # ties: the smallest count
  }
  k_node
}

# Bottom-up reconstruction of ancestral profiles with WGD-aware branches.
# `all_input_Xs` as in fit() (cells x segments per chromosome/allele),
# `k_node` from assign_wgd_counts(). For every chromosome/allele, the
# ancestor of two children at the parent's count is reconstructed with the
# existing G/B logic; a child with more WGDs than its parent is compared
# after doubling the ancestor (wgd_ancestor()). Returns per-edge event
# counts summed over the genome, the WGD-labelled edges, and the root
# k (WGD on the stem above the root if > 0).
reconstruct_wgd <- function(tree, all_input_Xs, chromosomes, alleles, k_node, w = 1) {
  n_tip <- length(tree$tip.label)
  n_node <- n_tip + tree$Nnode
  edge <- tree$edge
  post <- ape::reorder.phylo(tree, "postorder")$edge
  internal <- unique(post[, 1])
  edge_events <- numeric(nrow(edge))   # gain/loss/BFB events per edge
  root <- n_tip + 1

  for (chromosome in chromosomes) for (allele in alleles) {
    X <- all_input_Xs[[chromosome]][[allele]]
    if (is.null(X)) next
    tv <- if (allele == "CN") 2 else 1
    prof <- vector("list", n_node)
    for (i in seq_len(n_tip)) prof[[i]] <- as.numeric(X[tree$tip.label[i], ])
    for (node in internal) {
      kids <- post[post[, 1] == node, 2]
      if (length(kids) != 2) stop("reconstruct_wgd expects a binary tree")
      kn <- k_node[node]
      dk <- k_node[kids] - kn
      a <- prof[[kids[1]]]; b <- prof[[kids[2]]]
      ei <- match(kids, edge[, 2])
      if (all(dk == 0)) {
        g <- greedy_distance_cpp(a, b, target_val = tv * 2^kn)
        a_first <- mean(a) > mean(b)   # same ordering as reconstruct_tree()
        bb <- if (a_first) greedy_bfb_distance_cpp(a, b, 0) else greedy_bfb_distance_cpp(b, a, 0)
        if (bb$cost < g$cost) {
          prof[[node]] <- bb$ancestor
          edge_events[ei] <- edge_events[ei] +
            if (a_first) c(bb$a_cost, bb$b_cost) else c(bb$b_cost, bb$a_cost)
        } else {
          prof[[node]] <- g$ancestor
          edge_events[ei] <- edge_events[ei] + c(g$a_cost, g$b_cost)
        }
      } else if (sum(dk > 0) == 1) {
        hi <- which(dk > 0); lo <- 3 - hi
        r <- wgd_ancestor(list(a, b)[[lo]], list(a, b)[[hi]], 2^dk[hi])
        prof[[node]] <- r$x
        edge_events[ei[lo]] <- edge_events[ei[lo]] + r$steps
        edge_events[ei[hi]] <- edge_events[ei[hi]] + r$after
      } else {
        # both children doubled relative to the parent: compare halved profiles
        ha <- round(a / 2^dk[1]); hb <- round(b / 2^dk[2])
        g <- greedy_distance_cpp(ha, hb, target_val = tv * 2^kn)
        prof[[node]] <- g$ancestor
        edge_events[ei] <- edge_events[ei] +
          c(g$a_cost + greedy_distance_cpp(2^dk[1] * g$ancestor, a, tv)$cost,
            g$b_cost + greedy_distance_cpp(2^dk[2] * g$ancestor, b, tv)$cost)
      }
    }
  }

  dk_edge <- k_node[edge[, 2]] - k_node[edge[, 1]]
  wgd_edges <- which(dk_edge > 0)
  list(
    edge_events = edge_events,
    n_wgd_edge = dk_edge,
    wgd_clades = lapply(wgd_edges, function(e) {
      child <- edge[e, 2]
      if (child <= n_tip) tree$tip.label[child] else ape::extract.clade(tree, child)$tip.label
    }),
    root_k = k_node[root],
    total_events = sum(edge_events) + w * (sum(dk_edge) + k_node[root])
  )
}

# Tips descending from each WGD edge, including the root stem
wgd_clade_sets <- function(tree, rec) {
  sets <- rec$wgd_clades
  if (rec$root_k > 0) sets <- c(sets, list(tree$tip.label))
  sets
}

# For every true set of WGD descendants, the best Jaccard index with a
# predicted WGD clade; also the number of predicted WGD edges.
score_wgd_labels <- function(true_sets, pred_sets) {
  jac <- function(x, y) length(intersect(x, y)) / length(union(x, y))
  best <- vapply(true_sets, function(t) if (length(pred_sets)) max(vapply(pred_sets, jac, 0, y = t)) else 0, 0)
  list(best_jaccard = best, n_pred = length(pred_sets), n_true = length(true_sets))
}

# Clone coherence of a tree: for each clone, the best F1 between its cells
# and any clade (or single tip) of the tree. Returns per-clone F1 and whether
# the clone is monophyletic.
clone_coherence <- function(tree, clones) {
  clones <- clones[tree$tip.label]
  parts <- lapply(ape::prop.part(tree), function(ix) tree$tip.label[ix])
  parts <- c(parts, as.list(tree$tip.label))
  res <- lapply(split(names(clones), clones), function(cells) {
    f1 <- vapply(parts, function(p) 2 * length(intersect(p, cells)) / (length(p) + length(cells)), 0)
    data.frame(n = length(cells), best_f1 = max(f1),
               monophyletic = ape::is.monophyletic(tree, cells))
  })
  cbind(clone = names(res), do.call(rbind, res))
}
