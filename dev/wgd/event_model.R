# Event-model counterparts of the WGD prototype (see event_model.cpp).
# Requires: devtools::load_all(), Rcpp::sourceCpp("dev/wgd/event_model.cpp"),
# source("dev/wgd/wgd_proto.R") (for wgd_combine(), assign_wgd_counts(), ...).

# Genome-wide components under the event model, same layout as
# wgd_components(): S0 = sum of E(a, b), or of min(E, B) with use_bfb;
# W1, W2 = directional WGD costs.
event_components <- function(all_input_Xs, chromosomes, alleles, use_bfb = FALSE) {
  alleles <- alleles[alleles %in% names(all_input_Xs[[1]])]
  chromosomes <- chromosomes[chromosomes %in% names(all_input_Xs)]
  cells <- rownames(all_input_Xs[[chromosomes[1]]][[alleles[1]]])
  n <- length(cells)
  S0 <- W1 <- W2 <- matrix(0, n, n, dimnames = list(cells, cells))
  for (chromosome in chromosomes) for (allele in alleles) {
    X <- all_input_Xs[[chromosome]][[allele]][cells, , drop = FALSE]
    storage.mode(X) <- "double"
    m <- event_matrices(X, use_bfb)
    S0 <- S0 + m$S0
    W1 <- W1 + m$W1
    W2 <- W2 + m$W2
  }
  list(S0 = S0, W1 = W1, W2 = W2)
}

# Bottom-up ancestral reconstruction under the event model, with WGD counts
# per node from assign_wgd_counts(). For a node at count kn with children at
# counts k1, k2, the ancestor x minimises E(2^(k1-kn) x, child1) +
# E(2^(k2-kn) x, child2), which covers same-level children, one doubled
# child and two doubled children with one exact computation.
# Returns per-edge event counts, WGD edges and total events (+ w per WGD).
reconstruct_wgd_event <- function(tree, all_input_Xs, chromosomes, alleles, k_node, w = 1) {
  n_tip <- length(tree$tip.label)
  n_node <- n_tip + tree$Nnode
  edge <- tree$edge
  post <- ape::reorder.phylo(tree, "postorder")$edge
  internal <- unique(post[, 1])
  edge_events <- numeric(nrow(edge))
  root_events <- 0   # events from the diploid ancestor to the root profile
  root <- n_tip + 1

  for (chromosome in chromosomes) for (allele in alleles) {
    X <- all_input_Xs[[chromosome]][[allele]]
    if (is.null(X)) next
    tv <- if (allele == "CN") 2 else 1
    prof <- vector("list", n_node)
    for (i in seq_len(n_tip)) prof[[i]] <- as.numeric(X[tree$tip.label[i], ])
    for (node in internal) {
      kids <- post[post[, 1] == node, 2]
      if (length(kids) != 2) stop("reconstruct_wgd_event expects a binary tree")
      f <- 2^(k_node[kids] - k_node[node])
      r <- event_ancestor(prof[[kids[1]]], f[1], prof[[kids[2]]], f[2], t = tv * 2^k_node[node])
      prof[[node]] <- r$x
      ei <- match(kids, edge[, 2])
      edge_events[ei] <- edge_events[ei] + c(r$cost_p, r$cost_q)
    }
    # the root itself descends from the (possibly doubled) diploid ancestor
    root_events <- root_events +
      event_distance(rep(tv * 2^k_node[root], length(prof[[root]])), prof[[root]])
  }

  dk_edge <- k_node[edge[, 2]] - k_node[edge[, 1]]
  wgd_edges <- which(dk_edge > 0)
  list(
    edge_events = edge_events,
    root_events = root_events,
    n_wgd_edge = dk_edge,
    wgd_clades = lapply(wgd_edges, function(e) {
      child <- edge[e, 2]
      if (child <= n_tip) tree$tip.label[child] else ape::extract.clade(tree, child)$tip.label
    }),
    root_k = k_node[root],
    total_events = sum(edge_events) + root_events + w * (sum(dk_edge) + k_node[root])
  )
}
