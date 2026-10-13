# Whole-genome duplication (WGD) support for fit(wgd = TRUE).

# Number of WGDs per cell (0, 1 or 2) estimated from ploidy:
# round(log2(mean total copy number / 2)), clamped to 0..2. Total copy number
# is CN if present, otherwise A + B.
estimate_wgd_counts = function(data) {
  cn = if ("CN" %in% names(data)) data$CN else data$A + data$B
  ploidy = tapply(cn, data$cell_id, mean, na.rm = TRUE)
  k = round(log2(pmax(ploidy, 1e-9) / 2))
  stats::setNames(as.integer(pmin(pmax(k, 0), 2)), names(ploidy))
}

# Check user-supplied WGD counts against the cells of the fit
validate_wgd_counts = function(wgd_counts, cells) {
  if (is.null(names(wgd_counts))) {
    stop("wgd_counts must be a named vector (names = cell ids).", call. = FALSE)
  }
  missing = setdiff(cells, names(wgd_counts))
  if (length(missing)) {
    stop("wgd_counts is missing ", length(missing), " cell(s), e.g. ",
         paste(utils::head(missing, 3), collapse = ", "), call. = FALSE)
  }
  k = wgd_counts[cells]
  if (any(is.na(k)) || any(!k %in% 0:2)) {
    stop("wgd_counts must be 0, 1 or 2 for every cell.", call. = FALSE)
  }
  stats::setNames(as.integer(k), cells)
}

# Parsimony assignment of WGD counts (0..kmax) to every node of a rooted tree.
# Tips are fixed to `tip_counts` (named by tip label). Along an edge the count
# can only increase, at cost `penalty` per WGD; the root's parent has 0.
# Ties are resolved towards the smaller count (later WGD).
# Returns integer counts for nodes 1..(Ntip + Nnode), in ape numbering.
assign_wgd_counts = function(tree, tip_counts, penalty = 1, kmax = 2) {
  n_tip = length(tree$tip.label)
  n_node = n_tip + tree$Nnode
  ks = 0:kmax
  cost = matrix(Inf, n_node, kmax + 1)
  cost[cbind(seq_len(n_tip), tip_counts[tree$tip.label] + 1)] = 0
  child_cost = function(child, k_parent) {
    ok = ks >= k_parent
    min(cost[child, ok] + penalty * (ks[ok] - k_parent))
  }
  post = ape::reorder.phylo(tree, "postorder")$edge
  for (node in unique(post[, 1])) {
    kids = post[post[, 1] == node, 2]
    for (k in ks) cost[node, k + 1] = sum(vapply(kids, child_cost, numeric(1), k_parent = k))
  }
  root = n_tip + 1
  k_node = integer(n_node)
  k_node[root] = ks[which.min(cost[root, ] + penalty * ks)]
  pre = ape::reorder.phylo(tree, "cladewise")$edge
  for (i in seq_len(nrow(pre))) {
    parent = pre[i, 1]
    child = pre[i, 2]
    ok = ks >= k_node[parent]
    total = cost[child, ok] + penalty * (ks[ok] - k_node[parent])
    k_node[child] = ks[ok][which.min(total)]
  }
  k_node
}

# WGD-labelled edges of a tree: one row per edge where the WGD count
# increases, with the number of WGDs and the descendant cells.
wgd_edges = function(tree, k_node) {
  n_tip = length(tree$tip.label)
  dk = k_node[tree$edge[, 2]] - k_node[tree$edge[, 1]]
  idx = which(dk > 0)
  cells = lapply(idx, function(e) {
    child = tree$edge[e, 2]
    if (child <= n_tip) tree$tip.label[child] else ape::extract.clade(tree, child)$tip.label
  })
  tibble::tibble(
    edge = idx,
    parent = tree$edge[idx, 1],
    child = tree$edge[idx, 2],
    n_wgd = dk[idx],
    n_cells = lengths(cells),
    cells = cells
  )
}
