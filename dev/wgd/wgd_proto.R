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
