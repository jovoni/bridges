sum_across_alleles = function(Ds, alleles = c("A", "B")) {
  all_summed_Ds = lapply(names(Ds), function(chromosome) {
    Ds_chr = Ds[[which(names(Ds) == chromosome)]]
    Ds_to_sum = lapply(alleles, function(allele) {
      Ds_chr[[allele]]
    })
    Reduce('+', Ds_to_sum)
  })
  names(all_summed_Ds) = names(Ds)
  all_summed_Ds
}

get_target_val = function(allele) {
  if (allele == "CN") return(2)
  if (allele %in% c("A", "B")) return(1)
  stop()
}

find_minimal_distances = function(Ds_1, Ds_2) {
  alleles = names(Ds_1[[1]])
  chromosomes = names(Ds_1)

  all_min_Ds = lapply(chromosomes, function(chromosome) {
    D1_chr = Ds_1[[which(names(Ds_1) == chromosome)]]
    D2_chr = Ds_2[[which(names(Ds_2) == chromosome)]]
    min_Ds = lapply(alleles, function(allele) {
      D1 = D1_chr[[allele]]
      D2 = D2_chr[[allele]]
      pmin(D1, D2, na.rm = FALSE)
    })
    names(min_Ds) = alleles
    min_Ds
  })
  names(all_min_Ds) = chromosomes
  all_min_Ds
}

# Enhanced distance computation functions
compute_greedy_distances = function(
  all_input_Xs,
  chromosomes,
  alleles,
  g_dist_func
) {
  # Get the distance function from G_DISTS
  dist_func = get_g_dist(g_dist_func)
  # if (!g_dist_func %in% names(G_DISTS)) {
  #   stop(paste("Invalid G distance function:", g_dist_func,
  #              "\nAvailable options:", paste(names(G_DISTS), collapse = ", ")))
  # }
  # dist_func = G_DISTS[[g_dist_func]]

  alleles = alleles[alleles %in% names(all_input_Xs[[1]])]
  chromosomes = chromosomes[chromosomes %in% names(all_input_Xs)]

  all_Ds = lapply(chromosomes, function(chromosome) {
    Xchr = all_input_Xs[[which(names(all_input_Xs) == chromosome)]]
    Ds = lapply(alleles, function(allele) {
      input_X = Xchr[[which(names(Xchr) == allele)]]
      compute_distance_matrix(
        input_X,
        dist_func,
        target_val = get_target_val(allele)
      )
    })
    names(Ds) = alleles
    Ds
  })
  names(all_Ds) = chromosomes
  all_Ds
}

compute_avg_distances = function(
  all_input_Xs,
  chromosomes,
  alleles,
  bfb_penalty,
  b_dist_func
) {
  # Get the distance function from B_DISTS
  dist_func = get_b_dist(b_dist_func)

  alleles = alleles[alleles %in% names(all_input_Xs[[1]])]
  chromosomes = chromosomes[chromosomes %in% names(all_input_Xs)]

  all_Ds = lapply(chromosomes, function(chromosome) {
    Xchr = all_input_Xs[[which(names(all_input_Xs) == chromosome)]]
    Ds = lapply(alleles, function(allele) {
      input_X = Xchr[[which(names(Xchr) == allele)]]
      compute_distance_matrix(input_X, dist_func, penalty = bfb_penalty)
    })
    names(Ds) = alleles
    Ds
  })
  names(all_Ds) = chromosomes
  all_Ds
}


# Helper functions for users
# list_available_distances = function() {
#   list(
#     G_functions = names(G_DISTS),
#     B_functions = names(B_DISTS)
#   )
# }

# Convenience function to show all possible combinations
# show_distance_combinations = function() {
#   g_funcs = names(G_DISTS)
#   b_funcs = names(B_DISTS)
#
#   cat("Available G:B distance function combinations:\n")
#   cat("===========================================\n")
#   for (g in g_funcs) {
#     for (b in b_funcs) {
#       cat(sprintf("G: %-15s | B: %s\n", g, b))
#     }
#   }
#   cat("\nUsage example:\n")
#   cat('result = main(data, g_dist_func = "G_with_steps", b_dist_func = "A_each")\n')
# }

# Pairwise distances for one cells x segments matrix in a single C++ call
# (distance_matrices_cpp()), computed once per pair of distinct profiles.
# Gives the same values as compute_distance_matrix() with
# greedy_distance_cpp (G) and greedy_bfb_distance_cpp (B):
#   $min    = pmin(G, B), or G if use_bfb = FALSE
#   $greedy = G, $bfb = B    (only when keep_parts = TRUE)
# Diagonals are Inf.
fast_distance_matrices = function(input_X, use_bfb = TRUE, penalty = 0, keep_parts = FALSE) {
  storage.mode(input_X) = "double"
  res = distance_matrices_cpp(input_X, use_bfb, penalty, keep_parts)
  cells = rownames(input_X)
  for (k in c("min", "greedy", "bfb")) {
    if (!is.null(res[[k]])) dimnames(res[[k]]) = list(cells, cells)
  }
  res
}

# Total distance matrix for fit()'s default distances ("greedy_fast" and
# "bfb_fast" or NULL): sum over chromosomes and alleles of pmin(G, B),
# accumulated one chromosome/allele at a time so that only the total is
# kept in memory. Per-chromosome G and B matrices are returned only when
# keep_matrices = TRUE, in the same layout as compute_greedy_distances() /
# compute_avg_distances(). Rows are aligned by cell name across chromosomes.
#
# With wgd = TRUE, the directional WGD costs W1 (one doubling) and W2 (two
# doublings) are accumulated too (wgd_matrices_cpp()), and the distance is
# the cheapest explanation genome-wide:
#   D = min(D, w + W1, w + t(W1), 2w + W2, 2w + t(W2)),   w = wgd_penalty
# The WGD is charged once for the whole genome, while the pre-doubling
# ancestor is chosen per chromosome/allele. `D_no_wgd` keeps the distance
# without the WGD option.
#
# WGD costs are only computed from a cell to cells with more WGDs, according
# to `wgd_counts` (named by cell): W1 (one doubling) when the other cell has
# at least one more, W2 (two doublings) when it has two more. Other entries
# are Inf. Set
# wgd_all_pairs = TRUE to compute every ordered pair.
compute_fast_distances = function(all_input_Xs, chromosomes, alleles,
                                  bfb_penalty = 0, use_bfb = TRUE,
                                  keep_matrices = FALSE,
                                  wgd = FALSE, wgd_penalty = 1,
                                  wgd_counts = NULL, wgd_all_pairs = FALSE) {
  alleles = alleles[alleles %in% names(all_input_Xs[[1]])]
  chromosomes = chromosomes[chromosomes %in% names(all_input_Xs)]

  cells = rownames(all_input_Xs[[chromosomes[1]]][[alleles[1]]])
  D = matrix(0, length(cells), length(cells), dimnames = list(cells, cells))
  greedy_Ds = avg_Ds = NULL
  n_unique = list()
  if (wgd) {
    W1 = W2 = matrix(0, length(cells), length(cells), dimnames = list(cells, cells))
    if (is.null(wgd_counts)) {
      if (!wgd_all_pairs) stop("wgd_counts is required unless wgd_all_pairs = TRUE.", call. = FALSE)
      wgd_counts = stats::setNames(integer(length(cells)), cells)
    }
    wgd_counts = as.integer(wgd_counts[cells])
  }

  for (chromosome in chromosomes) {
    for (allele in alleles) {
      X = all_input_Xs[[chromosome]][[allele]]
      if (!setequal(rownames(X), cells) || nrow(X) != length(cells)) {
        stop("Chromosome ", chromosome, ", allele ", allele,
             " does not contain the same cells as the other chromosomes.", call. = FALSE)
      }
      X = X[cells, , drop = FALSE]
      m = fast_distance_matrices(X, use_bfb = use_bfb, penalty = bfb_penalty,
                                 keep_parts = keep_matrices)
      D = D + m$min
      if (wgd) {
        storage.mode(X) = "double"
        w = wgd_matrices_cpp(X, wgd_counts, wgd_all_pairs)
        W1 = W1 + w$W1
        W2 = W2 + w$W2
      }
      n_unique[[chromosome]][[allele]] = m$n_unique
      if (keep_matrices) {
        greedy_Ds[[chromosome]][[allele]] = m$greedy
        avg_Ds[[chromosome]][[allele]] = if (use_bfb) m$bfb else m$greedy
      }
    }
  }

  D_no_wgd = NULL
  if (wgd) {
    D_no_wgd = D
    D = pmin(D, wgd_penalty + W1, wgd_penalty + t(W1),
             2 * wgd_penalty + W2, 2 * wgd_penalty + t(W2))
  }

  list(D = D, D_no_wgd = D_no_wgd, greedy_Ds = greedy_Ds, avg_Ds = avg_Ds,
       n_unique = n_unique)
}

compute_distance_matrix = function(input_X, dist_func, ...) {
  sequences = lapply(rownames(input_X), function(r) {
    input_X[r, ]
  })
  names(sequences) = rownames(input_X)

  nodes <- sequences
  n <- length(nodes)

  # Assign unique names to leaf nodes if they don't have any
  if (is.null(names(nodes))) {
    names(nodes) <- paste0("seq", 1:n)
  }

  # Init D
  D <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      D[j, i] = D[i, j] = dist_func(nodes[[i]], nodes[[j]], ...)$cost
    }
  }
  diag(D) <- Inf
  rownames(D) = colnames(D) = names(nodes)
  D
}
