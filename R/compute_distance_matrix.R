
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
compute_greedy_distances = function(all_input_Xs,
                                    chromosomes = c(1:22, "X", "Y"),
                                    alleles = c("CN", "A", "B"),
                                    target_vals = c(2, 1, 1),
                                    g_dist_func = "G") {
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
      compute_distance_matrix(input_X, dist_func, target_val = target_vals[which(names(Xchr) == allele)])
    })
    names(Ds) = alleles
    Ds
  })
  names(all_Ds) = chromosomes
  all_Ds
}

compute_avg_distances = function(all_input_Xs,
                                 chromosomes = c(1:22, "X", "Y"),
                                 alleles = c("CN", "A", "B"),
                                 bfb_penalty = 0,
                                 b_dist_func = "avg") {
  # Get the distance function from B_DISTS
  dist_func = get_b_dist(b_dist_func)
  # if (!b_dist_func %in% names(B_DISTS)) {
  #   stop(paste("Invalid B distance function:", b_dist_func,
  #              "\nAvailable options:", paste(names(B_DISTS), collapse = ", ")))
  # }
  # dist_func = B_DISTS[[b_dist_func]]

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

compute_distance_matrix = function(input_X, dist_func, ...) {
  sequences = lapply(rownames(input_X), function(r) {input_X[r,]})
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
