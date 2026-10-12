
# Getter ####
get_b_dist = function(b_dist_func) {
  B_DISTS = list(
    "avg" = find_greedy_distance_with_avg,
    #"avg_v1" = find_greedy_distance_with_avg_v1,
    #"avg_v2" = find_greedy_distance_with_avg_v2,
    #"avg_v3" = find_greedy_distance_with_avg_v3,
    "A" = find_greedy_distance_A,
    #"A_each" = find_greedy_distance_A_each,
    "bfb" = find_greedy_distance_A_contig,
    #"B" = find_greedy_distance_B,
    #"B_each" = find_greedy_distance_B_each,
    "B_contig" = find_greedy_distance_B_contig,
    "bfb_fast" = greedy_bfb_distance_cpp
  )
  if (!is.character(b_dist_func) || length(b_dist_func) != 1 ||
      !b_dist_func %in% names(B_DISTS)) {
    stop("Unknown b_dist_func '", paste(b_dist_func, collapse = ", "),
         "'. Available: ", paste(names(B_DISTS), collapse = ", "), call. = FALSE)
  }
  B_DISTS[[b_dist_func]]
}

get_g_dist = function(g_dist_func) {
  G_DISTS = list(
    "G" = find_greedy_distance,
    "greedy" = find_greedy_distance_with_steps,
    "greedy_fast" = greedy_distance_cpp
  )
  if (!is.character(g_dist_func) || length(g_dist_func) != 1 ||
      !g_dist_func %in% names(G_DISTS)) {
    stop("Unknown g_dist_func '", paste(g_dist_func, collapse = ", "),
         "'. Available: ", paste(names(G_DISTS), collapse = ", "), call. = FALSE)
  }
  G_DISTS[[g_dist_func]]
}
