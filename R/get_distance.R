
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
  B_DISTS[[b_dist_func]]
}

get_g_dist = function(g_dist_func) {
  G_DISTS = list(
    "G" = find_greedy_distance,
    "greedy" = find_greedy_distance_with_steps,
    "greedy_fast" = greedy_distance_cpp
  )
  G_DISTS[[g_dist_func]]
}
